#include "Math/IFunction.h"
#include "Math/IFunctionfwd.h"
#include "Rtypes.h"
#include "TMath.h"
#include "TMinuitMinimizer.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TGraph.h"
#include "TAxis.h"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>


class LikelihoodFitter : public ROOT::Math::IMultiGenFunction{
private: 
  int _npar;
  double _total_rate;
  std::vector<int> _counts;
  std::vector<double> _time;
public:
  LikelihoodFitter() : ROOT::Math::IMultiGenFunction(), _npar(0) {}
  unsigned int NDim() const { return _npar-1;}
  LikelihoodFitter* Clone() const {
    LikelihoodFitter* f = new LikelihoodFitter;
    for(int i = 0; i<_npar; i++){
      f->AddSource(_counts[i],_time[i]);
    }
    f->SetTotalRate(_total_rate);
    return f;
  }
  void SetTotalRate(double r) { _total_rate = r; }
  double GetTotalRate(){ return _total_rate;}
  double GetExpectedRate(){
    double r = 0;
    for(int i=0;i<_npar; i++){
      r += _counts[i]/_time[i];
    }
    return r;
  }
  int AddSource(int counts, double time){
    _counts.push_back(counts);
    _time.push_back(time);
    return ++_npar;
  }
  int GetCounts(size_t i){ return _counts[i]; }
  int GetTime(size_t i){ return _time[i]; }
  double DoEval(const double* x) const;
};

double LikelihoodFitter::DoEval(const double* x) const
{
  //x are the parameters from MINUIT
  double likelihood = 0;
  //want the -loglikelyhood
  //each source contributes to likelihood from poisson distribution
  double sum=0;
  for(int i=0; i<_npar-1; i++){
    likelihood -= log(TMath::Poisson(_counts[i],x[i]*_time[i]));
    //cout<<x[i]<<"--------"<<i<<"------=--------"<<likelihood<<endl;
    sum += x[i];
  }
  //the last one must be equal to the difference 
  double last = _total_rate - sum;
  if(last < 0){
    //uh-oh, we're outside our domain! add some huge number...
    likelihood += 1000000;
  }
  else{
    likelihood -= log(TMath::PoissonI(_counts[_npar-1],last*_time[_npar-1]));
  }
  return likelihood;

}
  
TGraph* rate_likelihood(int npts=200, double rmin=0, double rmax = 0.2)
{
  vector<double> x;
  vector<double> y;
  
  LikelihoodFitter lf;
  lf.AddSource(9,1500);
  lf.AddSource(1,15./0.00008);
  //lf.AddSource(0,88);
  lf.AddSource(0,48);
  //lf.AddSource(20,44);

  const char* names[] = {"radio","nveto","mveto","external"};
 
  double maxlike = 0;
  double bestrate = 0;
  double integral = 0;
  
  //ROOT::Minuit2::Minuit2Minimizer minim("Simplex");
  TMinuitMinimizer minim("simplex");
  minim.SetFunction(lf);
  double minval = 0;
  for(int i=0; i<npts; i++){
    double rate = rmin + i*(rmax-rmin)/(npts-1);
    lf.SetTotalRate(rate);
    if(rate==0){
      //handle this specially
      std::vector<double> rates(lf.NDim());
      minval = lf.DoEval(&rates[0]);
    }
    else{
      for(size_t j=0; j<lf.NDim(); j++){
	minim.Clear();
	stringstream name("");
	name<<"rate"<<j;
	double base = (1.*lf.GetCounts(j))/lf.GetTime(j);
	//limit the rate to total rate or 10+sqrt(n)/t, whichiver is higher
	double limit = (10+5*sqrt(lf.GetCounts(j)))/lf.GetTime(j);
	if(limit > rate) limit = rate;
	minim.SetLimitedVariable(j,names[j],base,0.01,0,limit);
      }
      //minim.SetPrintLevel(3);
      minim.SetMaxFunctionCalls(1000);
      minim.SetMaxIterations(1000);
      if(!minim.Minimize()){
	minim.SetMaxFunctionCalls(5000);
	minim.SetMaxIterations(5000);
	if(!minim.Minimize())
	  continue;
      }
	
      //minim.PrintResults();
      minval = minim.MinValue();
    }
    x.push_back(rate);
    y.push_back( exp(-minval));
    //cout<<rate<<" "<<y[i]<<endl;
    if(y.back() > maxlike){
      maxlike = y.back();
      bestrate = rate;
    }
    if(x.size()>1){
      size_t last = x.size()-1;
      integral += (x[last]-x[last-1])*0.5*(y[last]+y[last-1]);
    }
  }
  cout<<"Maximum Likelihood is "<<maxlike<<" at a rate of "<<bestrate<<endl;
  cout<<"Expected best rate is "<<lf.GetExpectedRate()<<endl;
  //normalize everything
  double running=0;
  bool found90=0, found95=0;
  npts = x.size();
  for(int i=0; i<npts; i++){
    y[i] /= integral;
    if(i>0){
      running += (x[i]-x[i-1])*0.5*(y[i]+y[i-1]);
    }
    if(!found90 && running > 0.9){
      cout<<"90% Upper CL rate = "<<x[i]<<endl;
      found90 = true;
    }
    if(!found95 && running > 0.95){
      cout<<"95% Upper CL rate = "<<x[i]<<endl;
      found95 = true;
    }    
  }
  TGraph* g = new TGraph(npts,&x[0],&y[0]);
  g->Draw("alp");
  g->SetTitle("Normalized Likelihood for neutron rate");
  g->GetYaxis()->SetTitle("Normalized likelihood");
  g->GetXaxis()->SetTitle("Background rate in fiducial volume [cts/year/33 kg]");
  
  //now try to determine the prob for 1+ events
  //prob for 0 is convolution of e^-r p(r), where r is above + ar39 rate, const
  double ar39 = 800*365.25*33 * 1.E-9;
  ar39 = 0.1;
  double p0 = 0;
  for(int i=0; i<npts-1; i++){
    double ravg = (x[i] + ar39)*3;
    p0 += exp(-ravg)*y[i]*(x[i+1]-x[i]);
  }
  cout<<"The \"probability\" to receive 1 or more events in 3 years is "
      <<1-p0<<endl;
  return g;

}

#include <TMath.h>

double acceptance(double* x, double* par)
{
  double n = x[0];
  double fmean = par[0];
  double fcut = par[1];
  double k = n*fcut+1;
  
  double a = n-(k-1);
  double b = (k-1)+1;
  double c = 1-fmean;
  
  //return TMath::BinomialI(fmean,n,k);
  return TMath::BetaIncomplete(fmean,k,n-k+1);
  //return 1-TMath::BetaIncomplete(c,a,b);
}

double acceptance2(double* x, double* par)
{
  int n = x[0];
  double fmean = par[0];
  double fcut = par[1];
  int k = n*fcut+1;
  
  int a = n-k;
  int b = k+1;
  double c = 1-fmean;
  
  double sum=0;
  for(int i=k; i<=n; i++){
    sum += TMath::Binomial(n,i)*pow(fmean,i)*pow(1-fmean,n-i);
  }
  return sum;

}


TF1* betafrac(double fmean=0.27, double fcut = 0.7)
{
  double high=200, low=10;
  
  TF1* f1 = new TF1("fbeta",acceptance,low,high,2);
  f1->SetParNames("fmean","fcut");
  f1->SetParameters(fmean,fcut);
  
  double min=1.e-12, max = 1.e-1;
  if(fmean > fcut){
    min = 0.1; 
    max = 1.1;
  }
  
  TH2F* axes = new TH2F("axes","Ideal beta-like acceptance fraction",
			100,low,high,
			100,min,max);
  axes->Draw();
  axes->GetXaxis()->SetTitle("detected photoelectrons");
  axes->GetYaxis()->SetTitle("beta acceptance fraction");
  axes->GetYaxis()->CenterTitle();
  axes->GetYaxis()->SetTitleOffset(1.3);
  int nbins = high-low;
  
  //f1->SetNpx(nbins*20);
  //f1->SetLineColor(kRed);
  f1->SetTitle("fcut=0.7");
  f1->DrawCopy("same");
  f1->SetParameter("fcut",0.6);
  f1->SetLineStyle(9);
  f1->SetTitle("fcut=0.6");
  f1->DrawCopy("same");
  f1->SetParameter("fcut",0.5);
  f1->SetLineStyle(7);
  f1->SetTitle("fcut=0.5");
  f1->DrawCopy("same");
  gPad->SetLogy();
  gPad->BuildLegend();
  
  return f1;  

}

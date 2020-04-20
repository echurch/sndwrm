#include "TMath.h"
#include "TF1.h"
#include "TFeldmanCousins.h"
#include "TRolke.h"
#include "TGraph.h"

const double cm = 1;
const double m = 100*cm;
const double km = 1000*m;
const double fm = 1.e-15*m;

const double cm2 = cm*cm;
const double barn = 1.e-24*cm2;
const double pb = 1.e-12*barn;
const double cm3 = cm*cm*cm;

const double kg = 1;
const double ton = 1000*kg;

const double s = 1;
const double day = 3600*24*s;
const double year = 365.24*day;

const double keV = 1;
const double MeV = 1000 * keV;
const double GeV = 1000 * MeV;

const double c = 3.e8 * m/s;
const double pi = TMath::Pi();
const double N0 = 6.02e26 / kg;
const double hc = 197 * MeV * fm;
const double Mn = 932 * MeV;
const double ANa = 23;
const double ASi = 28.085;
const double AAr = 39.948;
const double AGe = 72.64;
const double AXe = 131.293;
const double ar39_atmo = 192./day/kg/keV;
const double skin = 0.9 * fm;

double HelmFactor(double Er, double A)
{
  double Mt = Mn* A;
  double q = sqrt(2*Mt*Er);
  //double r = sqrt( pow(1.23*pow(A,1./3.)-0.6,2)*fm*fm+7./3.*pi*pi*0.52*0.52*fm*fm-5*skin*skin);
  double r = 1.14*pow(A,1./3.)*fm;
  double qr = q*r/hc;
  double qs = q*skin/hc;
  return 3.*(sin(qr) - qr*cos(qr))/pow(qr,3) * exp(-qs * qs / 2.);
}

double VelocityFunction(double v, double vesc, double v0, double ve)
{
  double k = 1./(pow(pi,1.5)*v0*v0*v0 * (erf(vesc/v0) - 2*vesc/(sqrt(pi)*v0)*exp(-vesc*vesc/(v0*v0)) ) );
  double prefactor = pi*k*v0*v0*v/ve;
  double dist = prefactor;
  if(v <= vesc-ve)
    dist *= (exp(-(v-ve)*(v-ve)/(v0*v0)) - exp(-(v+ve)*(v+ve)/(v0*v0)) );
  else if(v <= vesc+ve)
    dist *= (exp(-(v-ve)*(v-ve)/(v0*v0)) - exp(-vesc*vesc/(v0*v0)));
  else dist = 0;
  return dist;
}

double VelocityIntegral(double vmin, double vesc, double v0, double ve)
{
  double k = 1./(pow(pi,1.5)*v0*v0*v0 * (erf(vesc/v0) - 2*vesc/(pow(pi,0.5)*v0)*exp(-vesc*vesc/(v0*v0)) ) );
  double prefactor = pow(pi, 1.5)* v0*v0*v0*k/(2*ve);
  
  double rate = prefactor;
  if(vmin <= vesc - ve)
    rate *= (erf((ve-vmin)/v0) + erf((ve+vmin)/v0) - 
	     4*ve/(sqrt(pi)*v0)*exp(-vesc*vesc/(v0*v0)));
  else if(vmin <= vesc+ve)
    rate *= (erf((ve-vmin)/v0) + erf((vesc)/v0) - 
	     2.*(ve+vesc-vmin)/(sqrt(pi)*v0)*exp(-vesc*vesc/(v0*v0)));
  else
    rate = 0;
  return rate;  
}

double DifferentialWimpRate(double Er, double sigma, double rho, 
			    double A, double Z, double fpfn, double Md,  
			    double vesc, double v0, double ve)
{
  double Mt = Mn*A;
  double mu = Mt*Md/(Mt+Md);
  double mun = Mn*Md/(Mn+Md);
  double F = HelmFactor(Er, A);
  double vmin = c*sqrt(Mt*Er/(2.*mu*mu));
  return N0*sigma*rho*Mt/(2.*A*mun*mun*Md) * pow(Z*fpfn+(A-Z),2)*F*F * 
    VelocityIntegral(vmin,vesc,v0,ve)*c*c;

}
/*
class IntegralAboveThreshold{
public:
  TF1* func;
  double max;
  double operator()(double* x, double*){ return func->Integral(*x,max); }
  IntegralAboveThreshold(TF1* f, double Max) : func(f), max(Max) {}
  TF1* MakeFunction(double Min, double Max){
    return new TF1("intabovethresh",*this,Min,Max,0);
  }
  ClassDef(IntegralAboveThreshold,1);
};
*/

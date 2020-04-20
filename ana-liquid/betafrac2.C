double acceptance(double n, double fmean, double fcut)
{
  //double n = x[0];
  //double fmean = par[0];
  //double fcut = par[1];
  double k = n*fcut; //+1 ?
  
  double a = n-(k-1);
  double b = (k-1)+1;
  double c = 1-fmean;
  
  //return TMath::BinomialI(fmean,n,k);
  return TMath::BetaIncomplete(fmean,k,n-k+1);
  //return 1-TMath::BetaIncomplete(c,a,b);
}

/*
void betafrac2()
{
  double fbeta = 0.27, frecoil = 0.71;
  TF1* f1 = new TF1("acc","acceptance([0],[1],x)",0,1);
  TH2F* axes = new TH2F("axes","axes",100,0,1,100,1.e-12,1.1);
  axes->Draw();
  for(int n=10; n<60; n+=20){
    f1->SetParameters(n,fbeta);
    f1->SetLineColor(kRed);
    f1->DrawCopy("same");
    f1->SetParameters(n,frecoil);
    f1->SetLineColor(kBlue);
    f1->DrawCopy("same");
  }
}
*/

double acceptance2(double accrecoil, double n, 
		   double fbeta=0.27, double frecoil=0.71)
{						
  TF1* f = (TF1*)(gROOT->GetFunction("acc"));
  if(!f)
    f = new TF1("acc","acceptance([0],[1],x)",0,1);
  f->SetParameters(n,frecoil);
  double fcut = f->GetX(accrecoil);
  return acceptance(n,fbeta,fcut);
}

void betafrac2()
{
  TH2F* axes = new TH2F("axes","axes",100,0.5,1,100,1.e-11,1);
  axes->Draw();
  axes->GetXaxis()->SetTitle("recoil acceptance");
  axes->GetYaxis()->SetTitle("beta acceptance");
  axes->GetYaxis()->CenterTitle();
  axes->GetYaxis()->SetTitleOffset(1.3);
  axes->SetTitle("beta leakage vs recoil acceptance");
  gPad->SetLogy();
  for(int n=10; n<60; n+=20){
    TF1* f1 = new TF1("f1","acceptance2(x,[0])",0.5,0.99);
    f1->SetParameter(0,n);
    char name[50];
    sprintf(name,"n=%d",n);
    f1->SetTitle(name);
    f1->DrawCopy("same");
  }
  gPad->BuildLegend();
}

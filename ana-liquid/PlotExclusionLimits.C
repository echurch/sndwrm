{
  gROOT->Reset();
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine(".L WimpSpectrum.C+");
  //  gROOT->ProcessLine(".x StandardVars.C");
  double sigma = 1.e-45*cm2;
  double rho = 0.3*GeV/cm3;
  double vesc = 544*km/s;
  double v0 = 250*km/s;
  double ve = 263*km/s;
  double Md = 100*GeV;

  gStyle->SetOptStat(000000);
  gStyle->SetTextSize(0.04);
  
  TH2F* axes = new TH2F("axes","DM 90% sensitivities",
			//			100,5,1000,100,1e-46,1e-39);
			//			100,0.01,10,100,1e-50,1e-41);
			100,0.01,10,100,1e-49,1e-42);

  axes->Draw("AXIS");
  gPad->SetLogx();
  gPad->SetLogy();
  axes->GetXaxis()->SetTitle("WIMP mass [TeV]");
  axes->GetXaxis()->CenterTitle();
  axes->GetXaxis()->SetTitleOffset(1.3);
  axes->GetYaxis()->SetTitle("Cross Section [cm^{2}]");
  axes->GetYaxis()->CenterTitle();
  axes->GetYaxis()->SetTitleOffset(1.3);
  /*

  axes->SetTitle("DarkSide-50 3 year Sensitivity");
  gStyle->SetTextSize(0.04);
  //draw mssm
  TGraph* g2 = new TGraph("buchmueller.png.dat");
  g2->SetFillColor(kGray+1);
  g2->Draw("f");
  TGraph* g3 = new TGraph("buchmueller.png.dat2");
  g3->SetFillColor(kGray+2);
  g3->Draw("f");
  TText* pmssm = new TText(450,2e-45,"MSSM");
  pmssm->SetTextColor(kGray+2);
  pmssm->Draw();


  
  //draw cogent
  TGraph* gcogent = new TGraph("cogent.dat");
  gcogent->SetLineColor(kGreen);
  gcogent->SetLineWidth(2);
  gcogent->Draw("l");
  TText* pcogent = new TText(13, 3e-41,"CoGeNT");
  pcogent->SetTextColor(kGreen);
  pcogent->Draw();

  //draw dama
  TGraph* gdama1 = new TGraph("dama1.dat");
  gdama1->SetLineWidth(2);
  gdama1->SetLineColor(kMagenta);
  gdama1->Draw("l");
  TText* pdama1 = new TText(14,2e-40,"DAMA/Na");
  pdama1->SetTextColor(kMagenta);
  pdama1->Draw();
  
  TGraph* gdama2 = new TGraph("dama2.dat");
  gdama2->SetLineWidth(2);
  gdama2->SetLineColor(kMagenta);
  gdama2->Draw("l");
  TText* pdama2 = new TText(100,1e-41,"DAMA/I");
  pdama2->SetTextColor(kMagenta);
  pdama2->Draw();
  
  //draw xenon
  TGraph* gxenon1 = new TGraph("xenon.dat");
  gxenon1->SetLineColor(kBlue);
  gxenon1->SetLineWidth(2);
  gxenon1->SetLineStyle(7);
  gxenon1->Draw("l");
  TText* pxenon = new TText(400,1.8e-44,"Xenon100");
  pxenon->SetTextColor(kBlue);
  pxenon->Draw();
  
  //draw cdms
  TGraph* gcdms1 = new TGraph("cdms.dat");
  gcdms1->SetLineColor(kRed);
  gcdms1->SetLineWidth(2);
  gcdms1->SetLineStyle(9);
  gcdms1->Draw("l");
  TText* pcdms = new TText(400,3e-43,"CDMSII");
  pcdms->SetTextColor(kRed);
  pcdms->Draw();
  */




  //draw LUX result
  /*
  TGraph* gdama1 = new TGraph("lux2017.csv");
  gdama1->SetLineWidth(2);
  gdama1->SetLineColor(kGreen);
  gdama1->Draw("l");
  TText* pdama1 = new TText(1.,2e-45,"LUX LIMIT");
  pdama1->SetTextColor(kGreen);
  pdama1->Draw();
  */

  // DS50 result
  //  TGraph* gdama2 = new TGraph("deep3600proj.csv");
  TGraph* gds50 = new TGraph("DS50_Limit.txt");
  gds50->SetLineWidth(2);
  gds50->SetLineColor(kOrange);
  gds50->Draw("l");
  TText* gds50t = new TText(1.9,1.3e-43,"Darkside50 LIMIT");
  gds50t->SetTextColor(kOrange);
  gds50t->SetTextAngle(13);
  gds50t->Draw();

  // DEAP result
  //  TGraph* gdama2 = new TGraph("deep3600proj.csv");
  TGraph* gdama2 = new TGraph("DEAP_3600_Limit.txt");
  gdama2->SetLineWidth(2);
  gdama2->SetLineColor(kMagenta);
  gdama2->Draw("l");
  TText* pdama2 = new TText(2.,1e-44,"DEAP3600 LIMIT");
  pdama2->SetTextColor(kMagenta);
  pdama2->SetTextAngle(13);
  pdama2->Draw();

  //  TGraph* gdama2 = new TGraph("deep3600proj.csv");
  TGraph* xenon1t = new TGraph("XENON1T_Limit.txt");
  xenon1t->SetLineWidth(2);
  xenon1t->SetLineColor(kGreen);
  xenon1t->Draw("l");
  TText* xenon1tt = new TText(1.0,6e-46,"XENON1T LIMIT");
  xenon1tt->SetTextColor(kGreen);
  xenon1tt->SetTextAngle(13);
  xenon1tt->Draw();
  
  //Darkside20k future
  TGraph* gxenon1 = new TGraph("darkside_200tyr.csv");
  gxenon1->SetLineColor(kBlue);
  gxenon1->SetLineWidth(2);
  gxenon1->SetLineStyle(7);
  gxenon1->Draw("l");
  TText* pxenon = new TText(0.0108,1e-42,"Darkside20k - 200 ton-yr");
  pxenon->SetTextColor(kBlue);
  pxenon->SetTextAngle(-78);
  pxenon->Draw();
  
  //GADMC future, Threshold 30 keV.
  TGraph* gcdms1 = new TGraph("gadmc_3ktyr.csv");
  gcdms1->SetLineColor(kRed);
  gcdms1->SetLineWidth(2);
  gcdms1->SetLineStyle(9);
  gcdms1->Draw("l");
  TText* pcdms = new TText(1.0,6e-49,"ARGO 3 kton-yr");
  pcdms->SetTextColor(kRed);
  pcdms->SetTextAngle(13);
  pcdms->Draw();

  //draw Nu Floor
  /*
  TGraph* gcdms2 = new TGraph("floor.csv");
  gcdms2->SetLineColor(kGray);
  gcdms2->SetLineWidth(2);
  gcdms2->SetLineStyle(9);
  gcdms2->Draw("l");
  TText* pcdms2 = new TText(20/1000.,2e-50,"Neutrino Floor");
  pcdms2->SetTextColor(kGray);
  pcdms2->Draw();
  */



  TRolke rolke;
  TFeldmanCousins fc;
  fc.SetMuMax(50000.); ///// For large bkgd
  fc.SetMuStep(1.);
  rolke.SetPoissonBkgKnownEff(0,0,20,0.5);  // N_obs, N_obs_bkgdregion, sig/bkgd (in PSD in MC), eff 
  double ul;
  ul = fc.CalculateUpperLimit(100,100);

  double counts = rolke.GetUpperLimit();
  ////counts = 4.6; // assumes a 50% efficiency from PSD
  ////counts = 7.8; // 2,2 FC and assumes a 50% efficiency from PSD
  counts = ul*2.0; // assumes a 50% efficiency from PSD
  std::cout << "ul is " << ul << std::endl;
  ////  double masstime = 100*kg*year; // Darkside
  double masstime = 3000000*kg*year ; // DUNE
  ////double masstime = 200000*kg*year ; // Darkside

  ul = fc.CalculateUpperLimit(1+10,1+10); // DUNE 100 keV 1 neutrons + 10 atmos nus
  counts100 = ul*2.0; // assumes a 50% efficiency from PSD

  ul = fc.CalculateUpperLimit((14.+1)+13+1,(14.+1)+13+1); // DUNE 75 keV neutrons (ss+rock) + atmos nus  + 1 leakage e
  counts75 = ul*2.0; // assumes a 50% efficiency from PSD

  ul = fc.CalculateUpperLimit(1,1); 
  counts5 = ul*2.0; // assumes a 50% efficiency from PSD


  ul = fc.CalculateUpperLimit((15+13)+25+10,(15+13)+25+10);   // DUNE 50 keV
  ul = fc.CalculateUpperLimit((2+13)+17+1,(2+13)+17+1);   // DUNE 50 keV
  counts9 = ul*2.0; // assumes a 50% efficiency from PSD


  //double counts = 2.3;
  rho = 0.3*GeV/cm3;
  v0 = 220*km/s;
  ve = 244*km/s;


  //vesc = 544*km/s;





  //  double emin = 30, emax=200; // Darkside50
  double emin = 100, emax=500; // DUNE!
  double A=AAr;
  ////emin = 30 ; // darkside20k

  cout<<"The 90% upper CL on mean counts is "<<counts<<endl;
  TF1* f1 = new TF1("f1","DifferentialWimpRate(x*[0],[1],[2],[3],[4],[5],[6],[7],[8],[9])",0,200);
  
  f1->SetParameters(keV,sigma,rho,A,18,1,Md,vesc,v0,ve);
  


  vector<double> x,y;
  /*
  for(double mass = 1; mass<10000; mass *=1.01){
    f1->SetParameter(6,mass*GeV);
    double rate = masstime*f1->Integral(emin,emax);
    if(rate==0) continue;
    double limit = sigma*counts5/rate;
    //cout<<mass<<" "<<rate<<" "<<limit<<endl;
    x.push_back(mass/1000);
    y.push_back(limit);
  }
  
  TGraph* g = new TGraph(x.size(), &x[0], &y[0]);
  g->Draw("l");
  g->SetLineWidth(4);
  //  TText* t = new TText(22,1.2e-45,"DarkSide-50");

  TText* t = new TText(0.035,4e-42,"DUNE Module 4: 3 kt-yrs, 50% eff. "); //+std::itoa(emin));
  TText* t2 = new TText(0.6,0.9e-48,"100 keV thresh: 1 bgkd evt "); //+std::itoa(emin));
  t2->SetTextAngle(13);
  t->Draw();
  t2->Draw();

  */

  x.clear(); y.clear();
  for(double mass = 1; mass<10000; mass *=1.01){
    f1->SetParameter(6,mass*GeV);
    double rate = masstime*f1->Integral(emin,emax);
    if(rate==0) continue;
    double limit = sigma*counts100/rate;
    //cout<<mass<<" "<<rate<<" "<<limit<<endl;
    x.push_back(mass/1000);
    y.push_back(limit);
  }
  
  TGraph* g100 = new TGraph(x.size(), &x[0], &y[0]);
  g100->Draw("l");
  g100->SetLineWidth(4);
  //  TText* t = new TText(22,1.2e-45,"DarkSide-50");

  TText* t100 = new TText(0.6,5.6e-48,"100 keV thresh: 11 bkd"); //+std::itoa(emin));
  t100->SetTextSize(0.03);
  t100->SetTextAngle(15);
  t100->Draw();



  /*
  emin = 100.;
  x.clear(); y.clear();
  for(double mass = 1; mass<10000; mass *=1.01){
    f1->SetParameter(6,mass*GeV);
    double rate = masstime*f1->Integral(emin,emax);
    if(rate==0) continue;
    double limit = sigma*counts/rate;
    //cout<<mass<<" "<<rate<<" "<<limit<<endl;
    x.push_back(mass/1000);
    y.push_back(limit);
  }
  
  TGraph* g3 = new TGraph(x.size(), &x[0], &y[0]);
  g3->Draw("l");
  g3->SetLineWidth(4);

  TText* t4 = new TText(0.6,1.2e-47,"100 keV thresh: 100 bkd evts"); //+std::itoa(emin));
  t4->SetTextAngle(13);
  t4->Draw();
  gPad->Update();
  */


  emin = 75.;
  x.clear(); y.clear();
  for(double mass = 1; mass<10000; mass *=1.01){
    f1->SetParameter(6,mass*GeV);
    double rate = masstime*f1->Integral(emin,emax);
    if(rate==0) continue;
    double limit = sigma*counts75/rate;
    //cout<<mass<<" "<<rate<<" "<<limit<<endl;
    x.push_back(mass/1000);
    y.push_back(limit);
  }
  
  TGraph* g2 = new TGraph(x.size(), &x[0], &y[0]);
  g2->Draw("l");
  g2->SetLineWidth(3);

  //  TText* t3 = new TText(0.018,4e-42,"75 keV thresh: 10 bkd evts"); //+std::itoa(emin));
  // t3->SetTextAngle(-85);
  TText* t3 = new TText(0.45,1.6e-48,"75 keV thresh: 29 bkd"); //+std::itoa(emin));
  t3->SetTextSize(0.03);
  t3->SetTextAngle(15);
  t3->Draw();
  gPad->Update();

  


  emin = 50.;
  x.clear(); y.clear();
  for(double mass = 1; mass<10000; mass *=1.01){
    f1->SetParameter(6,mass*GeV);
    double rate = masstime*f1->Integral(emin,emax);
    if(rate==0) continue;
    double limit = sigma*counts9/rate;
    //cout<<mass<<" "<<rate<<" "<<limit<<endl;
    x.push_back(mass/1000);
    y.push_back(limit);
  }
  
  TGraph* g751 = new TGraph(x.size(), &x[0], &y[0]);
  g751->Draw("l");
  g751->SetLineWidth(1);

  TText* t751 = new TText(0.3,6e-49,"50 keV thresh: 33 bkd"); //+std::itoa(emin));
  t751->SetTextSize(0.03);
  t751->SetTextAngle(15);
  t751->Draw();
  gPad->Update();



  //  cout<<"Sensitivity to 50 GeV WIMP is "<<g->Eval(50)<<" cm2"<<endl;
  // cout<<"Sensitivity to 100 GeV WIMP is "<<g->Eval(100)<<" cm2"<<endl;



  
  if(0){
    //draw xenon
    rolke.SetGaussBkgGaussEff(3,1.8,0.32, 0.03*0.32, 0.6);
    counts = rolke.GetUpperLimit();
    masstime = 4848*kg*day;
    emin = 8.4; 
    emax=44.6;
    A=AXe;  
    f1->SetParameters(keV,sigma,rho,A,1,1,Md,vesc,v0,ve);
    x.clear(); y.clear();
    for(double mass = 1; mass<10000; mass *=1.01){
      f1->SetParameter(6,mass*GeV);
      double rate = masstime*f1->Integral(emin,emax);
      if(rate==0) continue;
      double limit = sigma*counts/rate;
      //cout<<mass<<" "<<rate<<" "<<limit<<endl;
      x.push_back(mass/1000);
      y.push_back(limit);
    }
    TGraph* gxe = new TGraph(x.size(), &x[0], &y[0]);
    gxe->SetName("gxe");
    gxe->SetTitle("xenon");
    gxe->SetLineWidth(2);
    gxe->SetLineStyle(7);
    gxe->SetLineColor(kBlue);
    gxe->Draw("l");
    
    //draw cdms
    rolke.SetGaussBkgGaussEff(2,0.9, 0.32, 0.02, 0.2);
    counts = rolke.GetUpperLimit();
    masstime = 612*kg*day;
    //how to get it to scale properly???
    masstime *= 2.5;  
    emin = 10; 
    emax=100;
    A=AGe;  
    f1->SetParameters(keV,sigma,rho,A,18,1,Md,vesc,v0,ve);
    x.clear(); y.clear();
    for(double mass = 1; mass<10000; mass *=1.01){
      f1->SetParameter(6,mass*GeV);
      double rate = masstime*f1->Integral(emin,emax);
      if(rate==0) continue;
      double limit = sigma*counts/rate;
      //cout<<mass<<" "<<rate<<" "<<limit<<endl;
      x.push_back(mass);
      y.push_back(limit);
    }

    /*
    TGraph* gcdms = new TGraph(x.size(), &x[0], &y[0]);
    gcdms->SetName("gcdms");
    gcdms->SetTitle("cdms");
    gcdms->SetLineWidth(2);
    gcdms->SetLineStyle(9);
    gcdms->SetLineColor(kRed);
    gcdms->Draw("l");
    */

  }

 


}

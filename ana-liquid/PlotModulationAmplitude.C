{
  gROOT->Reset();
  gROOT->ProcessLine(".L WimpSpectrum.C+");
  gROOT->ProcessLine(".x StandardVars.C");

  double A[] = {AAr, AGe, AXe};
  const char* name[] = {"Argon","Germanium","Xenon"};
  const int ngraphs = +sizeof(A)/sizeof(double);
  TMultiGraph* mg = new TMultiGraph;
  TLegend* leg = new TLegend(0.6,0.7,0.9,0.85);
  const double maxint =300;
  TF1* f1 = new TF1("f1","DifferentialWimpRate(x*[0],[1],[2],[3],[4],[5],[6],[7],[8],[9])",0,maxint);
  bool integrated = true;
  bool fractional = true;

  f1->SetParameters(keV,sigma,rho,A[0],0,1,Md,vesc,v0,ve);
  for(int i=0; i<ngraphs; i++){
    const int npts = 200;
    const double emin = 1, emax = 201;
    
    double x[npts], y[npts];
    f1->SetParameter(3,A[i]);
    for(int pt = 0; pt<npts; pt++){
      x[pt] = emin + pt*(emax-emin)/(1.*npts);
      if(integrated){
	//this is the average rate
	f1->SetParameter(9,ve);
	double avg = kg*day*f1->Integral(x[pt],maxint);
	f1->SetParameter(9, ve+15*km/s);
	double max = kg*day*f1->Integral(x[pt],maxint);
	f1->SetParameter(9, ve-15*km/s);
	double min = kg*day*f1->Integral(x[pt],maxint);
	y[pt] = fabs(max-min)/(2. * (fractional ? avg : 1)) ;
      }  
      else{
	double avg = kg*day*keV*DifferentialWimpRate(x[pt]*keV,sigma,rho,A[i],
						     0,1,Md,vesc,v0,ve);
	double max = kg*day*keV*DifferentialWimpRate(x[pt]*keV,sigma,rho,A[i],
						     0,1,Md,vesc,v0,ve+15*km/s);
	double min = kg*day*keV*DifferentialWimpRate(x[pt]*keV,sigma,rho,A[i],
						     0,1,Md,vesc,v0,ve-15*km/s);
	y[pt] = fabs(max-min)/(2.*(fractional ? avg : 1));
      }
      
    }
    TGraph* g = new TGraph(npts,x,y);
    g->SetLineColor(i+1);
    g->SetMarkerColor(i+1);
    g->SetLineWidth(2);
    g->SetTitle(name[i]);
    leg->AddEntry(g,name[i], "l");
    mg->Add(g);
  }
  //if(gPad) gPad->Clear();
  mg->Draw("al");
  mg->GetXaxis()->SetRangeUser(0,200);
  if(!fractional)
    gPad->SetLogy();
  mg->GetYaxis()->UnZoom();
  leg->Draw();
  gPad->Update();
  if(integrated){
    mg->GetXaxis()->SetTitle("Threshold Energy [keV]");
    mg->GetYaxis()->SetTitle("Integrated Modulation [counts kg^{-1} day^{-1}]");
    if(fractional)
      mg->GetYaxis()->SetTitle("Fractional Integrated Modulation");
    mg->GetHistogram()->SetTitle("Modulation of WIMP-induced Nuclear Recoil Rate Above Threshold");
  }
  else{
    mg->GetXaxis()->SetTitle("Recoil Energy [keV]");
    mg->GetYaxis()->SetTitle("Modulation [counts kg^{-1} day^{-1} keV^{-1}]");
    if(fractional)
       mg->GetYaxis()->SetTitle("Fractional Modulation");
    mg->GetHistogram()->SetTitle("Modulation of WIMP-induced Nuclear Recoil Rate");

  }
  mg->GetYaxis()->CenterTitle(true);
  mg->GetYaxis()->SetTitleOffset(1.27);
  
  
  

}

{
  gROOT->Reset();
  gROOT->ProcessLine(".L WimpSpectrum.C+");
  gROOT->ProcessLine(".x StandardVars.C");

  //  double A[] = {ASi, AGe, AAr, AXe};
  std::vector<double>  A{ASi, AGe, AAr, AXe};
  const char* name[] = {"Silicon","Germanium", "Argon","Xenon"};
  const int ngraphs = sizeof(A)/sizeof(double);
  TMultiGraph* mg = new TMultiGraph;
  TLegend* leg = new TLegend(0.49,0.64,0.89,0.89);
  for(int i=0; i<ngraphs; i++){
    const int npts = 200;
    const double emin = 1, emax = 701;  
    //    const double emin = 50, emax = 800;  
   
    double x[npts], y[npts];
    for(int pt = 0; pt<npts; pt++){
      x[pt] = emin + pt*(emax-emin)/(1.*npts);
      y[pt] = kg*day*keV*DifferentialWimpRate(x[pt]*keV,sigma,rho,A[i],
					      0,1,Md,vesc,v0,ve);
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
  TCanvas* c1 = new TCanvas();
  c1->SetLogy();
  mg->Draw("al");
  mg->GetXaxis()->SetRangeUser(0,700);

  mg->GetYaxis()->UnZoom();

  leg->SetHeader(Form("m_{#chi}=%.0f GeV, #sigma_{n}=10^{%.0f} cm^{2}", 
		      Md/GeV, log10(sigma/cm2)));
  leg->SetNColumns(2);
  leg->SetBorderSize(1);
  leg->SetFillStyle(1001);
  leg->SetTextAlign(22);
  
  leg->Draw();

  gPad->Update();
  mg->GetXaxis()->SetTitle("Nuclear Recoil Energy [keV]");
  mg->GetYaxis()->SetTitle("Rate [counts day^{-1} kg^{-1} keV^{-1}]");
  mg->GetYaxis()->CenterTitle(true);
  mg->GetYaxis()->SetTitleOffset(1.27);
  mg->GetHistogram()->SetTitle("WIMP-induced Nuclear Recoil Spectrum");
  
  

}

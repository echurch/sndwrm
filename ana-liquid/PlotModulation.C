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
  
  f1->SetParameters(keV,sigma,rho,A[0],0,1,Md,vesc,v0,ve);
  for(int i=0; i<ngraphs; i++){
    const int npts = 200;
    const double emin = 70, emax = 201;
    
    double x[npts], y[npts];
    f1->SetParameter(3,A[i]);
    for(int pt = 0; pt<npts; pt++){
      x[pt] = 4.*pi*pt/(1.*npts);
      f1->SetParameter(9,ve+15*km/s*cos(x[pt]));
      y[pt] = kg*day*f1->Integral(emin,maxint);
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
  gPad->SetLogy();
  mg->GetYaxis()->UnZoom();
  leg->Draw();
  gPad->Update();
  mg->GetXaxis()->SetTitle("Threshold Energy [keV]");
  mg->GetYaxis()->SetTitle("Integrated Rate [counts kg^{-1} day^{-1}]");
  mg->GetYaxis()->CenterTitle(true);
  mg->GetYaxis()->SetTitleOffset(1.17);
  mg->GetHistogram()->SetTitle("Total WIMP-induced Nuclear Recoil Rate Above Threshold");
  
  

}

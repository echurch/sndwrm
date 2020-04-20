{
  gROOT->Reset();
  gROOT->ProcessLine(".L WimpSpectrum.C+");
  TF1* f1 = new TF1("f1","[0]*VelocityFunction([0]*x,[1],[2],[3])",0,900);
  f1->SetParameters(km/s, 544*km/s,250*km/s,263*km/s);
  f1->Draw();
  f1->SetTitle("WIMP Lab Speed Distribution");
  f1->GetXaxis()->SetTitle("WIMP Speed [km/s]");
  f1->GetYaxis()->SetTitle("dP/dV [(km/s)^{-1}]");
  gPad->Update();

}

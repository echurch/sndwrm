
import pdb
import numpy as np
from ROOT import TFile, TH1F

f = TFile("e_optphys.root")
h0 = TH1F("prim","Vertical distance of primary from mid-plane [m]",50,0.,6.)
h1 = TH1F("cphot","Count of opt photons/event penetrating plane at 5.5m vertical distance vs primary y location [m]",50,0.,6.)


Nent = f.Step.GetEntries()
Yprim = None

for step in range(Nent):

    entry = f.Step.GetEntry(step)

    if f.Step.PID==11 and f.Step.StepNum==1:
        Yprim = f.Step.Y/1000.0
        h0.Fill(Yprim)

        
    if f.Step.PID==0 and f.Step.TY/1000.0>5.5:
        h1.Fill(Yprim)

h1.Divide(h1,h0)

    
f2 = TFile("e_opthys_hists.root","RECREATE")
h0.Write()
h1.Write()
f2.Close()

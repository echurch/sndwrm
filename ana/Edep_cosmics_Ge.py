from ROOT import TFile, gPad, TH1F, TCanvas, gROOT, TGraph
import numpy as np
import pdb
#gROOT.SetBatch(1)



###### cuidado: this assumes Npq=Nps and I am properly transcribing 13*7*7
def exposeTime(fTGfromSS,fSim):
    tg = fTGfromSS.proton_spectra
    Npts = tg.GetN()
    x = np.zeros(Npts); y = np.zeros(Npts)
    for ii in range(Npts):
        y[ii] = tg.GetPointY(ii)
        x[ii] = tg.GetPointX(ii)
    dx = x[1:]-x[:-1]
    integral = (dx*y[1:]).sum() # 1/days*cm2
    integral = 13*7*7*2*2*2*1E8  ## These are the halfbox dimenstions in m, units now 1/days
    time = Npq/integral
    return time


if __name__=="__main__":

    fq = TFile("../build/proton_sqbbc_cosmic.root")
    fs = TFile("../build/proton_shielding_cosmic.root")

    Npq = fq.T1.GetEntries()
    Nps = fs.T1.GetEntries()
    #Ellen says units on y-axis are EU/day/cm2, with EU being MeV for protons
    fsource = TFile("tgraph-spectra_cosmic/job6b_spectra.root")
    Texp = exposeTime(fsource,fq)

    binsz = 10
    hs = TH1F ("hsn","",50,0.,binsz*100.)
    hq = TH1F ("hqn","",50,0.,binsz*100.)

    c1 = TCanvas()
    fs.T3.Draw("EnergyDepEvt>>hsn") ##,"","",1000000)
    fq.T3.Draw("EnergyDepEvt>>hqn") ##,"","",1000000)
    del c1
    print(f"fs and fq Nentries: {hs.Integral()} and {hq.Integral()}")

    c1 = TCanvas()
    hq.SetLineStyle(2)
    hs.SetLineStyle(1)

    hs.GetYaxis().SetTitle(f"OLTARIS protons - entries per day per {binsz} MeV")
    hs.GetXaxis().SetTitle("Energy deposited in pathfinder_psyche HPGe [MeV]")
    hs.Scale(1.0/Texp)
    hq.Scale(1.0/Texp)
    hs.SetMinimum(10);
    hs.Draw("Hist")
    hq.Draw("Hist,same")
    gPad.SetLogy(1)
    c1.SetLogy(1)
    #gPad.Update()
    hs.SetTitle("OLTARIS protons Shielding PhysicsList ")
    hq.SetTitle("SPQBBC PhysicsList ")
    gPad.BuildLegend(0.7,0.5,0.95,0.7)
    #gPad.Update()
    c1.Update()

    c1.SaveAs("ps_HPGe_OLTARIS6bhi_homog_SHIELDING_SPQBBC.png")
    del c1


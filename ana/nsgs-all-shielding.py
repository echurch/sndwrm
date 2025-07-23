
from ROOT import TFile, gStyle, gPad, gROOT, TH1D, TH2D
import pdb
from matplotlib import pyplot as plt
import numpy as np


# Don't pop up TCanvas
gROOT.SetBatch(1)

# First calculate total gamma flux rate from bottom for our no-shielding job

fnf0  = TFile("cavneutrs_FLR_7x7x31_auto_25mmLead_0cmBP.root")
fnf10 = TFile("cavneutrs_FLR_7x7x31_auto_25mmLead_10cmBP.root")
fnf20 = TFile("cavneutrs_FLR_7x7x31_auto_25mmLead_20cmBP.root")
fnf30 = TFile("cavneutrs_FLR_7x7x31_auto_25mmLead_30cmBP.root")
fnf40 = TFile("cavneutrs_FLR_7x7x31_auto_25mmLead_40cmBP.root")
fnf60 = TFile("cavneutrs_FLR_7x7x31_auto_25mmLead_60cmBP.root")

fnc0  = TFile("cavneutrs_CAV_7x7x31_auto_25mmLead_0cmBP.root")
fnc10 = TFile("cavneutrs_CAV_7x7x31_auto_25mmLead_10cmBP.root")
fnc20 = TFile("cavneutrs_CAV_7x7x31_auto_25mmLead_20cmBP.root")
fnc30 = TFile("cavneutrs_CAV_7x7x31_auto_25mmLead_30cmBP.root")
fnc40 = TFile("cavneutrs_CAV_7x7x31_auto_25mmLead_40cmBP.root")
fnc60 = TFile("cavneutrs_CAV_7x7x31_auto_25mmLead_60cmBP.root")

fgf0  = TFile("cavgam_FLR_7x7x31_auto_25mmLead_0cmBP.root")
fgf10 = TFile("cavgam_FLR_7x7x31_auto_25mmLead_10cmBP.root")
fgf20 = TFile("cavgam_FLR_7x7x31_auto_25mmLead_20cmBP.root")
fgf30 = TFile("cavgam_FLR_7x7x31_auto_25mmLead_30cmBP.root")
fgf40 = TFile("cavgam_FLR_7x7x31_auto_25mmLead_40cmBP.root")
fgf60 = TFile("cavgam_FLR_7x7x31_auto_25mmLead_60cmBP.root")

fgc0  = TFile("cavgam_CAV_7x7x31_auto_25mmLead_0cmBP.root")
fgc10 = TFile("cavgam_CAV_7x7x31_auto_25mmLead_10cmBP.root")
fgc20 = TFile("cavgam_CAV_7x7x31_auto_25mmLead_20cmBP.root")
fgc30 = TFile("cavgam_CAV_7x7x31_auto_25mmLead_30cmBP.root")
fgc40 = TFile("cavgam_CAV_7x7x31_auto_25mmLead_40cmBP.root")
fgc60 = TFile("cavgam_CAV_7x7x31_auto_25mmLead_60cmBP.root")


def gs_fl():
        ## Gamma floor flux

        Nevts = fgf0.T3.GetMaximum("Event")
        h10 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
        fgf0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/20/70*12*56)","colz") # 12.6E4 gs/m^2 * surface size launched, then divided by surf size counted
        Ntot0 = h10.Integral()/Nevts
        
        NtotPoly = []
        filelist = [fgf0, fgf10, fgf20, fgf30, fgf40, fgf60]
        for file in filelist:

            file.T3.Draw("Event")
            Nevts = file.T3.GetMaximum("Event")
            
            h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
            file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/20/70*12*56)","colz") 
            NtotPoly.append(h1.Integral()/Nevts)
            del h1
    


        fig, ax = plt.subplots()
        ax.plot([0,10,20,30,40,60],NtotPoly/np.repeat(Ntot0,6),label="B-Poly ",linestyle='-',marker='.',color="black")


        ay2 = ax.twinx()
        ay2.plot([0,10,20,30,40,60],NtotPoly,linestyle='-',marker='*',color='green',alpha=0.5)
        ay2.set_yticks(np.round(np.arange(0.,40.,4.),decimals=3))
        ay2.set_yticklabels(np.round(np.arange(0,40.,4.),decimals=3),color='green')
        ay2.set_ylabel('gamma floor flux with BP+Pb shielding [Hz/m2]',color='green')
        #ay2.set_ylim(0,40)
        ax.legend()
        ax.set_title("gamma floor [12.6E4 gs/cm2/sec] attenuation for Edep>1MeV")
        ax.set_xlabel("Floor B-Poly thickness, along w 2.5cm floor Pb, and 23cm BPE in lower 2/3 of walls")
        ax.set_ylabel("Flux fraction compared to no-shielding")
        ax.set_yscale('log')
        ay2.set_yscale('log')
        #plt.show()
        plt.tight_layout()

        plt.savefig("./gamma-attenutation_floor-PolyandPbshielding.png")


def ns_fl():
# First calculate total neutron flux rate from bottom for our no-shielding job

    Nevts = fnf0.T3.GetMaximum("Event")
    h10 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    fnf0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&EnergyDepEvt>1.0&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/2.94E-2/20/70)","colz") # 2.94E-2 ns/m^2 * surface size launched
    Ntot0 = h10.Integral()/Nevts


    # Now calculate total neutron flux rate from bottom for our non-zero shielding jobs
    NtotPoly = []
    filelist = [fnf0, fnf10, fnf20, fnf30, fnf40, fnf60]
    for file in filelist:

        file.T3.Draw("Event")
        Nevts = file.T3.GetMaximum("Event")
    
        h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
        file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.0&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/2.94E-2/20/70)","colz") 
        NtotPoly.append(h1.Integral()/Nevts)
        del h1
    


    fig, ax = plt.subplots()
    ax.plot([0,10,20,30,40,60],NtotPoly/np.repeat(Ntot0,6),label="B-Poly ",linestyle='-',marker='.',color="black")

    ay2 = ax.twinx()
    ay2.plot([0,10,20,30,40,60],NtotPoly,linestyle='-',marker='*',color='green',alpha=0.5)
    ay2.set_yticks(np.round(np.arange(0.,1.0,0.1),decimals=3))
    ay2.set_yticklabels(np.round(np.arange(0,1.0,0.1),decimals=3),color='green')
    ay2.set_ylabel('neutron floor flux with BP+Pb shielding [Hz]',color='green')
#    ay2.set_ylim(0,1.2)
    ax.legend()
    ax.set_title("neutron floor [2.94E-2 ns/cm2/sec] attenuation for Edep>1MeV")
    ax.set_xlabel("Floor B-Poly thickness [cm], along w 2.5cm floor Pb, and 23cm BP in lower 2/3 of walls")
    ax.set_ylabel("Rate fraction compared to no-shielding")
    ax.set_yscale('log')
    ay2.set_yscale('log')
    plt.tight_layout()
    #plt.show()
    plt.savefig("./neutron-attenutation_floor-PolyandPbshielding.png")


def gs_fl_cav():
    ## Gamma+CAV gamma rate

    Nevts = fgf0.T3.GetMaximum("Event")
    h10 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    fgf0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/20/70)","colz") # 12.6E4 gs/m^2 * surface size launched, then NOT divided by surf size counted
    Ntot0 = h10.Integral()/Nevts

    Nevts2 = fgc0.T3.GetMaximum("Event")
    h11 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    fgc0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/(3.14159*16*70+3.14159*16*16))","colz") # 12.6E4 gs/m^2 * surface size launched, then NOT divided by surf size counted
    Ntot2 = h11.Integral()/Nevts2

    NtotPoly = []
    filelist = [fgf0, fgf10, fgf20, fgf30, fgf40, fgf60]
    for file in filelist:

        file.T3.Draw("Event")
        Nevts = file.T3.GetMaximum("Event")
    
        h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
        file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/20/70)","colz") 
        NtotPoly.append(h1.Integral()/Nevts)
        del h1

    NtotPoly2 = []
    filelist2 = [fgc0, fgc10, fgc20, fgc30, fgc40, fgc60]
    for file in filelist2:
        
        file.T3.Draw("Event")
        Nevts2 = file.T3.GetMaximum("Event")
    
        h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
        file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/(3.14159*16*70+3.14159*16*16))","colz") 
        NtotPoly2.append(h1.Integral()/Nevts2)
        del h1


    fig, ax = plt.subplots()
    ax.plot([0,10,20,30,40,60],[x+y for x,y in zip(NtotPoly,NtotPoly2)]/np.repeat(Ntot0+Ntot2,6),label="B-Poly ",linestyle='-',marker='.',color="black")

    ay2 = ax.twinx()
    ay2.plot([0,10,20,30,40,60],[x+y for x,y in zip(NtotPoly,NtotPoly2)],linestyle='-',marker='*',color='green',alpha=0.5)
    ay2.set_yticks(np.round(np.arange(0.,40000.,5000.),decimals=1))
    ay2.set_yticklabels(np.round(np.arange(0,40000.,5000.),decimals=1),color='green')
    ay2.set_ylabel('gamma rate with BP+Pb shielding [Hz]',color='green')
    #ay2.set_ylim(0,40)
    ax.legend()
    ax.set_title("gamma cavern+floor [12.6E4 gs/cm2/sec]  for Edep>1MeV")
    ax.set_xlabel("Floor B-Poly thickness [cm], along w 2.5cm floor Pb, and 23cm BP in lower 2/3 of walls")
    ax.set_ylabel("Rate attenuation compared to no-shielding")
    ax.set_yticks(np.round(np.arange(0.,1.01,0.1),decimals=1))
    ax.set_yticklabels(np.round(np.arange(0,1.01,0.1),decimals=1),color='black')
    #plt.show()
    plt.tight_layout()
    #ax.set_yscale('log')
    plt.savefig("./gamma-attenutation_floorandcyl-PolyandPbshielding.png")


def ns_fl_cav():
    ## neutrons+CAV gamma rate

    Nevts = fnf0.T3.GetMaximum("Event")
    h10 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    fnf0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<30E3&&abs(X)<6.5E3)*1./(1.0/2.94E-2/20/70)","colz") # 2.94 ns/m^2 * surface size launched, then NOT divided by surf size counted
    Ntot0 = h10.Integral()/Nevts

    Nevts2 = fnc0.T3.GetMaximum("Event")
    h11 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    fnc0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<30E3&&abs(X)<6.5E3)*1./(1.0/2.94E-2/(3.14159*16*70+3.14159*16*16))","colz") # 2.94 ns/m^2 * surface size launched, then NOT divided by surf size counted
    Ntot2 = h11.Integral()/Nevts2

    NtotPoly = []
    filelist = [fnf0, fnf10, fnf20, fnf30, fnf40, fnf60]
    for file in filelist:

        file.T3.Draw("Event")
        Nevts = file.T3.GetMaximum("Event")
    
        h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
        file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<30E3&&abs(X)<6.5E3)*1./(1.0/2.94E-2/20/70)","colz") 
        NtotPoly.append(h1.Integral()/Nevts)
        del h1

    NtotPoly2 = []
    filelist2 = [fnc0, fnc10, fnc20, fnc30, fnc40, fnc60]
    for file in filelist2:
        
        file.T3.Draw("Event")
        Nevts2 = file.T3.GetMaximum("Event")
    
        h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
        file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<30E3&&abs(X)<6.5E3)*1./(1.0/2.94E-2/(3.14159*16*70+3.14159*16*16))","colz") 
        NtotPoly2.append(h1.Integral()/Nevts2)
        del h1


    fig, ax = plt.subplots()
    ax.plot([0,10,20,30,40,60],[x+y for x,y in zip(NtotPoly,NtotPoly2)]/np.repeat(Ntot0+Ntot2,6),label="B-Poly ",linestyle='-',marker='.',color="black")

    ay2 = ax.twinx()

    ay2.plot([0,10,20,30,40,60],[x+y for x,y in zip(NtotPoly,NtotPoly2)],linestyle='-',marker='*',color='green',alpha=0.5)
    ay2.set_yticks(np.round(np.arange(0.,1.4,0.1),decimals=1))
    ay2.set_yticklabels(np.round(np.arange(0,1.4,0.1),decimals=1),color='green')
    ay2.set_ylabel('neutron rate with BP+Pb shielding [Hz]',color='green')
#    ay2.set_ylim(0,1.4)
    ax.legend()
    ax.set_title("neutron cavern+floor [2.94E-2 ns/cm2/sec]  for Edep>1MeV")
    ax.set_xlabel("Floor B-Poly thickness [cm], along w 2.5cm floor Pb, and 23 cm PB in lower 2/3 of walls")
    ax.set_ylabel("Rate attenuation compared to no-shielding")
    ax.set_yticks(np.round(np.arange(0.,1.01,0.1),decimals=1))
    ax.set_yticklabels(np.round(np.arange(0,1.01,0.1),decimals=1),color='black')
    #plt.show()
    plt.tight_layout()
    #ax.set_yscale('log')
    plt.savefig("./neutron-attenutation_floorandcyl-PolyandPbshielding.png")


def gs_fl_cav_low():
    ## Gamma+CAV gamma rate

    Nevts = fgf0.T3.GetMaximum("Event")
    h10 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    fgf0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&Y<0&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/20/70)","colz") # 12.6E4 gs/m^2 * surface size launched, then NOT divided by surf size counted
    Ntot0 = h10.Integral()/Nevts

    Nevts2 = fgc0.T3.GetMaximum("Event")
    h11 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    fgc0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&Y<0&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/(3.14159*16*70+3.14159*16*16))","colz") # 12.6E4 gs/m^2 * surface size launched, then NOT divided by surf size counted
    Ntot2 = h11.Integral()/Nevts2

    NtotPoly = []
    filelist = [fgf0, fgf10, fgf20, fgf30, fgf40, fgf60]
    for file in filelist:

        file.T3.Draw("Event")
        Nevts = file.T3.GetMaximum("Event")
    
        h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
        file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&Y<0&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/20/70)","colz") 
        NtotPoly.append(h1.Integral()/Nevts)
        del h1

    NtotPoly2 = []
    filelist2 = [fgc0, fgc10, fgc20, fgc30, fgc40, fgc60]
    for file in filelist2:
        
        file.T3.Draw("Event")
        Nevts2 = file.T3.GetMaximum("Event")
    
        h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
        file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&Y<0&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/(3.14159*16*70+3.14159*16*16))","colz") 
        NtotPoly2.append(h1.Integral()/Nevts2)
        del h1


    fig, ax = plt.subplots()
    ax.plot([0,10,20,30,40,60],[x+y for x,y in zip(NtotPoly,NtotPoly2)]/np.repeat(Ntot0+Ntot2,6),label="B-Poly ",linestyle='-',marker='.',color="black")

    ay2 = ax.twinx()
    ay2.plot([0,10,20,30,40,60],[x+y for x,y in zip(NtotPoly,NtotPoly2)],linestyle='-',marker='*',color='green',alpha=0.5)
    ay2.set_yticks(np.round(np.arange(0.,40000.,5000.),decimals=1))
    ay2.set_yticklabels(np.round(np.arange(0,40000.,5000.),decimals=1),color='green')
    ay2.set_ylabel('gamma rate with BP+Pb shielding [Hz]',color='green')
    #ay2.set_ylim(0,40)
    ax.legend()
    ax.set_title("gamma cavern+floor [12.6E4 gs/cm2/sec] for Edep>1MeV -- Lower Half")
    ax.set_xlabel("Floor B-Poly thickness [cm], along w 2.5cm floor Pb, and 23cm BP in lower 2/3 of walls")
    ax.set_ylabel("Rate attenuation compared to no-shielding")
    ax.set_yticks(np.round(np.arange(0.,1.01,0.1),decimals=1))
    ax.set_yticklabels(np.round(np.arange(0,1.01,0.1),decimals=1),color='black')

    ax.set_yscale('log')
    ay2.set_yscale('log')
    #plt.show()
    plt.tight_layout()

    plt.savefig("./gamma-attenutation_floorandcyl-PolyandPbshielding_low.png")


def ns_fl_cav_low():
    ## neutrons+CAV gamma rate

    Nevts = fnf0.T3.GetMaximum("Event")
    h10 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    fnf0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&Y<0&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/2.94E-2/20/70)","colz") # 2.94 ns/m^2 * surface size launched, then NOT divided by surf size counted
    Ntot0 = h10.Integral()/Nevts

    Nevts2 = fnc0.T3.GetMaximum("Event")
    h11 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    fnc0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&Y<0&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/2.94E-2/(3.14159*16*70+3.14159*16*16))","colz") # 2.94 ns/m^2 * surface size launched, then NOT divided by surf size counted
    Ntot2 = h11.Integral()/Nevts2

    NtotPoly = []
    filelist = [fnf0, fnf10, fnf20, fnf30, fnf40, fnf60]
    for file in filelist:

        file.T3.Draw("Event")
        Nevts = file.T3.GetMaximum("Event")
    
        h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
        file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&Y<0&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/2.94E-2/20/70)","colz") 
        NtotPoly.append(h1.Integral()/Nevts)
        del h1

    NtotPoly2 = []
    filelist2 = [fnc0, fnc10, fnc20, fnc30, fnc40, fnc60]
    for file in filelist2:
        
        file.T3.Draw("Event")
        Nevts2 = file.T3.GetMaximum("Event")
    
        h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
        file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&Y<0&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/2.94E-2/(3.14159*16*70+3.14159*16*16))","colz") 
        NtotPoly2.append(h1.Integral()/Nevts2)
        del h1


    fig, ax = plt.subplots()
    ax.plot([0,10,20,30,40,60],[x+y for x,y in zip(NtotPoly,NtotPoly2)]/np.repeat(Ntot0+Ntot2,6),label="B-Poly ",linestyle='-',marker='.',color="black")

    ay2 = ax.twinx()

    ay2.plot([0,10,20,30,40,60],[x+y for x,y in zip(NtotPoly,NtotPoly2)],linestyle='-',marker='*',color='green',alpha=0.5)
    ay2.set_yticks(np.round(np.arange(0.,1.4,0.1),decimals=1))
    ay2.set_yticklabels(np.round(np.arange(0,1.4,0.1),decimals=1),color='green')
    ay2.set_ylabel('neutron rate with BP+Pb shielding [Hz]',color='green')
#    ay2.set_ylim(0,1.4)
    ax.legend()
    ax.set_title("neutron cavern+floor [2.94E-2 ns/cm2/sec]  for Edep>1MeV -- Lower Half")
    ax.set_xlabel("Floor B-Poly thickness [cm], along w 2.5cm floor Pb, and 23 cm PB in lower 2/3 of walls")
    ax.set_ylabel("Rate attenuation compared to no-shielding")
    ax.set_yticks(np.round(np.arange(0.,1.01,0.1),decimals=1))
    ax.set_yticklabels(np.round(np.arange(0,1.01,0.1),decimals=1),color='black')

    ax.set_yscale('log')
    ay2.set_yscale('log')
    #plt.show()
    plt.tight_layout()
    #ax.set_yscale('log')
    plt.savefig("./neutron-attenutation_floorandcyl-PolyandPbshielding_low.png")



        
if __name__ == "__main__":
    ns_fl_cav()
    #gs_fl_cav()        
    #ns_fl()
    #gs_fl() # The only one that shows rate in per m2 per sec.
        
    # ns_fl_cav_low()
    #gs_fl_cav_low()        

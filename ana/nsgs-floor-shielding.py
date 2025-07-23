
from ROOT import TFile, gStyle, gPad, gROOT, TH1D, TH2D
import pdb
from matplotlib import pyplot as plt
import numpy as np


# Don't pop up TCanvas
gROOT.SetBatch(1)

# First calculate total gamma flux rate from bottom for our no-shielding job

f0 = TFile("cavneutr_bot_7x7x31_auto_None_0.root")
fP20 = TFile("cavneutr_bot_7x7x31_auto_Poly_20.root")
fP40 = TFile("cavneutr_bot_7x7x31_auto_Poly_40.root")
fP60 = TFile("cavneutr_bot_7x7x31_auto_Poly_60.root")
fP80 = TFile("cavneutr_bot_7x7x31_auto_Poly_80.root")
fP100 = TFile("cavneutr_bot_7x7x31_auto_Poly_100.root")
fH20 = TFile("cavneutr_bot_7x7x31_auto_H2O_20.root")
fH40 = TFile("cavneutr_bot_7x7x31_auto_H2O_40.root")
fH60 = TFile("cavneutr_bot_7x7x31_auto_H2O_60.root")
fH80 = TFile("cavneutr_bot_7x7x31_auto_H2O_80.root")
fH100 = TFile("cavneutr_bot_7x7x31_auto_H2O_100.root")
fPb20 = TFile("cavneutr_bot_7x7x31_auto_Pb_2.root")
fPb40 = TFile("cavneutr_bot_7x7x31_auto_Pb_4.root")
fPb60 = TFile("cavneutr_bot_7x7x31_auto_Pb_6.root")
fPb80 = TFile("cavneutr_bot_7x7x31_auto_Pb_8.root")
fPb100 = TFile("cavneutr_bot_7x7x31_auto_Pb_10.root")


# First calculate total neutron flux rate from bottom for our no-shielding job

Nevts = f0.T3.GetMaximum("Event")
h10 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
f0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&EnergyDepEvt>1.0&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/2.94E-2/20/70)","colz") # 2.94E-2 ns/m^2 * surface size launched, then divided by surf size counted
Ntot0 = h10.Integral()/Nevts


# Now calculate total neutron flux rate from bottom for our non-zero shielding jobs
NtotPoly = []
filelist = [fP20, fP40, fP60, fP80, fP100]
for file in filelist:

    file.T3.Draw("Event")
    Nevts = f0.T3.GetMaximum("Event")
    
    h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.0&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/2.94E-2/20/70)","colz") 
    NtotPoly.append(h1.Integral()/Nevts)
    del h1
    
NtotH2O = []
filelist = [fH20, fH40, fH60, fH80, fH100]
for file in filelist:

    file.T3.Draw("Event")
    Nevts = f0.T3.GetMaximum("Event")
    
    h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.0&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/2.94E-2/20/70)","colz")
    NtotH2O.append(h1.Integral()/Nevts)
    del h1

NtotPb = []
filelist = [fPb20, fPb40, fPb60, fPb80, fPb100]
for file in filelist:

    file.T3.Draw("Event")
    Nevts = f0.T3.GetMaximum("Event")
    
    h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.0&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/2.94E-2/20/70)","colz")
    NtotPb.append(h1.Integral()/Nevts)
    del h1


fig, ax = plt.subplots()
ax.plot([20,40,60,80,100],NtotPoly/np.repeat(Ntot0,5),label="B-Poly ",linestyle='-',marker='.',color="black")
ax.plot([20,40,60,80,100],NtotH2O/np.repeat(Ntot0,5),label="H2O ",linestyle='-',marker='.',color='blue')
ax.plot([20,40,60,80,100],NtotPb/np.repeat(Ntot0,5),label="Lead",linestyle='-',marker='.',color='red')    

ax2 = ax.secondary_xaxis('bottom')
ax2.set_xlabel('Shield Lead thickness [cm]')
# Set the color of the secondary x-axis
ax2.spines['bottom'].set_color('red')
ax2.spines["bottom"].set_position(("outward", 35.))
ax2.set_xticks([20, 40, 60, 80, 100])
ax2.set_xticklabels([2, 4, 6, 8, 10])
ax2.tick_params(axis='x', colors='red')
ax2.xaxis.label.set_color('red')

ay2 = ax.twinx()
ay2.plot([20,40,60,80,100],NtotH2O,linestyle='-',marker='*',color='green',alpha=0.5)
ay2.set_yticks(np.round(np.arange(0.,0.06,0.01),decimals=3))
ay2.set_yticklabels(np.round(np.arange(0,0.06,0.01),decimals=3),color='green')
ay2.set_ylabel('n floor rate with H2O shielding [Hz]',color='green')
##ay2.set_yscale('log')

ax.legend()
ax.set_title("neutron floor [2.94E-6 ns/cm2/sec] attenuation for Edep>1MeV")
ax.set_xlabel("Shield B-Poly, H2O thickness [cm]")
ax.set_ylabel("Rate fraction compared to no-shielding")
ax.set_yscale('log')

plt.tight_layout()
#plt.show()
plt.savefig("./neutron-attenutation_floor-shielding.png")



f0 = TFile("cavgam_bot_7x7x31_auto_None.root")
fP20 = TFile("cavgam_bot_7x7x31_auto_Poly_20.root")
fP40 = TFile("cavgam_bot_7x7x31_auto_Poly_40.root")
fP60 = TFile("cavgam_bot_7x7x31_auto_Poly_60.root")
fP80 = TFile("cavgam_bot_7x7x31_auto_Poly_80.root")
fP100 = TFile("cavgam_bot_7x7x31_auto_Poly_100.root")
fH20 = TFile("cavgam_bot_7x7x31_auto_20_H2O.root")
fH40 = TFile("cavgam_bot_7x7x31_auto_40_H2O.root")
fH60 = TFile("cavgam_bot_7x7x31_auto_60_H2O.root")
fH80 = TFile("cavgam_bot_7x7x31_auto_80_H2O.root")
fH100 = TFile("cavgam_bot_7x7x31_auto_100_H2O.root")
fPb20 = TFile("cavgam_bot_7x7x31_auto_Lead_2.root")
fPb40 = TFile("cavgam_bot_7x7x31_auto_Lead_4.root")
fPb60 = TFile("cavgam_bot_7x7x31_auto_Lead_6.root")
fPb80 = TFile("cavgam_bot_7x7x31_auto_Lead_8.root")
fPb100 = TFile("cavgam_bot_7x7x31_auto_Lead_10.root")


Nevts = f0.T3.GetMaximum("Event")
h10 = TH2D("hname0","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
f0.T3.Draw("Z:X>>hname0","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/20/70*12*56)","colz") # 12.6E4 gs/m^2 * surface size launched, then divided by surf size counted
Ntot0 = h10.Integral()/Nevts

NtotPoly = []
filelist = [fP20, fP40, fP60, fP80, fP100]
for file in filelist:

    file.T3.Draw("Event")
    Nevts = f0.T3.GetMaximum("Event")
    
    h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/20/70*12*56)","colz") 
    NtotPoly.append(h1.Integral()/Nevts)
    del h1
    
NtotH2O = []
filelist = [fH20, fH40, fH60, fH80, fH100]
for file in filelist:

    file.T3.Draw("Event")
    Nevts = f0.T3.GetMaximum("Event")
    
    h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/20/70*12*56)","colz") 
    NtotH2O.append(h1.Integral()/Nevts)
    del h1

NtotPb = []
filelist = [fPb20, fPb40, fPb60, fPb80, fPb100]
for file in filelist:

    file.T3.Draw("Event")
    Nevts = f0.T3.GetMaximum("Event")
    
    h1 = TH2D("hname","htitle",40,-6E3,6.E3,40,-28E3,+28E3)
    file.T3.Draw("Z:X>>hname","(Y>-6.972E3&&EnergyDepEvt>1.&&abs(Z)<28E3&&abs(X)<6E3)*1./(1.0/12.6E4/20/70*12*56)","colz") 
    NtotPb.append(h1.Integral()/Nevts)
    del h1


fig, ax = plt.subplots()
ax.plot([20,40,60,80,100],NtotPoly/np.repeat(Ntot0,5),label="B-Poly ",linestyle='-',marker='.',color="black")
ax.plot([20,40,60,80,100],NtotH2O/np.repeat(Ntot0,5),label="H2O ",linestyle='-',marker='*',color='blue')
ax.plot([20,40,60,80,100],NtotPb/np.repeat(Ntot0,5),label="Lead",linestyle='-',marker='o',color='red')    

ax2 = ax.secondary_xaxis('bottom')
ax2.set_xlabel('Shield Lead thickness [cm]')
# Set the color of the secondary x-axis
ax2.spines['bottom'].set_color('red')
ax2.spines["bottom"].set_position(("outward", 35.))
ax2.set_xticks([20, 40, 60, 80, 100])
ax2.set_xticklabels([2, 4, 6, 8, 10])
ax2.tick_params(axis='x', colors='red')
ax2.xaxis.label.set_color('red')

ay2 = ax.twinx()
ay2.plot([20,40,60,80,100],NtotH2O,linestyle='-',marker='*',color='green',alpha=0.5)
ay2.set_yticks(np.round(np.arange(0.,20.,3.),decimals=3))
ay2.set_yticklabels(np.round(np.arange(0,20.,3.),decimals=3),color='green')
ay2.set_ylabel('g floor flux with H2O shielding [Hz/m2]',color='green')

ax.legend()
ax.set_title("gamma floor [12.6 gs/cm2/sec] attenuation for Edep>1MeV")
ax.set_xlabel("Shield B-Poly, H2O thickness [cm]")
ax.set_ylabel("Flux fraction compared to no-shielding")
#plt.show()
plt.tight_layout()
#ax.set_yscale('log')
plt.savefig("./gamma-attenutation_floor-shielding.png")




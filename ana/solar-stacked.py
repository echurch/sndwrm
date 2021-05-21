import numpy as np
from matplotlib import pyplot as plt
from ROOT import TFile, TH1, THStack
from skhep.visual import MplPlotter as skh_plt

import pdb

#neutrons
#f0n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_0.root")
#f0n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_recal_0.root")
f0n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_shinyg10_0.root")
#f0a39 = TFile("/Volumes/EC2TB/G4/mod4/Ar39_fidv_recal_0.root")
f0a39 = TFile("/Volumes/EC2TB/G4/mod4/Ar39_fidv_shinyg10_0.root")
#f0a42 = TFile("/Volumes/EC2TB/G4/mod4/Ar42_fidv_recal_0.root")
f0a42 = TFile("/Volumes/EC2TB/G4/mod4/Ar42_fidv_shinyg10_0.root")
f0r = TFile("/Volumes/EC2TB/G4/mod4/Rn222_fidv_shinyg10_0.root")
#f0h = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_recal_0.root")
f0h = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_shinyg10_0.root")
f08 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_b8_shinyg10_0.root")
#f0c = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_CNO_recal_0.root")
f0c = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_CNO_shinyg10_0.root")
h0nCC = f0n.H23
h0hES = f0h.H22
h0hCC = f0h.H23
h08ES = f08.H22
h08CC = f08.H23
h0cES = f0c.H22
h0cCC = f0c.H23
h0a39 = f0a39.H23
h0a42 = f0a42.H23
h0r = f0r.H23

#f1n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_1.root")
#f1n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_recal_1.root")
f1n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_shinyg10_1.root")
#f1a39 = TFile("/Volumes/EC2TB/G4/mod4/Ar39_fidv_recal_1.root")
f1a39 = TFile("/Volumes/EC2TB/G4/mod4/Ar39_fidv_shinyg10_1.root")
#f1a42 = TFile("/Volumes/EC2TB/G4/mod4/Ar42_fidv_recal_1.root")
f1a42 = TFile("/Volumes/EC2TB/G4/mod4/Ar42_fidv_shinyg10_1.root")
f1r = TFile("/Volumes/EC2TB/G4/mod4/Rn222_fidv_shinyg10_1.root")
f1h = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_recal_1.root")
f1h = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_shinyg10_1.root")
f18 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_b8_shinyg10_1.root")
#f1c = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_CNO_recal_1.root")
f1c = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_CNO_shinyg10_1.root")
h1nCC = f1n.H23
h1hES = f1h.H22
h1hCC = f1h.H23
h18ES = f18.H22
h18CC = f18.H23
h1cES = f1c.H22
h1cCC = f1c.H23
h1a39 = f1a39.H23
h1a42 = f1a42.H23
h1r = f1r.H23

#f2n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_2.root")
#f2n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_recal_2.root")
f2n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_shinyg10_2.root")
#f2a39 = TFile("/Volumes/EC2TB/G4/mod4/Ar39_fidv_2.root")
#f2a42 = TFile("/Volumes/EC2TB/G4/mod4/Ar42_fidv_shinyg10_2.root")
f2r = TFile("/Volumes/EC2TB/G4/mod4/Rn222_fidv_shinyg10_2.root")
#f2h = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_recal_2.root")
f2h = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_shinyg10_2.root")
f28 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_b8_shinyg10_2.root")
#f2c = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_CNO_recal_2.root")
f2c = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_CNO_shinyg10_2.root")
h2nCC = f2n.H23
h2hES = f2h.H22
h2hCC = f2h.H23
h28ES = f28.H22
h28CC = f28.H23
h2cES = f2c.H22
h2cCC = f2c.H23
#h2a39 = f2a39.H23
#h2a42 = f2a42.H23
h2r = f2r.H23

#f3n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_3.root")
#f3n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_recal_3.root")
f3n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_shinyg10_3.root")
#f3a39 = TFile("/Volumes/EC2TB/G4/mod4/Ar39_fidv_3.root")
#f3a42 = TFile("/Volumes/EC2TB/G4/mod4/Ar42_fidv_3.root")
f3r = TFile("/Volumes/EC2TB/G4/mod4/Rn222_fidv_shinyg10_3.root")
#f3h = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_3.root")
f3h = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_shinyg10_3.root")
f38 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_b8_shinyg10_3.root")
#f3c = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_CNO_recal_3.root")
f3c = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_CNO_shinyg10_3.root")
h3nCC = f3n.H23
h3hES = f3h.H22
h3hCC = f3h.H23
h38ES = f38.H22
h38CC = f38.H23
h3cES = f3c.H22
h3cCC = f3c.H23
#h3a39 = f3a39.H23
#h3a42 = f3a42.H23
h3r = f3r.H23

# Calculate normalizaions, ala xs_flux-avgd * N_Ar * Flux_int

N_Ar = 1.4 * 6*12*20 * 1E6 /40. *6.022E23
kgAr = (6*12*20)*1E6 * 1.4 * 1E-3
# In my smaller (less tall in y) fidV I would have just as many histo entries, now weighted thusly.
N_Ar = 1.4 * 6*9*20 * 1E6 /40. *6.022E23
kgAr = (6*9*20)*1E6 * 1.4 * 1E-3
# Below from Marley runs' log files

xctot = 0.00101851 * 1E-40 #cm^2
x8tot = 0.0406549 * 1E-40
xhtot = 0.111884  * 1E-40

# From https://arxiv.org/pdf/astro-ph/0402114.pdf, Bahcall/Pinsonneault
# Tacking on rough oscillation factors as from Fig 5 of https://www.sciencedirect.com/science/article/pii/S2212686414000211
Fluxh = 7.88E3 * 0.3 # /cm2/sec/atom
Flux8 = 5.79E6 * 0.3
Fluxc = (5.71+5.03+0.0591)*1E8 * 0.5

#Normalize to 1 year. I threw 1500*4 evts for each solar run, 20k*4 for coldcryo neutrons.
Nthrown = 1500*4
tsim8 = Nthrown/(N_Ar*x8tot*Flux8 )
wt8 = 3.14E7/tsim8
Nthrown = 1500*4
tsimh = Nthrown/(N_Ar*xhtot*Fluxh )
wth = 3.14E7/tsimh
Nthrown = 1500*4
tsimc = Nthrown/(N_Ar*xctot*Fluxc )
wtc = 3.14E7/tsimc


# neutrons
# SS has neutrons at level of ~1E-11 /cm3/sec from https://book
Nthrown = 20000*4
actnSS = 1E-11
VolSS = (12*12*0.02*2 + 12*58*0.02*4)*1E6 
tsimn = Nthrown/(VolSS*actnSS )
wtn = 3.14E7/tsimn

# Ar39
# UAr has 1Bq/kg/1500 of Ar39
Nthrown = 20000*2
actAr39 = 1.0/1500
tsimn = Nthrown/(kgAr*actAr39)
wta39 = 3.14E7/tsimn
# Ar42
# UAr has 1Bq/kg/1500 of Ar39, 10^-5 of that of Ar42
Nthrown = 20000*2
actAr42 = 9.2*1E-5/1500  # https://iopscience.iop.org/article/10.1088/1742-6596/718/6/062004/pdf
tsimn = Nthrown/(kgAr*actAr42)
wta42 = 3.14E7/tsimn
# Rn222
# 
Nthrown = 10000*4
actRn222 = 2.0E-6   # 10x worse than world-leading Darkside/kg
tsimn = Nthrown/(kgAr*actRn222)
wtr = 3.14E7/tsimn

h08ES.Add(h18ES)
h08ES.Add(h28ES)
h08ES.Add(h38ES)
h08CC.Add(h18CC)
h08CC.Add(h28CC)
h08CC.Add(h38CC)

h0hES.Add(h18ES)
h0hES.Add(h28ES)
h0hES.Add(h38ES)
h0hCC.Add(h18CC)
h0hCC.Add(h28CC)
h0hCC.Add(h38CC)

h0cES.Add(h1cES)
h0cES.Add(h2cES)
h0cES.Add(h3cES)
h0cCC.Add(h1cCC)
h0cCC.Add(h2cCC)
h0cCC.Add(h3cCC)

h0nCC.Add(h1nCC)
h0nCC.Add(h2nCC)
h0nCC.Add(h3nCC)

h0a39.Add(h1a39)
#h0a39.Add(h2a39)
#h0a39.Add(h3a39)

h0a42.Add(h1a42)
#h0a42.Add(h2a42)
#h0a42.Add(h3a42)


h0r.Add(h1r)
h0r.Add(h2r)
h0r.Add(h3r)

'''
hsES = THStack("hsES","stackedES");
hsES.Add(h0hES,wth)
hsES.Add(h08ES,wt8)
hsES.Add(h0cES,wtc)
hsES.Add(h18ES,wt8)
hsES.Add(h1hES,wtc)
hsES.Add(h1cES,wtc)
hsES.Add(h28ES,wt8)
hsES.Add(h2hES,wth)
hsES.Add(h2cES,wtc)
hsES.Add(h38ES,wt8)
hsES.Add(h3hES,wth)
hsES.Add(h3cES,wtc)
'''

Nb = h0hES.GetNbinsX()
n0ES = []; n1ES =[]; n2ES = []
n0CC = []; n1CC =[]; n2CC = []; n3CC = []
na39 = []; na42 = []
nr = []
wts = []

n0ESe = np.zeros((3,Nb))
n0CCe = np.zeros((4,Nb))


for ii in range(Nb):
    if h0hES.GetBinContent(ii):
        n0ES = np.concatenate((n0ES,np.repeat(h0hES.GetBinCenter(ii),int(h0hES.GetBinContent(ii)))),axis=0)
    if h08ES.GetBinContent(ii):
        n1ES = np.concatenate((n1ES,np.repeat(h08ES.GetBinCenter(ii),int(h08ES.GetBinContent(ii)))),axis=0)
    if h0cES.GetBinContent(ii):
        n2ES = np.concatenate((n2ES,np.repeat(h0cES.GetBinCenter(ii),int(h0cES.GetBinContent(ii)))),axis=0)
    n0ESe[0,ii]  = h0hES.GetBinError(ii)
    n0ESe[1,ii]  = h08ES.GetBinError(ii)
    n0ESe[2,ii]  = h0cES.GetBinError(ii)

    if h0hCC.GetBinContent(ii):
        n0CC = np.concatenate((n0CC,np.repeat(h0hCC.GetBinCenter(ii),int(h0hCC.GetBinContent(ii)))),axis=0)
    if h08CC.GetBinContent(ii):
        n1CC = np.concatenate((n1CC,np.repeat(h08CC.GetBinCenter(ii),int(h08CC.GetBinContent(ii)))),axis=0)
    if h0cCC.GetBinContent(ii):
        n2CC = np.concatenate((n2CC,np.repeat(h0cCC.GetBinCenter(ii),int(h0cCC.GetBinContent(ii)))),axis=0)
    if h0nCC.GetBinContent(ii):
        n3CC = np.concatenate((n3CC,np.repeat(h0nCC.GetBinCenter(ii),int(h0nCC.GetBinContent(ii)))),axis=0)
    if h0a39.GetBinContent(ii):
        na39 = np.concatenate((na39,np.repeat(h0a39.GetBinCenter(ii),int(h0a39.GetBinContent(ii)))),axis=0)
    if h0a42.GetBinContent(ii):
        na42 = np.concatenate((na42,np.repeat(h0a42.GetBinCenter(ii),int(h0a42.GetBinContent(ii)))),axis=0)
    if h0r.GetBinContent(ii):
        nr = np.concatenate((nr,np.repeat(h0r.GetBinCenter(ii),int(h0r.GetBinContent(ii)))),axis=0)

    n0CCe[0,ii]  = h0hCC.GetBinError(ii)
    n0CCe[1,ii]  = h08CC.GetBinError(ii)
    n0CCe[2,ii]  = h0cCC.GetBinError(ii)
    n0CCe[3,ii]  = h0nCC.GetBinError(ii)
#    n0ES[0,ii] = h0nES.GetBinContent(ii)


wts.append( np.ones((len(n3CC)))*wtn)
wts.append( np.ones((len(na42)))*wta42)
wts.append( np.ones((len(na39)))*wta39)
wts.append( np.ones((len(nr)))*wtr)
wts.append( np.ones((len(n0ES)))*wth)
wts.append( np.ones((len(n1ES)))*wt8)
wts.append( np.ones((len(n2ES)))*wtc)
wts.append( np.ones((len(n0CC)))*wth)
wts.append( np.ones((len(n1CC)))*wt8)
wts.append( np.ones((len(n2CC)))*wtc)


nESCCn = []
nESCCn.append(n3CC)
nESCCn.append(na42)
nESCCn.append(na39)
nESCCn.append(nr)
nESCCn.append(n0ES)
nESCCn.append(n1ES)
nESCCn.append(n2ES)
nESCCn.append(n0CC)
nESCCn.append(n1CC)
nESCCn.append(n2CC)


# Just up to 10 MeV
binsz = (h0hES.GetBinCenter(1) - h0hES.GetBinCenter(0)) ##* 5 # makes 20 bins in 20 MeV, down from 100 bins in original histos.

uedge = h0hES.GetBinLowEdge(Nb)+binsz
uedge = h0hES.GetBinLowEdge(Nb)/2+binsz # just to 10 MeV
#labelv = np.array(["neutrons","Ar42","Ar39","Rn","hepES","B8ES","CNOES","hepCC","8BCC","CNOCC"])
labelv = np.array(["neutrons " + str(actnSS) + "/cm3/s","Ar42 "+"Bq/L/1500/9E5","Ar39 " + "Bq/L/1500","Rn "+str(actRn222)+"Bq/kg","hepES","8BES","CNOES","hepCC","8BCC","CNOCC"])

# This little snippet swaps out the Rn for suppressed Ar42
suppress42 =  "Ar42 Bq/L/1500/5E8"
labelv[3] = suppress42
wts[3] = wts[1]/500.
nESCCn[3] = nESCCn[1]
linewidth = 4*np.ones(10)
linewidth[3] = 1
color = np.array(["blue","orange","green","navajowhite","salmon","lightgray","black","red","gray","black"])
cts, be, er = skh_plt.hist(nESCCn,weights=wts,errorbars=False, histtype='step',label=labelv,bins=np.arange(0.,uedge,binsz),color=color) #,stacked='true'
hdict = dict()
hdict["counts"]=cts
hdict["binedges"]=be
hdict["err"]=er

#fileout = "solar_neutrinos_9mhi-shinyg10-50mattn_10MeVmax"
fileout = "solar_neutrinos_9mhi-10MeVmax-shinyg10-50mattn"
plt.legend()
plt.title(fileout + ' Spectra in 6x9x20 m3')
plt.xlabel('Energy [MeV]')
plt.ylabel('Events / 1.6 kTonne-yr / ' + str(binsz) + ' MeV')
plt.yscale('log')              
#plt.show()

                                                                                                                                
plt.savefig(fileout+'.png')
plt.close()



# full 20 MeV
binsz = (h0hES.GetBinCenter(1) - h0hES.GetBinCenter(0)) * 5 # makes 20 bins in 20 MeV, down from 100 bins in original histos.
uedge = h0hES.GetBinLowEdge(Nb)+binsz
uedge = h0hES.GetBinLowEdge(Nb)+binsz 
#labelv = np.array(["neutrons","Ar42","Ar39","Rn","hepES","8BES","CNOES","hepCC","8BCC","CNOCC"])
labelv = np.array(["neutrons " + str(actnSS) + "/cm3/s","Ar42 "+"Bq/L/1500/9.2E5","Ar39 " + "Bq/L/1500","Rn "+str(actRn222)+"Bq/kg","hepES","8BES","CNOES","hepCC","8BCC","CNOCC"])
# snippet
labelv[3] = suppress42

cts, be, er = skh_plt.hist(nESCCn,weights=wts,errorbars=False, histtype='step',label=labelv,bins=np.arange(0.,uedge,binsz),color=color) #,stacked='true'
hdict = dict()
hdict["counts"]=cts
hdict["binedges"]=be
hdict["err"]=er

#fileout = "solar_neutrinos_9mhi-shinyg10-50mattn"
fileout = "solar_neutrinos_9mhi-shinyg10-50mattn"
plt.legend()
plt.title(fileout + ' Spectra in 6x9x20 m3')
plt.xlabel('Energy [MeV]')
plt.ylabel('Events / 1.6 kTonne-yr / ' + str(binsz) + ' MeV')
plt.yscale('log')              
#plt.show()
#pdb.set_trace()
                                                                                                                                
plt.savefig(fileout+'.png')
plt.close()


import numpy as np
from matplotlib import pyplot as plt
from ROOT import TFile, TH1, THStack
from skhep.visual import MplPlotter as skh_plt

import pdb

height = 4.5
height = 3.0
Nyr = 5
me = 0.511
ResConst = 0.06 ## presume ResConst/sqrt(E) charge resolution, so this is the sig/Q at Q.
Q = 2.458
fXe = 0.03
Nlife = 1.0E28 #  I use 1.0E28 or 29

Q208 = 2.614
Q214 = 2.4477
ResConst = 0.03
nu2bb = True
Nlife2nu = 2.16E22 * 0.69


def specGam(hin,Qgam):
    Nbins = hin.GetNbinsX()*40
    Maxx = hin.GetNbinsX()*hin.GetBinWidth(1)
    spec = np.zeros(Nbins)

    Nevts = hin.Integral()

    # spec = Int dE' N(E') S(E'-E), with N a delta function and S(E'-E) = N2 exp(-((E'-E)/sqrt2 sig)^2), sig = ResConst/sqrt(E')
    for ii in range(Nbins):
        E = Maxx*ii/Nbins
        sig = ResConst*E*np.sqrt(Q)/np.sqrt(E)
        # Int exp(-ax^2/a) = sqrt(2pi/a)
        spec[ii] = Nevts * np.sqrt(sig/2/np.pi) * np.exp( -1/2*((E-Qgam)/sig)**2 )

    spec = np.nan_to_num(spec)
    spec *= Nevts/spec.sum()*40. # *40 to keep the y-axis scaled to the per 40keV bin size along w other bkgds.

    return spec


def spec0nbb(hin):
    Nbins = hin.GetNbinsX()*40
    Maxx = hin.GetNbinsX()*hin.GetBinWidth(1)
    spec = np.zeros(Nbins)

    Ndec = Nyr/Nlife ## approx from 1.0 - exp(-Nyr/Nlife)
    Nevts = 8* (3*height*20) * 1E6*3.3/136*6.022E23 * fXe * Ndec

    # spec = Int dE' N(E') S(E'-E), with N a delta function and S(E'-E) = N2 exp(-((E'-E)/sqrt2 sig)^2), sig = ResConst/sqrt(E')
    for ii in range(Nbins):
        E = Maxx*ii/Nbins
        sig = ResConst * E*np.sqrt(Q)/np.sqrt(E)  #/np.sqrt(Q) * Q
        # Int exp(-ax^2/a) = sqrt(2pi/a)
        spec[ii] = Nevts * np.sqrt(sig/2/np.pi) * np.exp( -1/2*((E-Q)/sig)**2 )

    spec = np.nan_to_num(spec)
    spec *= Nevts/spec.sum()*40. # *40 to keep the y-axis scaled to the per 40keV bin size along w other bkgds.

    return spec


def spec2nbb(hin):
    Nbins = hin.GetNbinsX()*40
    Maxx = hin.GetNbinsX()*hin.GetBinWidth(1)
    spec = np.zeros(Nbins)
    specth = np.zeros(Nbins)

    Ndec = Nyr/Nlife2nu ## approx from 1.0 - exp(-Nyr/Nlife)
    Nevts = 8* (3*height*20) * 1E6*3.3/136*6.022E23 * fXe * Ndec

    T0 = Q # - 2*me
    
    for ii in range(Nbins):
        K = Maxx*ii/Nbins
        # Int exp(-ax^2/a) = sqrt(2pi/a)

        if T0>K:
            # https://www.annualreviews.org/doi/10.1146/annurev-nucl-102711-094904
            specth[ii] =  K * (T0-K)**5 * (1 + 2*K + 4*K**2/3. + K**3/3. + K**4/30.)

        # convolution for energy-dependent energy smearing
        sig =  ResConst * K*np.sqrt(Q)/np.sqrt(K)
        for jj in range(-200*40,200*40,1):
            tau = (ii-jj)*Maxx/Nbins
            gauss = 1/np.sqrt(2*np.pi)/sig * np.exp( -tau**2/2/sig**2 )
            spec[ii] += specth[ii] * gauss

        if not ii%400:
            print ("spec2nbb is on " + str(ii) + "th iteration of " + str(Nbins))

    spec = np.nan_to_num(spec)            
    spec *= Nevts/spec.sum()*40. # *40 to keep the y-axis scaled to the per 40keV bin size along w other bkgds.
    return spec


    
f0n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_0.root")
#f0a42 = TFile("/Volumes/EC2TB/G4/mod4/Ar42_fidv_2.root")
f08 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_8b_shinyg10_0.root")
f0g = TFile("/Volumes/EC2TB/G4/mod4/tl208_G10_0.root")
#f0r = TFile("/Volumes/EC2TB/G4/mod4/Rn222_fidv_shinyg10_0.root")

f1n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_1.root")
f1a42 = TFile("/Volumes/EC2TB/G4/mod4/Ar42_fidv_1.root")
f18 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_8b_shinyg10_1.root")
f1g = TFile("/Volumes/EC2TB/G4/mod4/tl208_G10_1.root")

f2n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_2.root")
f2a42 = TFile("/Volumes/EC2TB/G4/mod4/Ar42_fidv_2.root")
f28 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_8b_shinyg10_2.root")
f2g = TFile("/Volumes/EC2TB/G4/mod4/tl208_G10_2.root")

f3n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_3.root")
f3a42 = TFile("/Volumes/EC2TB/G4/mod4/Ar42_fidv_3.root")
f38 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_8b_shinyg10_3.root")
f3g = TFile("/Volumes/EC2TB/G4/mod4/tl208_G10_3.root")

f4n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_4.root")
f5n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_5.root")
f6n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_6.root")
f7n = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_7.root")

# call a function to create with ~8% res, say.
h08ES = f08.H22
h08CC = f08.H23
h0g = f0g.H23
h0n = f0n.H23

h18ES = f18.H22
h18CC = f18.H23
h1a42 = f1a42.H23
h1g = f1g.H23
h1n = f1n.H23

h28ES = f28.H22
h28CC = f28.H23
h2a42 = f2a42.H23
h2g = f2g.H23
h2n = f2n.H23

h38ES = f38.H22
h38CC = f38.H23
h3a42 = f3a42.H23
h3g = f3g.H23
h3n = f3n.H23

h4n = f4n.H23
h5n = f5n.H23
h6n = f6n.H23
h7n = f7n.H23



hbb = spec0nbb(h0n) # construct this one from a guassian of width Eres. Only pass and use h0n for binning properties.
h2bb = spec2nbb(h0n)



# Calculate normalizaions, ala xs_flux-avgd * N_Ar * Flux_int
N_Ar = 1.4 * 6*2*height*40 * 1E6 /40. *6.022E23
kgAr = (6*2*height*40)*1E6 * 1.4 * 1E-3

# xs from Marley runs' log files
x8tot = 0.0406549 * 1E-40

# From https://arxiv.org/pdf/astro-ph/0402114.pdf, Bahcall/Pinsonneault
# Tacking on rough oscillation factors as from Fig 5 of https://www.sciencedirect.com/science/article/pii/S2212686414000211
Flux8 = 5.79E6 * 0.3

#Normalize 8B events
Nthrown = 5000*4
tsim8 = Nthrown/(N_Ar*x8tot*Flux8 )
wt8 = Nyr*3.14E7/tsim8
if height==3.0:
    wt8 /= 2. # Ckv/Scint single e vs ee, Brodsky et al say 3x w 75% efficiency 


# neutrons
# SS has neutrons at level of ~1E-11 /cm3/sec from https://book
Nthrown = 10000*8
actnSS = 1E-11
VolSS = (12*12*0.02*2 + 12*58*0.02*4)*1E6  # 2cm thick
tsimn = Nthrown/(VolSS*actnSS )
wtn = Nyr*3.14E7/tsimn

# Ar42
# UAr has 1Bq/kg/1500 of Ar39, 10^-5 of that of Ar42
Nthrown = 20000*3
actAr42 = 9.2*1E-5/1500/500.  # https://iopscience.iop.org/article/10.1088/1742-6596/718/6/062004/pdf and Henning's approved UAr42 reduction
if height==3.0:
    actAr42 /= 10. # aggro, not signed off yet by Henning
    actAr42 /= 2. # Ckv/Scint single e vs ee 

tsima42 = Nthrown/(kgAr*actAr42)
wta42 = Nyr*3.14E7/tsima42

# Tl208: 36% of Th232 at 50 mBq/kq
# 0.050 Bq/kg *0.36* 6*60*0.003*1E6*1.65 *1E-3 * piE7 is 1.007E9 decays.
# This presumes 3mm G10 covers full top/bottom of cryostat.


Nthrown = 423333*4
kgG10 = 2*6*60*0.003*1E6*1.65*1E-3
act208 = 50.0E-3 * 0.36
tsim208 = Nthrown/(kgG10*act208)
wt208 = Nyr*3.14E7/tsim208


# Rn222
# 
Nthrown = 10000*4
actRn222 = 2.0E-6   # 10x worse than world-leading Darkside/kg
tsimn = Nthrown/(kgAr*actRn222)
wtr = Nyr*3.14E7/tsimn


h08ES.Add(h18ES)
h08ES.Add(h28ES)
h08ES.Add(h38ES)
h08CC.Add(h18CC)
h08CC.Add(h28CC)
h08CC.Add(h38CC)

h0n.Add(h1n)
h0n.Add(h2n)
h0n.Add(h3n)
h0n.Add(h4n)
h0n.Add(h5n)
h0n.Add(h6n)
h0n.Add(h7n)

h1a42.Add(h2a42)
h1a42.Add(h3a42)

h0g.Add(h1g)
h0g.Add(h2g)
h0g.Add(h3g)


#  smoothe the sparse gamma SiPM counts bins above 1.2 MeV, then calculate from h0g the number expected for y=3m fidv. 
# suppressed by exp(-150cm/33cm). 150cm is the extra length needed to travel. 33cm is the compton int length in Ar.
h0gsum = h0g.Integral() # sum in a 1.2 to 2.6 MeV:
fidscale = np.exp(-150./33.)
nbinsgt12 = (2.6-1.2)/h0g.GetBinWidth(1)
countspb  = h0gsum/nbinsgt12
'''
for ii in range(h0g.GetNbinsX()):
    if h0g.GetBinCenter(ii)>1.2:
        if height==3.0:
            h0g.SetBinContent(ii,countspb*fidscale)
'''

## We will treat Tl208 as a Gaussian with Resolution Res around 2.614 MeV, and not use the simulation SiPM histograms other than for normalization.
h0g.Scale(wt208/h0g.Integral()*fidscale)
h208 = specGam(h0g,Q208) # This function actually uses h0g.Integral()


Nb = h08ES.GetNbinsX()
n1ES = []; n1CC =[];
na42 = []
nn = []
ng208 = []
nr = []
nbb = []
n2bb = []


# Stuff the Root Histos into  arrays, so that I can histo them with skhep_hist()
for ii in range(Nb):
    if h08ES.GetBinContent(ii):
        n1ES = np.concatenate((n1ES,np.repeat(h08ES.GetBinCenter(ii),int(h08ES.GetBinContent(ii)))),axis=0)
    if h08CC.GetBinContent(ii):
        n1CC = np.concatenate((n1CC,np.repeat(h08CC.GetBinCenter(ii),int(h08CC.GetBinContent(ii)))),axis=0)
    if h0n.GetBinContent(ii):
        nn = np.concatenate((nn,np.repeat(h0n.GetBinCenter(ii),int(h0n.GetBinContent(ii)))),axis=0)
    if h1a42.GetBinContent(ii):
        na42 = np.concatenate((na42,np.repeat(h1a42.GetBinCenter(ii),int(h1a42.GetBinContent(ii)))),axis=0)
    if h0g.GetBinContent(ii):
        ng208 = np.concatenate((ng208,np.repeat(h0g.GetBinCenter(ii),int(h0g.GetBinContent(ii)))),axis=0)
#    nbb = np.concatenate((nbb,np.repeat(h0g.GetBinCenter(ii),int(hbb[ii]))),axis=0)
#    n2bb = np.concatenate((n2bb,np.repeat(h0g.GetBinCenter(ii),int(h2bb[ii]))),axis=0)

    
wts = []
#wts.append( np.ones((len(nbb)))*1.0)
#wts.append( np.ones((len(n2bb)))*1.0)  ## added this
wts.append( np.ones((len(nn)))*wtn)
wts.append( np.ones((len(na42)))*wta42)
#wts.append( np.ones((len(nr)))*wtr)
wts.append( np.ones((len(n1ES)))*wt8)
wts.append( np.ones((len(n1CC)))*wt8)
##wts.append( np.ones((len(ng208)))*wt208)
if nu2bb:
    wts = []
#    wts.append( np.ones((len(nbb)))*1.0)
#    wts.append( np.ones((len(n2bb)))*1.0)  ## added this
    wts.append( np.ones((len(na42)))*wta42)
    wts.append( np.ones((len(n1ES)))*wt8)
##    wts.append( np.ones((len(ng208)))*wt208) ## lose this and use sepc208 instead


# Above does not work for the rescaled h0g histogram, so do it by itself here:
# making the rescaling code above un-needed!
if height==3.0:
    ng208 = []
    for ii in range(Nb):
        if h0g.GetBinCenter(ii)>=1.2 and h0g.GetBinCenter(ii)<=2.6:
            ng208 = np.concatenate((ng208,np.repeat(h0g.GetBinCenter(ii),1.)),axis=0)

## using spec208 instead
##wts[-1] = np.ones(len(ng208))*wt208*countspb*fidscale
##wt208 = wts[-1][0]

nESCCn = []
#nESCCn.append(nbb)
#nESCCn.append(n2bb)  ## added this
nESCCn.append(nn)
nESCCn.append(na42)
#nESCCn.append(nr)
nESCCn.append(n1ES)
nESCCn.append(n1CC)
nESCCn.append(ng208)
if nu2bb:
    nESCCn = []
#    nESCCn.append(nbb)
#    nESCCn.append(n2bb)  ## added this
    nESCCn.append(na42)
    nESCCn.append(n1ES)
##    nESCCn.append(ng208)
    


# Just up to 10 MeV
binsz = (h08ES.GetBinCenter(1) - h08ES.GetBinCenter(0)) 

uedge = h08ES.GetBinLowEdge(Nb)+binsz

exposure = "tau="+str(Nlife) + " yrs, " + str(int(fXe*100))+ "%Xe136"
labelv = np.array(["0nubb " + exposure,"neutrons " + str(actnSS) + "/cm3/s SS","Ar42 "+"Bq/L/1500/5E8","8BES","8BCC","Tl208  50 mBq/kg G10"])
if height == 3.0:
    labelv = np.array(["0nubb " + exposure,"neutrons " + str(actnSS) + "/cm3/s SS","Ar42 "+"Bq/L/1500/5E8/10","8BES","8BCC"])##,"Tl208  5 mBq/kg G10"])
#linewidth = 4*np.ones(10)
#linewidth[3] = 1
color = np.array(["black","blue","orange","navajowhite","salmon"])##,"green"])

exp2nubb = "tau="+str(Nlife2nu) + " yrs, " + str(int(fXe*100))+ "%Xe136"
if nu2bb:
    labelv = np.array(["Ar42 "+"Bq/L/1500/5E8/10","8BES"])##,"Tl208  5 mBq/kg G10"])
    color = np.array(["orange","navajowhite"])##,"green"])


cts, be, er = skh_plt.hist(nESCCn,weights=wts,errorbars=False, histtype='step',label=labelv,bins=np.arange(0.,uedge,binsz),color=color)#,stacked='true')

hdict = dict()
hdict["counts"]=cts
hdict["binedges"]=be
hdict["err"]=er

binpk = int(Nb/uedge*Q)+1

sig = 0
bkgd = 0
plt.plot(np.arange(0,4.000,0.001),hbb,label='0nubb: '+exposure,color='red')
plt.plot(np.arange(0,4.000,0.001),h2bb,label='2nubb',color='blue')
plt.plot(np.arange(0,4.000,0.001),h208,label='Tl208 50 mBq/kg',color='green')

fileout = "Xe1360nubb_res"+str(100*ResConst)+"_"+str(height)+"m"+"_0nbbtau_"+str(Nlife)

#plt.title(fileout + ' Spectra in 6x'+str(2*height)+'x40 m3')
plt.xlabel('Energy [MeV]')
if height==4.5:
    plt.ylabel('Events / 0.284 kTonneXe136-'+str(Nyr)+'yr / ' + str(binsz) + ' MeV')
if height==3.0:
    plt.ylabel('Events / 0.189 kTonneXe136-'+str(Nyr)+'yr / ' + str(binsz) + ' MeV')
plt.yscale('log')              




if height==3.0:
    for ii in range(binpk-1,binpk+2,1):
        sig += hbb[ii]
        if h0g.GetBinCenter(ii)>=1.2 and h0g.GetBinCenter(ii)<=2.6:
##            bkgd += wt208 + wta42*h1a42.GetBinContent(ii)  + wt8*h08ES.GetBinContent(ii)
            bkgd +=  wta42*h1a42.GetBinContent(ii)  + wt8*h08ES.GetBinContent(ii)
        else:
            bkgd += wta42*h1a42.GetBinContent(ii)  + wt8*h08ES.GetBinContent(ii)
#        print("sig/bkd bin: " + str(ii))
    plt.xlim((2.0,3.0))

    if nu2bb:

#    plt.text(2.7,1.0,"s/sqrt(b) = " + str(round(sig/np.sqrt(bkgd),3)), fontsize=11)
        plt.ylim((0.1,300.0))
        plt.yscale('linear')
        plt.axvspan(Q-2*ResConst*Q, Q+2*ResConst*Q, color='y', alpha=0.5, lw=0)


if height==3.0:
    hbkgd = []
    for ii in range(binpk-6,binpk+5,1):
        if h0g.GetBinCenter(ii)>=1.2 and h0g.GetBinCenter(ii)<=2.6 and not ii==binpk+5-1:
##            hbkgd.append(wt208 + wta42*h1a42.GetBinContent(ii)  + wt8*h08ES.GetBinContent(ii))
            hbkgd.append(wta42*h1a42.GetBinContent(ii)  + wt8*h08ES.GetBinContent(ii))
        elif not ii==binpk+5-1:
##            hbkgd.append(wta42*h1a42.GetBinContent(ii)  + wt8*h08ES.GetBinContent(ii))
            hbkgd.append(wta42*h1a42.GetBinContent(ii)  )
    hc2bb = []; hc0bb = []
    hc208 = [];
    for ii in range(int(binpk*40-5.5*40),int(binpk*40+5.5*40),1):

        if not (ii-20)%40:
            sum2bb = 0; sum0bb = 0; sum208 = 0;
        sum2bb+=h2bb[ii]
        sum0bb+=hbb[ii]
        sum208+=h208[ii]
        if not ii==int(binpk*40-5.5*40) and not (ii-20)%40:
            hc2bb.append(sum2bb)
            hc0bb.append(sum0bb)
            hc208.append(sum208)


    nbkgdat = np.random.choice(np.arange(uedge/Nb*(binpk-6),uedge/Nb*(binpk+5),uedge/Nb)[:-1]+uedge/Nb/2.,int(np.array(hbkgd).sum()),p=np.array(hbkgd)/np.array(hbkgd).sum())
    nbkg208  = np.random.choice(np.arange(uedge/Nb*(binpk-6),uedge/Nb*(binpk+5),uedge/Nb)[:-1]+uedge/Nb/2.,int(np.array(hc208).sum()),p=np.array(hc208)/np.array(hc208).sum())
    bkg2nu  = np.random.choice(np.arange(uedge/Nb*(binpk-5.5),uedge/Nb*(binpk+5.5),uedge/Nb)[:-1]+uedge/Nb/2.,int(np.array(hc2bb).sum()),p=np.array(hc2bb)/np.array(hc2bb).sum())
    bkg0nu  = np.random.choice(np.arange(uedge/Nb*(binpk-5.5),uedge/Nb*(binpk+5.5),uedge/Nb)[:-1]+uedge/Nb/2.,int(np.array(hc0bb).sum()),p=np.array(hc0bb)/np.array(hc0bb).sum())


    ctbgd,_ = np.histogram(nbkgdat,bins=np.arange(uedge/Nb*(binpk-6),uedge/Nb*(binpk+5),uedge/Nb))
    ct208,_ = np.histogram(nbkg208,bins=np.arange(uedge/Nb*(binpk-6),uedge/Nb*(binpk+5),uedge/Nb))
    ct2nu,_ = np.histogram(bkg2nu,bins=np.arange(uedge/Nb*(binpk-6),uedge/Nb*(binpk+5),uedge/Nb))
    ct0nu,_ = np.histogram(bkg0nu,bins=np.arange(uedge/Nb*(binpk-6),uedge/Nb*(binpk+5),uedge/Nb))
    

    plt.errorbar(np.arange(uedge/Nb*(binpk-6),uedge/Nb*(binpk+5),uedge/Nb)[:-1]+uedge/Nb/2.,ctbgd+ct2nu+ct0nu,xerr=None,yerr=np.sqrt(ctbgd+ct208+ct2nu+ct0nu),color='black',ls='none',marker='.',label='pseudo-data')
    if nu2bb:
        plt.xlim((2.3,2.64))
    
plt.legend()#loc='lower left')    
# Now pick random #s from each histogram.

plt.savefig(fileout+'-notitle.png')

plt.close()



#### Make this plot again with 1/2 bin size. #####

cts, be, er = skh_plt.hist(nESCCn,weights=[x/2. for x in wts],errorbars=False, histtype='step',label=labelv,bins=np.arange(0.,uedge,binsz),color=color)#,stacked='true')

hdict = dict()
hdict["counts"]=cts
hdict["binedges"]=be
hdict["err"]=er

binpk = int(Nb/uedge*Q)+1

sig = 0
bkgd = 0
plt.plot(np.arange(0,4.000,0.001),hbb/2.,label='0nubb',color='red')
plt.plot(np.arange(0,4.000,0.001),h2bb/2.,label='2nubb',color='blue')
plt.plot(np.arange(0,4.000,0.001),h208/2.,label='Tl208 50 mBq/kg',color='green')
    
fileout += "_20keVbins"

#plt.title(fileout + ' Spectra in 6x'+str(2*height)+'x40 m3')
plt.xlabel('Energy [MeV]')
if height==4.5:
    plt.ylabel('Events / 0.284 kTonneXe136-'+str(Nyr)+'yr / ' + str(binsz/2.) + ' MeV')
if height==3.0:
    plt.ylabel('Events / 0.189 kTonneXe136-'+str(Nyr)+'yr / ' + str(binsz/2.) + ' MeV')


        
if height==3.0:
    
    plt.ylim((0.1,150.0))
    plt.yscale('linear')
    plt.axvspan(Q-2*ResConst*Q, Q+2*ResConst*Q, color='y', alpha=0.5, lw=0)

    hbkgd = []
    for ii in range(binpk-6,binpk+5,1):
        if h0g.GetBinCenter(ii)>=1.2 and h0g.GetBinCenter(ii)<=2.6 and not ii==binpk+10-1:
##            hbkgd.append(wt208/2. + wta42*h1a42.GetBinContent(ii)/2.  + wt8*h08ES.GetBinContent(ii)/2.)
##            hbkgd.append(wt208/2. + wta42*h1a42.GetBinContent(ii)/2.  + wt8*h08ES.GetBinContent(ii)/2.)
            hbkgd.append(wta42*h1a42.GetBinContent(ii)/2.  + wt8*h08ES.GetBinContent(ii)/2.)
            hbkgd.append(wta42*h1a42.GetBinContent(ii)/2.  + wt8*h08ES.GetBinContent(ii)/2.)
        elif not ii==binpk+10-1:
##            hbkgd.append(wta42*h1a42.GetBinContent(ii)/2. + wt8*h08ES.GetBinContent(ii)/2.)
##            hbkgd.append(wta42*h1a42.GetBinContent(ii)/2. + wt8*h08ES.GetBinContent(ii)/2.)
            hbkgd.append(wta42*h1a42.GetBinContent(ii)/2.)
            hbkgd.append(wta42*h1a42.GetBinContent(ii)/2.)

    binpk *= 2
    hc2bb = []; hc0bb = []
    hc208 = []
    for ii in range(int(binpk*20-11*20),int(binpk*20+11*20),1):

        if not (ii-10)%20:
            sum2bb = 0; sum0bb = 0; sum208 = 0
        sum2bb+=h2bb[ii]/2. # scale down from 40 keV bins to 20 keV
        sum0bb+=hbb[ii]/2.
        sum208+=h208[ii]/2.
        if not ii==int(binpk*20-11*20) and not (ii-10)%20:
            hc2bb.append(sum2bb)
            hc0bb.append(sum0bb)
            hc208.append(sum208)


    Nb *= 2



    nbkgdat = np.random.choice(np.arange(uedge/Nb*(binpk-12),uedge/Nb*(binpk+11),uedge/Nb)[:-1]+uedge/Nb/2.,int(np.array(hbkgd).sum()),p=np.array(hbkgd)/np.array(hbkgd).sum())

    bkg2nu  = np.random.choice(np.arange(uedge/Nb*(binpk-12),uedge/Nb*(binpk+11),uedge/Nb)[:-1]+uedge/Nb/2.,int(np.array(hc2bb).sum()),p=np.array(hc2bb)/np.array(hc2bb).sum())
    bkg0nu  = np.random.choice(np.arange(uedge/Nb*(binpk-12),uedge/Nb*(binpk+11),uedge/Nb)[:-1]+uedge/Nb/2.,int(np.array(hc0bb).sum()),p=np.array(hc0bb)/np.array(hc0bb).sum())
    bkg208  = np.random.choice(np.arange(uedge/Nb*(binpk-12),uedge/Nb*(binpk+11),uedge/Nb)[:-1]+uedge/Nb/2.,int(np.array(hc208).sum()),p=np.array(hc208)/np.array(hc208).sum())                      

#    binpk /= 2
    ctbgd,_ = np.histogram(nbkgdat,bins=np.arange(uedge/Nb*(binpk-12),uedge/Nb*(binpk+11),uedge/Nb))
#    binpk *= 2
    ct208,_ = np.histogram(nbkg208,bins=np.arange(uedge/Nb*(binpk-12),uedge/Nb*(binpk+11),uedge/Nb))
    ct2nu,_ = np.histogram(bkg2nu,bins=np.arange(uedge/Nb*(binpk-12),uedge/Nb*(binpk+11),uedge/Nb))
    ct0nu,_ = np.histogram(bkg0nu,bins=np.arange(uedge/Nb*(binpk-12),uedge/Nb*(binpk+11),uedge/Nb))
    


    plt.errorbar(np.arange(uedge/Nb*(binpk-11),uedge/Nb*(binpk+12),uedge/Nb)[:-1]+uedge/Nb/2.,ctbgd+ct208+ct2nu+ct0nu,xerr=None,yerr=np.sqrt(ctbgd+ct208+ct2nu+ct0nu),color='black',ls='none',marker='.',label='pseudo-data')

    if nu2bb:
        plt.xlim((2.3,2.64))

plt.legend()#loc='lower left')    

plt.savefig(fileout+'-notitle.png')

plt.close()

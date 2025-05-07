import numpy as np
from matplotlib import pyplot as plt
from ROOT import TFile, TH1, THStack

## must have done: pip install scikit_hep==3.1. Later versions don't understand skhep.visual. EC, 6-Oct-2021
#from skhep.visual import MplPlotter as skh_plt

import pdb


fnsnc = TFile("NeutExt7M_nsnc.root") ## 700k
fns = TFile("NeutExt7M_ns.root")
f32H2O = TFile("NeutExt7M_32cmH2O.root")
f32BP = TFile("NeutExt7M_32cmBP.root")
f20BPSE = TFile("NeutExt7M_20cmBPSE.root")
fIB = TFile("NeutExt700_ib.root") ## 700k
froof = TFile("Neut700_roof.root")
f8B = TFile("marley_optphys_b8.root")

## "s" for spectrum
snsnc = []
sns = []
s32H2O = []
s32BP = []
s20BPSE = []
sIB= []
s8B = []
sn = []
wts = []

fvcuts = "abs(entry.X)<3.0E3 and abs(entry.Y)<3.0E3"


ttnsnc = fnsnc["T3"]
print("Looping on NSNC")
for entry in ttnsnc:

    if abs(entry.Z)<24E3 and eval(fvcuts) and entry.EnergyDepEvt>0.:
        snsnc.append(entry.EnergyDepEvt)

sn.append(snsnc)
Nthrown = ttnsnc.GetEntries()
actnn = 2.94E-6 # /cm2 /sec
Voln = (18*18*2 + 18*72*4)*1E4 ## scaling all 6 sides to the one (possibly shielded) long y-z side. Makes sense for ns, nsnc runs
tsimn = Nthrown/(Voln*actnn )
wtnsnc = 3.14E7/tsimn

wts.append(np.ones(len(snsnc)) * wtnsnc)

ttns = fns["T3"]
print("Looping on NS")
for entry in ttns:

    if abs(entry.Z)<24E3 and eval(fvcuts) and  entry.EnergyDepEvt>0.:
        sns.append(entry.EnergyDepEvt)
        
sn.append(sns)

Nthrown = ttns.GetEntries()
Voln = (18*18*2 + 18*72*4)*1E4 ## if 4 (3)scaling all 6 (just 5, will add the roof separately) sides to the one (possibly shielded) long y-z side. Makes sense for ns, nsnc runs
tsimn = Nthrown/(Voln*actnn )
wtns = 3.14E7/tsimn
wts.append(np.ones(len(sns)) * wtns)

''' now read in and calculate wts for roof. Add the spectrum to all the shielding results, properly wtd.'''
ttroof = froof["T3"]
sroof = []
wtroof = 0.0
### comment this section out to exclude roof. Change 3s to 4s below
'''
print("Looping on Roof")
for entry in ttroof:

    if abs(entry.Z)<24.E3 and eval(fvcuts) and  entry.EnergyDepEvt>0.:
        sroof.append(entry.EnergyDepEvt)

Nthrown = ttroof.GetEntries()
actnn = 2.94E-6 # /cm2 /sec ## Same for roof as sides,floor
Voln = (18*72*1)*1E4 ## Just one side
tsimn = Nthrown/(Voln*actnn )
wtroof = 3.14E7/tsimn
'''


tt = f32H2O["T3"]
print("Looping on 32H2O")
for entry in tt:

    if abs(entry.Z)<24E3 and eval(fvcuts) and entry.EnergyDepEvt>0.:
        s32H2O.append(entry.EnergyDepEvt)


Nthrown = tt.GetEntries()
Voln = (18*18*2 + 18*72*4)*1E4 ## scaling  5 sides to the one (shielded) long y-z side.
tsimn = Nthrown/(Voln*actnn )
wt32H2O = 3.14E7/tsimn
wts.append(np.concatenate((np.ones(len(s32H2O)) * wt32H2O, np.ones(len(sroof)) * wtroof)))
s32H2O.extend(sroof)
sn.append(s32H2O)


tt = f32BP["T3"]
print("Looping on 32BP")
for entry in tt:

    if abs(entry.Z)<24E3 and eval(fvcuts) and entry.EnergyDepEvt>0.:
        s32BP.append(entry.EnergyDepEvt)
Nthrown = tt.GetEntries()
Voln = (18*18*2 + 18*72*4)*1E4 ## scaling  5 sides to the one (shielded) long y-z side.           
tsimn = Nthrown/(Voln*actnn )
wt32BP = 3.14E7/tsimn
wts.append(np.concatenate((np.ones(len(s32BP)) * wt32BP, np.ones(len(sroof)) * wtroof)))
s32BP.extend(sroof)
sn.append(s32BP)


tt = f20BPSE["T3"]
print("Looping on 20BPSE")
for entry in tt:

    if abs(entry.Z)<24E3 and eval(fvcuts) and entry.EnergyDepEvt>0.:
        s20BPSE.append(entry.EnergyDepEvt)
Nthrown = tt.GetEntries()
Voln = (18*18*2 + 18*72*4)*1E4 ## scaling  5 sides to the one (shielded) long y-z side.           
tsimn = Nthrown/(Voln*actnn )
wt20BPSE = 3.14E7/tsimn
wts.append(np.concatenate((np.ones(len(s20BPSE)) * wt20BPSE, np.ones(len(sroof)) * wtroof)))           
s20BPSE.extend(sroof)
sn.append(s20BPSE)


tt = f8B["T3"]
print("Looping on 8B")
for entry in tt:

    if abs(entry.Z)<24E3 and eval(fvcuts) and entry.EnergyDepEvt>0.:
        s8B.append(entry.EnergyDepEvt)

Nthrown = tt.GetEntries()
NAr = 1.4 * 10*10*56 * 1E6 /40. * 6.022E23
NAr = 1.4 * 6*6*48 * 1E6 /40. * 6.022E23
xs8Btot = 0.0406549 * 1E-40 # cm^2
fl8B = 5.79E6 * 0.3 # /1cm2/sec
tsimn = Nthrown/(xs8Btot*NAr*fl8B )
wt8B = 3.14E7/tsimn
wts.append(np.ones(len(s8B)) * wt8B)
sn.append(s8B)

'''
tt = fIB["T3"]
print("Looping on IBeams")
for entry in tt:

    if abs(entry.Z)<28E3 and eval(fvcuts) and entry.EnergyDepEvt>0.:
        sIB.append(entry.EnergyDepEvt)

Nthrown = tt.GetEntries()
actIB = 1E-11 # /cm3 / sec
VolIB = (1.*15*0.005*39)*5*1E6 ## about 5 of these long walls of IBeams, equivalently
tsimn = Nthrown/(VolIB*actIB )
wtIB = 3.14E7/tsimn
wts.append(np.ones(len(sIB)) * wtIB)
sn.append(sIB)
'''

kT = 12*12*56*1.4/1E3
kT = 8*8*56*1.4/1E3
kT = 10*10*56*1.4/1E3
kT = 10*10*48*1.4/1E3
kT = 6*6*48*1.4/1E3
bins = np.arange(0,11.0,0.25)
bsz = bins[1]-bins[0]

plt.hist(sn,bins=bins,histtype='step',label=("NSNC","NS","32cmH2O","32cmBP","20cmBPSE","8B"),weights=wts)
plt.xlabel('Energy [MeV]')
plt.ylabel('Events / ' +str(round(kT,4))+ ' kTonne-yr / ' + str(bsz) + ' MeV')
plt.yscale('log')  
plt.legend()
#plt.show()
#pdb.set_trace()

if fvcuts:
    cuts = "_xy10"
    # cuts = "_xy8"
    cuts = "_xy6"
else:
    cuts = ""
    
plt.savefig("./Shielding-6sides_z48"+str(cuts)+".png")


pass

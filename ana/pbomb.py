from ROOT import TFile,TH1F
import numpy as np
from matplotlib import pyplot as plt
import uproot
from skhep.visual import MplPlotter as skh_plt

import pdb



files = ["../build/pbomb_optphys_40.root","../build/pbomb_optphys_60.root","../build/pbomb_optphys_80.root","../build/pbomb_optphys_40_nopc.root","../build/pbomb_optphys_60_nopc.root","../build/pbomb_optphys_80_nopc.root"]
files = ["../build/pbomb_optphys.root"]  #~3100 evts 80%pcc

qe = 1.

xf = []; yf = []; zf = []
Nevts = None
sumgvv = []
weightsv = []
weightsvf = []
sumgvvf = []
labelv = []
#h3 = TH3F("pr","photresp",3,0.,3.,3,0.,6.,3,0.,20.) # bins of 3 in each half of x,y,z

for file in files:

    print ("Reading " + str(file))

    f = uproot.open(file)
    steps = f['Steps']
    steps = steps.arrays(namedecode='utf-8')
    Nsteps = len(steps['Time'])
    trks = f['Tracks']
    trks = trks.arrays(namedecode='utf-8')
    Nphot = 1250 # len(trks['Event'])
    Nevts = trks['Event'].max()

    eventO = -12
    sumg = 0
    sumgv = []
    sumgf = 0
    sumgvf = []
    fidV = 0
    cntf = 0
    newEvt = 0
    xprim = yprim = zprim = None

    for step in range(Nsteps):
#    for step in range(400000):

        # below are  the Primaries (SProcess==null).  Also, insist that we're onto a new event.

        if  not fidV and b'null' in steps['SProcess'][step] :
#            fidV = abs(steps['X'][step]) < 2200. and ( ( steps['Y'][step] > 1000. and steps['Y'][step] < 3000.) or (steps['Y'][step] < -1000. and steps['Y'][step] > -3000. ) ) and abs(steps['Z'][step]) < 27000.  
            fidV = abs(steps['X'][step]) < 3000. and abs(steps['Y'][step]) < 12000. and abs(steps['Z'][step]) < 20000.  
            if fidV:
                cntf += 1
#                print ("Event/trkID/stepNum: " + str(steps['Event'][step]) +"/" + str(steps['TrkID'][step]) + "/" + str(steps['StepNum'][step]))
## Need to loop over all steps in an event and get the sumg for each
                xprim = steps['X'][step]
                yprim = steps['Y'][step]
                zprim = steps['Z'][step]

        if b'SiPM' in steps['TVolume'][step]: 
            sumg +=  1

        if b'SiPM' in steps['TVolume'][step]  and fidV: 
            sumgf +=  1
#           h3.Fill(abs(steps['X'][step]),steps['Y'][step],steps['Z'][step],1./0.1/Nevts) # detected photons/0.1 MeV/Nevents
            xf.append(abs(xprim)/1000.) # to meters
            yf.append(abs(yprim)/1000.)
            zf.append(abs(zprim)/1000.)


        if ((eventO != steps['Event'][step] and steps['Event'][step]>0)  or  step>=(Nsteps-1) ):
            newEvt = 1
            sumgv.append(sumg)
            sumg = 0
            if fidV:
                sumgvf.append(sumgf)
            sumgf = 0
            fidV = 0

            
        eventO = steps['Event'][step] 
        newEvt = 0

    weights =  np.repeat(1.0/Nevts,len(sumgv))
    weightsv.append(weights)
    weightsf =  np.repeat(1.0/Nevts,len(sumgvf))
    weightsvf.append(weightsf)
    sumgvv.append(sumgv)
    sumgvvf.append(sumgvf)
    sumgv = []
    sumgvf = []
    labelv.append("-".join(file.split(".root")[0].split("_")[2:]))
    print ("total in fv: " + str(cntf))
    cntf



# The photon bombs we're using here are 1290 photons, as comes from 100 keV n.r. This is equivalent to x=0.0716 MeV ee.
# 1290 photons = x . 24000 . 0.75  ... where the last factor is Birk's
wt = 1.0/0.0716/((Nevts+1.)/27) # per MeV per event, presuming equal distribution of launched pbombs over 27 xyz boxes
pperMeV = np.histogramdd(np.array([xf,yf,zf]).T,bins=(np.arange(0.,3.1,1.),np.arange(0.,6.1,2.),np.arange(0.,20.1,6.6)),weights=np.ones(np.array([xf,yf,zf]).T.shape[0])*wt)
pdb.set_trace()
np.save("Photons-per-MeV-histxyz",pperMeV[0])
np.save("Photons-per-MeV-binsxyz",pperMeV[1])
exit()


binsz = 1250*qe/20.

cts, be, er = skh_plt.hist(sumgvv,weights=weightsv,bins=np.arange(0,Nphot*qe,binsz),errorbars=False, histtype='step',label=labelv)#,stacked='true')
hdict = dict()
hdict["counts"]=cts
hdict["binedges"]=be
hdict["err"]=er

fileout = "pbomb_ngamma"
plt.legend()
plt.title(fileout+ ' 40,60,80 fpcc ' + 'w/wo cathode SiPMs, 1250 photons/100keV')
plt.xlabel('Counted photons - no q.e.')
plt.ylabel('n gammas - normalized')
#plt.yscale('log')
plt.savefig(fileout+'.png')



if  len(sumgvvf):

    cts, be, er = skh_plt.hist(sumgvvf,weights=weightsvf,bins=np.arange(0,Nphot*qe,binsz),errorbars=False, histtype='step',label=labelv)#,stacked='false')
    hdict = dict()
    hdict["counts"]=cts
    hdict["binedges"]=be
    hdict["err"]=er

    fileout = "pbomb_ngamma_fid"
    plt.legend()
    plt.title(fileout+ ' 40,60,80 fpcc ' + 'w/wo cathode SiPMs, 1250 photons/100keV from fidv')
    plt.xlabel('Counted photons - no q.e.')
    plt.ylabel('n gammas - normalized')
#plt.yscale('log')
    plt.savefig(fileout+'.png')

else:
    print ("Not creating pbomb_ngamma_fid, since no events satisfied fidV cut.")


npa = np.array(sumgvv)
meanp = np.mean(npa,axis=1)
stdp = np.std(npa,axis=1)

meanpf =  []
stdpf = []
for ii in len(sumgvvf):
    meanpf.append( np.array(sumgvvf[ii])[np.array(sumgvvf[ii])!=0].mean())
    stdpf.append( np.array(sumgvvf[ii])[np.array(sumgvvf[ii])!=0].std())
    
print("Means/StdDevs for all pbombs inside acrylic: " + str(meanp) + " -- " + str(stdp))
print("Means/StdDevs for all pbombs inside fidV: " + str(meanpf) + " -- " + str(stdpf))

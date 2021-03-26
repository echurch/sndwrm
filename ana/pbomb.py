from ROOT import TFile,TH1F
import numpy as np
from matplotlib import pyplot as plt
import uproot
from skhep.visual import MplPlotter as skh_plt

import pdb



files = ["../build/pbomb_optphys_40.root","../build/pbomb_optphys_60.root","../build/pbomb_optphys_80.root","../build/pbomb_optphys_40_nopc.root","../build/pbomb_optphys_60_nopc.root","../build/pbomb_optphys_80_nopc.root"]
qe = 1.

sumgvv = []
weightsv = []
weightsvf = []
sumgvvf = []
labelv = []
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

    for step in range(Nsteps):
#    for step in range(400000):

        # below are  the Primaries (SProcess==null).  Also, insist that we're onto a new event.

        if  not fidV and b'null' in steps['SProcess'][step] :
            fidV = abs(steps['X'][step]) < 2200. and ( ( steps['Y'][step] > 1000. and steps['Y'][step] < 3000.) or (steps['Y'][step] < -1000. and steps['Y'][step] > -3000. ) ) and abs(steps['Z'][step]) < 27000.  
            if fidV:
                cntf += 1
#                print ("Event/trkID/stepNum: " + str(steps['Event'][step]) +"/" + str(steps['TrkID'][step]) + "/" + str(steps['StepNum'][step]))
## Need to loop over all steps in an event and get the sumg for each


        if b'SiPM' in steps['TVolume'][step]: 
            sumg +=  1

        if b'SiPM' in steps['TVolume'][step]  and fidV: 
            sumgf +=  1

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



if  range(len(sumgvv)):

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
for ii in range(len(sumgvvf)):
    meanpf.append( np.array(sumgvvf[ii])[np.array(sumgvvf[ii])!=0].mean())
    stdpf.append( np.array(sumgvvf[ii])[np.array(sumgvvf[ii])!=0].std())
    
print("Means/StdDevs for all pbombs inside acrylic: " + str(meanp) + " -- " + str(stdp))
print("Means/StdDevs for all pbombs inside fidV: " + str(meanpf) + " -- " + str(stdpf))

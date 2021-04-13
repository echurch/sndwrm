from ROOT import TFile,TH1F
import numpy as np
from matplotlib import pyplot as plt
import uproot
from skhep.visual import MplPlotter as skh_plt

import pdb



#files = ["../build/marley_optphys_cno.root"]
#files = ["../build/marley_optphys_8B.root"]
#files = ["../build/marley_optphys_SN-liv.root"]
files = ["../build/marley_optphys_SN-liv-CEvNS.root"]
qe = 1.
birk = 0.7
sumgvv = []
weightsv = []
weightsvf = []
weightsvftype = []
sumgvvf = []
sumgvvftype = []
labelv = []
labelvtype = []

for file in files:

    print ("Reading " + str(file))

    f = uproot.open(file)
    steps = f['Steps']
    steps = steps.arrays(namedecode='utf-8')
    Nsteps = len(steps['Time'])
    trks = f['Tracks']
    trks = trks.arrays(namedecode='utf-8')
    maxE = 0.010 # len(trks['Event']) ~ 10MeV (3 MeV) [12 MeV] {100 keV} max scale for 8B (CNO) [SN-liv] {CEvNS}
    phpE = 12500 # 24000 (12500) photons/MeV ionizing e's (n.r.)
    Nphot = phpE*maxE
    Nevts = trks['Event'].max()+1

    ES = np.zeros(int(Nevts))  # nue-e scattering false by default


    for trk in  range(len(trks['PID'])):
        if trks['PID'][trk]==12:
            ES[int(trks['Event'][trk])] = 1


    eventO = -12
    sumg = 0
    sumgv = []
    sumgf = 0
    sumgfes = 0
    sumgfcc = 0
    sumgvf = []
    sumgvfes = []
    sumgvfcc = []
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

        if b'SiPM' in steps['TVolume'][step] and steps['PID'][step]==0: 
            sumg +=  1

        if b'SiPM' in steps['TVolume'][step]  and fidV and steps['PID'][step]==0: 
            sumgf +=  1

        if ES[int(steps['Event'][step])] and b'SiPM' in steps['TVolume'][step] and steps['PID'][step]==0:
            sumgfes +=1/phpE/birk
        elif (not ES[int(steps['Event'][step])]) and b'SiPM' in steps['TVolume'][step] and steps['PID'][step]==0:
            sumgfcc +=1/phpE/birk

        if ((eventO != steps['Event'][step] and steps['Event'][step]>0)  or  step>=(Nsteps-1) ):
            newEvt = 1
            sumgv.append(sumg)
            sumg = 0

            sumgvf.append(sumgf)

            if ES[int(steps['Event'][step-1])]:
                sumgvfes.append(sumgfes)
            else:
                sumgvfcc.append(sumgfcc)
            sumgf = 0
            sumgfes = 0
            sumgfcc = 0
            fidV = 0
            if not eventO%10:
                print("Onto event: " + str(steps['Event'][step]))

            
        eventO = steps['Event'][step] 
        newEvt = 0

    weights =  np.repeat(1.0/Nevts,len(sumgv))
    weightsv.append(weights)
    weightsf =  np.repeat(1.0/Nevts,len(sumgvf))
    weightsvf.append(weightsf)
    sumgvv.append(sumgv)
    sumgvvf.append(sumgvf)
    sumgvvftype.append(sumgvfes)
    sumgvvftype.append(sumgvfcc)

    weightsvftype.append(np.repeat(1.0/Nevts,len(sumgvfes)))
    weightsvftype.append(np.repeat(1.0/Nevts,len(sumgvfcc)))

    sumgv = []
    sumgvf = []
    sumgvfes = []
    sumgvfcc = []
    labelv.append("-".join(file.split(".root")[0]))

    labelvtype.append(file.split(".root")[0]+"_ES")
    labelvtype.append(file.split(".root")[0]+"_CC")
    print ("total in fv: " + str(cntf))
    cntf

binsz = Nphot*qe/20.

cts, be, er = skh_plt.hist(sumgvv,weights=weightsv,bins=np.arange(0,Nphot*qe,binsz),errorbars=False, histtype='step',label=labelv)#,stacked='true')
hdict = dict()
hdict["counts"]=cts
hdict["binedges"]=be
hdict["err"]=er

fileout = "marley_ngamma"
plt.legend()
plt.title(fileout + '8B w/wo cathode SiPMs, Marley photons')
plt.xlabel('Counted photons - no q.e.')
plt.ylabel('n gammas - normalized')
#plt.yscale('log')
plt.savefig(fileout+'.png')
plt.close()



if  range(len(sumgvv)):

    Nphot /= phpE
    binsz /= phpE

    cts, be, er = skh_plt.hist(sumgvvftype,weights=weightsvftype,bins=np.arange(0,Nphot*qe,binsz),errorbars=False, histtype='step',label=labelvtype,stacked='true')
    hdict = dict()
    hdict["counts"]=cts
    hdict["binedges"]=be
    hdict["err"]=er

#    fileout = "marley_nue_solar_8B"
#    fileout = "marley_nue_solar_CNO"
    fileout = "marley_nue_SN_liv_COH"
    plt.legend()
    plt.title(fileout+ ' ES/CC evts ' + 'w cathode SiPMs, MeV/'+str(phpE)+' photons')
    plt.xlabel('Summed Energy [MeV] - no q.e.,birks=0.7')
    plt.ylabel('Events - normalized')
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

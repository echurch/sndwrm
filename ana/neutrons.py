from ROOT import TFile,TH1F
import numpy as np
from matplotlib import pyplot as plt
import uproot
from skhep.visual import MplPlotter as skh_plt
import glob
import pdb



#files = ["../build/marley_optphys_cno.root"]
#files = ["../build/marley_optphys_8B.root"]
#files = ["../build/marley_optphys_SN-liv.root"]
#files = ["../build/marley_optphys_SN-liv-CEvNS.root"]
files = glob.glob("../build/neutron_coldcryoskin_*.root")
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

    maxE = 12. # len(trks['Event']) ~ 10MeV (3 MeV) [12 MeV] {100 keV} max scale for 8B (CNO) [SN-liv] {CEvNS}
    phpE = 24000 # 24000 (12500) photons/MeV ionizing e's (n.r.)
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
#        if  not fidV and b'null' in steps['SProcess'][step] :
#            fidV = abs(steps['X'][step]) < 2200. and ( ( steps['Y'][step] > 1000. and steps['Y'][step] < 3000.) or (steps['Y'][step] < -1000. and steps['Y'][step] > -3000. ) ) and abs(steps['Z'][step]) < 27000.  

        if  not fidV and  steps['PID'][step]==22 : ## take presence of gamma to mean n capture or robust inelastic thingy happened.
            fidV = abs(steps['X'][step]) < 3000. and abs( steps['Y'][step]) < 6000 and abs(steps['Z'][step]) < 20000.  
            if fidV:
                cntf += 1
#                print ("Event/trkID/stepNum: " + str(steps['Event'][step]) +"/" + str(steps['TrkID'][step]) + "/" + str(steps['StepNum'][step]))
## Need to loop over all steps in an event and get the sumg for each

        if b'SiPM' in steps['TVolume'][step] and steps['PID'][step]==0: 
            sumg +=  1

        if b'SiPM' in steps['TVolume'][step]  and fidV and steps['PID'][step]==0: 
            sumgf +=  1/phpE/birk

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
    labelv.append(file.split(".root")[0]+"_tot")

    labelvtype.append(file.split(".root")[0]+"_ES")
    labelvtype.append(file.split(".root")[0]+"_CC")
    print ("total in fv: " + str(cntf))
    cntf

binsz = Nphot/phpE*qe/20.

pdb.set_trace()
cts, be, er = skh_plt.hist([item for sublist in sumgvvf for item in sublist],weights=[item for sublist in weightsvf for item in sublist],bins=np.arange(0,Nphot/phpE*qe,binsz),errorbars=False, histtype='step',label=labelv)#,stacked='true')
hdict = dict()
hdict["counts"]=cts
hdict["binedges"]=be
hdict["err"]=er

fileout = "neutron_capture"
plt.legend()
plt.title(fileout + 'summed SiPM energy from coldcryoskin photons - capture inside acrylic volume')
plt.xlabel('Counted photons - no q.e.')
plt.ylabel('n gammas - normalized')
#plt.yscale('log')
plt.savefig(fileout+'.png')
plt.close()


npa = np.array(sumgvv)

from ROOT import TFile
from matplotlib import pyplot as plt
import numpy as np
import math
import pdb
import math

from skhep.visual import MplPlotter as skh_plt



#thresh = 0.020 # MeV 
#thresh = 0.050 # MeV 
sep = 20 # mm
thresh = 0.075 # MeV 
#sep = 50 # mm
sepfed = 1000 # mm

fname="../build/tl208_G10.root"
f = TFile(fname)
#f = TFile("../build/k40_liquid.root")
fout = fname.split("/")[-1].split("_")[0]


def fdist (x,y,z,tke):
    d = []
    for ii in range(len(x)-1):  # have checked that vector is at least 2 elements long
        for jj in range(ii+1,len(x)-1,1):  # have checked that vector is at least 2 elements long
            if  tke[ii]>thresh and tke[jj]>thresh:
                d.append( math.sqrt((x[ii] - x[jj])**2+(y[ii] - y[jj])**2+(z[ii] - z[jj])**2) )

    return d


# add all n.r. vertices from this event now that we have one above threshold inside the fid volume.
def allvertices (x,y,z,t,ID,trkin):
    f.Tracks.GetEntry(trkin)
    evt = f.Tracks.Event

    # Scan all trks near this one, requiring that they're in the event
    for trk in range(-100,100,1):
        if trkin+trk<0:
            continue
        f.Tracks.GetEntry(trkin+trk)
        if f.Tracks.Event<evt:
            continue
        if f.Tracks.Event>evt:
            break


        # below appends all e- vtxs/trks over threshold and within sepfed [mm] of the edge of our 6x6x20 m^3 box or just inside the box.
        if (f.Tracks.KEnergy > thresh) and (f.Tracks.PID==11 ) and ( abs(f.Tracks.Startx)<(3000+sepfed) and abs(f.Tracks.Startz)<(10000+sepfed)  ):  # y is drift direction
            x.append(f.Tracks.Startx)
            y.append(f.Tracks.Starty)
            z.append(f.Tracks.Startz)
            t.append(f.Tracks.KEnergy)

    f.Tracks.GetEntry(trkin)
    return x,y,z,t




Nent = f.Tracks.GetEntries()
event = 0
#Nent = 50000


trkke0 = []
nke0 = []
nke0_sep = []
nke0_bare = []
nke0_nonbare = []

trkke = []
trkid = []
xke = []
yke = []
zke = []
dists = []
neutroncap = False
cnttrks = 0
trkloop = False


for trk in range(Nent):
    fidv = True
    f.Tracks.GetEntry(trk)
    dist = []
    if not f.Tracks.Event == event:  # We are onto a new event, so let's do our checks, then reinitialize for next evt.
        event = f.Tracks.Event

        if not f.Tracks.Event%100000:
            print ("Event " + str(f.Tracks.Event))

        nke0.append(f.Tracks.KEnergy)

    
        if len(trkke)>1:
            dist = fdist(xke,yke,zke,trkke)
#            pdb.set_trace()
            dists.append(dist) # dist is a list of distances between vtxes for just finished event

        trkke = []
        xke = []
        yke = []
        zke = []
        cnttrks = 0
        trkloop = False

    fidv = abs(f.Tracks.Startx) < 3000 and abs(f.Tracks.Starty) < 3000 and abs(f.Tracks.Startz) < 10000

    #  Every track every event
    if (f.Tracks.KEnergy > thresh) and (f.Tracks.PID==11 ) and fidv and not trkloop:
        '''
        trkke.append(f.Tracks.KEnergy)
        xke.append(f.Tracks.Startx)
        yke.append(f.Tracks.Startx)
        zke.append(f.Tracks.Startx)
        '''
        xke,yke,zke,trkke = allvertices(xke,yke,zke,trkke,f.Tracks.TrkID,trk)
        # Now that we've found 1 vtx inside fidvolume, we've looped over all candidate vtxes in this evt and added them if they pass.
        # Set trkloop=True so we  don't do this again for this evt. Will later see if any of 'em are close enough together.
        trkloop = True


    if (f.Tracks.PID==11 ) and fidv:
        trkke0.append(f.Tracks.KEnergy)
        
    


fig = plt.figure()
H,__,__ = skh_plt.hist(trkke0,bins=np.arange(0,0.2,0.01), errorbars=True, histtype='step')
plt.title('Electron KE')
plt.xlabel('MeV')
plt.ylabel('Entriess/10keV')
plt.savefig('KE_'+str(thresh)+'_e_'+fout+'.png')



fig = plt.figure()
flatdists = [y for x in dists for y in x]
H2,__,__ = skh_plt.hist(flatdists,bins=np.arange(0.,100.,4.), errorbars=True, histtype='step')
plt.title('Distances between vertices w  KE recoil> ' + str(thresh) + ' MeV in an event')
plt.xlabel('mm')
plt.ylabel('Entries/mm')
plt.savefig('dists_'+str(thresh)+'_'+fout+'_e_.png')






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
thresh = 0.010 # MeV 
#thresh = 0.100 # MeV 

#sep = 50 # mm
sepfed = 1000 # mm

#f = TFile("../build/neutron_shield.root")
#f = TFile("../build/neutron.root")

#f = TFile("../build/neutron_plastic20cm.root")
##f = TFile("../build/neutron_outsideshield_dunedistn_liquid.root")
##f = TFile("../build/neutron_outsideshield25cm_dunedistn_liquid.root")

fname = "../build/neutron_outsidefoamwoodss_dunedistn_liquid.root"
normtoktyr = 3.064 # See macros/neutron_fw.mac with 0-water shield assumption. With 40 cm H20 I get #1000 reduction
normtoktyr = 3.064/80.*2  # With 20cm water shield I get ~80. The extra 2 is because I generated n's isotropically, whereas they come out of walls at 1E-5.

####fname = "../build/neutron_outsidefoamwoodss_cs_10cm_dunedistn_liquid.root"
##fname = "../build/neutron_outsidefoamwoodss_cs_5cm_dunedistn_22x8acryl_liquid.root"
##normtoktyr = 238854/15500000. *10. # For cryoskin. See the macros/neutrons_fw_cs.mac ### MY *10 to get it to 1E-9 n/cm3/sec 4-Feb-2020

f = TFile(fname)

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
    ## get the x,y,z coords of the primary neutron
    xyzp = None

    # Scan all trks near this one, requiring that they're in the event
    for trk in range(-100,100,1):
        if trkin+trk<0:
            continue
        f.Tracks.GetEntry(trkin+trk)
        if f.Tracks.Event<evt:
            continue
        if f.Tracks.Event>evt:
            break

        if f.Tracks.PID==2112 and f.Tracks.TrkID==1:
            xyzp = (f.Tracks.Startx,f.Tracks.Starty,f.Tracks.Startz)

        '''
        if f.Tracks.trkID==trkin:
            continue
        '''

        # below appends all n.r. vtxs/trks over threshold and within sepfed [mm] of the edge of our 6x6x20 m^3 box or just inside the box.
        if (f.Tracks.KEnergy > thresh) and (f.Tracks.PID>100E6) and ( abs(f.Tracks.Startx)<(3000+sepfed) and abs(f.Tracks.Startz)<(10000+sepfed)  ):  # y is drift direction
            x.append(f.Tracks.Startx)
            y.append(f.Tracks.Starty)
            z.append(f.Tracks.Startz)
            t.append(f.Tracks.KEnergy)


    f.Tracks.GetEntry(trkin)
    return x,y,z,t,xyzp




Nent = f.Tracks.GetEntries()
Nent = 1000000

event = 0


xyzprimeevt = np.empty((0,3),dtype=float)
trkke0 = []
nke0 = []
nke0_sep = []
nke0_sep_gam = []
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

        if f.Tracks.TrkID==1 and len(trkke)==0:
            nke0_bare.append(KE_n_primary)

        if f.Tracks.TrkID==1 and len(trkke)>0:
            nke0_nonbare.append(KE_n_primary)
    
        if len(trkke)>1:
            dist = fdist(xke,yke,zke,trkke)
#            pdb.set_trace()
            dists.append(dist) # dist is a list of distances between vtxes for just finished event
            oldlen = len(nke0_sep_gam)
            if len(dist):
                if any(d>sep for d in dist):
                    nke0_sep.append(KE_n_primary)
                    nke0_sep_gam.append(KE_n_primary)
                    oldlen = -12

            # here we keep trk of events where even though there aren't 2 separate vtxs, there is an evident neutroncapture
            if neutroncap and len(nke0_sep_gam)==oldlen: 
                nke0_sep_gam.append(f.Tracks.KEnergy)

        trkke = []
        xke = []
        yke = []
        zke = []
        neutroncap = False
        cnttrks = 0
        trkloop = False


    fidv = abs(f.Tracks.Startx) < 3000 and abs(f.Tracks.Starty) < 3000 and abs(f.Tracks.Startz) < 10000

    #  Every track every event
    if (f.Tracks.KEnergy > thresh) and (f.Tracks.PID>100E6) and fidv and not trkloop:
        '''
        trkke.append(f.Tracks.KEnergy)
        xke.append(f.Tracks.Startx)
        yke.append(f.Tracks.Startx)
        zke.append(f.Tracks.Startx)
        '''
        xke,yke,zke,trkke,xyzprime = allvertices(xke,yke,zke,trkke,f.Tracks.TrkID,trk)
        # Now that we've found 1 vtx inside fidvolume, we've looped over all candidate vtxes in this evt and added them if they pass.
        # Set trkloop=True so we  don't do this again for this evt. Will later see if any of 'em are close enough together.
        trkloop = True
        if xyzprime is not None:
            xyzprimeevt = np.vstack((xyzprimeevt,np.array(xyzprime)))


## 17-Feb-2019, EC. New indent here. Force gamma captures to be identified only if we also have an n.r. in fidv overthresh.
## The problem with identifying ncaptures generally, as I had it, is I then end up subtracting events from the 1vtx-over-thresh sample.
        if f.Tracks.PID==22 and f.Tracks.KEnergy>1.0: # MeV
            neutroncap = True
#           print ("High energy gamma.")

    if f.Tracks.PID>100E6 and fidv:
        trkke0.append(f.Tracks.KEnergy)
    if f.Tracks.PID==2112 and f.Tracks.TrkID==1:
        KE_n_primary = f.Tracks.KEnergy
        
    if  (f.Tracks.PID==2112 and f.Tracks.TrkID==1) or (f.Tracks.PID>100E6 and f.Tracks.KEnergy>thresh and fidv):
        cnttrks+=1
    
fout = fname.split("/")[-1].split(".")[0]

fig = plt.figure()
H,__,__ = skh_plt.hist(trkke0,bins=np.arange(0,0.2,0.01), errorbars=True, histtype='step')
plt.title('Nuclear recoil KE')
plt.xlabel('MeV')
plt.ylabel('Entries/10keV')
plt.yscale('log')
plt.savefig('NucleusKE_'+str(thresh)+fout+'.png')

fig = plt.figure()
H2,__,__ = skh_plt.hist(nke0,bins=np.arange(0.5,8.0,0.05), errorbars=True, histtype='step',label='Original n')
plt.title('Initial Neutron KE')
plt.xlabel('MeV')
plt.ylabel('Entries/10keV')
plt.savefig('NeutronKE0_'+str(thresh)+'MeV-'+str(sep)+'mm_dunedistn.png')


fig = plt.figure()
flatdists = [y for x in dists for y in x]
H2,__,__ = skh_plt.hist(flatdists,bins=np.arange(0.,100.,4.), errorbars=True, histtype='step')
plt.title('Distances between vertices w  KE recoil> ' + str(thresh) + ' MeV in an event')
plt.xlabel('mm')
plt.ylabel('Entries/mm')
plt.savefig('dists_'+str(thresh)+fout+'.png')

'''
fig = plt.figure()
H2,__,__ = skh_plt.hist(nke0,bins=np.arange(0.5,3.0,0.05), errorbars=True, histtype='step',label='Primary neutron spectrum used')
H0,__,__ = skh_plt.hist(nke0_bare,bins=np.arange(0.5,3.0,0.05), errorbars=True, histtype='step',label='Primary neutron and sub-thresh scatt-nuclei only')

plt.title('Initial Neutron and non-interacting spectra')
plt.xlabel('MeV')
plt.ylabel('Entries/10keV')
plt.legend()
plt.savefig('NeutronKE_interact'+str(thresh)+'_'+str(sep)+fout+'.png')

'''
fig = plt.figure()

H0,bins,__ = skh_plt.hist(nke0,bins=np.arange(0.5,8.0,0.05), errorbars=True, histtype='step',label='Primary neutron  only',Fill=True,weights=np.repeat(normtoktyr,len(nke0)))
H4,bins,__ = skh_plt.hist(nke0_bare,bins=np.arange(0.5,8.0,0.05), errorbars=True, histtype='step',label='primary neutron below thresh',stacked=True,Fill=True,weights=np.repeat(normtoktyr,len(nke0_bare)))

## Now make the histos for point of origin of primary neutron for evts in which at least one n.r. over thresh happened inside fidv.
fig = plt.figure()
##fig, ax = plt.subplots(2, 2)
labels = []
##fig.add_subplot()

eps = 1.
xymax = 6000.
zmax = 31000.

pdb.set_trace()
np.save('neutron_xyz.npy',xyzprimeevt)
xyzprimeevt[:,0] = np.clip(xyzprimeevt[:,0],-xymax+eps,xymax-eps)
xyzprimeevt[:,1] = np.clip(xyzprimeevt[:,1],-xymax+eps,xymax-eps)
xyzprimeevt[:,2] = np.clip(xyzprimeevt[:,2],-zmax+eps,zmax-eps)


### Two things: Keep the xymax,zmax away from below bin edges, and (2) note that the bizarre indexing for imshow, ala
### https://stackoverflow.com/questions/11367683/imshow-and-histogram2d-cant-get-them-to-work

h2a = np.histogram2d(xyzprimeevt[:,0],xyzprimeevt[:,1],bins=(np.arange(-6200.,+6200.,100.),np.arange(-6200.,+6200.,100.)),weights=np.repeat(normtoktyr,len(xyzprimeevt[:,0])))
plt.imshow(h2a[0][::-1,:],aspect='auto',extent=[-6200,6200,-6200,6200])
#plt.pcolor(np.arange(-6200.,+6200.,100.),np.arange(-6200.,+6200.,100.),h2a[0])
#plt.xlim(-6200.,6200.)
#plt.ylim(-6200.,6200.)
plt.title('Primary neutron xy')
plt.savefig('NeutronKE_xy_'+str(thresh)+'_'+str(sep)+fout+'.png')


fig = plt.figure()
h2b = np.histogram2d(xyzprimeevt[:,0],xyzprimeevt[:,2],bins=(np.arange(-6200.,+6200.,100.),np.arange(-32000.,+32000.,100.)), weights=np.repeat(normtoktyr,len(xyzprimeevt[:,0])))
plt.imshow(h2b[0][::-1,:],aspect='auto',extent=[-32000,32000,-6200,6200])
#plt.pcolor(np.arange(-6200.,+6200.,100.),np.arange(-32000.,+32000.,100.),h2b[0])
#plt.xlim(-6200.,6200.)
#plt.ylim(-32000.,32000.)
plt.title('Primary neutron xz')
plt.savefig('NeutronKE_xz_'+str(thresh)+'_'+str(sep)+fout+'.png')

fig = plt.figure()
h2c = np.histogram2d(xyzprimeevt[:,1],xyzprimeevt[:,2],bins=(np.arange(-6200.,+6200.,100.),np.arange(-32000.,+32000.,100.)), weights=np.repeat(normtoktyr,len(xyzprimeevt[:,0])))
plt.imshow(h2c[0][::-1,:],aspect='auto',extent=[-32000,32000,-6200,6200])
#plt.pcolor(np.arange(-6200.,+6200.,100.),np.arange(-32000.,+32000.,100.),h2c[0])
#plt.xlim(-6200.,6200.)
#plt.ylim(-32000.,32000.)
plt.title('Primary neutron yz')
#plt.legend([h2a,h2b,h2c],labels)
plt.savefig('NeutronKE_yz_'+str(thresh)+'_'+str(sep)+fout+'.png')



fig = plt.figure()
## Try to stack from smallest to biggest, else the tiny contributions on top are not seen. EC, 5-Feb-2020.
labels = list(['multisite-cut: '+str(round(len(nke0_sep_gam)*normtoktyr,3))+". Net: "+str(round((len(nke0)-len(nke0_bare)-len(nke0_sep_gam))*normtoktyr,3)),'detected (above thresh): '+str(round((len(nke0)-len(nke0_bare))*normtoktyr,3))])
H2,bins,__ = skh_plt.hist(nke0,bins=np.arange(0.5,8.0,0.05), histtype='step',label='Primary neutron spectrum: '+str(round(normtoktyr*len(nke0)/1.E6,2))+"*E6", weights=np.repeat(normtoktyr,len(nke0)),errorbars=True)

#plt.errorbar((bins[0:-1]+bins[1:])/2.,H0-H4,yerr=np.sqrt(H0**2+H4**2),label=labels[1])
plt.plot((bins[0:-1]+bins[1:])/2.,H0-H4,label=labels[1],marker="_",linestyle='None')
H5,__,__ = skh_plt.hist(nke0_sep_gam,bins=bins,Fill=True,label=labels[0],weights=np.repeat(normtoktyr,len(nke0_sep_gam)),errorbars=True)   #### stacked=True,

plt.suptitle('Initial Neutron and detectable spectra')
plt.title('theshold: '+str(thresh)+' MeV, Pos resolution: '+str(sep)+' mm')
plt.xlabel('MeV')
plt.ylabel('Entries/50keV/kt-yr')
plt.ylim(1.E-3, 1.E+5)
plt.yscale('log')
plt.legend()
plt.savefig('NeutronKE_interact'+str(thresh)+'_'+str(sep)+fout+'.png')
np.save('neutron_spec_orig_H2.npy',H2)
np.save('neutron_spec_tagged_H5.npy',H5)
np.save('neutron_spec_orig-thresh_H2-H4.npy',H0-H4)

fig = plt.figure()
H0,__,__ = skh_plt.hist(np.repeat(nke0_nonbare,normtoktyr),bins=np.arange(0.5,8.0,0.05), errorbars=True, histtype='step',label='Primary neutron interacts in the volume')
H3,__,__ = skh_plt.hist(np.repeat(nke0_sep,normtoktyr),bins=np.arange(0.5,8.0,0.05), errorbars=True, histtype='step',label='2 vtxes '+str(thresh)+' MeV and '+str(sep) +' mm separation')
H4,__,__ = skh_plt.hist(np.repeat(nke0_sep_gam,normtoktyr),bins=np.arange(0.5,8.0,0.05), errorbars=True, histtype='step',label='2 vtxes '+str(thresh)+' MeV and '+str(sep) +' mm separation OR a 1+ MeV gamma')

plt.title('Interacting Neutron KE -- and post cuts')
plt.xlabel('MeV')
plt.ylabel('Entries/50keV/kt-yr')
plt.legend()
plt.savefig('NeutronKE_cuts_'+str(thresh)+'_'+str(sep)+fout+'.png')


fig = plt.figure()
H2,bins = np.histogram(np.repeat(nke0,normtoktyr),bins=np.arange(0.5,8.0,0.05))
yerr=abs(H4-H0)/H2*np.sqrt(1/H4+1/H0+1/H2)
yerr[np.isinf(abs(yerr))] = 0.0
plt.errorbar(bins[:-1]+np.diff(bins)/2,(H0-H4)/H2,yerr)
plt.title('Interacting n KE post cuts as fraction of original distribution')
plt.xlabel('MeV')
plt.ylabel('Entries/50keV/kt-yr')

plt.savefig('NeutronKE0_cutfrac_'+str(thresh)+'MeV_'+str(sep)+'mm_'+fout+'.png')


fig = plt.figure()
H2,bins = np.histogram(np.repeat(nke0,normtoktyr),bins=np.arange(0.5,8.0,0.05))
yerr=abs(H4-H0)/H0*np.sqrt(1/H4+1/H0)
yerr[np.isinf(abs(yerr))] = 0.0
plt.errorbar(bins[:-1]+np.diff(bins)/2,(H0-H4)/H2,yerr)
plt.title('Interacting n KE post cuts as fraction of visible n distribution')
plt.xlabel('MeV')
plt.ylabel('Entries/50keV/kt-yr')

plt.savefig('NeutronKE_cutfrac_'+str(thresh)+'MeV_'+str(sep)+'mm_'+fout+'.png')


from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import pyplot as plt
import pdb

fig = plt.figure()
ax = fig.gca(projection='3d')
#ax.set_aspect("equal")
ax.set_aspect("auto")

# draw cube
def rect_prismz(x_range, y_range, z_range,col):
    # TODO: refactor this to use an iterator
    xx, yy = np.meshgrid(x_range, y_range)
    ax.plot_wireframe(xx, yy, np.ones(xx.shape)*z_range[0], color="r", alpha=0.1)
    ax.plot_surface(xx, yy, np.ones(xx.shape)*z_range[0], color=col, alpha=0.1)
    ax.plot_wireframe(xx, yy, np.ones(xx.shape)*z_range[1], color="r", alpha=0.1)
    ax.plot_surface(xx, yy, np.ones(xx.shape)*z_range[1], color=col, alpha=0.1)

def rect_prismx(x_range, y_range, z_range,col):
    yy, zz = np.meshgrid(y_range, z_range)
    ax.plot_wireframe(np.ones(yy.shape)*x_range[0], yy, zz, color="r", alpha=0.1)
    ax.plot_surface(np.ones(yy.shape)*x_range[0], yy, zz, color=col, alpha=0.1)
    ax.plot_wireframe(np.ones(yy.shape)*x_range[1], yy, zz, color="r", alpha=0.1)
    ax.plot_surface(np.ones(yy.shape)*x_range[1], yy, zz, color=col, alpha=0.1)

def rect_prismy(x_range, y_range, z_range,col):
    xx, zz = np.meshgrid(x_range, z_range)
    ax.plot_wireframe(xx, np.ones(xx.shape)*y_range[0], zz, color="r", alpha=0.1)
    ax.plot_surface(xx, np.ones(xx.shape)*y_range[0], zz, color=col, alpha=0.1)
    ax.plot_wireframe(xx, np.ones(xx.shape)*y_range[1], zz, color="r", alpha=0.1)
    ax.plot_surface(xx, np.ones(xx.shape)*y_range[1], zz, color=col, alpha=0.1)

arr = np.load("/Users/chur558/geant4.10.05/examples/extended/radioactivedecay/rdecay02/ana/Photons-per-MeV-3MeV-e-shinyg10-histxyz.npy")
arr = np.load("/Users/chur558/geant4.10.05/examples/extended/radioactivedecay/rdecay02/ana/Photons-per-MeV-3MeV-e-shinyg10-histxyz.npy")

rect_prismz(np.linspace(-6,6,100), np.linspace(-6,6,100), np.array([-28, 28]),'r')
rect_prismx(np.array([-6, 6]), np.linspace(-6,6,100), np.linspace(-28,28,100),'r')
rect_prismy(np.linspace(-6,6,100), np.array([-6,6]), np.linspace(-28,28,100),'r')
rect_prismz(np.linspace(-3,3,100), np.linspace(-6,6,100), np.array([-20, 20]),'b')
rect_prismx(np.array([-3, 3]), np.linspace(-6,6,100), np.linspace(-20,20,100),'b')
rect_prismy(np.linspace(-3,3,100), np.array([-6,6]), np.linspace(-20,20,100),'b')
plt.ylim(-10,10)
plt.xlim(-10,10)



zmin, zmax = ax.get_zlim()
ax.text(0, 0, 40., "Beam axis")
ax.text(0, 20, 0., "Up")

pdb.set_trace()
ax.text(2.5,1,0,str(round(arr[2,1,1],0)) + " Photons/MeV")
#ax.text(2.5,12,0,str(round(arr[2,2,1],0)) )
ax.text(1.,11,0,str(round(arr[1,2,1],0)) )
#ax.text(2.5,1,19,str(round(arr[2,0,2],0)) )
ax.text(0.,1,18,str(round(arr[1,1,2],0)) )
ax.scatter(2.5,1,0)
#ax.scatter(2.5,12,0)
ax.scatter(1.,11,0)
#ax.scatter(2.5,1,19)
ax.scatter(0.,1,18)


ax.view_init(20, 45)



plt.savefig("PhotonCountMeV-e3MeV-shinyg10.png")


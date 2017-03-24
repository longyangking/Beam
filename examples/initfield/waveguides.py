import sys
sys.path.append("../..")

import Beam
import numpy as np
import matplotlib.pyplot as plt

Nx = 400
Ny = 400

def dn(z):
    x = np.linspace(-Nx/2,Nx/2,Nx)
    y = np.linspace(-Ny/2,Ny/2,Ny)
    x,y = np.meshgrid(x,y)

    R = Nx/4
    delta = Nx/13

    dn =  np.exp(-((x - R)**2 + (y)**2)/(delta**2)) \
        + np.exp(-((x + R)**2 + (y)**2)/(delta**2)) \
        + np.exp(-((x - R*np.cos(np.pi/3))**2 + (y - R*np.sin(np.pi/3))**2)/(delta**2)) \
        + np.exp(-((x + R*np.cos(np.pi/3))**2 + (y - R*np.sin(np.pi/3))**2)/(delta**2)) \
        + np.exp(-((x - R*np.cos(np.pi/3))**2 + (y + R*np.sin(np.pi/3))**2)/(delta**2)) \
        + np.exp(-((x + R*np.cos(np.pi/3))**2 + (y + R*np.sin(np.pi/3))**2)/(delta**2)) 
    return 7.0*10**-4*dn

if __name__ == '__main__':
    dx = 1.0*10**-6
    dy = 1.0*10**-6
    Lx = dx*Nx
    Ly = dy*Ny
    delta = 100.0*10**-6
    n0 = 1.45
    k0 = 2.0*np.pi/(633.0*10**-9)

    Z = 30.0*10**-3

    #plt.imshow(dn(0))
    #plt.show()
    
    x = np.linspace(0,Lx,Nx)
    y = np.linspace(0,Ly,Ny)
    x,y = np.meshgrid(x,y)
    
    initfield = np.exp(-((x-Lx/2)**2+(y-Ly/2)**2)/(delta**2))

    bpm = Beam.BPM(n0=n0,k0=k0,dx=dx,dy=dy,Nx=Nx,Ny=Ny,Z=Z,dn=dn,initfield=initfield,verbose=False)
    bpm.solve() 
    viewer = Beam.Viewer(bpm)
    viewer.viewintensity()
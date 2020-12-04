import sys
sys.path.append("../..")

import Beam
import numpy as np

def dn(z):
    Nx = 100
    Ny = 100
    return np.zeros((Nx,Ny))

if __name__ == '__main__':
    dx = 1.0*10**-6
    dy = 1.0*10**-6
    Nx = 100
    Ny = 100
    Lx = dx*Nx
    Ly = dy*Ny
    delta = 15.0*10**-6
    n0 = 1.45
    k0 = 2.0*np.pi/(633.0*10**-9)

    Z = 1.0*10**-3
    
    x = np.linspace(0,Lx,Nx)
    y = np.linspace(0,Ly,Ny)
    x,y = np.meshgrid(x,y)
    
    initfield = 1.0*(x>Lx/4)*(x<Lx*3/4)*(y>Ly/4)*(y<Ly*3/4)

    bpm = Beam.BPM(n0=n0,k0=k0,dx=dx,dy=dy,Nx=Nx,Ny=Ny,Z=Z,dn=dn,initfield=initfield,verbose=True)
    bpm.solve() 
    viewer = Beam.Viewer(bpm)
    #viewer.viewrecord(0)
    viewer.viewintensity()

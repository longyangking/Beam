import sys
sys.path.append("../..")

import Beam
import numpy as np
import matplotlib.pyplot as plt

Nx = 300
Ny = 300

def dn(z):
    return np.zeros((Nx,Ny))

if __name__ == '__main__':
    dx = 1.0*10**-6
    dy = 1.0*10**-6
    Lx = dx*Nx
    Ly = dy*Ny
    delta = 15.0*10**-6
    n0 = 1.45
    k0 = 2.0*np.pi/(633.0*10**-9)

    Z = 30.0*10**-3
    
    x = np.linspace(0,Lx,Nx)
    y = np.linspace(0,Ly,Ny)
    x,y = np.meshgrid(x,y)
    
    initfield = np.exp(-((x-Lx/2)**2+(y-Ly/2)**2)/(delta**2))

    bpm = Beam.BPM3rd(n0=n0,k0=k0,dx=dx,dy=dy,Nx=Nx,Ny=Ny,dn=dn,verbose=True)
    fields = bpm.calculate(
            zs=[0, Z],
            field0=initfield
        ) 

    field = fields[0]
    plt.imshow(np.real(field))
    plt.show()

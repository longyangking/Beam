import sys
sys.path.append("../..")

import Beam
import numpy as np
import matplotlib.pyplot as plt

Nx = 50
Ny = 50

def dn(z):
    epsilon = np.zeros((Nx,Ny))
    for i in range(Nx):
        for j in range(Ny):
                epsilon[i,j] = 7.0**-4*(np.exp(-(((i-Nx/4.0)/(Nx/10.0))**2 + ((j-Ny*3/4.0)/(Ny/10.0))**2)**3)\
                    + np.exp(-(((i-Nx/4.0)/(Nx/10.0))**2 + ((j-Ny/4.0)/(Ny/10.0))**2)**3)\
                    + np.exp(-(((i-Nx*3/4.0)/(Nx/10.0))**2 + ((j-Ny/4.0)/(Ny/10.0))**2)**3))
    return epsilon

if __name__ == '__main__':
    dx = 1.0*10**-6
    dy = 1.0*10**-6
    ax = dx*Nx
    ay = dy*Ny
    n0 = 1.45
    k0 = 2.0*np.pi/(633.0*10**-9)

    fig1 = plt.figure()
    plt.imshow(dn(0),interpolation='none')   

    bandstructure = Beam.Bandstructure(ax,ay,n0,k0,Nx,Ny,dn,verbose=True)
    bands=3
    kx,ky,bulkband =  bandstructure.bulkband([-np.pi,0],[0,0],bands=bands)
    kxs = kx
    kys = ky
    bulkbands = bulkband

    kx,ky,bulkband =  bandstructure.bulkband([0,0],[np.pi,np.pi],bands=bands)
    kxs = np.concatenate((kxs,kx))
    kys = np.concatenate((kys,ky))
    bulkbands = np.concatenate((bulkbands,bulkband))

    kx,ky,bulkband =  bandstructure.bulkband([np.pi,np.pi],[2*np.pi,np.pi],bands=bands)
    kxs = np.concatenate((kxs,kx))
    kys = np.concatenate((kys,ky))
    bulkbands = np.concatenate((bulkbands,bulkband))

    fig2 = plt.figure()
    points = np.size(bulkbands,axis=0)
    for i in range(np.size(bulkbands,axis=1)):
        plt.scatter(range(points),np.real(bulkbands[:,i]))
    plt.title('Band Structure')
    plt.xlabel('Wavevector')
    plt.ylabel(r'Propagation Factor')
    plt.show()

    plt.show()
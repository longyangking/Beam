# Author: Yang Long<longyang_123@yeah.net>
#
# License: LGPL-2.1
import numpy as np

PA,PR,FU = 0,1,2

class BPM:
    def __init__(self,n0,k0,dx,dy,Nx,Ny,Z,dn=None,source=None,initfield=None,verbose=False):
        self.n0 = n0
        self.k0 = k0
        self.dx = dx
        self.dy = dy
        self.dz = 0.5*n0*k0*(0.5*dx**2+0.5*dy**2)
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = int(Z/dz)+1
        self.dn = dn                # dn should be a function about space
        self.source = source        # source should be a function about source
        self.initfield = initfield  # initfield should be a matrix

        self.fieldr = np.zeros([3,Nx,Ny],)
        self.fieldi = np.zeros([3,Nx,Ny])

        self.verbose = verbose

        if self.verbose:
            self.pin = 0
            self.record = np.zeros([50,Nx,Ny],dtype=complex)

    def laplacian(self,Z):
        Ztop = Z[0:-2,1:-1]
        Zleft = Z[1:-1,0:-2]
        Zbottom = Z[2:,1:-1]
        Zright = Z[1:-1,2:]
        Zcenter = Z[1:-1,1:-1]
        return (Ztop + Zleft + Zbottom + Zright - 4 * Zcenter)

    def solve(self):
        a =  1.0/(2.0*self.n0*self.k0)
        for t in range(self.Nz):
            sourcer = np.real(self.source(dz*t))
            deltafieldr = laplacian(self.fieldr[PR,:,:])
            #self.fieldi[FU,1:-1,1:-1] = self.fieldi[PA,1:-1,1:-1]  \
            #    + a*dz/dx**2 * deltafieldr \
            #    - dz * Vz(dz*t)* psi_r[PR,1:-1,1:-1] \
            #    + dz * sourcer[1:-1,1:-1]
                
            #deltapsi = laplacian(psi_i[PR,:,:])
            #psi_r[FU,1:-1,1:-1] = psi_r[PA,1:-1,1:-1]  \
            #    - a*dz/dx**2 * deltapsi \
            #    + dz * Vz(dz*(t+0.5))* psi_i[PR,1:-1,1:-1]\
            #    - dz * sourcei[1:-1,1:-1]

            if self.verbose:
                if t%(self.Nz/20) == 0:
                    print 'Finished = {num}%'.format(num=100.0*t/self.Nz)
                if t%(self.Nz/50) == 0:
                    self.record[self.pin,:,:] = self.fieldr[FR,:,:] + 1j*self.fieldi[FR,:,:]
                    self.pin += 1
    
    def output(self):
        return self.fieldr[FR,:,:] + 1j*self.fieldi[FR,:,:]

    def info(self):
        return self.record
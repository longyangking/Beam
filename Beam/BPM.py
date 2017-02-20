# Author: Yang Long<longyang_123@yeah.net>
#
# License: LGPL-2.1
import numpy as np

PR,FU = 0,1

class BPM:
    def __init__(self,n0,k0,dx,dy,Nx,Ny,Z,dn=None,source=None,initfield=None,verbose=False):
        self.n0 = n0
        self.k0 = k0
        self.dx = dx
        self.dy = dy
        self.dz = 0.5*n0*k0*(0.5*dx**2+0.5*dy**2)
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = int(Z/self.dz)+1
        self.dn = dn                # dn should be a function about space
        self.source = source        # source should be a function about source
        self.initfield = initfield  # initfield should be a matrix

        self.fieldr = np.zeros([2,Nx,Ny])
        self.fieldi = np.zeros([2,Nx,Ny])

        self.verbose = verbose

        if self.verbose:
            print 'dz = {dz}, Nz = {Nz}'.format(dz=self.dz,Nz=self.Nz)
            self.pin = 0
            self.record = np.zeros([50,Nx,Ny],dtype=complex)

    def laplacian(self,Z):
        Ztop = Z[0:-2,1:-1]
        Zleft = Z[1:-1,0:-2]
        Zbottom = Z[2:,1:-1]
        Zright = Z[1:-1,2:]
        Zcenter = Z[1:-1,1:-1]
        return (Ztop + Zleft + Zbottom + Zright - 4 * Zcenter)

    def solve(self,initfield=None,source=None):
        a =  1.0/(2.0*self.n0*self.k0)
        if initfield is not None:
            self.fieldr[PR] = np.real(initfield)
            self.fieldi[PR] = np.imag(initfield)
        elif self.initfield is not None:
            self.fieldr[PR] = np.real(self.initfield)
            self.fieldi[PR] = np.imag(self.initfield)

        if source is not None:
            source = self.source
        
        dz = self.dz
        dx = self.dx # TODO this will be removed in the future

        for t in range(self.Nz):
            # Imaginary part of field in t+1/2
            if source is None:
                sourcer = np.zeros([self.Nx,self.Ny])
            else:
                sourcer = np.real(source(dz*t))

            deltafieldr = self.laplacian(self.fieldr[PR,:,:])
            self.fieldi[FU,1:-1,1:-1] = self.fieldi[PR,1:-1,1:-1]  \
                + a*dz/dx**2 * deltafieldr \
                - dz*(-self.k0*self.dn(dz*t)[1:-1,1:-1])* self.fieldr[PR,1:-1,1:-1] \
                + dz*sourcer[1:-1,1:-1]

            # Boundary (will be replaced by PML)
            self.fieldi[FU,0,:] = 2*self.fieldi[FU,1,:] - self.fieldi[FU,2,:]
            self.fieldi[FU,-1,:] = 2*self.fieldi[FU,-2,:] - self.fieldi[FU,-3,:]
            self.fieldi[FU,:,0] = 2*self.fieldi[FU,:,1] - self.fieldi[FU,:,2]
            self.fieldi[FU,:,-1] = 2*self.fieldi[FU,:,-2] - self.fieldi[FU,:,-3]
            
            # Real part of field in t+1
            if source is None:
                sourcei = np.zeros([self.Nx,self.Ny])
            else:
                sourcei = np.imag(source(dz*(t+0.5)))

            deltapsi = self.laplacian(self.fieldi[FU,:,:])
            self.fieldr[FU,1:-1,1:-1] = self.fieldr[PR,1:-1,1:-1]  \
                - a*dz/dx**2 * deltapsi \
                + dz * (-self.k0*self.dn(dz*(t+0.5))[1:-1,1:-1])* self.fieldi[PR,1:-1,1:-1]\
                - dz * sourcei[1:-1,1:-1]

            # Boundary (will be replaced by PML)
            self.fieldr[FU,0,:] = 2*self.fieldr[FU,1,:] - self.fieldr[FU,2,:]
            self.fieldr[FU,-1,:] = 2*self.fieldr[FU,-2,:] - self.fieldr[FU,-3,:]
            self.fieldr[FU,:,0] = 2*self.fieldr[FU,:,1] - self.fieldr[FU,:,2]
            self.fieldr[FU,:,-1] = 2*self.fieldr[FU,:,-2] - self.fieldr[FU,:,-3]

            # Increment
            self.fieldr[PR] = self.fieldr[FU]
            self.fieldi[PR] = self.fieldi[FU]

            if self.verbose:
                if t%(self.Nz/20) == 0:
                    print 'Finished = {num}%'.format(num=100.0*t/self.Nz)
                if t%(self.Nz/50) == 0 and self.pin < 50:
                    self.record[self.pin,:,:] = self.fieldr[PR,:,:] + 1j*self.fieldi[PR,:,:]
                    self.pin += 1
        if self.verbose:
            print 'Calculation complete!'
    
    def output(self):
        return self.fieldr[PR,:,:] + 1j*self.fieldi[PR,:,:]

    def info(self):
        return self.record

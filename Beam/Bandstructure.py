# Author: Yang Long <longyang_123@yeah.net>
#
# License: LGPL-2.1

import numpy as np
import scipy as sp
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

class Bandstructure:
    def __init__(self,ax,ay,n0,k0,Nx,Ny,dn,Z0=None,Nz=None,tolerance=10.0**-12,verbose=False):
        self.ax = 1.0*ax
        self.ay = 1.0*ay
        self.Z0 = Z0

        self.n0 = n0
        self.k0 = k0
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz

        self.dx = self.ax/self.Nx
        self.dy = self.ay/self.Ny
        if (Z0 is not None) and (Nz is not None):
            self.dz = Z0/Nz
        else:
            self.dz = self.dx
        
        self.dn = dn

        self.tolerance = tolerance
        self.verbose = verbose
       
    def bulkband(self,kstart,kend,kpoints=15,bands=3):
        kxstart, kystart = kstart
        kxend, kyend = kend
        kx = np.linspace(kxstart,kxend,kpoints)
        ky = np.linspace(kystart,kyend,kpoints)

        bulkband = np.zeros((kpoints,bands),dtype=complex)
        for i in range(kpoints):
            bulkband[i] = self.eigenvalue(kx=kx[i],ky=ky[i],bands=bands)
        return kx,ky,bulkband
    
    #def eigenvalue(self,kx,ky,bands):
    #    eigenValues,eigenVectors = self.eigensystem(kx,ky,bands)
    #    return eigenValues[:bands]
        
    def eigenvalue(self,kx,ky,bands):
        C1 = self.dz/2/self.n0/self.k0/self.dx**2
        C2 = self.dz/2/self.n0/self.k0/self.dy**2
        C0 = self.dz*self.k0

        H = sparse.lil_matrix((self.Nx*self.Ny,self.Nx*self.Ny),dtype=complex)

        if self.Nz is not None:
            vectorZ = np.zeros((self.Nz,self.Nx*self.Ny,bands),dtype=complex)
            betaZ = np.zeros((self.Nz,bands),dtype=complex)

        if self.verbose:
            print 'Calculation with wave-vector ({kx},{ky})'.format(kx=kx,ky=ky)

        for n in range(self.Ny):
            for m in range(self.Nx):
                if m==0:
                    H[m + n*self.Nx, self.Nx-1 + n*self.Nx] = C1*np.exp(-1j*kx)
                else:
                    H[m + n*self.Nx, m - 1 + n*self.Nx] = C1

                if m==self.Nx-1:
                    H[m + n*self.Nx, n*self.Nx] = C1*np.exp(1j*kx)
                else:
                    H[m + n*self.Nx, m + 1 + n*self.Nx] = C1
                
                if n==0:
                    H[m + n*self.Nx, m + (self.Ny - 2)*self.Nx] = C2*np.exp(-1j*ky) # TODO Debug
                else:
                    H[m + n*self.Nx, m + (n-1)*self.Nx] = C2
                
                if n==self.Ny-1:
                    H[m + n*self.Nx, m] = C2*np.exp(1j*ky)
                else:
                    H[m + n*self.Nx, m + (n+1)*self.Nx] = C2

        if self.Nz is not None:
            for t in range(self.Nz):
                dn = self.dn(self.dz*t)
                for n in range(self.Ny):
                    for m in range(self.Nx):
                        H[m + n*self.Nx, m + n*self.Nx] = dn[m,n]*C0 - 2*C1 - 2*C2

                if self.verbose:
                    print 'Subprocess: step {step}/{totalstep} ...'.format(step=t+1,totalstep=self.Nz)
                eigenValues,eigenVectors = linalg.eigs(H,k=bands,which='LR',tol=self.tolerance)
                self.betaZ[t,:] = eigenValues
                self.vectorZ[t,:,:] = eigenVectors
        else:
            dn = self.dn(0)
            for n in range(self.Ny):
                for m in range(self.Nx):
                    H[m + n*self.Nx, m + n*self.Nx] = dn[m,n]*C0 - 2*C1 - 2*C2
            eigenValues,eigenVectors = linalg.eigs(H,k=bands,which='LR',tol=self.tolerance)

        if self.Nz is None:
            return eigenValues[:bands] #,eigenVectors[bands]
        
        
        U = np.zeros((bands,bands),dtype=complex)
        Uz = np.eye(bands)
        for t in range(self.Nz):
            if t != self.Nz-1:
                for ri in range(bands):
                    for rj in range(bands):
                        U[ri,rj] = vectorZ[:,ri,t+1].dot(vectorZ[:,rj,t])
                Uz = U*np.diag(np*exp(i*betaZ[:,t]))*Uz
            else:
                for ri in range(bands):
                    for rj in range(bands):
                        U[ri,rj] = vectorZ[:,ri,0].dot(vectorZ[:,rj,t])
                Uz = U*np.diag(np*exp(i*betaZ[:,t]))*Uz
        eigenValues,eigenVectors = np.imag(np.log(np.linalg.eig(Uz)))
        return eigenValues[:bands]
        
    def edgestate(self):
        pass

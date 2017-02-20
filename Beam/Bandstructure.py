# Author: Yang Long <longyang_123@yeah.net>
#
# License: LGPL-2.1

import numpy as np
import scipy as sp

class Bandstructure:
    def __init__(self,ax,ay,n0,k0,Nx,Ny,Z0=0,Nz=10,dn):
        self.ax = ax
        self.ay = ay
        self.Z0 = Z0
        self.n0 = n0
        self.k0 = k0
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        
        self.H = sp.sparse.lil_matrix((Nx*Ny,Nx*Ny))
        self.dn = dn
       
    def solve(self,Nx=None,Ny=None,Nz=None,dn=None,kpoints=15):
        if Nx is None:
            Nx = self.Nx
        if Ny is None:
            Ny = self.Ny
        if Nz is None:
            Nz = self.Nz
        if dn is None:
            dn = self.dn

        # TODO
        

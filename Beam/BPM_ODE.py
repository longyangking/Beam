import numpy as np
import scipy as sp 
from scipy.sparse import lil_matrix
from scipy.integrate import solve_ivp
import datetime

class BPM:
    '''
    Calculate beam propagation based on the third-party ODE solver from Scipy
    '''
    def __init__(self, dn, n0, k0, Nx, Ny, dx, dy, verbose=False):
        self.dn = dn
        self.n0 = n0
        self.k0 = k0
        self.Nx = Nx
        self.Ny = Ny
        self.dx = dx
        self.dy = dy

        self.verbose = verbose

        if self.verbose:
            print("Initiating the H0...", end="")

        self.H0 = self.__init_H0()

        if self.verbose:
            print("Done!")
            # stiffness = self.__stiffness()
            # if stiffness > 0.75:
            #     print("It may be stiff. The value: [{0}]".format(stiffness))
            # else:
            #     print("It may be non-stiff. The value: [{0}]".format(stiffness))

    def __init_H0(self):
        N0 = self.Nx*self.Ny
        H0 = lil_matrix((N0, N0), dtype=complex)

        # Parameter in BPM
        eps0 = 2*np.pi  # TODO TEST
        kappa0 = 1.0    # TODO TEST

        for i in range(self.Nx):
            for j in range(self.Ny):
                index0 = j*self.Nx + i
                H0[index0, index0] = eps0
                if i-1 >= 0:
                    index = j*self.Nx + i - 1
                    H0[index0, index] = kappa0
                if i+1 < self.Nx:
                    index = j*self.Nx + i + 1 
                    H0[index0, index] = kappa0
                if j-1 >= 0:
                    index = (j-1)*self.Nx + i
                    H0[index0, index] = kappa0
                if j+1 < self.Ny:
                    index = (j+1)*self.Nx + i
                    H0[index0, index] = kappa0

        return H0

    def hamiltonian(self, z):
        N0 = self.Nx*self.Ny
        H = lil_matrix((N0, N0), dtype=complex)

        dn0 = 1.0   # TODO TEST
        for i in range(self.Nx):
            for j in range(self.Ny):
                dn = self.dn(z)
                index = j*self.Nx + i
                H[index, index] = dn0*dn[i, j]

        return self.H0 + H

    def __stiffness(self):
        H = self.hamiltonian(0)
        #ldas = sp.sparse.linalg.eigsh(H)
        # ldas_sum = np.trace(H)
        # ldas_product = np.linalg.det(H)
        
        # if ldas_sum != 0:
        #     stiffness = ldas_product/ldas_sum
        # else:
        #     stiffness = float("inf")

        # return stiffness

    def calculate(self, zs, field0, method="RK45", vectorized=False, verbose=False):
        support_methods = ["RK45", "RK23", "DOP853", "BDF"]
        if method not in support_methods:
            raise Exception("Unsupported Method [{method}]!".format(method=method))

        if verbose:
            starttime = datetime.datetime.now()
            print("Perform the ODE solver with the method [{method}]...".format(method=method), end="")

        def dpsi(t, psi):
            H = self.hamiltonian(t)
            dpsi = -1j*H.dot(psi)

            return dpsi

        psi0 = np.array(field0, dtype=complex).reshape((self.Nx*self.Ny, ))

        z0, zf = np.min(zs), np.max(zs)
        sol = solve_ivp(
            fun=dpsi,
            t_span=[z0, zf], 
            y0=psi0,
            method=method,
            t_eval=zs,
            vectorized=vectorized
        )


        if verbose:
            endtime = datetime.datetime.now()
            second = (endtime - starttime).seconds
            minute = int(np.floor(second/60))
            second = second - 60*minute
            print("Done! Spend time: [{minute} min {second} sec]".format(
                minute=minute,
                second=second
            ))

        psis = sol.y
        fields = np.zeros((len(zs), self.Nx, self.Ny), dtype=complex)
        for i in range(len(zs)):
            fields[i, :, :] = psis[:,i].reshape(self.Nx, self.Ny)

        return fields
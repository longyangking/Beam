import numpy as np
import scipy as sp 
from scipy.integrate import solve_ivp

class BPM:
    '''
    Calculate beam propagation based on the third-party ODE solver from Scipy
    '''
    def __init__(self, dnfunc, n0, k0, Nx, Ny, dx, dy, verbose=False):
        self.dnfunc = dnfunc
        self.n0 = n0
        self.k0 = k0
        self.Nx = Nx
        self.Ny = Ny
        self.dx = dx
        self.dy = dy

        self.verbose = verbose

    def hamiltonian(self, z):
        H = 

    def calculate(self, zs, psi0, method="RK4"):
        support_methods = ["RK4", "RK23", "DOP853", "BDF"]
        if method not in support_methods:
            raise Exception("Unsupported Method [{method}]!".format(method=method))

        if self.verbose:
            starttime = datetime.datetime.now()
            print("Perform the ODE solver with the method [{method}]...".format(method=method), end="")

        def dpsi(t, psi):
            dpsi = 

            return dpsi

        z0, zf = np.min(zs), np.max(zs)
        sol = solve_ivp(
            fun=dpsi,
            t_span=[z0, zf], 
            y0=psi0,
            method=method,
            t_eval=zs
        )


        if self.verbose:
            endtime = datetime.datetime.now()
            second = (endtime - starttime).seconds
            minute = int(np.round(second/60))
            second = second - 60*minute
            print("Done! Spend time: [{minute} min {second} sec]".format(
                minute=minute,
                second=second
            ))

        return sol.y
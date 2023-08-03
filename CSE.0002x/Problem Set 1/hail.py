################################################################################
# CSE.0002x
# Problem Set 1: hail
# Name:

import math
import matplotlib.pyplot as plt
from IVPlib_pset1 import IVP
import IVPlib_pset1 as IVPlib


################################################################################
## Hail class definition (a subclass of IVP)
################################################################################


class HailIVP(IVP):

    def evalf(self, u, t):
        """
        Calculates forcing for hail trajectory given the current state and time.
        IMPORTANT: Some of the provided helper functions assume u[0] is the velocity (V)
        and u[1] = the altitude (z). And thus f[0] = dV/dt and f[1] = dz/dt.

        Args:
            u (float list): current state
            t (float): current time

        Returns:
            f (float list): current forcing
        """
        #### BEGIN SOLUTION ####
        #### REMEMBER: use self.get_p to get values of any parameters needed
        
        # Assign values of V and z.
        V = u[0]
        z = u[1]
        
        # get parameter values
        rhoa = self.get_p["rhoa"] # kg/m^3
        rhop = self.get_p["rhop"] # kg/m^3
        dp = self.get_p["dp"] # m
        CD = self.get_p["CD"]
        g = self.get_p["g"] # m/s^2
        
        # calculate mass
        m = (math.pi/6) * rhop * dp**3
        
        # calculate cross-sectional area of the hail
        Ap = (math.pi/4) * dp**2
        
        # calculate D (drag force) as a function of V.
        D = 0.5 * rhoa * V**2 * Ap * CD
        
        # calculate f(V, z) = (dV/dt, dz/dt)
        f = []
        dV_dt = g - (D/m)  # calculate dV/dt
        dz_dt = -V  # calculate dz/dt
        f.append(dV_dt)
        f.append(dz_dt)
        return f
        
        #### END SOLUTION ####


################################################################################
## Functions to solve hail trajectory IVP
################################################################################


def hail_Verror(hail_IVP, t, V):
    """
    Calculates the exact (analytic solution) hail velocity Vex and the absolute error e.

    Args:
        hail_IVP (IVP object): IVP for hail trajectory
        t (float list): time values t[n]
        V (float list): velocity values from numerical method, V[n]

    Returns:
        e (float list): velocity error e[n] = abs(V[n]-Vex[n])
        Vex (float list): exact (analytic solution) velocity values, Vex[n]
    """
    #### BEGIN SOLUTION ####
    # get parameter values
    rhoa = hail_IVP.get_p["rhoa"] # kg/m^3
    rhop = hail_IVP.get_p["rhop"] # kg/m^3
    dp = hail_IVP.get_p["dp"] # m
    CD = hail_IVP.get_p["CD"]
    g = hail_IVP.get_p["g"] # m/s^2
    
    # calculate mass
    m = (math.pi/6) * rhop * dp**3
    
    # calculate cross-sectional area of the hail
    Ap = (math.pi/4) * dp**2
    
    # calculate Vterm (terminal velocity)
    Vterm = (2*m*g)/(rhoa * Ap * CD)
    
    # calculate Vex (exact velocity from analytic method)
    Vex = []
    C = 0  # C = 0, if V[0] = 0. Assume V[0] = 0.
    for velocity in V:
        Vex.append(Vterm * math.tanh((g*t)/Vterm + C))
    
    # calculate error: e[n] = abs(V[n]-Vex[n])
    if len(V) != len(Vex):
        raise ValueError("Something went terribly wrong.")
    e = []
    for i in range(len(V)):
        e.append(abs(V[i] - Vex[i]))
    
    return e, Vex
    #### END SOLUTION ####


def hail_Veplot(t, V, Vex, e, method='numerical'):
    """
    Produces a single figure with subplots of:
        * the numerical and exact velocity (overlayed) vs. time
        * the error vs. time
    The figure will be labeled with the method name.

    Args:
        t (float list): time values t[n]
        V (float list): velocity values from numerical method, V[n]
        Vex (float list): exact (analytic solution) velocity values, Vex[n]
        e (float list): velocity error e[n]
        method (string, optional): name of the numerical method used to generate results

    Returns:
        axs (array of Axes): handle to the Axes objects that comprise the figure's subplots
    """
    #### BEGIN SOLUTION ####
    raise NotImplementedError("Implement hail_Veplot")
    #### END SOLUTION ####


def hail_Vzplot(t, V, z, method='numerical'):
    """
    Produces a single figure with subplots of:
        * velocity vs. time
        * z vs. time
    The figure will be labeled with the method name.

    Args:
        t (float list): time values t[n]
        V (float list): velocity values from numerical method, V[n]
        z (float list): z values from numerical method, z[n]
        method (string, optional): name of the numerical method used to generate results

    Returns:
        axs (array of Axes): handle to the Axes objects that comprise the figure's subplots
    """
    #### BEGIN SOLUTION ####
    raise NotImplementedError("Implement hail_Vzplot")
    #### END SOLUTION ####


def hail_run_case(hail_IVP, dt, mlist):
    """
    Solves the hail trajectory problem described in hail_IVP with a timestep dt for
    the methods in mlist and plots the solution and error as functions of time.

    Args:
        hail_IVP (IVP object): describes IVP case to be simulated
        dt (float): timestep
        mlist (function list): list of numerical methods to run on case

    Returns:
        results (dictionary): maps each method name to a tuple of two elements:
            * the Axes objects in the figure that plots the V and error for that method
            * the Axes objects in the figure that plots the V and z for that method
    """
    results = {}

    for method in mlist:
        t, u = IVPlib.solve(hail_IVP, dt, method)
        V = []
        z = []
        for i in range(len(t)):
            V.append(u[i][0])
            z.append(u[i][1])
        e, Vex = hail_Verror(hail_IVP, t, V)
        emax = max(e)

        mname = method.__name__
        print(f"Max error = {emax:.3e}  for  {mname}")
        axs_ve = hail_Veplot(t, V, Vex, e, mname)
        axs_vz = hail_Vzplot(t, V, z, mname)
        results[mname] = (axs_ve, axs_vz)

    return results


def hail_run_conv_study(hail_IVP, dtv, mlist, nflist):
    """
    Solves the hail trajectory problem described in hail_IVP for all of the timesteps
    in dtv and using all of the methods in mlist. Then plots the error as a function
    of dt and the number of function evaluations.

    Args:
        hail_IVP (IVP object): describes IVP case to be simulated
        dtv (float list): timesteps
        mlist (function list): list of numerical methods to run on case
        nflist (float list): for each method in mlist, gives the number of function
            evaluations per timestep the method uses
    """
    fig_evsdt, axs_evsdt = plt.subplots()
    fig_evsnf, axs_evsnf = plt.subplots()

    for m in range(len(mlist)):
        emax = []
        nf = []
        method = mlist[m]
        for dt in dtv:
            t, u = IVPlib.solve(hail_IVP, dt, method)
            V = []
            for i in range(len(t)):
                V.append(u[i][0])
            e, Vex = hail_Verror(hail_IVP, t, V)
            emax.append(max(e))
            nf.append(nflist[m]*len(t))

        axs_evsdt.loglog(dtv, emax, label=method.__name__)
        axs_evsnf.loglog(nf, emax, label=method.__name__)

    axs_evsdt.grid(True, which='both', ls='-')
    axs_evsdt.set_xlabel(r'$\Delta t$')
    axs_evsdt.set_ylabel('Max error')
    axs_evsdt.legend(loc='lower right')

    axs_evsnf.grid(True, which='both', ls='-')
    axs_evsnf.set_xlabel('Number of f evaluations')
    axs_evsnf.set_ylabel('Max error')
    axs_evsnf.legend(loc='lower left')


if __name__ == '__main__':
    uI = [0, 5000]
    tI = 0
    tF = 20 # s

    p = {}
    p['rhoa'] = 1.0 # kg/m^3
    p['rhop'] = 700.0 # kg/m^3
    p['dp'] = 0.02 # m
    p['CD'] = 0.5
    p['g'] = 9.81 # m/s^2

    hail_IVP = HailIVP(uI, tI, tF, p)

    # Set up which methods to run
    mlist = [IVPlib.step_FE, IVPlib.step_RK4] # change this to run only one of the methods
    nflist = [1, 4] # number of f evaluations used per timestep each method

    # Simulate a single case to demonstrate basic behavior of methods with dt=1.0
    hail_run_case(hail_IVP, 1.0, mlist)

    # Run varying dt for all methods and plot emax vs dt and emax vs function calls
    dtv = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28]
    hail_run_conv_study(hail_IVP, dtv, mlist, nflist)

    plt.show()

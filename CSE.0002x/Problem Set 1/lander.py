################################################################################
# CSE.0002x
# Problem Set 1: lander
# Name: Atakan Devrent

import math
import matplotlib.pyplot as plt
from IVPlib_pset1 import IVP
import IVPlib_pset1 as IVPlib


################################################################################
## Lander class definition (a subclass of IVP)
################################################################################


class LanderIVP(IVP):

    def evalf(self, u, t):
        """
        Calculates forcing for a Martian lander given the current state and time.

        Args:
            u (float list): current state
            t (float): current time

        Returns:
            f (float list): current forcing
        """
        #### BEGIN SOLUTION ####
        
        # Unpack u to V, z
        V = u[0]
        z = u[1]
        #print("u[1] =", u[1])#debug
        
        # load parameter values
        
        # assign atmospheric parameters
        #print("z =", z)#debug
        param = self.atmosphere(z)  # dict
        Ta = param["Ta"]
        pa = param["pa"]
        rhoa = param["rhoa"]
        g = param["g"]

        m_l = self.get_p("m_l")  # kg
        A_l = self.get_p("A_l")  # m^2
        A_p = self.get_p("A_p")  # m^2
        CD_l = self.get_p("CD_l")  # constant
        CD_p = self.get_p("CD_p")  # constant
        theta_e = math.radians(self.get_p("theta_e"))  # rad
        theta_p = math.radians(self.get_p("theta_p"))  # rad
        
        D_l = 0.5 * rhoa * V**2 * A_l * CD_l
        D_p = 0.5 * rhoa * V**2 * A_p * CD_p
        
        # params depending on if parachute is open or not.
        V_p = self.get_p("V_p")  # (m/s) Parachute is deployed if current velocity <= V_p.
        if V >= V_p:  # parachute is not deployed.
            theta = theta_e
            D = D_l
        else:  # parachute is deployed.
            theta = theta_p
            D = D_l + D_p
        
        # calculate f(u, t) where u = (V, z)
        f = []
        dV_dt = g * math.cos(theta) - D/m_l  # f[0] is dV/dt
        dz_dt = -V * math.cos(theta)  # f[1] is dz_dt
        f.append(dV_dt)
        f.append(dz_dt)
        return f
        #### END SOLUTION ####

    def atmosphere(self, z):
        """
        Returns in a dictionary the properties of the Martian atmosphere at the
        altitude z (assumed in meters). Currently, the properties returned are:
            'Ta': temperature (K)
            'pa': pressure (kPa)
            'rhoa': density (kg/m^3)
            'g': gravity (m/s^2)

        Args:
            z (float): altitude in meters

        Returns:
            props (dictionary): property keys and values
        """
        #print("Type of z:", type(z))#debug
        if z > 7000:
            Ta = 235 - 0.00108*(z-7000) # K
        else:
            Ta = 235 - 0.00222*(z-7000) # K
        pa = 0.699*math.exp(-0.00009*z) # kPa
        rhoa = pa / (0.1921*Ta) # kg/m^3

        G = 6.674e-11 # gravitational constant (m^3/ kg/s^2)
        Mm = 6.41693e23 # mass of Mars (kg)
        Rm = 3389.5*1000 # radius of Mars (m)
        g = G*Mm / (Rm+z)**2 # m/s^2

        props = {}
        props['Ta'] = Ta
        props['pa'] = pa
        props['rhoa'] = rhoa
        props['g'] = g
        return props


################################################################################
## Functions to solve Martian lander trajectory IVP
################################################################################


def lander_Vzaplot(t, V, z, a, method='numerical'):
    """
    Produces a single figure with subplots of:
        * velocity vs. time
        * z vs. time
        * acceleration vs. time
    The figure will be labeled with the method name.

    Args:
        t (float list): time values t[n]
        V (float list): velocity values from numerical method, V[n] (in m/s)
        z (float list): z values from numerical method, z[n] (in km)
        a (float list): acceleration values from numerical method, a[n] (in m/s**2)
        method (string, optional): name of the numerical method used to generate results

    Returns:
        axs (array of Axes): handle to the Axes objects that comprise the figure's subplots
    """
    #### BEGIN SOLUTION ####
    fig, axs = plt.subplots(3, 1, sharex=True)
    axs[0].plot(t, V, color = "b", label = str(method))
    axs[0].set_ylabel("$V$ (m/s)")
    axs[0].legend(loc = "upper right")
    axs[1].plot(t, z, color = "b")
    axs[1].set_ylabel("$z$ (km)")
    axs[2].plot(t, a, color = "b")
    axs[2].set_ylabel("$a$ (m/s^2)")
    axs[2].set_xlabel("$t$ (sec)")
    
    return axs
    
    
    #### END SOLUTION ####

def lander_run_case(lander_IVP, dt, mlist):
    """
    Solves the Martian lander problem described in lander_IVP with a timestep dt for
    the methods in mlist. For each method in mlist, this function will:
        * plot the solution by calling lander_Vzaplot
        * print the time (in sec) and altitude (in km) at which the parachute deployed
        * print the altitude (in km) and velocity (in m/s) at the final time

    Args:
        lander_IVP (IVP object): describes IVP case to be simulated
        dt (float): timestep
        mlist (function list): list of numerical methods to run on case

    Returns:
        results (dictionary): maps each method name to a tuple of three elements:
            * the Axes objects in the figure corresponding to that method
            * the time of parachute deployment
            * the altitude (z) at which deployment occurred
    """
    #### BEGIN SOLUTION ####
    results = {}
    
    for method in mlist:
        parachute_deployed = False
        t, u = IVPlib.solve(lander_IVP, dt, method)
        V = []
        z = []
        a = []
        for i in range(len(t)):
            V.append(u[i][0])
            z.append(u[i][1])
            fa = lander_IVP.evalf(u[i], t[i])
            a.append(fa[0])  # f[0] is dV/dT = a
            
            deployment_condition = (u[i][0] <= lander_IVP.get_p("V_p")) and (parachute_deployed == False)
            if deployment_condition:
                parachute_deployed = True
                #print(u[i][0], "<=", lander_IVP.get_p("V_p"))#debug
                #print(parachute_deployed)#debug
                deployment_time = t[i]
                deployment_z = u[i][1]
        
        mname = method.__name__
        z_km = []
        for altitude in z:
            z_km.append(altitude/1000)  # convert to km
        axs_object = lander_Vzaplot(t, V, z_km, a, method = mname)
        deployment_z_km = deployment_z/1000
        results[mname] = (axs_object, deployment_time, deployment_z_km)
        final_z = z_km[len(z_km)-1]
        final_V = V[len(V)-1]
        
        print("Method:", mname)
        print("------------------------------------------")
        print("Parachute deployed at t, z:", "{0:.2e}".format(deployment_time), "s,", "{0:.2e}".format(deployment_z_km), "km")
        print("Final z, V:", "{0:.2e}".format(final_z), "km,", "{0:.2e}".format(final_V), "m/s")
        
    return results
            
    
    #### END SOLUTION ####


if __name__ == '__main__':
    uI = [5800.0, 125*1000.0] # Initial velocity (m/s) and altitude (m)
    tI = 0
    tF = 300.0 # s

    p = {}
    p['m_l'] = 3300.0 # kg
    p['V_p'] = 470.0 # parachute deploy velocity (m/s)
    p['A_l'] = 15.9 # m^2
    p['CD_l'] = 1.7
    p['A_p'] = 201.0 # m^2
    p['CD_p'] = 1.2
    p['theta_e'] = 83.0 # flight path angle during entry, in degrees
    p['theta_p'] = 70.0 # flight path angle during parachute, in degrees

    lander_IVP = LanderIVP(uI, tI, tF, p)

    # Set up which method(s) to run
    mlist = [IVPlib.step_RK4]

    # Simulate a single case to demonstrate basic behavior of methods with dt=0.1 s
    lander_run_case(lander_IVP, 0.1, mlist)

    plt.show()

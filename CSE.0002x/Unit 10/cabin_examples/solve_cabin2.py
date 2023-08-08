import matplotlib.pyplot as plt
import IVPlib_rev4 as IVPlib
from IVPlib_rev4 import IVP

class Cabin2IVP(IVP):
    def evalf(self, u, t):
        """        
        Args:
            u (float list): current temperature of cabin.
            t (float): current time
    
        Returns:
            f (float list): returns du/dt
        """

        m0cc = self.get_p('m0')*self.get_p('cc')
        m1cc = self.get_p('m1')*self.get_p('cc')

        hA_01    = self.get_p('h_01')   * self.get_p('A_01')

        hA_0r    = self.get_p('h_0r')   * self.get_p('A_0r')
        hA_0b    = self.get_p('h_0b')   * self.get_p('A_0b')
        hA_0f    = self.get_p('h_0f')   * self.get_p('A_0f')
        hA_0top  = self.get_p('h_0top') * self.get_p('A_0top')
        hA_0bot  = self.get_p('h_0bot') * self.get_p('A_0bot')

        hA_1l    = self.get_p('h_1l')   * self.get_p('A_1l')
        hA_1b    = self.get_p('h_1b')   * self.get_p('A_1b')
        hA_1f    = self.get_p('h_1f')   * self.get_p('A_1f')
        hA_1top  = self.get_p('h_1top') * self.get_p('A_1top')
        hA_1bot  = self.get_p('h_1bot') * self.get_p('A_1bot')
        
        u_out    = self.get_p('u_out')
        u_ground = self.get_p('u_ground')
        
        Q0_stove  = self.get_p('Q_stove')

        A00 = -(hA_01 + hA_0r + hA_0b + hA_0f + hA_0top + hA_0bot)/m0cc
        A01 = hA_01/m0cc
        
        A10 = hA_01/m1cc
        A11 = -(hA_01 + hA_1l + hA_1b + hA_1f + hA_1top + hA_1bot)/m1cc
        
        Q0_out    = (hA_0r + hA_0b + hA_0f + hA_0top)*u_out(t)
        Q0_ground = hA_0bot*u_ground(t)
        b0 = (Q0_out + Q0_ground + Q0_stove(t))/m0cc
        
        Q1_out    = (hA_1l + hA_1b + hA_1f + hA_1top)*u_out(t)
        Q1_ground = hA_1bot*u_ground(t)
        b1 = (Q1_out + Q1_ground)/m1cc
        
        f0 = A00*u[0] + A01*u[1] + b0
        f1 = A10*u[0] + A11*u[1] + b1

        return [f0, f1]

def u_ground(t):
    """
        Args:
            t (float): current time in seconds
    
        Returns:
            u_ground (float): ground temperature at time t in deg C
    """
    return 0.0

def u_out(t):
    """
        Args:
            t (float): current time in seconds
    
        Returns:
            u_out (float): outside temperature at time t in deg C
    """
    return 0.0

def Q_stove(t):
    """
        Args:
            t (float): current time in seconds
    
        Returns:
            Q_stove (float): heat generated by stove at time t in Watts
    """
    return 5.0e3


A_0r     = 25.  # room0 area of right (m^2)
A_0fb    = 25.  # room0 area of front, back (m^2)
A_0bt    = 50.  # room0 area of bottom, top (m^2)
A_01     = 25.  # area of shared wall between room0 and room1 (m^2)
A_1l     = 25.  # room1 area of left (m^2)
A_1fb    = 12.5 # room1 area of front, back (m^2)
A_1bt    = 25.  # room1 area of bottom, top (m^2)
h_ground = 0.04 # ground heat transfer coefficient (W/(m^2 C))
h_roof   = 0.4 # roof heat transfer coefficient (W/(m^2 C))
h_int    = 2.0 # interior walls heat transfer coefficient (W/(m^2 C))
h_ext    = 2.0 # exterior walls heat transfer coefficient (W/(m^2 C))
m0       = 200. # mass of air in main room
m1       = 100. # mass of air in bath room
cc       = 700.0 # J / (kg C)
        
u0I   = [0.0, 0.0] # Initial temperature of cabin (C)
tFmin = 300.0 # final time to simulate to (min)
dtmin = 5e-1 # time increment to give solutions at (min)

# Convert times to seconds
tF = tFmin*60
dt = dtmin*60

# Initialize CabinIVP object
p = {}
p['A_0r']     = A_0r
p['A_0f']     = A_0fb
p['A_0b']     = A_0fb
p['A_0bot']   = A_0bt
p['A_0top']   = A_0bt

p['A_01']     = A_01

p['A_1l']     = A_1l
p['A_1f']     = A_1fb
p['A_1b']     = A_1fb
p['A_1bot']   = A_1bt
p['A_1top']   = A_1bt

p['h_0r']   = h_ext
p['h_0f']   = h_ext
p['h_0b']   = h_ext
p['h_0bot'] = h_ground
p['h_0top'] = h_roof

p['h_01']   = h_int

p['h_1l']   = h_ext
p['h_1f']   = h_ext
p['h_1b']   = h_ext
p['h_1bot'] = h_ground
p['h_1top'] = h_roof

p['m0'] = m0
p['m1'] = m1
p['cc'] = cc

p['u_out']    = u_out
p['u_ground'] = u_ground
p['Q_stove']  = Q_stove

cabin2IVP = Cabin2IVP(u0I, 0.0, tF, p)

# Solve cabin IVP
t, v = IVPlib.solve(cabin2IVP, dt, IVPlib.step_FE)

# Extract cabin and bathroom temperature from v and scale time to minutes
u0 = []
u1 = []
for n in range(len(t)):
    t[n] = t[n]/60.0 # convert to minutes
    u0.append(v[n][0])
    u1.append(v[n][1])

# Plot   
fig, ax = plt.subplots()
ax.scatter(t,u0,marker='o',label='main')
ax.scatter(t,u1,marker='x',label='bath')
ax.set_xlabel('t (min)')
ax.set_ylabel('temperature (C)')
ax.grid(True)
ax.legend()

import numpy as np
import matplotlib.pyplot as plt

def calc_xf(V0, th0, g):
    """
    Calculates the impact location (xf) of a projectile
    given the initial speed, angle, and gravity.

    Args:
        V0 (float or numpy array of floats): Initial speed
        th0 (loat or numpy array of floats): Initial angle (degrees)
        g (float): gravity.

    Returns:
        xf as either a float or numpy array (if V0 and th0 are numpy arrays)

    """
    xf = (V0**2)/g*np.sin(2*th0*np.pi/180.)
    return xf

def run_xfMC(V0lim, th0lim, g, N):
    """
    Perform a Monte Carlo simulation of sample size N calculating xf
    the point of impact of a projectile on the ground with  initial
    speed V0 and initial angle th0 are uniformally distributed.

    Args:
        V0lim (tuple of floats): V0lim[0] gives minimum V0 and
        V0lim[1] gives maximum V0
        th0lim (tuple of floats): th0lim[0] gives minimum th0 and
        th0lim[1] gives maximum th0.  Note: th0lim will be in degrees.
        g (float): gravity
        N (integer): size of sample to run Monte Carlo on

    Returns:
        xf (numpy array of floats): point of impact for all instances in sample

    """
    rng = np.random.default_rng()

    V0  = rng.uniform(V0lim[0],   V0lim[1], N)
    th0 = rng.uniform(th0lim[0], th0lim[1], N)

    xf = calc_xf(V0,th0,g)

    return xf


def calc_ptarget(xf, xftarget):
    """
    Calculate probability (fraction) of instances of xf which
    fall between xftarget[0] < xf <  xftarget[1]

    Args:
        xf (numpy array of floats): impact locations
        xftarget (tuple of floats): gives desired impact location range

    Returns:
        Ptarget (float): estimated probability

    """
    Ntarget = np.count_nonzero(np.logical_and(xf>xftarget[0],xf<xftarget[1]))
    Ptarget = Ntarget/len(xf)
    return Ptarget

# Main routine

V0lim = (28.,32.)
th0lim = (25.,35.)
xftarget = (78.,81.)

N   = int(1e8)

xf = run_xfMC(V0lim, th0lim, 9.81, N)
plt.hist(xf,min(200,int(N**0.5)),density=True)
plt.xlabel('$x_f$ (m)')
plt.ylabel('Density')

Ptarget = calc_ptarget(xf, xftarget)
xfmean = xf.mean()
print("N = {:.1e}: Ptarget = {:.4f} xfmean = {:.2f} m".format(N,Ptarget,xfmean))
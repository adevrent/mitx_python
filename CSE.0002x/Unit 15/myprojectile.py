import numpy as np
import matplotlib.pyplot as plt

def calc_xf(V0, th0, g):
    """Calculates impact location xf

    Args:
        V0 (ndarray of floats): Values of initial velocity, randomly drawn from a distribution
        th0 (ndarray of floats): Values of initial angle in DEGREES, randomly drawn from a distribution
        g (float): gravitational constant

    Returns:
        xf (ndarray of floats): impact location
    """
    xf = (V0**2/g)*np.sin(2*th0*np.pi/180)
    return xf

def uniform_V0_th0(V0lim, th0lim, N):
    """Creates length N arrays of V0 and th0 according to the uniform distribution between limits.

    Args:
        V0lim (tuple of floats): limits of the uniform distribution of V0
        th0lim  (tuple of floats): limits of the uniform distribution of th0 in DEGREES
        N (int): number of values to be generated.
        
    Returns:
        V0 (ndarray of floats): An array of initial velocities to be used in a Monte Carlo simulation.
        th0 (ndarray of floats): An array of initial angles to be used in a Monte Carlo simulation.
    """
    rng = np.random.default_rng()
    V0 = rng.uniform(V0lim[0], V0lim[1], N)
    th0 = rng.uniform(th0lim[0], th0lim[1], N)
    
    return V0, th0

def triangular_V0_th0(V0lim, th0lim, N):
    """Creates length N arrays of V0 and th0 according to the triangular distribution between limits.

    Args:
        V0lim (tuple of floats): limits of the triangular distribution of V0
        th0lim  (tuple of floats): limits of the triangular distribution of th0 in DEGREES
        N (int): number of values to be generated.
        
    Returns:
        V0 (ndarray of floats): An array of initial velocities to be used in a Monte Carlo simulation.
        th0 (ndarray of floats): An array of initial angles to be used in a Monte Carlo simulation.
    """
    rng = np.random.default_rng()
    V0 = rng.triangular(V0lim[0], 0.5*(V0lim[0] + V0lim[1]), V0lim[1], N)
    th0 = rng.triangular(th0lim[0], 0.5*(th0lim[0] + th0lim[1]), th0lim[1], N)
    
    return V0, th0

N = 1000000
V0, th0 = triangular_V0_th0((28, 32), (25, 35), N)
g = 9.81  # m/s^2

# Calculate xf array and statistics
xf = calc_xf(V0, th0, g)
mean = xf.mean()

# Print statistics
print("Expected value (mean) of xf =", mean)

# Plot histogram
fig, ax = plt.subplots()
fig.set_size_inches(10, 8)
fig.text(0.1, 0, f"Distribution of $x_f$ for triangular distribution of $V_0$ and $\Theta_0$ using a Monte Carlo simulation with an $N=${N:.2E} sample size.", fontsize="large")
ax.hist(xf, bins=200, density=True)
ax.set_title("Probability density of impact location $x_f$")
ax.set_xlabel("$x_f$")
ax.set_ylabel("Probability density")
ax.vlines(mean, 0, 0.07, "r", label="mean")
ax.legend()
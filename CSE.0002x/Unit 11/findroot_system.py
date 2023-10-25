################################################################################
# Intro to CSE
#
# Newton method implemented to find the root of a 
# system of equations.

import numpy as np
import matplotlib.pyplot as plt

def calcr(x):   
    return calcr_all(x, calc_r=True)
  
def calcdrdx(x):   
    return calcr_all(x, calc_drdx=True)

def calcr_all(x, calc_r=False, calc_drdx=False):
    """
    Calculate r and dr/dx

    Args:

    x: NumPy 1D array

    calc_r: boolean.
            If True, then calculate r(x)

    calc_drdx: boolean.
            If True, then calculate dr/dx(x)

    Returns:
        if calc_r and calc_drdx: return tuple (r, r_x)
        else return just r or just r_x
    """
    
    # Calculate r if requested
    if calc_r:
        r0 = x[0]**2 + x[1]**2 - 1
        r1 = x[1] - np.sin(x[0])
        r = np.array([r0, r1])

    # Calculate dr/dx if requested
    if calc_drdx:
        r0_x0 = 2*x[0]
        r0_x1 = 2*x[1]
        r1_x0 = -np.cos(x[0])
        r1_x1 = 1.
        r_x = np.array([[r0_x0, r0_x1],[r1_x0, r1_x1]])

    if calc_r and calc_drdx:
        return (r, r_x)
    elif calc_r:
        return r
    elif calc_drdx:
        return r_x
    else:
        raise ValueError("calc_r and calc_drdx both False")

def findroot_Newton(calc_r, calc_drdx, x0, n, ax=None):
    """
        Runs the Newton-Raphson method for finding an approximation to a root of r

        Args:
            calc_r (scalar function reference): function which returns r(x) 
                when called as calc_r(x)
                
            calc_drdx (scalar function reference): function which returns dr/dx(x) 
                when called as calc_drdx(x)

            x0 (float): initial guess of x

            n (int): number of iterations to take
                       
            ax (matplotlib pyplot Axes object reference): if provided,
            the plot the locations in the xk, r(xk) plane for each
            Newton iterate
                
        Returns:
            x: approximation to a root of f, i.e., such that r(x) = 0

        """

    xk = x0
    print(f"Iteration {0:2d}: x = {xk}")
    for i in range(n):
        rk = calc_r(xk)

        if (ax is not None): # Plot location of current iterate 
            # WARNING: This only plots state 0 and 1 locations!
            ax.plot(xk[0], xk[1], 'r*')
            plt.pause(1.0)

        drdxk = calc_drdx(xk)
        dx = np.linalg.solve(drdxk, -rk)
        xk = xk + dx
        
        print(f"Iteration {i+1:2d}: x = {xk}")

    if (ax is not None): # Plot location of final iterate
        ax.plot(xk[0], xk[1], 'g*')

    return xk


################################################################################
## Main body
################################################################################

if __name__ == "__main__":
    x0circ = np.linspace(-1.,1.,1001)
    x0sin  = np.linspace(-2.,2.,1001)
    fig, ax = plt.subplots()
    ax.plot(x0circ,  (1-x0circ**2)**0.5,'b-')
    ax.plot(x0circ, -(1-x0circ**2)**0.5,'b-')
    ax.plot(x0sin, np.sin(x0sin),'r-')
    ax.grid(True)
    ax.axis('Equal')
    ax.set_xlabel('$x_0$')
    ax.set_ylabel('$x_1$')
    

    print()
    print("Running Newton")
    findroot_Newton(calcr, calcdrdx, [1.5, 1.5], 6, ax=ax)
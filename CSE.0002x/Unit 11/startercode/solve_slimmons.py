################################################################################
# CSE.0002x
# Problem Set 2: solve_slimmons
# DO NOT MODIFY THIS FILE!

"""
This sets up and solves the equilibrium and then the temperature evolution
after Intro to CSE students reduce their heating in Slimmons residence hall.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from HVACIVP import HVACIVP
from mylinsolver import myGE, myGE_vec


def printnpfloat(u):
    """
    Args: u (float)

    Returns: string with u in desired format
    """
    return f"{u:.6e}"


def uground(t):
    """
    Args:
        t (float): current time in seconds

    Returns:
        uground (float): ground temperature at time t in deg C
    """
    return 5.0


def uair(t):
    """
    Args:
        t (float): current time in seconds

    Returns:
        uout (float): outside air temperature at time t in deg C
    """
    return -10.0


def Qroom(t):
    """
    Args:
        t (float): current time in seconds

    Returns:
        Qroom (float): heat generated by heating system in room at time t in Watts
    """
    return 1.5e3


def Qdin(t):
    """
    Args:
        t (float): current time in seconds

    Returns:
        Qdin (float): heat generated by heating system in dining room at time t in Watts
    """
    return 1.5e4


def Qc(t):
    """
    Args:
        t (float): current time in seconds

    Returns:
        Qc (float): heat generated by a heater in C rooms at time t in Watts
    """
    if t is None:
        return 1.5e3
    else:
        return 0.9e3


def Qs(t):
    """
    Args:
        t (float): current time in seconds

    Returns:
        Qs (float): heat generated by a heater in S rooms at time t in Watts
    """
    if t is None:
        return 1.5e3
    else:
        return 0.9e3


def Qe(t):
    """
    Args:
        t (float): current time in seconds

    Returns:
        Qe (float): heat generated by a heater in E rooms at time t in Watts
    """
    if t is None:
        return 1.5e3
    else:
        return 0.9e3


def get_HVACIVP_pdict():
    """
    Construct and return the dictionary of parameters for an HVACIVP object
    """
    p = {}
    p['Afloor'] = 20. # floor (and roof) area (m^2)
    p['Afb'] = 10. # front, back area (m^2)
    p['Awall'] = 10. # area of side walls (m^2)
    p['Adinfloor'] = 5*p['Afloor'] # dining room floor (and ceiling) area (m^2)
    p['Adinfb'] = 10*p['Afb'] # dining room front, back area (m^2)

    p['hint'] = 2.0 # interior walls heat transfer coefficient (W/(m^2 C))
    p['hext'] = 2.0 # exterior walls heat transfer coefficient (W/(m^2 C))
    p['hground'] = 0.04 # ground heat transfer coefficient (W/(m^2 C))
    p['hfloor'] = 1.0 # floor (internal) heat transfer coefficient (W/(m^2 C))
    p['hroof'] = 0.4 # roof heat transfer coefficient (W/(m^2 C))

    p['mroom'] = 50. # mass of air in main room
    p['mdin'] = 10*p['mroom'] # mass of air in dining room
    p['cc'] = 700.0 # J / (kg C)

    p['uair'] = uair
    p['uground'] = uground
    p['Qroom'] = Qroom
    p['Qdin'] = Qdin
    p['Qc'] = Qc
    p['Qs'] = Qs
    p['Qe'] = Qe

    return p


if __name__ == '__main__':
    np.set_printoptions(formatter={'float': printnpfloat})

    uI = np.zeros(271) # initial temperature of cabin (C)
    tI = 0.0
    tFmin = 50.0 # final time to simulate to (min)
    dtmin = 2e-1 # time increment to give solutions at (min)

    # Convert times to seconds
    tF = tFmin*60
    dt = dtmin*60

    # Initialize HVACIVP object
    p = get_HVACIVP_pdict()

    # Set some plotting variables
    norm = colors.Normalize(vmin=10.0, vmax=30.0, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.jet)
    fig_eq, ax_eq = plt.subplots()
    ax_eq.set_xlim(-5, 40)
    ax_eq.set_ylim(-5, 15)
    ax_eq.set_aspect(1.0)
    fig_eq.colorbar(mapper, ax=ax_eq)
    plotdict = {}
    plotdict['fig'] = fig_eq
    plotdict['ax'] = ax_eq
    plotdict['norm'] = norm
    plotdict['mapper'] = mapper

    # Construct IVP to solve equilibrium problem
    slimmonsIVP_eq = HVACIVP('slimmons.csv', uI, tI, tF, p)

    # Solve for equilibrium condition using default linear solver (np.linalg.solve)
    print("Using default linear solver to find equilibrium condition")
    u_eq_def = slimmonsIVP_eq.solve_eq(plotdict=plotdict)

    # Solve for equilibrium condition using Gaussian elimination (in myGE)
    print()
    print("Using myGE to find equilibrium condition")
    u_eq_myGE = slimmonsIVP_eq.solve_eq(linsolve=myGE)
    print(f"Does solution = default solver solution: {np.allclose(u_eq_def, u_eq_myGE)}")

    # Solve for equilibrium condition using Gaussian elimination with vectorization (in myGE_vec)
    print()
    print("Using myGE_vec to find equilibrium condition")
    u_eq_myGE_vec = slimmonsIVP_eq.solve_eq(linsolve=myGE_vec)
    print(f"Does solution = default solver solution: {np.allclose(u_eq_def, u_eq_myGE_vec)}")

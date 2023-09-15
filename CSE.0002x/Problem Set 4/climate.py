################################################################################
# CSE.0002x
# Problem Set 4: climate
# edX Username: adevrent

import numpy as np
import matplotlib.pyplot as plt
from IVPlib_pset4 import IVP
import IVPlib_pset4 as IVPlib


################################################################################
## Climate class definition (a subclass of IVP)
################################################################################


class ClimateIVP(IVP):

    def evalf(self, u, t):
        """
        Calculates right-hand side in the climate model for atmospheric
        temperature, and concentrations of H2O and CO2, given the current
        state and time

        Args:
            u (NumPy ndarray): current state with:
                u[0] = atmospheric temperature (K)
                u[1] = H2O concentration (ppm)
                u[2] = CO2 concentration (ppm)
            t (float): current time

        Returns:
            f (NumPy ndarray): current right-hand side
        """
        day_per_yr = 365.
        sec_per_yr = 60.*60.*24.*day_per_yr

        heatcapC = self.get_p('heatcapC')
        S0       = self.get_p('S0')
        sigma    = self.get_p('sigma')
        tauH     = self.get_p('tauH')/day_per_yr
        tauC     = self.get_p('tauC')
        Href     = self.get_p('Href')
        Cref     = self.get_p('Cref')
        alpha    = self.get_p('alpha')

        T, H, C = u

        eps = self.calceps(u)

        incfT = (1-alpha)*S0/4
        outfT = (1-0.5*eps)*sigma*T**4

        PHC = self.get_p('PHC')(self, t)

        fT = (1/heatcapC) * (incfT - outfT) * sec_per_yr
        fH = (1/tauH) * (Href - H) + PHC[0]
        fC = (1/tauC) * (Cref - C) + PHC[1]
        f = np.array([fT, fH, fC])

        return f

    def evalf_u(self, u, t):
        """
        Calculates gradient of f with respect to u

        Args:
            u (NumPy ndarray): current state with:
                u[0] = atmospheric temperature (K)
                u[1] = H2O concentration (ppm)
                u[2] = CO2 concentration (ppm)
            t (float): current time

        Returns:
            f_u (NumPy ndarray)
        """
        day_per_yr = 365.
        sec_per_yr = 60.*60.*24.*day_per_yr

        heatcapC = self.get_p('heatcapC')
        sigma    = self.get_p('sigma')
        tauH     = self.get_p('tauH')/day_per_yr
        tauC     = self.get_p('tauC')

        T, H, C = u

        eps = self.calceps(u)
        eps_u = self.calceps_u(u)

        outfT_eps = -0.5*sigma*T**4
        outfT_T = 4*(1-0.5*eps)*sigma*T**3
        outfT_u = outfT_eps*eps_u
        outfT_u[0] += outfT_T

        f_u = np.zeros((3,3))
        f_u[0,:] = (1/heatcapC) * (-outfT_u) * sec_per_yr
        f_u[1,1] = -(1/tauH)
        f_u[2,2] = -(1/tauC)

        return f_u

    def calceps(self, u):
        """
        Calculate Earth's emissivity

        This assumes T is in Kelvin and [H20] and [CO2] are ppm.

        Args:
            u (NumPy ndarray): current state

        Returns:
            eps: emissivity
        """
        T, H, C = u

        Cref = self.get_p('Cref')

        Patm = 1013.25 # hPa
        eps = 1.24*(H*Patm*(1e-6)/T)**(1/7) + 0.06/(np.log(2))*np.log(C/Cref)

        return eps

    def calceps_u(self, u):
        """
        Calculate gradient of Earth's emissivity with respect the state

        This assumes T is in Kelvin and [H20] and [CO2] are ppm.

        Args:
            u (NumPy ndarray): current state

        Returns:
            eps_u: deps/du
        """
        T, H, C = u

        Patm = 1013.25  # hPa
        eps_T = -(1.24/7)*((H*Patm*(1e-6)  )**(1/7))*T**(-8/7)
        eps_H =  (1.24/7)*((  Patm*(1e-6)/T)**(1/7))*H**(-6/7)
        eps_C = 0.06/(np.log(2))/C

        return np.array([eps_T, eps_H, eps_C])


################################################################################
## Helper functions for plotting and simulation (DO NOT MODIFY)
################################################################################


def plot_THC(t, u, title="Climate evolution", figaxs=None):
    """
    Plot T, [H2O], and [CO2] versus time using a 3 by 1 plot
    with shared x axes.

    Args:
        t (NumPy ndarray): time values
        u (NumPy ndarray): states where
                           u[n,0] is T at t=t[n]
                           u[n,1] is [H2O] at t=t[n]
                           u[n,2] is [CO2] at t=t[n]
        title (string): title to use on plot
        figaxs (tuple): figaxs[0] is a figure handle to plot to
            and figaxs[1] is an array of axs handles to plot to

    Returns:
        fig, axs of the plots
    """
    if figaxs is None:
        fig, axs = plt.subplots(3, 1, sharex=True)
    else:
        fig, axs = figaxs

    axs_T = axs[0]
    axs_H = axs[1]
    axs_C = axs[2]
    axs_T.plot(t, u[:, 0] - u[0,0], 'g-')
    axs_T.set_title(title)
    axs_T.set_ylabel('$\Delta T$ (K)')
    axs_H.plot(t, u[:, 1], 'r-')
    axs_H.set_ylabel('H2O (ppm)')
    axs_C.plot(t, u[:, 2], 'b-')
    axs_C.set_ylabel('CO2 (ppm)')
    axs_C.set_xlabel('t (yr)')

    return fig, axs


def plot_PHC(t, PHC, title="Climate evolution", figaxs=None):
    """
    Plot production rates of [H2O] and [CO2] versus time using a 2 by 1 plot
    with shared x axes.

    Args:
        t (NumPy ndarray): time history
        PHC (NumPy ndarray): PHC[n,0] = production rate of H2O at t=t[n]
                             PHC[n,1] = production rate of CO2 at t=t[n]
        title (string): title to use on plot
        figaxs (tuple): figaxs[0] is a figure handle to plot to
            and figaxs[1] is an array of axs handles to plot to

    Returns:
        fig, axs of the plots
    """
    if figaxs is None:
        fig, axs = plt.subplots(2, 1, sharex=True)
    else:
        fig, axs = figaxs

    axs_PH2O = axs[0]
    axs_PCO2 = axs[1]
    axs_PH2O.plot(t, PHC[:, 0], 'r-')
    axs_PH2O.set_title(title)
    axs_PH2O.set_ylabel('Production H2O (ppm/yr)')
    axs_PCO2.plot(t, PHC[:, 1], 'b-')
    axs_PCO2.set_ylabel('Production CO2 (ppm/yr)')
    axs_PCO2.set_xlabel('t (yr)')

    return fig, axs


def run_climate_case(C, tauH, tauC, PHC_scenario):
    """
    Run a single climate case for the PHC_scenario
    with the values of C, tauH, and tauC provided.
    All other parameters are nominal.

    Args:
        C (float): effective specific heat capacity (J/(m^2 K))
        tauH (float): H2O adjustment timescale (days)
        tauC (float): CO2 adjustment timescale (yrs)
        PHC_scenario: reference to a function that gives production rates
            of H2O and CO2 at time t (t is in years)

    Returns:
        t (NumPy ndarray of floats): time history
        u (NumPy ndarray): states where
                           u[n,0] is T at t=t[n]
                           u[n,1] is [H2O] at t=t[n]
                           u[n,2] is [CO2] at t=t[n]
        PHC (NumPy ndarray): PHC[n,0] = production rate of H2O at t=t[n]
                             PHC[n,1] = production rate of CO2 at t=t[n]
    """
    # Set-up climate IVP dictionary
    p = {}
    p['heatcapC'] = C # J/(m^2*K)
    p['S0'] = 1370 # W/m^2
    p['sigma'] = 5.67e-8 # W/(m^2*K^4)
    p['tauH'] = tauH # days
    p['tauC'] = tauC # yr
    p['Href'] = 0.85e4 # ppm
    p['Cref'] = 280 # ppm
    p['alpha'] = 0.3

    # Production rate info
    p['PH0'] = 10.0 # initial production of H2O (ppm/yr)
    p['PC0'] = 8.0  # initial production of CO2 (ppm/yr)
    p['PHC'] = PHC_scenario # function to calculate H2O and CO2 production

    # Initial condition
    TI = 288.7 # K
    HI = 8500.25 # ppm
    CI = 420 # ppm
    uI = np.array([TI, HI, CI])

    # Time range to consider
    tI = 0
    tF = 2e2 # yr

    # Instantiate ClimateIVP object
    climate_IVP = ClimateIVP(uI, tI, tF, p)

    # Solve climate IVP
    dt = 1.0e-0
    method_dict = {}
    method_dict['use_evalf_u'] = True
    t, u = IVPlib.solve(climate_IVP, dt, IVPlib.step_BE, method_dict)

    nt = len(t)
    PHC = np.zeros((nt, 2))
    for i in range(nt):
        PHC[i, :] = PHC_scenario(climate_IVP, t[i])

    return t, u, PHC


######################################################################
# Time-dependent scenarios of H2O and CO2 production
######################################################################


def PHC_growth(cIVP, t, tgrowth=60):
    """
    Calculate the human-caused production of H2O and CO2
    in ppm/yr as a function of time for the growth scenario.
    In this scenario, the production rates are initially 'PH0'
    and 'PC0' (stored in the cIVP parameter dictionary) and then:

    from 0 < t < tgrowth:
        the production rates increase linearly from their initial values
        such that at t=tgrowth, the production rates are 3 times their
        initial rates
    for t > tgrowth:
        the production rates remain at 3 times their initial rates

    Args:
        cIVP (ClimateIVP object): contains the initial production rates
            in its parameter dictionary
        t (float): current time in years

    Returns:
        NumPy ndarray containing production rate of [H20, CO2]
    """
    #### BEGIN SOLUTION #####
    
    PH0      = cIVP.get_p('PH0')
    PC0      = cIVP.get_p('PC0')
    
    if t <= tgrowth:
        H20_growth = t * (2*PH0)/tgrowth
        H20 = PH0 + H20_growth
        
        CO2_growth = t * (2*PC0)/tgrowth
        CO2 = PC0 + CO2_growth
    else:
        H20 = 3 * PH0
        CO2 = 3 * PC0
    
    prod_rate = np.array((H20, CO2))  # array of the production rate at t
    
    return prod_rate
    
    
    #### END SOLUTION #####


def PHC_decline(cIVP, t, tdecline=10):
    """
    Calculate the human-caused production of H2O and CO2
    in ppm/yr as a function of time for the decline scenario.
    In this scenario, the production rates are initially 'PH0'
    and 'PC0' (stored in the cIVP parameter dictionary) and then:

    from 0 < t < tdecline:
        the production rates decrease linearly from their initial values
        such that at t=tdecline, the production rates are half their
        initial rates
    for t > tdecline:
        the production rates remain at half their initial rates

    Args:
        cIVP (ClimateIVP object): contains the initial production rates
            in its parameter dictionary
        t (float): current time in years

    Returns:
        NumPy ndarray containing production rate of [H20, CO2]
    """
    #### BEGIN SOLUTION #####
    PH0      = cIVP.get_p('PH0')
    PC0      = cIVP.get_p('PC0')
    
    if t <= tdecline:
        H20_decline = t * (0.5*PH0)/tdecline
        H20 = PH0 - H20_decline
        
        CO2_decline = t * (0.5*PC0)/tdecline
        CO2 = PC0 - CO2_decline
    else:
        H20 = 0.5 * PH0
        CO2 = 0.5 * PC0
    
    prod_rate = np.array((H20, CO2))  # array of the production rate at t
    
    return prod_rate
    #### END SOLUTION #####


######################################################################
# Deterministic climate simulation with nominal parameters
######################################################################


def run_climate_nominal(PHC_scenario, scenario_title):
    """
    Run climate model with nominal values of parameters

    Args:
        PHC_scenario: reference to a function that gives production rates
            of H2O and CO2 at time t (t is in years)
        scenario_title (string): title that will be used to label plots, etc.

    Returns:
        dTmax (float):
    """
    print(f"Scenario = {scenario_title}")
    print("Nominal case:")

    C = 4.5e8
    tauH = 9.0
    tauC = 60.
    t, u, PHC = run_climate_case(C, tauH, tauC, PHC_scenario)

    plot_PHC(t, PHC, title=scenario_title)
    plot_THC(t, u, title=scenario_title)

    #### BEGIN SOLUTION #####
    dTmax = u[:, 0].max() - u[0, 0]
    print(f"dTmax = {dTmax.round(2)} degrees K")
    
    return dTmax
    #### END SOLUTION #####


######################################################################
# Monte Carlo climate simulation with sampled parameters
######################################################################


def sample_climate_results(PHC_scenario, scenario_title, Nsample):
    """
    Run a Monte Carlo sample of the climate model, and plot the results of
    max temperature rise in a histogram.

    Args:
        Same args as those for run_climate_MC

    Returns:
        C (NumPy array of floats): sampled values of the effective heat capacity
        tauH (NumPy array of floats): sampled values of H2O's adjustment time
        tauC (NumPy array of floats): sampled values of CO2's adjustment time
        dTmax (NumPy array of floats): outcomes of the maximum temperature rise
            according to the simulation with corresponding parameters
        fighist: the Figure object containing the histogram
        axshist: the Axes object in which the histogram is drawn
        histbins (tuple): the ENTIRE tuple output of the call to axshist.hist()

        The values across the four arrays at the same index correspond to a single outcome
        within the entire sample.
    """
    #### BEGIN SOLUTION #####
    rng = np.random.default_rng()
    C = rng.triangular(1.5 * 10**8, 4.5 * 10**8, 7.5 * 10**8, Nsample)  # J/(K m^2)
    tauH = rng.triangular(7.0, 9.0, 11.0, Nsample)  # days
    tauC = rng.triangular(30.0, 60.0, 90.0, Nsample)  # years
    
    dTmax = np.zeros(Nsample)
    for i in range(Nsample):  # loop over sample size, NOT time
        C_temp = C[i]
        tauH_temp = tauH[i]
        tauC_temp = tauC[i]

        t, u, PHC = run_climate_case(C_temp, tauH_temp, tauC_temp, PHC_scenario)
        # PH_temp = PHC[:, 0]
        # PC_temp = PHC[:, 1]
        dTmax[i] = u[:, 0].max() - u[0, 0]
    
    fighist, axhist = plt.subplots()
    histbins = axhist.hist(dTmax, density=True)
    axhist.set_title(f"{scenario_title}, Nsample = {Nsample}")
    axhist.set_ylabel("PDF")
    axhist.set_xlabel("$\Delta T$ max (K)")
    
    return C, tauH, tauC, dTmax, fighist, axhist, histbins
    #### END SOLUTION #####


def calc_percentiles(dTmax):
    """
    Args:
        dTmax (NumPy array of floats): a sample of maximum temperature rise outcomes
            from a Monte Carlo climate simulation

    Returns:
        * the 5th-percentile value of the dTmax data
        * the 50th-percentile value of the dTmax data
        * the 95th-percentile value of the dTmax data
    """
    #### BEGIN SOLUTION #####
    perc_5th, perc_50th, perc_95th = np.percentile(dTmax, [5, 50, 95])
    
    print(f"dTmax   5%: {perc_5th.round(2)} K")
    print(f"dTmax  50%: {perc_50th.round(2)} K")
    print(f"dTmax  95%: {perc_95th.round(2)} K")
    
    return perc_5th, perc_50th, perc_95th
    #### END SOLUTION #####


def calc_mean_temp_rise_CI(dTmax):
    """
    Args:
        dTmax (NumPy array of floats): a sample of maximum temperature rise outcomes
            from a Monte Carlo climate simulation

    Returns:
        * the mean maximum temperature rise
        * the lower end of the 95% confidence interval on the mean
        * the upper end of the 95% confidence interval on the mean
    """
    #### BEGIN SOLUTION #####
    sample_mean = dTmax.mean()
    # Let Theta a random variable which is an estimation of the true sample mean.
    Theta_stdev = dTmax.std() / np.sqrt(len(dTmax))
    
    lower_end_95 = sample_mean - (1.96 * Theta_stdev)
    higher_end_95 = sample_mean + (1.96 * Theta_stdev)
    
    print(f"dTmax mean: {sample_mean.round(2)} with 95% C.I. = [{lower_end_95.round(2)}, {higher_end_95.round(2)}]")
    
    return sample_mean, lower_end_95, higher_end_95
    
    #### END SOLUTION #####


def calc_threshold_probability_CI(dTmax):
    """
    Args:
        dTmax (NumPy array of floats): a sample of maximum temperature rise outcomes
            from a Monte Carlo climate simulation

    Returns:
        * the estimated probability that the max temperature rise falls below 0.5 K
        * the lower end of the 95% confidence interval on the probability
        * the upper end of the 95% confidence interval on the probability
    """
    #### BEGIN SOLUTION #####
    lt_05 = np.count_nonzero(np.logical_and(dTmax < 0.5, True))
    sample_probability = lt_05 / len(dTmax)
    Theta_stdev = np.sqrt(sample_probability * (1 - sample_probability)) / np.sqrt(len(dTmax))
    
    lower_95 = sample_probability - (1.96 * Theta_stdev)
    upper_95 = sample_probability + (1.96 * Theta_stdev)
    
    print(f"P(dTmax < 0.5 K) = {sample_probability:.3f} with 95% C.I. = [{lower_95:.3f}, {upper_95:.3f}]")
    
    return sample_probability, lower_95, upper_95
    #### END SOLUTION #####


def plot_threshold_variability(C, tauH, tauC, dTmax, scenario_title):
    """
    Plot the sampled set of parameter values across 2D subplots that are
    pairwise combinations of the three parameters. Highlight in blue
    those outcomes whose simulated dTmax was below 0.5 K, and use red markers
    for all other data points.

    Args:
        C (NumPy array of floats): as returned by sample_climate_results
        tauH (NumPy array of floats): as returned by sample_climate_results
        tauC (NumPy array of floats): as returned by sample_climate_results
        dTmax (NumPy array of floats): as returned by sample_climate_results
        scenario_title (string): as given to run_climate_MC

    Returns:
        figscat: the Figure object containing the scatter plots
        axsscat: the array of Axes objects in which the scatter plots are drawn
    """
    #### BEGIN SOLUTION #####
    lt_05_indices = np.where(dTmax < 0.5)
    gt_05_indices = np.where(dTmax >= 0.5)
    
    # print("lt_05_indices:", lt_05_indices)
    # print("gt_05_indices:", gt_05_indices)
    
    figscat, axsscat = plt.subplots(3, 1, figsize=(10, 6))
    
    axsscat[0].scatter(C[lt_05_indices], tauH[lt_05_indices], c="blue", marker="*")
    axsscat[0].scatter(C[gt_05_indices], tauH[gt_05_indices], c="red", marker="*")
    axsscat[0].set_xlabel("C (J/(m^2 K))")
    axsscat[0].set_ylabel("tauH (days)")
    axsscat[0].set_title(scenario_title)
    
    axsscat[1].scatter(C[lt_05_indices], tauC[lt_05_indices], c="blue", marker="*")
    axsscat[1].scatter(C[gt_05_indices], tauC[gt_05_indices], c="red", marker="*")
    axsscat[1].set_xlabel("C (J/(m^2 K))")
    axsscat[1].set_ylabel("tauC (yrs)")
    
    axsscat[2].scatter(tauC[lt_05_indices], tauH[lt_05_indices], c="blue", marker="*")
    axsscat[2].scatter(tauC[gt_05_indices], tauH[gt_05_indices], c="red", marker="*")
    axsscat[2].set_xlabel("tauC (yrs)")
    axsscat[2].set_ylabel("tauH (days)")
    
    return figscat, axsscat
    #### END SOLUTION #####


def run_climate_MC(PHC_scenario, scenario_title, Nsample):
    """
    Run a Monte Carlo sample of the climate model, and perform percentile,
    confidence interval, and parameter variability analysis on results.

    Args:
        PHC_scenario: reference to a function that gives production rates
            of H2O and CO2 at time t (t is in years)
        scenario_title (string): title that will be used to label plots, etc.
        Nsample (int): sample size for Monte Carlo simulation
    """
    print(f"Scenario = {scenario_title}")
    print(f"Running Monte Carlo with sample size = {Nsample}")

    C, tauH, tauC, dTmax, _, _, _ = sample_climate_results(PHC_scenario, scenario_title, Nsample)
    calc_percentiles(dTmax)
    calc_mean_temp_rise_CI(dTmax)
    calc_threshold_probability_CI(dTmax)
    plot_threshold_variability(C, tauH, tauC, dTmax, scenario_title)


######################################################################
# Main block
######################################################################


if __name__ == '__main__':
    try:
        run_climate_nominal(PHC_growth, 'Growth scenario')
    except NotImplementedError:
        pass
    print()

    try:
        run_climate_nominal(PHC_decline, 'Decline scenario')
    except NotImplementedError:
        pass
    print()

    try:
        run_climate_MC(PHC_growth, 'Growth scenario', 100)
    except NotImplementedError:
        pass
    print()

    try:
        run_climate_MC(PHC_decline, 'Decline scenario', 100)
    except NotImplementedError:
        pass
    print()

    plt.show()

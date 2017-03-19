"""
Ruth's version of Burke's test_comp_grid.py as a function.
"""

import numpy as np
import matplotlib.pyplot as plt
import KeplerPORTs_utils as kpu


def calc_comp(kepid, period, radius):
    """
    Calculate the completeness at a given radius and period for a KIC star.
    This includes the probability of transiting.
    parameters:
    ----------
    kepid: (int)
        The KIC id.
    period: (float)
        The target period.
    radius: (float)
        The target radius.
    returns:
    --------
    The Completeness.
    FIXME: Interpolate instead of finding nearest.
    """

    # Instantiate pipeline completeness class structure
    doit = kpu.kepler_single_comp_data()
    doit.id = kepid
    doit.period_want = np.linspace(10.0, 700.0, 6000)
    doit.rp_want = np.linspace(0.5, 2.5, 3000)
    doit.rstar = 0.98
    doit.logg = 4.44
    doit.deteffver = 2
    doit.ecc = 0.0
    doit.dataspan = 1426.7
    doit.dutycycle = 0.879
    doit.pulsedurations = [1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 6.0, 7.5, 9.0,
                           10.5, 12.0, 12.5, 15.0]
    doit.cdpps = [36.2, 33.2, 31.0, 29.4, 28.0, 26.1, 25.4, 24.2, 23.1, 22.4,
                  21.9, 21.8, 21.7, 21.5]
    doit.mesthresh = np.full_like(doit.pulsedurations,7.1)

    # Calculate completeness over the grid of periods and radii.
    x1, x2 = kpu.kepler_single_comp(doit)
    X, Y = np.meshgrid(doit.period_want, doit.rp_want)

    # Select completeness for a given period and radius
    target_p = find_nearest(doit.period_want, period)
    target_r = find_nearest(doit.rp_want, radius)
    mp = doit.period_want == target_p  # period size = (6000,)
    mr = doit.rp_want == target_r  # radius size = (3000,)
    ind_p = np.arange(len(doit.period_want))[mp]
    ind_r = np.arange(len(doit.rp_want))[mr]
    probdet, probtot = x1, x2  # shapes = (3000, 6000), (3000, 6000)
    return probtot[ind_r, ind_p]


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx]


if __name__ == "__main__":
    print(calc_comp(10593626, 365.25, 1))

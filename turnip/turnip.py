# Create a plot of N planets vs age.

import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as sst

import gyro as g

plotpar = {'axes.labelsize': 18,
           'font.size': 10,
           'legend.fontsize': 15,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)


def merge_data(DATA_DIR="/Users/ruthangus/projects/turnip/turnip/data"):
    """
    Merge the McQuillan periods with the list of KOIs.
    """

    # Load planet data
    planets = pd.read_csv(os.path.join(DATA_DIR, "cumulative.csv"),
                          comment="#")
    planets = planets[planets.koi_pdisposition == "CANDIDATE"]
    print(planets.keys())

    # load McQuillan data
    koi_p = pd.read_csv(os.path.join(DATA_DIR, "Table_1.txt"), comment="#")
    kep_p = pd.read_csv(os.path.join(DATA_DIR, "Table_1_Periodic.txt"),
                        comment="#")
    kep_non_p = pd.read_csv(os.path.join(DATA_DIR,
                                         "Table_2_Non_Periodic.txt"),
                            comment="#")

    # Merge the period data
    periods = pd.concat([koi_p, kep_p], join="outer")
    periods.to_csv(os.path.join(DATA_DIR, "mcquillan_periods.csv"))

    # Merge the planet data with the rotation period data.
    df = pd.merge(planets, periods, on="kepid")
    df.to_csv(os.path.join(DATA_DIR, "koi_periods.csv"))


def planet_frequency(periods, planets):
    """
    Make distribution plots
    """
    # Don't count planet hosts more than once.
    unique, m = np.unique(planets.kepid.values, return_index=True)
    planets = planets.iloc[m]

    # Calculate ages
    all_gyro = g.gyro_age(periods.P_rot.values, periods.teff.values,
                          np.zeros(len(periods.P_rot.values)),
                          periods.logg.values)
    planet_gyro = g.gyro_age(planets.P_rot.values, planets.teff.values,
                             np.zeros(len(planets.P_rot.values)),
                             planets.logg.values)
    all_ages = all_gyro.barnes07(version="mh")
    planet_ages = planet_gyro.barnes07(version="mh")

    m0 = all_ages < 13.7
    m1 = planet_ages < 13.7
    make_hist(all_ages[m0], planet_ages[m1], "Age")
    make_hist(periods.P_rot.values, planets.Prot.values, "Period")
    make_hist(periods.teff.values, planets.teff.values, "Teff")
    make_hist(periods.R_var.values, planets.R_var.values, "Activity")


def make_hist(all_x, planet_x, name):
    """
    Make histogram plots.
    """

    # Remove NaNs
    m0, m1 = np.isfinite(all_x), np.isfinite(planet_x)

    print(name, sst.ks_2samp(all_x[m0], planet_x[m1]))

    # Make histograms of age, period, teff, etc
    all_x_hist, bins0 = np.histogram(all_x[m0])
    planet_x_hist, bins1 = np.histogram(planet_x[m1])
    all_x_hist_n, bins0 = np.histogram(all_x[m0], normed=True)
    planet_x_hist_n, bins1 = np.histogram(planet_x[m1], normed=True)

    # Re-format the shapes of the histograms for plotting (prepend with 0)
    all_x_hist_n = extend(all_x_hist_n)
    planet_x_hist_n = extend(planet_x_hist_n)
    all_x_hist, planet_x_hist = extend(all_x_hist), extend(planet_x_hist)

    # Poissonian uncertainties.
    all_yerr = 1./(all_x_hist**.5) * (all_x_hist_n/all_x_hist)
    planet_yerr = 1./(planet_x_hist**.5) * (planet_x_hist_n/planet_x_hist)

    plt.clf()
    plt.step(bins0, all_x_hist_n, color="k", label="All stars")
    plt.errorbar(bins0-(bins0[1]-bins0[0])*.5, all_x_hist_n,
                 yerr=all_yerr, fmt="k.", capsize=0, ms=.1)
    plt.step(bins1, planet_x_hist_n, color="r", label="Planet hosts")
    plt.errorbar(bins1-(bins1[1]-bins1[0])*.5, planet_x_hist_n,
                 yerr=planet_yerr, fmt="r.", capsize=0, ms=.1)
    plt.axvline(np.percentile(all_x[m0], 50), color="k", ls="--")
    plt.axvline(np.percentile(planet_x[m1], 50), color="r", ls="--")
    plt.xlabel("{}".format(name))
    plt.ylabel("$N$")
    plt.legend()
    xmax = bins0[-1]
    plt.xlim(0, xmax)
    plt.savefig("{}_hist".format(name))


def extend(x):
    x = list(x)
    x.insert(0, 0)
    return np.array(x)


if __name__ == "__main__":

    DATA_DIR = "/Users/ruthangus/projects/turnip/turnip/data"
    periods = pd.read_csv(os.path.join(DATA_DIR, "mcquillan_periods.csv"))
    planets = pd.read_csv(os.path.join(DATA_DIR, "koi_periods.csv"))
    planet_frequency(periods, planets)

# Create a plot of N planets vs age.

import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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

    all_gyro = g.gyro_age(periods.P_rot.values, periods.teff.values,
                          np.zeros(len(periods.P_rot.values)),
                          periods.logg.values)
    planet_gyro = g.gyro_age(planets.P_rot.values, planets.teff.values,
                             np.zeros(len(planets.P_rot.values)),
                             planets.logg.values)
    all_ages = all_gyro.barnes07(version="barnes")
    planet_ages = planet_gyro.barnes07(version="barnes")

    m0 = np.isfinite(all_ages) * (all_ages < 13.7)
    m1 = np.isfinite(planet_ages) * (planet_ages < 13.7)
    m2 = np.isfinite(periods.P_rot.values)
    m3 = np.isfinite(planets.Prot.values)
    m4 = np.isfinite(periods.teff.values)
    m5 = np.isfinite(planets.teff.values)

    plt.clf()
    plt.hist(all_ages[m0], normed=True, histtype="stepfilled", color="b",
             alpha=.5)
    plt.hist(planet_ages[m1], normed=True, histtype="stepfilled", color="r",
             alpha=.5)
    plt.xlabel("$\mathrm{Age~(Gyr)}$")
    plt.ylabel("$N$")
    plt.savefig("age_hist")

    plt.clf()
    plt.hist(periods.P_rot.values[m2], normed=True, histtype="stepfilled",
             color="b", alpha=.5)
    plt.hist(planets.Prot.values[m3], normed=True, histtype="stepfilled",
             color="r", alpha=.5)
    plt.xlabel("$P_{\mathrm{rot}}~\mathrm{(Days)}$")
    plt.ylabel("$N$")
    plt.savefig("period_hist")

    plt.clf()
    plt.hist(periods.teff.values[m4], normed=True, histtype="stepfilled",
             color="b", alpha=.5)
    plt.hist(planets.teff.values[m5], normed=True, histtype="stepfilled",
             color="r", alpha=.5)
    plt.xlabel("$T_{\mathrm{eff}}$")
    plt.ylabel("$N$")
    plt.savefig("teff_hist")



if __name__ == "__main__":

    DATA_DIR = "/Users/ruthangus/projects/turnip/turnip/data"
    periods = pd.read_csv(os.path.join(DATA_DIR, "mcquillan_periods.csv"))
    planets = pd.read_csv(os.path.join(DATA_DIR, "koi_periods.csv"))
    planet_frequency(periods, planets)

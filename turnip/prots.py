# plot rotation period vs orbital period

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import re
from gyro import gyro_age
import teff_bv as tbv
import scipy.stats as sps
from calc_completeness import calc_comp
# np.set_printoptions(threshold=np.nan, linewidth=9)

plotpar = {'axes.labelsize': 18,
           'text.fontsize': 10,
           'legend.fontsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)

DATA_DIR = "/Users/ruthangus/projects/turnip/turnip/data/"


def save_data(nbins):

    fnames = glob.glob(os.path.join(DATA_DIR, "/koi_results/*h5"))

    koi, period, errp, errm, lnerrp, lnerrm = [], [], [], [], [], []
    for i, fn in enumerate(fnames):
        df = pd.read_hdf(fn, key="samples")
        phist, bins = np.histogram(df.ln_period.values, nbins)
        ln_p = bins[phist == max(phist)][0]
        period.append(np.exp(ln_p))
        lnerrp.append(np.percentile(df.ln_period.values, 84) - ln_p)
        lnerrm.append(ln_p - np.percentile(df.ln_period.values, 16))
        errp.append(np.exp(lnerrp[i]/ln_p))
        errm.append(np.exp(lnerrm[i]/ln_p))
        koi.append(re.findall('\d+', fn)[0])
    table = pd.DataFrame({"koi": np.array(koi), "period": np.array(period),
                          "errp": np.array(errp), "errm": np.array(errm)})
    table.to_csv("period_point_estimates.csv")


def make_histogram():
    df = pd.read_csv("period_point_estimates.csv")
    plt.clf()
    plt.hist(df.period, 20)
    plt.savefig("gp_period_hist")


def make_df():
    df = pd.read_csv("period_point_estimates.csv")
    planets = pd.read_csv(os.path.join(DATA_DIR, "cumulative.csv"),
                          skiprows=155)

    kois = []
    for i, k in enumerate(planets.kepoi_name.values):
        # print(planets.kepoi_name.values[i])
        # print(type(planets.kepoi_name.values[i]))
        koi_str = re.findall('\d+', planets.kepoi_name.values[i])[0]
        kois.append(int(koi_str))
    planets["koi"] = kois
    joint = pd.merge(planets, df, on="koi")
    joint.to_csv("planet_periods.csv")


def plot_periods():
    df = pd.read_csv("planet_periods.csv")

    m = np.log(df.period.values) > 1
    lnporb = np.log(df.koi_period.values[m])
    lnprot = np.log(df.period.values[m])
    porb = df.koi_period.values[m]
    prot = df.period.values[m]
    radius = np.log(df.koi_prad.values[m])
    teff = df.koi_steff.values[m]

    plt.clf()
    plt.scatter(porb, prot, s=5*radius, c=teff, vmin=4400, vmax=7000)
    plt.loglog()
    plt.colorbar()
    plt.xlabel("$\ln(\mathrm{Orbital~period})$")
    plt.ylabel("$\ln(\mathrm{Rotation~period})$")
    plt.subplots_adjust(bottom=.15)
    plt.savefig("period_period")

    # find the short rotators
    m = np.log(df.period.values) < 1
    print(df.koi.values[m])
    # import kplr
    # client = kplr.API()
    # for i, k in enumerate(df.koi.values[m]):
    #     print(k)
    #     star = client.koi("{}.01".format(k))
    #     star.get_light_curves(fetch=True)


def plot_radii():
    df = pd.read_csv("planet_periods.csv")

    m = np.log(df.period.values) > 1
    prot = df.period.values[m]
    radius = np.log(df.koi_prad.values[m])
    teff = df.koi_steff.values[m]
    logg = df.koi_slogg.values[m]
    feh = np.zeros(len(logg))
    gyro = gyro_age(prot, teff, feh, logg)
    age = gyro.barnes07("mh")

    m = np.isfinite(age)

    plt.clf()
    plt.scatter(np.log(age[m]), np.log(radius[m]), c=teff[m], s=10, vmin=4400,
                vmax=7000)
    plt.colorbar()
    plt.xlabel("$\ln(\mathrm{Age,~Gyr})$")
    plt.ylabel("$\ln(\mathrm{Radius}, R_J)$")
    plt.subplots_adjust(bottom=.15)
    plt.savefig("age_radius")

    l = age[m] < 3.295
    print(len(radius[m][l]))
    print(len(radius[m][~l]))
    plt.clf()
    plt.hist(radius[m][l], 50, normed=True, alpha=.5, label="young")
    plt.hist(radius[m][~l], 40, normed=True, alpha=.5, label="old")
    plt.legend()
    plt.xlabel("Radius")
    plt.savefig("radius_hist")
    print(sps.ks_2samp(radius[m][l], radius[m][~l]))

    cum_young = np.cumsum(radius[m][l]) / sum(radius[m][l])
    cum_old = np.cumsum(radius[m][~l]) / sum(radius[m][~l])
    plt.clf()
    plt.plot(cum_young, label="young")
    plt.plot(cum_old, label="old")
    plt.savefig("radius_cdf")

    # # print(np.unique(df.kepid.values[m]))
    # for i in np.unique(df.kepid.values[m]):
    #     print("KIC", str(int(i)).zfill(9))

    n = radius[m][l] < .5
    n2 = radius[m][~l] < .5
    print(len(radius[m][l][n]))
    print(len(radius[m][~l][n2]))
    plt.clf()
    plt.hist(radius[m][l][n], 50, normed=True, alpha=.5, label="young")
    plt.hist(radius[m][~l][n2], 40, normed=True, alpha=.5, label="old")
    plt.legend()
    plt.xlabel("Radius")
    plt.savefig("radius_hist_hj")
    print(sps.ks_2samp(radius[m][l][n], radius[m][~l][n2]))

    n = radius[m] < .5
    plt.clf()
    plt.scatter(np.log(age[m][n]), np.log(radius[m][n]), c=teff[m][n], s=10,
                vmin=4400, vmax=7000)
    plt.colorbar()
    plt.xlabel("$\ln(\mathrm{Age,~Gyr})$")
    plt.ylabel("$\ln(\mathrm{Radius}, R_J)$")
    plt.subplots_adjust(bottom=.15)
    plt.savefig("age_radius_hj")


def plot_completeness():
    df = pd.read_csv("planet_periods.csv")
    comp = np.zeros((len(df.kepid.values)))
    print(df.kepid.values[:10])
    for i, kepid in enumerate(df.kepid.values[:10]):
        print(i, "of", len(df.kepid.values))
        print("id = ", kepid)
        comp[i] = calc_comp(kepid, 365.25, 1.)
        print(comp[i])
    df["probtot"] = comp

    plt.clf()
    plt.plot(comp[:10], df.period.values[:10], "k.")
    plt.savefig("comp_vs_period")


if __name__ == "__main__":
    # save_data(100)
    # make_histogram()
    # make_df()
    # plot_periods()
    # plot_radii()
    plot_completeness()

# plot rotation period vs orbital period

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import re

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

    find the short rotators
    m = np.log(df.period.values) < 1
    print(df.koi.values[m])
    # import kplr
    # client = kplr.API()
    # for i, k in enumerate(df.koi.values[m]):
    #     print(k)
    #     star = client.koi("{}.01".format(k))
    #     star.get_light_curves(fetch=True)


if __name__ == "__main__":
    # save_data(100)
    # make_histogram()
    make_df()
    plot_periods()

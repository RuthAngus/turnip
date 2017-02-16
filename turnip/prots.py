# plot rotation period vs orbital period

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import re

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
        kois.append(re.findall('\d+', planets.kepoi_name.values[i])[0])
    planets["koi"] = kois
    joint = pd.merge(planets, df, on="koi")
    joint.to_csv("planet_periods.csv")

def plot_periods():
    df = pd.read_csv("planet_periods.csv")
    plt.clf()
    print(df.period)
    plt.plot(df.koi_period.values, df.period.values)
    plt.savefig("period_period")


if __name__ == "__main__":
    # save_data(100)
    # make_histogram()
    make_df()
    plot_periods()

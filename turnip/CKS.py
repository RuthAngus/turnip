"""
Plot the CKS ages against the gyro ages.
Underpredicts the ages of M dwarfs/ the M dwarfs rotate too rapidly?
Overpredicts the ages of hot stars /the hot stars rotate too slowly?
I should get the M right and underpredict the hot stars in the van saders
model.
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import teff_bv as tbv

plotpar = {'axes.labelsize': 18,
           'text.fontsize': 10,
           'legend.fontsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)


def gyro_age(par, period, bv):
    bvc = a*(bv - c)**b
    return (period/bvc)**(1./n) * 1e-3


def uncertainty_hist(df):
    plt.clf()
    mean_uncert = .5*(df.iso_sage_err1.values/df.iso_sage.values -
                      df.iso_sage_err2.values/df.iso_sage.values)
    plt.hist(mean_uncert, 30)
    plt.axvline(np.median(mean_uncert), color="k")
    print(np.median(mean_uncert))
    plt.savefig("cks_uncertainty_hist")


if __name__ == "__main__":

    my = pd.read_csv("data/koi_periods_0712.csv")
    # my = pd.read_csv("data/Table_1_Periodic_mod.csv")
    df = pd.read_csv("data/cks_physical_merged.csv")

    uncertainty_hist(df)

    # a, b, n, c = .7725, .601, .5189, .4  # Barnes
    a, b, n, c = .407, .325, .566, .495  # MH
    par = [a, b, n, c]

    id_starname = []
    for koi in my.KOI.values:
        id_starname.append("K{}".format(str(koi).zfill(5)))
    my["id_starname"] = id_starname

    new = pd.merge(my, df, on="id_starname", how="inner")
    # new = pd.merge(my, df, on="id_kic", how="inner")

    bvs = tbv.teff2bv(new.teff.values, new.logg.values, new.iso_smet.values)
    new["bv"] = bvs
    m = bvs > c
    new = new.iloc[m]

    cks_age = new.iso_sage.values
    cks_age_errp = new.iso_sage_err1.values
    cks_age_errm = new.iso_sage_err2.values
    lncks_age = new.iso_slogage.values
    lncks_age_errp = new.iso_slogage_err1.values
    lncks_age_errm = new.iso_slogage_err2.values

    my_age = gyro_age(par, new.period.values, new.bv.values)

    plt.clf()
    plt.errorbar(cks_age, my_age, xerr=[-cks_age_errm, cks_age_errp],
                 fmt="k.", alpha=.1, zorder=0)
    plt.scatter(cks_age, my_age, c=new.teff.values, s=10, zorder=1)
    plt.colorbar()
    xs = np.linspace(0, 14, 100)
    # plt.plot(cks_age, my_age, "k.")
    plt.plot(xs, xs, "--")
    plt.xlabel("CKS age (Gyr)")
    plt.ylabel("Gyro age (Gyr)")
    plt.ylim(0,  14)
    plt.savefig("cks_vs_gyro")
    # plt.savefig("cks_vs_gyro_barnes")
    # plt.savefig("cks_vs_gyro_mcquillan")


    plt.clf()
    plt.errorbar(lncks_age, np.log10(my_age*1e9), xerr=[-lncks_age_errm,
                                                        lncks_age_errp],
                 fmt="k.", alpha=.05, zorder=0)
    plt.scatter(lncks_age, np.log10(my_age*1e9),
                c=new.teff.values, s=10, zorder=1)
    plt.colorbar()
    xs = np.linspace(min(lncks_age), max(lncks_age), 100)
    plt.plot(xs, xs, "w--", lw=2)
    plt.xlabel("$\log_{10}(\mathrm{CKS~age})$")
    plt.ylabel("$\log_{10}(\mathrm{Gyro~age})$")
    plt.subplots_adjust(left=.15, bottom=.15)
    plt.ylim(8, 10.6)
    plt.xlim(9, 10.2)
    plt.savefig("log_cks_vs_gyro.pdf")
    # plt.savefig("log_cks_vs_gyro_barnes")
    # plt.savefig("log_cks_vs_gyro_mcquillan")

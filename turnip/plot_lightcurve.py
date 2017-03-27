# Plot Kepler light curve

import os
import numpy as np
import matplotlib.pyplot as plt
import sys

from kepler_data import load_kepler_data

plotpar = {'axes.labelsize': 18,
           'text.fontsize': 10,
           'legend.fontsize': 18,
           'xtick.labelsize': 13,
           'ytick.labelsize': 13,
           'text.usetex': True}
plt.rcParams.update(plotpar)


def plot(kepid):
    path = "/Users/ruthangus/.kplr/data/lightcurves/"
    LC_DIR = os.path.join(path, "{}".format(str(int(kepid)).zfill(9)))
    print(LC_DIR)
    x, y, yerr = load_kepler_data(LC_DIR)

    plt.clf()
    plt.plot(x-x[0], y, "k.")
    plt.xlabel("$\mathrm{Time~(Days)}$")
    plt.ylabel("$\mathrm{Normalised~flux}$")
    plt.xlim(0, 20)
    plt.subplots_adjust(left=.18, bottom=.15)
    plt.savefig("{}_lc".format(kepid))


if __name__ == "__main__":
    # kepid = sys.argv[1]
    # plot(kepid)

    # Load list of all KIDs in ~/.kplr/data/lightcurves/
    kics = np.genfromtxt("kics.txt", dtype=str).T
    N = len(kics)
    for i, kic in enumerate(kics):
        print(kic, i, "of", N)
        plot(kic)

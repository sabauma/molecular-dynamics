
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colo
from matplotlib import rc
from CollisionData import CollisionData
from scipy.stats import gaussian_kde
from scipy.stats.mstats import mquantiles
from pylab import *
rc('text', usetex=True)


import os
import sys

import numpy as np

def main(fname):
    collision_data = CollisionData.data_from_file(fname)

    if collision_data is None:
        return

    energies  = collision_data.energies * 6.24150934e15
    distances = collision_data.distances

    if len(energies) < 100:
        return

    for d, e in zip(distances, energies):
        print d, ',', e

    #energies.sort()

    low  = np.min(energies)
    high = np.max(energies)

    print fname, np.max(energies)

    hist = gaussian_kde(energies)
    n, bins = np.histogram(energies, bins=100)

    fig = plt.figure()
    ax  = fig.add_subplot(1, 1, 1)

    # print np.mean(xe_dist)

    # update the view limits
    #ax.plot(diff_xs, diff_hist(diff_xs), c='r', marker='.', label='He-Xe')
    #ax.plot(he_xs, he_hist(he_xs), c='g', marker='.', label='He-He')
    #ax.plot(bins[:-1], n, c='b', marker='.')
    ax.scatter(energies, distances)
    ax.set_xlim(low, high)
    ax.set_ylim(0.0, np.max(distances))
    ax.set_title("Collision Energy Histogram")
    ax.set_xlabel("Collision Energy (keV)")
    ax.set_ylabel("Relative Density")
    ax.legend()

    fig.savefig(os.path.splitext(fname)[0] + ".png", dpi=250)
    plt.close()
    fig = None
    ax = None

if __name__ == '__main__':
    for fname in sys.argv[1:]:
        main(fname)

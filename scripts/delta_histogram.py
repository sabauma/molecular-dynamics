
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

    t1 = collision_data.type1s
    t2 = collision_data.type2s
    deltaVs = collision_data.deltaVs
    distances = collision_data.distances

    low  = np.min(deltaVs)
    high = np.max(deltaVs)

    #diff_indices = np.where(t1 != t2)
    #diff_deltaVs = deltaVs[diff_indices]
    #diff_hist = gaussian_kde(diff_deltaVs)
    #diff_xs   = np.linspace(low, high, 200)

    #he_indices = np.where(np.logical_and(t1 == 0, t2 == 0))
    #he_deltaVs = deltaVs[he_indices]
    #he_hist = gaussian_kde(he_deltaVs)
    #he_xs   = np.linspace(low, high, 200)
    #he_dist = distances[he_indices]

    xe_indices = np.where(np.logical_and(t1 == 7, t2 == 7))
    xe_deltaVs = deltaVs[xe_indices]

    xe_xs   = np.linspace(low, high, 200)
    xe_dist = distances[xe_indices]
    #print np.mean(xe_dist), np.median(xe_dist), np.mean(distances), np.median(distances)
    print mquantiles(distances, prob=[0.8, 0.85, 0.9, 0.95, 0.975, 0.99])

    if len(xe_deltaVs) <= 1:
        return

    xe_hist = gaussian_kde(xe_deltaVs)

    fig = plt.figure()
    ax  = fig.add_subplot(1, 1, 1)

    # print np.mean(xe_dist)

    # update the view limits
    #ax.plot(diff_xs, diff_hist(diff_xs), c='r', marker='.', label='He-Xe')
    #ax.plot(he_xs, he_hist(he_xs), c='g', marker='.', label='He-He')
    ax.plot(xe_xs, xe_hist(xe_xs), c='b', marker='.', label='Xe-Xe')
    ax.set_xlim(0, high)
    ax.set_title("Collision Radius vs Difference in Velocity")
    ax.set_xlabel("Delta V (m/s)")
    ax.set_ylabel("Relative Density")
    ax.legend()

    fig.savefig(os.path.splitext(fname)[0] + ".png", dpi=250)
    plt.close()
    fig = None
    ax = None

if __name__ == '__main__':
    for fname in sys.argv[1:]:
        main(fname)

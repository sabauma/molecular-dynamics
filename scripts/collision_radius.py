
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colo
from matplotlib import rc
from CollisionData import CollisionData
rc('text', usetex=True)

from pylab import *

import sys

import numpy as np

def main(fname):
    collision_data = CollisionData.data_from_file(fname)

    fig, ax = plt.subplots()

    t1 = collision_data.type1s
    t2 = collision_data.type2s
    deltaVs = collision_data.deltaVs
    distances = collision_data.distances

    del collision_data

    diff_indices = np.where(t1 != t2)[::10]
    diff_deltaVs = deltaVs[diff_indices]
    diff_distances = distances[diff_indices]

    he_indices = np.where(np.logical_and(t1 == 0, t2 == 0))[::10]
    he_deltaVs = deltaVs[he_indices]
    he_distances = distances[he_indices]

    xe_indices = np.where(np.logical_and(t1 == 7, t2 == 7))[::10]
    xe_deltaVs = deltaVs[xe_indices]
    xe_distances = distances[xe_indices]

    print np.mean(xe_distances)

    ax.scatter(diff_deltaVs, diff_distances, c='r', marker='.', label='He-Xe')
    ax.scatter(he_deltaVs, he_distances, c='g', marker='.', label='He-He')
    ax.scatter(xe_deltaVs, xe_distances, c='b', marker='.', label='Xe-Xe')
    ax.set_xlim(0, np.max(deltaVs))
    ax.set_ylim(np.min(distances), np.max(distances))
    ax.set_title("Collision Radius vs Difference in Velocity")
    ax.set_xlabel("Delta V (m/s)")
    ax.set_ylabel("Collision Radius (m)")
    ax.legend()
    fig.savefig(fname + ".png", dpi=250)
    plt.clf()

if __name__ == '__main__':
    for fname in sys.argv[1:]:
        main(fname)

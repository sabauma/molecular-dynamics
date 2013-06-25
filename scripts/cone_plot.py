
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
rc('text', usetex=True)

from pylab import *

import numpy as np
import sys

GRID_WIDTH  = 1000
GRID_HEIGHT = 1000

def data_from_file(fname):
    jdata = None
    positions = None
    velocities = None
    types = None
    with open(fname, 'r') as handle:
        data = np.genfromtxt(handle, delimiter=',')
        positions = data[:,0:3]
        velocities = data[:,3:6]
        if np.shape(data)[1] == 7:
            types = data[:,6]
    return (positions, velocities, types)

def main(fname):
    pos, vel, types = data_from_file(fname)

    x_min = np.min(pos[:,1])
    x_max = np.max(pos[:,1])
    z_min = np.min(pos[:,2])
    z_max = np.max(pos[:,2])

    x_range = x_max - x_min
    z_range = z_max - z_min

    X = np.ones((GRID_WIDTH, GRID_HEIGHT, 3))

    colors = { 0 : np.array([0.8, 0.0, 0.0]),
               1 : np.array([0.0, 0.8, 0.0]),
               2 : np.array([0.7, 0.7, 0.0]),
               3 : np.array([0.5, 0.0, 0.5]),
               4 : np.array([0.0, 0.5, 0.5]),
               5 : np.array([0.2, 0.8, 0.2]),
               6 : np.array([0.8, 0.2, 0.8]),
               7 : np.array([0.0, 0.0, 0.8]) }

    gridxs = ((pos[:,1] - x_min) / x_range * (GRID_WIDTH - 1)).astype(int)
    gridys = ((pos[:,2] - z_min) / z_range * (GRID_HEIGHT - 1)).astype(int)

    X[gridys, gridxs, :] = [colors[i] for i in types]

    plt.imshow(X, origin='lower', aspect=0.5, cmap=cm.jet)
    plt.savefig(fname + ".png", dpi=200)
    plt.clf()

if __name__ == '__main__':
    for fname in sys.argv[1:]:
        print fname
        main(fname)

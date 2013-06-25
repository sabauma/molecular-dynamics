
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colo
from matplotlib import rc
from CollisionData import CollisionData
rc('text', usetex=True)

from pylab import *

import sys

import numpy as np

GRID_WIDTH  = 1000
GRID_HEIGHT = 1000

def main(fname):
    collision_data = CollisionData.data_from_file(fname)
    time = collision_data.times
    pos  = collision_data.positions
    ene  = collision_data.energies

    ene_min = np.min(ene)
    ene_max = np.max(ene)

    jet = plt.get_cmap('jet')
    cNorm = colo.Normalize(vmin=0, vmax=ene_max)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    ep = np.zeros((len(ene), 3))
    for i in range(len(ene)):
        ep[i,:] = scalarMap.to_rgba(ene[i])[0:3]

    x_min = np.min(pos[:,0])
    x_max = np.max(pos[:,0])
    z_min = np.min(pos[:,2])
    z_max = np.max(pos[:,2])

    x_range = x_max - x_min
    z_range = z_max - z_min

    X = np.zeros((GRID_WIDTH, GRID_HEIGHT, 3))

    grid_x = ((pos[:,0] - x_min) / x_range * GRID_WIDTH).astype(int) - 1
    grid_y = ((pos[:,2] - z_min) / z_range * GRID_HEIGHT).astype(int) - 1

    X[grid_y, grid_x] = ep

    xs = np.arange(0, GRID_WIDTH) * x_range / GRID_WIDTH + x_min
    ys = np.arange(0, GRID_HEIGHT) * z_range / GRID_HEIGHT + z_min
    extents = (x_min, x_max, z_min, z_max)

    plt.imshow(X, origin='lower', cmap=cmx.jet, aspect='auto', extent=extents)
    plt.savefig(fname + ".png", dpi=150)
    plt.clf()

if __name__ == '__main__':
    for fname in sys.argv[1:]:
        main(fname)


from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
rc('text', usetex=True)

from pylab import *

import json
import glob
import numpy as np
import sys

GRID_WIDTH  = 1000
GRID_HEIGHT = 1000

POOL_SIZE = 8

def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

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

def main(args):
    for fname, snapshot in args:
        snap = None
        with open(snapshot, 'r') as infile:
            snap = json.load(infile)

        time = None
        if snap is not None:
            time = snap['t']
        snap = None

        print(fname, snapshot)
        pos, vel, types = data_from_file(fname)

        x_min = np.min(pos[:,1])
        x_max = np.max(pos[:,1])
        z_min = np.min(pos[:,2])
        z_max = np.max(pos[:,2])

        x_range = x_max - x_min
        z_range = z_max - z_min

        X = np.ones((GRID_WIDTH, GRID_HEIGHT, 3))

        colors = { 0 : np.array([0.8, 0.0, 0.0]),
                   1 : np.array([0.8, 0.0, 0.0]),
                   2 : np.array([0.7, 0.7, 0.0]),
                   3 : np.array([0.0, 0.8, 0.0]),
                   4 : np.array([0.0, 0.5, 0.5]),
                   5 : np.array([0.2, 0.8, 0.2]),
                   6 : np.array([0.8, 0.2, 0.8]),
                   7 : np.array([0.0, 0.0, 0.8]) }

        gridxs = ((pos[:,1] - x_min) / x_range * (GRID_WIDTH - 1)).astype(int)
        gridys = ((pos[:,2] - z_min) / z_range * (GRID_HEIGHT - 1)).astype(int)

        fig = plt.figure()
        ax  = fig.add_subplot(1, 1, 1)

        X[gridys, gridxs, :] = [colors[i] for i in types]

        im = ax.imshow(X, origin='lower', aspect=1.0, cmap=cm.jet)

        ax.set_title("t= %0.5f ns" % (time * 1.0e9))
        ax.set_xlabel("Red = D \n Green = D2 \n Blue = Xe")
        im.axes.get_xaxis().set_visible(False)
        im.axes.get_yaxis().set_visible(False)

        plt.savefig(fname + ".png", dpi=200)
        fig.clf()
        plt.clf()

if __name__ == '__main__':
    saves = glob.glob("SavePoint*.txt")
    snaps = glob.glob("Snapshot*.dat")
    saves.sort()
    snaps.sort()
    pairs = list(zip(saves, snaps))
    pool = Pool(processes=POOL_SIZE)

    pool.map_async(main, chunks(pairs, POOL_SIZE)).get(99999999)

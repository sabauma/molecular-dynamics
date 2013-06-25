
from matplotlib import rc
rc('text', usetex=True)

from pylab import *

import glob
import itertools
import math
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
import os
import sys

try:
    import simplejson as json
except ImportError:
    print "No simplejson for you!"
    import json

Length = 5.65e-10
Time   = 6.48148e-12
Velocity = Length / Time

BINS = 200

def data_from_file(fname):
    jdata = None
    positions = None
    velocities = None
    types = None
    with open(fname, 'r') as handle:
        data = np.loadtxt(handle, delimiter=',')
        positions = data[:,0:3]
        velocities = data[:,3:6]
        if np.shape(data)[1] == 7:
            types = data[:,6]
    return (positions, velocities, types)

def make_histogram(fname):
    positions, velocities, types = data_from_file(fname)

    # Extract the radii from the positions
    data = np.apply_along_axis(numpy.linalg.norm, 1, velocities) * Velocity

    data_he = data[np.where(types == 0)]
    data_xe = data[np.where(types == 3)]

    fig = plt.figure()
    ax1  = fig.add_subplot(2, 1, 1)
    ax2  = fig.add_subplot(2, 1, 2)

    lower = 0.0
    upper = np.max(data)

    n, bins = np.histogram(data_he, BINS, range=(lower, upper))
    n = n.astype(float) / len(data)

    left = np.array(bins[:-1])
    right = np.array(bins[1:])
    bottom = np.zeros(len(left))
    top = bottom + n

    # we need a (numrects x numsides x 2) numpy array for the path helper
    # function to build a compound path
    XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T

    # get the Path object
    barpath = path.Path.make_compound_path_from_polys(XY)

    # make a patch out of it
    patch = patches.PathPatch(barpath, facecolor='blue', edgecolor='gray', alpha=0.8)
    ax1.add_patch(patch)

    # update the view limits
    ax1.set_xlim(left[0], right[-1])
    ax1.set_ylim(bottom.min(), top.max())
    ax1.set_title("Helium")

    n, bins = np.histogram(data_xe, BINS, range=(lower, upper))
    n = n.astype(float) / len(data)

    left  = np.array(bins[:-1])
    right = np.array(bins[1:])
    bottom = np.zeros(len(left))
    top = bottom + n

    # we need a (numrects x numsides x 2) numpy array for the path helper
    # function to build a compound path
    XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T

    # get the Path object
    barpath = path.Path.make_compound_path_from_polys(XY)

    # make a patch out of it
    patch = patches.PathPatch(barpath, facecolor='blue', edgecolor='gray', alpha=0.8)
    ax2.add_patch(patch)

    # update the view limits
    ax2.set_xlim(left[0], right[-1])
    ax2.set_ylim(bottom.min(), top.max())
    ax2.set_title("Xenon")

    fig.savefig(os.path.splitext(fname)[0] + ".png", dpi=125)
    plt.close()

for fname in sys.argv[1:]:
    print fname
    make_histogram(fname)

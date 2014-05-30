
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
import matplotlib.cm as cm
import numpy as np
import numpy.linalg
import os
import sys
from CollisionData import CollisionData

def histogram_from_data(data):

    radii     = np.sqrt((data.positions * data.positions).sum(axis=1))
    ene, bins = np.histogram(radii, 50, weights=data.energies)
    n, _      = np.histogram(radii, 50)

    ind = np.where(n == 0.0)
    n[ind] = 1.0

    left  = np.array(bins[:-1])
    right = np.array(bins[1:])

    middle = 0.5 * (right + left)

    return (middle, ene / n)

def make_histogram(data, fname):

    middle, ene = histogram_from_data(data)

    fig = plt.figure()
    ax  = fig.add_subplot(1, 1, 1)

    # update the view limits
    ax.set_xlim(middle[0], middle[-1])
    #ax.set_ylim(0.0, 10.0)
    ax.set_xlabel("Radius")
    ax.set_ylabel("Density")

    ax.plot(middle, ene, c='r')

    fig.savefig(os.path.splitext(fname)[0] + ".png", dpi=150)
    plt.close()
    fig = None
    ax = None

for fname in sys.argv[1:]:
    print fname
    data = CollisionData.data_from_file(fname)
    make_histogram(data, fname)



from matplotlib import rc
rc('text', usetex=True)

from pylab import *

import glob
import itertools
import math
import matplotlib.patches as patches
import matplotlib.path as path
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

def json_data(pattern):
    files = glob.glob(pattern)
    for fname in sorted(files):
        jdata = None
        with open(fname, 'r') as handle:
            jdata = json.load(handle)
        yield jdata

def json_from_file(fname):
    jdata = None
    with open(fname, 'r') as handle:
        jdata = json.load(handle)
    return jdata

def make_histogram(fnames, labels=None):

    if labels is None:
        labels = ["" for _ in fname]

    fig = plt.figure()
    ax  = fig.add_subplot(1, 1, 1)

    data = json_from_file(fname)

    types     = np.array([p['Type'] for p in data['Particles']])
    positions = np.array([np.array(p['Position']) for p in data['Particles']])

    hydro_p = positions[np.where(types == 0)]
    xenon_p = positions[np.where(types == 3)]

    data = np.array([numpy.linalg.norm(p) for p in positions])
    hydro_data = data[np.where(types == 0)]
    xenon_data = data[np.where(types == 3)]

    lower = 0.0
    upper = np.max(data)

    n, bins = np.histogram(data, 50, range=(lower, upper))
    h_n, _  = np.histogram(hydro_data, 50, range=(lower, upper))
    x_n, _  = np.histogram(xenon_data, 50, range=(lower, upper))

    left   = np.array(bins[:-1])
    right  = np.array(bins[1:])
    bottom = np.zeros(len(left))
    top    = bottom + n

    shell_volume = 1.0 / 3.0  * (right ** 3.0 - left ** 3.0)

    bubble_density = len(positions) / (1.0 / 3.0 * np.max(data) ** 3.0)

    total_density = n / shell_volume / bubble_density

    middle = (bins[1:] + bins[:-1]) / 2.0
    ax.scatter(middle, total_density, c='r')
    ax.scatter(middle, x_n / shell_volume / bubble_density, c='g')
    ax.scatter(middle, h_n / shell_volume / bubble_density, c='b')

    # update the view limits
    ax.set_xlim(middle[0], middle[-1])
    ax.set_ylim(0.0, 2.0)
    ax.set_xlabel("Radius")
    ax.set_ylabel("Frequency")

    fig.savefig(os.path.splitext(fname)[0] + ".png", dpi=500)
    plt.close()
    fig = None
    ax = None

for fname in sys.argv[1:]:
    print fname
    make_histogram(fname)

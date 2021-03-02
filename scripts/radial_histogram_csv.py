
from matplotlib import rc
rc('text', usetex=True)

from pylab import *

import glob
import itertools
import math
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import numpy.linalg
import os
import sys

try:
    import simplejson as json
except ImportError:
    print("No simplejson for you!")
    import json

COLORS = [ np.array([0.8, 0.0, 0.0]),
           np.array([0.0, 0.8, 0.0]),
           np.array([0.0, 0.0, 0.8]),
           np.array([0.5, 0.0, 0.5]),
           np.array([0.0, 0.5, 0.5]),
           np.array([0.2, 0.8, 0.2]),
           np.array([0.8, 0.2, 0.8]),
           np.array([0.7, 0.7, 0.0]) ]

def setAxLinesBW(ax):
    """
    Take each Line2D in the axes, ax, and convert the line style to be 
    suitable for black and white viewing.
    """
    MARKERSIZE = 3

    COLORMAP = {
        'b': {'marker': None, 'dash': (None,None)},
        'g': {'marker': None, 'dash': [5,5]},
        'r': {'marker': None, 'dash': [5,3,1,3]},
        'c': {'marker': None, 'dash': [1,3]},
        'm': {'marker': None, 'dash': [5,2,5,2,5,10]},
        'y': {'marker': None, 'dash': [5,3,1,2,1,10]},
        'k': {'marker': 'o', 'dash': (None,None)} #[1,2,1,10]}
        }

    for line in ax.get_lines():
        origColor = line.get_color()
        line.set_color('black')
        #line.set_dashes(COLORMAP[origColor]['dash'])
        #line.set_marker(COLORMAP[origColor]['marker'])
        line.set_markersize(MARKERSIZE)

def setFigLinesBW(fig):
    """
    Take each axes in the figure, and for each line in the axes, make the
    line viewable in black and white.
    """
    for ax in fig.get_axes():
        setAxLinesBW(ax)

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

def histogram_from_data(data, bubble_radius, total_particles, lower=0.0, upper=None):
    if upper is None:
        upper = bubble_radius

    n, bins = np.histogram(data, 50, range=(lower, upper))

    left  = np.array(bins[:-1])
    right = np.array(bins[1:])

    shell_volume = 1.0 / 3.0 * (right ** 3.0 - left ** 3.0)

    bubble_density = total_particles / ((1.0 / 3.0) * bubble_radius ** 3.0)

    total_density = n / (shell_volume * bubble_density)

    middle = 0.5 * (right + left)

    return (middle, total_density)

def make_histogram(fname):
    positions, _, types = data_from_file(fname)

    unique_types    = set(types)
    total_particles = len(positions)

    # Extract the radii from the positions
    data = np.apply_along_axis(numpy.linalg.norm, 1, positions)
    data *= 5.65e-10
    positions = None
    bubble_radius = np.max(data)

    middle, total_density = histogram_from_data(data,
                                                bubble_radius,
                                                total_particles,
                                                0.0,
                                                np.max(data))

    fig = plt.figure()
    ax  = fig.add_subplot(1, 1, 1)

    # update the view limits
    #ax.set_xlim(middle[0], middle[-1])
    #ax.set_ylim(0.0, 10.0)
    ax.set_xlabel("Radius ($\\mu m$)")
    ax.set_ylabel("Density")

    # ax.scatter(middle, total_density, c='r')

    for typ in sorted(unique_types):

        indices = np.where(types == typ)
        species_radii = data[indices]

        middle, total_density = histogram_from_data(species_radii,
                                                    bubble_radius,
                                                    total_particles,
                                                    0.0,
                                                    np.max(data))

        ax.plot(middle / 1.0e-6, total_density)

    setFigLinesBW(fig)
    fig.savefig(os.path.splitext(fname)[0] + ".png")
    plt.close()
    fig = None
    ax = None

for fname in sys.argv[1:]:
    print(fname)
    make_histogram(fname)

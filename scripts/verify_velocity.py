
from matplotlib import rc

rc('text', usetex=True)

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from pylab import *

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
        line.set_dashes(COLORMAP[origColor]['dash'])
        line.set_marker(COLORMAP[origColor]['marker'])
        line.set_markersize(MARKERSIZE)

def setFigLinesBW(fig):
    """
    Take each axes in the figure, and for each line in the axes, make the
    line viewable in black and white.
    """
    for ax in fig.get_axes():
        setAxLinesBW(ax)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for fname in sys.argv[1:]:
    with open(fname, 'r') as infile:
        data = np.genfromtxt(infile, delimiter=',')

    time     = data[:,0] / 1.0e-6
    velocity = data[:,2]

    label_name = os.path.splitext(fname)[0]
    label_name = label_name.replace('_', ' ')
    label_name = label_name.replace('%', '\\%')

    ax.plot(time, velocity, label=label_name)

    print fname, ": ", np.max(np.abs(velocity))

ax.set_xlabel("Time ($\\mu$s)")
ax.set_ylabel("Velocity (m/s)")

setFigLinesBW(fig)
ax.legend(bbox_to_anchor=(0.75, 0.25), loc=2, borderaxespad=0.0)

plt.show()

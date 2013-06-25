
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
import scipy.signal

windows = [1.0e-10, 1.0e-11, 2.0e-12]

def simulate_window(time, force, window_size):
    newdat = []
    for i in xrange(1000, len(time), 100):
        indices = np.where(time >= time[i] - window_size)
        td = time[indices]
        fd = force[indices]
        t_bar = np.mean(td)
        f_bar = np.mean(fd)
        dt = td[-1] - td[0]
        newdat.append([t_bar, f_bar / dt])
    return np.array(newdat)

def chunky_window(time, force, size):
    newdat = []
    for i in xrange(size, len(time), size / 2):
        td = time[i : i + size]
        fd = force[i : i + size]
        t_bar = np.mean(td)
        f_bar = np.mean(fd)
        dt = td[-1] - td[0]
        newdat.append([t_bar, f_bar])

    return np.array(newdat)

handle = open(sys.argv[1], 'r')

data = np.genfromtxt(handle, delimiter=',')[::50]


time  = data[:,0]
force = data[:,1]

avg = chunky_window(time, force, 200)

fig = plt.figure()
ax  = fig.add_subplot(2, 1, 1)
ax_gauge = fig.add_subplot(2, 1, 2)

ax.set_yscale('log')
ax_gauge.set_yscale('log')

ax.set_xlabel("Time")
ax.set_ylabel(r"$\frac{\Delta m}{A}$")
ax.set_xlim(-np.max(time) * 0.05, np.max(time) * 1.05)
ax.scatter(time, force)
ax.plot(avg[:,0], avg[:,1], c='r')

for size in windows:
    gauge = simulate_window(time, force, size)
    ax_gauge.plot(gauge[:,0], gauge[:,1], label=str(size))

ax_gauge.legend(bbox_to_anchor=(1.00, 1.00), loc=2, borderaxespad=0.0)

show()

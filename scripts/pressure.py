# -*- coding: utf-8 -*-
# TODO: Overlay radius onto the plot.
from matplotlib import rc

rc('text', usetex=True)

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from pylab import *

GAMMA = 5.0 / 3.0
a3 = 2.179270629335e-20
R0 = 2.2e-6
SurfaceTension = 7.5e-2
P0 = 1.0e5 + 2.0 * SurfaceTension / R0

RADI = 0.95 * R0
RISO = 1.00 * R0

infile = open(sys.argv[1], "r")
data = np.genfromtxt(infile, delimiter=',')
infile.close()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

time     = data[:,0]
radius   = data[:,1]
velocity = data[:,2]
pressure = data[:,3]
dpdt     = data[:,4]

k = np.array([ 1.0 if r > RISO else GAMMA for r in radius ])

volume = radius ** 3

ideal_pressure = P0 * ((R0 ** 3.0) / (volume)) ** k

ax.set_yscale('log')
ax.plot(time, pressure, label="pressure")
ax.plot(time, ideal_pressure, label="ideal\_pressure")
ax.legend(bbox_to_anchor=(1.00, 1.00), loc=2, borderaxespad=0.0)

show()

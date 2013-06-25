
from matplotlib import rc
rc('text', usetex=True)

from pylab import *

import math
import glob
import itertools
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

CONE_SOLID_ANGLE = 0.00061246
R_H2 = 1.09e-10
R_XE = 1.51e-10

RADIUS = {'He' : 1.09e-10, 'Xe' : 1.51e-10}

def json_data(pattern):
    files = glob.glob(pattern)
    for fname in sorted(files):
        jdata = None
        with open(fname, 'r') as handle:
            jdata = json.load(handle)
        yield jdata

real_time    = []
sim_time     = []
collisions   = []
wall_hits    = []
cone_hits    = []
cell_crosses = []
radii        = []
particle_vol = []

for data in json_data("Snapshot*.dat"):
    real_time.append(data['real_time'])
    sim_time.append(data['t'])

    collisions.append(data['TotalCollisions'])
    wall_hits.append(data['TotalWallCollisions'])
    cone_hits.append(data['TotalConeBoundaryCollisions'])
    cell_crosses.append(data['TotalCellCrossings'])

    radii.append(data['r'])

    total_volume = 0.0

    for b in data['bins']:
        for elem, percent in b['percent'].iteritems():
            total_volume += (4.0 / 3.0 * math.pi * RADIUS[elem] ** 3.0) * b['n'] * percent

    particle_vol.append(total_volume)

print radii

# Convert to arrays for data processing
real_time    = np.array(real_time)
sim_time     = np.array(sim_time)
collisions   = np.array(collisions)
wall_hits    = np.array(wall_hits)
cone_hits    = np.array(cone_hits)
cell_crosses = np.array(cell_crosses)
radii        = np.array(radii)
particle_vol = np.array(particle_vol)

dt = sim_time - np.concatenate((np.zeros(1), sim_time[:-1]))

collision_rate = collisions / dt
wall_rate = wall_hits / dt
bubble_volume = 1.0 / 3.0 * CONE_SOLID_ANGLE * (radii ** 3.0)
packing_factor = particle_vol / bubble_volume

print np.max(collision_rate)

fig = plt.figure()
ax1 = fig.add_subplot(3, 1, 1)
ax2 = fig.add_subplot(3, 1, 2)
ax3 = fig.add_subplot(3, 1, 3)

ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
#plot(x, np.log(volume), x, np.log(constant))
#ax1.scatter(sim_time, collision_rate, label="collision rate")
#ax1.scatter(sim_time, wall_rate, label="wall rate")
ax1.plot(sim_time, real_time, label="clock vs sim")
ax1.legend(bbox_to_anchor=(1.00, 1.00), loc=2, borderaxespad=0.0)

##indices = np.where(sim_time >= 0.6e-14 + 4.1676960e-6)
##sim_time = sim_time[indices]
#dt       = dt[indices]
ax2.plot(sim_time, dt, label="real time per simulation time")
ax2.legend(bbox_to_anchor=(1.00, 1.00), loc=2, borderaxespad=0.0)

ax3.plot(sim_time, packing_factor, label="packing factor")
ax3.legend(bbox_to_anchor=(1.00, 1.00), loc=2, borderaxespad=0.0)

print np.max(packing_factor)

show()


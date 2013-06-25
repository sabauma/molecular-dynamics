
from matplotlib import rc

rc('text', usetex=True)

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from pylab import *
import math

infile = open(sys.argv[1], "r")
data = np.genfromtxt(infile, delimiter=',')
infile.close()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

FREQUENCY = 30.0e3
OMEGA = 2.0 * 3.14159265358979323 * FREQUENCY
T_INI = 1.675146020165000e-5

GAMMA = 5.0 / 3.0
a3 = 2.179270629335e-20
R0 = 2.2e-6
SurfaceTension = 7.5e-2
P0 = 1.0e5 + 2.0 * SurfaceTension / R0

AMP = 1.65 * P0

def wave_function(x):
    return -AMP * sin(OMEGA * (time + T_INI))

time     = data[:,0]
radius   = data[:,1]
velocity = data[:,2]
pressure = data[:,3]

dv = velocity[1:] - velocity[:-1]
dt = time[1:] - time[:-1]

# pressure = pressure ** (5.0 / 3.0)
y = np.repeat(radius[-1], len(radius))
amp = np.apply_along_axis(wave_function, 0, time)

#ax.set_yscale('log')
ax.plot(time, radius, label=r"$R$")
#ax.plot(time, amp, label=r"$P_a$")
#ax.plot(time, pressure, label=r"$P_g$")
#ax.plot(time[1:], dv / dt, label=r"a")
#ax.plot(time, y, label=r"$\bar{V}$")
ax.legend(bbox_to_anchor=(1.00, 1.00), loc=2, borderaxespad=0.0)

fig.savefig(os.path.splitext(sys.argv[1])[0] + ".png", dpi=500)
show()


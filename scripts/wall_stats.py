
from matplotlib import rc

rc('text', usetex=True)

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from pylab import *

infile = open(sys.argv[1], "r")
data = np.genfromtxt(infile, delimiter=',')
infile.close()

time     = data[:,0]
radius   = data[:,1]
velocity = data[:,2]
pressure = data[:,3]
dpdt     = data[:,4]

fig = plt.figure()
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2)

ax1.set_yscale('log')

neg_indices = velocity < 0.0

ax1.set_xlim(0.0, np.max(time) * 1.125)

ax1.plot(time, pressure, label="Pressure (Pa)")
ax1.legend(bbox_to_anchor=(1.00, 1.00), loc=2, borderaxespad=0.0)

ax2.set_yscale('log')
ax2.set_xlim(ax1.get_xlim())
ax2.scatter(time, np.abs(velocity), c=neg_indices, label="Velocity $(\\frac{m}{s})$")
ax2.legend(bbox_to_anchor=(1.00, 1.00), loc=2, borderaxespad=0.0)

fig.savefig(os.path.splitext(sys.argv[1])[0] + ".png", dpi=500)
show()


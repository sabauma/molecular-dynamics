
from matplotlib import rc

rc('text', usetex=True)

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm

fig = plt.figure()

infile = open(sys.argv[1], "r")
data = np.genfromtxt(infile, delimiter=',', skiprows=1)
infile.close()

ambient   = data[:,0]
max_r     = data[:,1]
wall_temp = data[:,2]
min_r     = data[:,3]
max_v     = data[:,4]
ratio     = data[:,5]
fusions   = data[:,6]
expansion = max_r / ambient

indices = np.where(fusions <= 0.0)
#ambient = ambient[indices]
#max_r     = max_r[indices]
#wall_temp = wall_temp[indices]
#min_r     = min_r[indices]
#max_v     = max_v[indices]
#ratio     = ratio[indices]
#fusions   = fusions[indices]
#expansion = max_r / ambient

x = np.column_stack((ambient, expansion, wall_temp, ambient * expansion))  #stack explanatory variables into an array
x = sm.add_constant(x, prepend=True) #add a constant
y = np.log(fusions)
y[indices] = -1.0

res = sm.OLS(y,x).fit() #create a model and fit it
print res.params
print res.bse
print res.summary()

max_r_unique     = np.array([i for i in sorted(set(max_r))])
wall_temp_unique = np.array([i for i in sorted(set(wall_temp))])
ambient_r_unique = np.array([i for i in sorted(set(ambient))])

max_r_to_ind     = dict((j, i) for i, j in enumerate(max_r_unique))
wall_temp_to_ind = dict((j, i) for i, j in enumerate(wall_temp_unique))
ambient_to_ind = dict((j, i) for i, j in enumerate(ambient_r_unique))

###################################################################################
# Bounding frames, ambient x expansion
###################################################################################
#exp[ -20.9791 + 5.322e6 * ambient + 0.3437 * expansion + 0.0043 * wall_temp ]

example_ambient   = np.arange(0.5e-6, 5.0e-6, 0.5e-7)
example_expansion = np.arange(5.0, 50.0, 0.5)

X, Y = np.meshgrid(example_ambient, example_expansion)
Z1 = np.exp(-6.8669 + -1.832e6 * X + 0.0289 * Y + 0.0039 * 1000 + 1.624e5 * X * Y)
Z2 = np.exp(-6.8669 + -1.832e6 * X + 0.0289 * Y + 0.0039 * 500  + 1.624e5 * X * Y)

#undesired = np.where(np.logical_or(Z1 < 0.0, Z1 > 50.0))
#Z1[undesired] = np.nan

#ax = Axes3D(fig)
#ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, cmap=cm.coolwarm,
        #linewidth=0, antialiased=False )

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
c1 = ax1.contour(X, Y, Z1, levels=[1.0, 3.0, 5.0, 10.0, 25.0, 50.0, 75.0, 100.0], cmap=cm.jet)
c2 = ax2.contour(X, Y, Z2, levels=[1.0, 3.0, 5.0, 10.0, 25.0, 50.0, 75.0, 100.0], cmap=cm.jet)
#ax.plot_surface(X, Y, Z2, rstride=1, cstride=1, cmap=cm.coolwarm,
        #linewidth=0, antialiased=False)
#sc = ax.plot_surface(max_r, wall_temp, fusions, c=ambient, label="Ambient Radius")
ax1.set_title("Fusion rate at 1000C")
ax1.set_xlabel("Ambient Radius (m)")
ax1.set_ylabel("Expansion Factor")
ax1.clabel(c1, inline=1, fontsize=10)

ax2.set_title("Fusion rate at 500C")
ax2.set_xlabel("Ambient Radius (m)")
ax2.set_ylabel("Expansion Factor")
ax2.clabel(c2, inline=1, fontsize=10)
#ax.set_zlabel("Fusion Events")
#fig.colorbar(sc)

plt.show()
fig.clf()


####################################################################################
# Wall velocity vs radius
####################################################################################

fus = fusions

ax1 = fig.add_subplot(3,1,1)
ax2 = fig.add_subplot(3,1,2)
ax3 = fig.add_subplot(3,1,3)
ax1.scatter(ambient, fus, label="Wall Temp")
ax1.set_xlim((np.min(ambient) / 1.1, np.max(ambient) * 1.1))
ax1.set_xlabel("Ambient Radius")
ax1.set_ylabel("Fusions")

ax2.scatter(expansion, fus, label="Wall Temp")
ax2.set_xlim((np.min(expansion) / 1.1, np.max(expansion) * 1.1))
ax2.set_xlabel("Expansion Ratio")
ax2.set_ylabel("Fusions")


ax3.scatter(wall_temp, fus, label="Wall Temp")
ax3.set_xlim((np.min(wall_temp) / 1.1, np.max(wall_temp) * 1.1))
ax3.set_xlabel("Wall Temperature")
ax3.set_ylabel("Fusions")

fig.savefig("plots.png", dpi=250)
fig.clf()

####################################################################################
# Wall velocity vs radius
####################################################################################

ax = fig.add_subplot(1,1,1)
ax.set_xlim((np.min(max_r) / 1.1, np.max(max_r) * 1.1))
ax.set_ylim((np.min(max_v) / 1.1, np.max(max_v) * 1.1))
sc = ax.scatter(max_r, max_v, c=wall_temp, label="Wall Temp")
fig.colorbar(sc)
ax.set_xlabel("Maximum Radius (m)")
ax.set_ylabel("Maximum Wall Velocity (m/s)")
ax.set_title("Wall Velocity vs Initial Radius")

fig.savefig("wall_velocity_vs_max_radius.png", dpi=250)
fig.clf()

####################################################################################
# Wall velocity vs expansion ratio
####################################################################################

ax = fig.add_subplot(1,1,1)
ax.set_xlim((np.min(expansion) / 1.1, np.max(expansion) * 1.1))
ax.set_ylim((np.min(max_v) / 1.1, np.max(max_v) * 1.1))
sc = ax.scatter(expansion, max_v, c=ambient, label="Wall Temp")
fig.colorbar(sc)
ax.set_xlabel("Expansion Ratio")
ax.set_ylabel("Maximum Wall Velocity (m/s)")
ax.set_title("Wall Velocity vs Expansion Ratio")

fig.savefig("wall_velocity_vs_expansion_ratio.png", dpi=250)
fig.clf()

####################################################################################
# Fusion rate vs radius
####################################################################################

ax = fig.add_subplot(1,1,1)
ax.set_xlim((np.min(max_r) / 1.1, np.max(max_r) * 1.1))
ax.set_ylim((np.min(fusions) / 1.1, np.max(fusions) * 1.1))
sc = ax.scatter(max_r, fusions, c=wall_temp, label="Wall Temp")
fig.colorbar(sc)
ax.set_xlabel("Maximum Radius (m)")
ax.set_ylabel("Fusion Rate")
ax.set_title("Fusion Rate vs Initial Radius")

fig.savefig("fusions_vs_max_radius.png", dpi=250)
fig.clf()

####################################################################################
# Fusion rate vs wall velocity
####################################################################################

ax = fig.add_subplot(1,1,1)
ax.set_xlim((np.min(max_v) / 1.1, np.max(max_v) * 1.1))
ax.set_ylim((np.min(fusions) / 1.1, np.max(fusions) * 1.1))
sc = ax.scatter(max_v, fusions, c=wall_temp, label="Wall Temp")
fig.colorbar(sc)
ax.set_xlabel("Wall Velocity (m/s)")
ax.set_ylabel("Fusion Rate")
ax.set_title("Fusion Rate vs Wall Velocity")

fig.savefig("fusions_vs_wall_velocity.png", dpi=250)
fig.clf()

####################################################################################
# Fusion rate vs radius * temperature
####################################################################################

mesh = np.zeros((len(max_r_unique), len(wall_temp_unique)))

for r, t, f in zip(max_r, wall_temp, fusions):
    mesh[max_r_to_ind[r], wall_temp_to_ind[t]] = f

X, Y = np.meshgrid(max_r_unique, wall_temp_unique)

ax = Axes3D(fig)
ax.plot_wireframe(X, Y, mesh)
#sc = ax.plot_surface(max_r, wall_temp, fusions, c=ambient, label="Ambient Radius")
ax.legend(bbox_to_anchor=(1.00, 1.00), loc=2, borderaxespad=0.0)
ax.set_title("Fusion rate vs Temperature x Max radius")
ax.set_xlabel("Maximum Radius (m)")
ax.set_ylabel("Temperature (C)")
ax.set_zlabel("Fusion Events")
#fig.colorbar(sc)

fig.savefig("fusion_rate_vs_max_radius_wall_temp.png", dpi=250)
fig.clf()

####################################################################################
# Fusion rate vs radius * wall velocity
####################################################################################

fix_fusions = fusions
fix_fusions[np.where(fix_fusions == 0.0)] = 1.0

ax = fig.add_subplot(1,1,1)
ax.set_xlim((np.min(max_r) / 1.1, np.max(max_r) * 1.1))
ax.set_ylim((np.min(max_v) / 1.1, np.max(max_v) * 1.1))
sc = ax.scatter(max_r, max_v, c=np.log(fix_fusions), label="Wall Temp")
fig.colorbar(sc)
ax.set_xlabel("Maximum Radius (m)")
ax.set_ylabel("Wall Velocity (m/s)")
ax.set_title("Fusion Rate vs Radius x Wall Velocity")

fig.savefig("fusions_vs_max_radius_wall_velocity.png", dpi=250)
fig.clf()

####################################################################################
# Fusion rate vs Expansion ratio * ambient radius
####################################################################################


mesh = np.zeros((len(ambient_r_unique), len(wall_temp_unique)))

for r, t, f in zip(ambient, wall_temp, fusions):
    mesh[ambient_to_ind[r], wall_temp_to_ind[t]] = f

X, Y = np.meshgrid(ambient_r_unique, wall_temp_unique)

ax = Axes3D(fig)
ax.plot_surface(X, Y, mesh, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0)
# sc = ax.scatter(ambient, expansion, fusions, c=ambient, label="Ambient Radius")
ax.legend(bbox_to_anchor=(1.00, 1.00), loc=2, borderaxespad=0.0)
ax.set_title("Fusion rate vs Expansion ratio x Ambient Radius")
ax.set_xlabel("Ambient Radius (m)")
ax.set_ylabel("Expansion Ratio")
ax.set_zlabel("Fusion Events")

fig.savefig("fusion_rate_vs_expansion_ratio_ambient_radius.png", dpi=250)
fig.clf()

####################################################################################
# Fusion rate vs radius * wall velocity
####################################################################################

fix_fusions = fusions
fix_fusions[np.where(fix_fusions == 0.0)] = 1.0

ax = fig.add_subplot(1,1,1)
ax.set_xlim((np.min(max_r) / 1.1, np.max(max_r) * 1.1))
ax.set_ylim((np.min(wall_temp) / 1.1, np.max(wall_temp) * 1.1))
sc = ax.scatter(max_r, wall_temp, c=np.log(fix_fusions), label="Wall Temp")
fig.colorbar(sc)
ax.set_xlabel("Maximum Radius (m)")
ax.set_ylabel("Wall Temperature (C)")
ax.set_title("Fusion Rate vs Radius x Wall Temperature")

fig.savefig("fusions_vs_max_radius_wall_temperature.png", dpi=250)
fig.clf()


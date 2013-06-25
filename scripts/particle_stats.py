
import glob
import itertools
import math
import numpy as np
import numpy.linalg
import os
import sys

def data_from_file(fname):
    positions = None
    velocities = None
    types = None
    with open(fname, 'r') as handle:
        data = np.loadtxt(handle, delimiter=',')
        positions = data[:,0:3]
        velocities = data[:,3:6]
        if np.shape(data)[1] == 7:
            types = data[:,6]
    return (positions, velocities, types)

def run_statistics(fname):
    positions, velocities, types = data_from_file(fname)
    norm = numpy.linalg.norm

    unique_types = set(types)
    for typ in sorted(unique_types):
        indices = np.where(types == typ)
        species_positions = positions[indices]
        species_velocities = np.apply_along_axis(norm, 1, velocities[indices])

        species_radii = np.apply_along_axis(norm, 1, species_positions)

        print "Species %f" % typ
        print "Radius Mean: ", np.mean(species_radii)
        print "Radius Variance: ", np.var(species_radii)
        print "Velocity Mean: ", np.mean(species_velocities)
        print "Velocity Variance: ", np.var(species_velocities)
        print
        print

run_statistics(sys.argv[1])

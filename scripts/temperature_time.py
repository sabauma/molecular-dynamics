
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

def compute_temps(fname):
    jdata = json_from_file(fname)

    n          = np.array([p['n'] for p in jdata['bins']])
    radii      = np.array([0] + [r['r'] for r in jdata['bins']])
    temperatures   = [t['temperature'] for t in jdata['bins']]
    temperature_xe = np.array([x['Xe'] for x in temperatures])
    temperature_he = np.array([x['He'] for x in temperatures])

    avgs = [np.mean(tempterature_he), np.mean(temperature_xe)]
    errs = [np.var(temperature_he), np.var(temperature_xe)]

    return (avgs, errs)


def main():
    for fname in sys.argv[1:]:
        make_histogram(fname)

if __name__ == '__main__': main()

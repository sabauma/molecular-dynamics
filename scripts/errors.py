from matplotlib import rc
rc('text', usetex=True)

from pylab import *

import glob
import itertools
try:
    import simplejson as json
except ImportError:
    import json
import math
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
import os
import sys

def json_data(pattern):
    files = glob.glob(pattern)
    for fname in sorted(files):
        print fname
        jdata = None
        with open(fname, 'r') as handle:
            jdata = json.load(handle)
        yield jdata

def json_from_file(fname):
    jdata = None
    with open(fname, 'r') as handle:
        jdata = json.load(handle)
    return jdata

def take(n, iterable):
    return itertools.islice(iterable, n)

errors = np.array([data['RPKEquation']['Error'] for data in take(250, json_data('SavePoint*.txt'))])

fig = plt.figure()
ax  = fig.add_subplot(1, 1, 1)

ax.plot(errors)

plt.show()

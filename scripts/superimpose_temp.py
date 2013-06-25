
from matplotlib import rc

rc('text', usetex=True)

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from pylab import *
from snapshot import Snapshot
from make_histogram import make_histograms
import string

def getnum(fname):
    return "".join([i for i in fname if not str.isalpha(i)])[:-1]


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

datas = [ Snapshot.read_snapshot(fname) for fname in sys.argv[1:] ]

make_histograms(datas, "temperature", map(getnum, sys.argv[1:]))

print datas


import numpy as np
import snapshot
import sys

def compute_average(field, snap):
    count = np.array(snap.get_bin_statistic("n"))
    perc = snap.get_bin_statistic("percentage")
    data = snap.get_bin_statistic(field)

    for elt in data[0].keys():
        element_perc = np.array([ i[elt] for i in perc])
        element_data = np.array([ i[elt] for i in data ])
        avg = sum(element_data * element_perc * count) / sum(element_perc * count)
        print "Average ", field, " for ", elt, " is ", avg

for fname in sys.argv[2:]:
    snap = snapshot.Snapshot.read_snapshot(fname)
    compute_average(sys.argv[1], snap)

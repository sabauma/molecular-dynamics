
import sys

from snapshot import Snapshot
from make_histogram import make_histogram

for fname in sys.argv[1:]:
    snap = Snapshot.read_snapshot(fname)
    make_histogram(snap, "temperature")
    make_histogram(snap, "max_temperature")
    make_histogram(snap, "ionization")
    make_histogram(snap, "kinetic_energy", scalar=6.24150974e18)


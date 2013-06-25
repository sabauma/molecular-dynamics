
import sys

from snapshot import Snapshot
from make_histogram import make_histogram

for fname in sys.argv[1:]:
    snap = Snapshot.read_snapshot(fname)
    make_histogram(snap, "temperature", fname)
    make_histogram(snap, "max_temperature", fname)
    make_histogram(snap, "ionization", fname)
    make_histogram(snap, "kinetic_energy", fname, scalar=6.24150974e18)


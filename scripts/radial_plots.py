

from snapshot import Snapshot
from make_histogram import make_histogram
from multiprocessing import Pool

import errno
import os
import sys
import glob

POOL_SIZE = 8

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def process_files(fnames):
    for fname in fnames:
        snap = Snapshot.read_snapshot(fname)

        count = snap.get_bin_statistic("n")
        percent = snap.get_species_statistic("percentage")
        count[count == 0] = 1.0
        energy_scale = 6.24150974e18 / count

        make_histogram(snap, "temperature")
        make_histogram(snap, "max_temperature")
        make_histogram(snap, "ionization")
        make_histogram(snap, "kinetic_energy", scalar=energy_scale)
        make_histogram(snap, "percentage", scalar=snap.get_bin_statistic("density"))

def rename(dir, pattern, titlePattern):
    for pathAndFilename in glob.iglob(os.path.join(dir, pattern)):
        title, ext = os.path.splitext(os.path.basename(pathAndFilename))
        os.rename(pathAndFilename, os.path.join(dir, titlePattern % title + ext))

def main():
    files = sys.argv[1:]
    pool = Pool(processes=POOL_SIZE)

    pool.map(process_files, chunks(files, POOL_SIZE))

    # Create classification directories
    make_sure_path_exists("temperature")
    make_sure_path_exists("max_temperature")
    make_sure_path_exists("ionization")
    make_sure_path_exists("kinetic_energy")
    make_sure_path_exists("percentage")

    rename('.', 'temperature*.png', './temperature/%s')
    rename('.', 'max_temperature*.png', './max_temperature/%s')
    rename('.', 'ionization*.png', './ionization/%s')
    rename('.', 'kinetic_energy*.png', './kinetic_energy/%s')
    rename('.', 'percentage*.png', './percentage/%s')

if __name__ == '__main__':
    main()

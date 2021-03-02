
import numpy as np
import collections

try:
    import simplejson as json
except ImportError:
    print("Was unable to import simplejson!")
    import json

class Snapshot(object):

    def __init__(self):
        self.__data = None

    def __getitem__(self, index):
        return self.__data[index]

    def get_species_statistic(self, field_name):
        cache = collections.defaultdict(list)
        bins  = self.__data["bins"]

        for bin in bins:
            for k, v in bin[field_name].items():
                cache[k].append(v)

        retval = dict((k, np.array(v)) for k, v in cache.items())

        return retval

    def get_bin_statistic(self, field_name):
        bins  = self.__data["bins"]
        return np.array([bin[field_name] for bin in bins])


    @staticmethod
    def read_snapshot(fname):
        retval = Snapshot()

        with open(fname, 'r') as handle:
            retval.__data = json.load(handle)

        return retval

def main():
    snap = Snapshot.read_snapshot("example.dat")
    print(snap.get_species_statistic("temperature"))
    print(snap.get_bin_statistic("density"))

if __name__ == '__main__': main()


import numpy

class CollisionData(object):

    def __init__(self):
        self.times     = None
        self.positions = None
        self.deltaVs   = None
        self.energies  = None
        self.distances = None
        self.type1s    = None
        self.type2s    = None

    @staticmethod
    def data_from_file(fname):
        retval = CollisionData()
        with open(fname, 'r') as handle:
            data = numpy.genfromtxt(handle, delimiter=',', skip_header=1)

            if len(numpy.shape(data)) == 1:
                return None

            retval.times     = data[:,0]
            retval.positions = data[:,1:4]
            retval.deltaVs   = data[:,4]
            retval.energies  = data[:,5]
            retval.distances = data[:,6]
            retval.type1s    = data[:,7]
            retval.type2s    = data[:,8]
        return retval

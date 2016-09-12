import numpy as np
import inspect  # Used for storing the input
from aquifer_parameters import param_maq


class AquiferData:
    def __init__(self, model, kaq, Haq, c, z, npor, ltype, aqnumber):
        # All input variables except model should be numpy arrays
        # That should be checked outside this function
        self.model = model
        # Needed for heads
        self.kaq = kaq
        self.Naq = len(kaq)
        self.Haq = Haq
        self.T = self.kaq * self.Haq
        self.Tcol = self.T[:, np.newaxis]
        self.c = c
        # Needed for tracing
        self.z = np.atleast_1d(z)
        self.Nlayers = len(self.z) - 1
        self.npor = npor
        self.ltype = ltype
        self.aqnumber = aqnumber
        # tag indicating whether an aquifer is Laplace (confined on top)
        if self.ltype[0] == 'a':
            self.ilap = 1
        else:
            self.ilap = 0
        #
        self.area = 1e200  # Smaller than default of ml.aq so that inhom is found

    def initialize(self):
        self.elementlist = []  # Elementlist of aquifer
        d0 = 1.0 / (self.c * self.T)
        d0[:-1] += 1.0 / (self.c[1:] * self.T[:-1])
        dp1 = -1.0 / (self.c[1:] * self.T[1:])
        dm1 = -1.0 / (self.c[1:] * self.T[:-1])
        A = np.diag(dm1, -1) + np.diag(d0, 0) + np.diag(dp1, 1)
        w, v = np.linalg.eig(A)
        # sort lab in decending order, hence w in ascending order
        index = np.argsort(abs(w))
        w = w[index]
        v = v[:, index]
        if self.ilap:
            self.lab = np.zeros(self.Naq)
            self.lab[1:] = 1.0 / np.sqrt(w[1:])
            self.zeropluslab = self.lab  # to be deprecated when new lambda is fully implemented
            v[:, 0] = self.T / np.sum(self.T)  # first column is normalized T
        else:
            self.lab = 1.0 / np.sqrt(w)
        self.eigvec = v
        self.coef = np.linalg.solve(v, np.diag(np.ones(self.Naq))).T

    def add_element(self, e):
        self.elementlist.append(e)

    def isinside(self, x, y):
        raise Exception('Must overload AquiferData.isinside()')

    def storeinput(self, frame):
        self.inputargs, _, _, self.inputvalues = inspect.getargvalues(frame)


class Aquifer(AquiferData):
    def __init__(self, model, kaq, Haq, c, z, npor, ltype, aqnumber):
        AquiferData.__init__(self, model, kaq, Haq, c, z, npor, ltype, aqnumber)
        self.inhomlist = []
        self.area = 1e300  # Needed to find smallest inhom

    def initialize(self):
        # cause we are going to call initialize for inhoms
        AquiferData.initialize(self)
        for inhom in self.inhomlist:
            inhom.initialize()
        for inhom in self.inhomlist:
            inhom.create_elements()

    def add_inhom(self, inhom):
        self.inhomlist.append(inhom)
        return len(self.inhomlist) - 1  # returns number in the list

    def find_aquifer_data(self, x, y):
        rv = self
        for inhom in self.inhomlist:
            if inhom.isinside(x, y):
                if inhom.area < rv.area:
                    rv = inhom
        return rv
        # Not used anymore I think 5 Nov 2015
        # def find_aquifer_number(self, x, y):
        #    rv = -1
        #    for i,inhom in enumerate(self.inhomlist):
        #        if inhom.isinside(x, y):
        #            if inhom.area < rv.area:
        #                rv = i
        #    return rv

    def find_layer(self, z):
        '''
        If z in aquifer layer, returns number of aquifer layer
        If z in leaky layer, returns number of aquifer layer below it '''
        for i in range(self.Nlayers):
            if z > self.z[i + 1]:
                return self.aqnumber[i]

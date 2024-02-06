import inspect  # Used for storing the input

import numpy as np

from .aquifer_parameters import param_maq
from .constant import ConstantStar


class AquiferData:
    def __init__(self, model, kaq, c, z, npor, ltype):
        # All input variables except model should be numpy arrays
        # That should be checked outside this function
        self.model = model
        # Needed for heads
        self.kaq = np.atleast_1d(kaq)
        self.naq = len(kaq)
        self.c = np.atleast_1d(c)
        self.hstar = None
        # Needed for tracing
        self.z = np.atleast_1d(z)
        self.Hlayer = self.z[:-1] - self.z[1:]  # thickness of all layers
        self.nlayers = len(self.z) - 1
        self.npor = np.atleast_1d(npor)
        self.ltype = np.atleast_1d(ltype)
        # tag indicating whether an aquifer is Laplace (confined on top)
        if self.ltype[0] == "a":
            self.ilap = 1
        else:
            self.ilap = 0
        #
        self.area = 1e200  # Smaller than default of ml.aq so that inhom is found
        self.layernumber = np.zeros(self.nlayers, dtype="int")
        self.layernumber[self.ltype == "a"] = np.arange(self.naq)
        self.layernumber[self.ltype == "l"] = np.arange(self.nlayers - self.naq)
        if self.ltype[0] == "a":
            self.layernumber[
                self.ltype == "l"
            ] += 1  # first leaky layer below first aquifer layer
        self.zaqtop = self.z[:-1][self.ltype == "a"]
        self.zaqbot = self.z[1:][self.ltype == "a"]
        self.Haq = self.zaqtop - self.zaqbot
        self.T = self.kaq * self.Haq
        self.Tcol = self.T[:, np.newaxis]
        self.zlltop = self.z[:-1][self.ltype == "l"]
        self.zllbot = self.z[1:][self.ltype == "l"]
        if self.ltype[0] == "a":
            self.zlltop = np.hstack((self.z[0], self.zlltop))
            self.zllbot = np.hstack((self.z[0], self.zllbot))
        self.Hll = self.zlltop - self.zllbot
        self.nporaq = self.npor[self.ltype == "a"]
        if self.ltype[0] == "a":
            self.nporll = np.ones(len(self.npor[self.ltype == "l"]) + 1)
            self.nporll[1:] = self.npor[self.ltype == "l"]
        else:
            self.nporll = self.npor[self.ltype == "l"]

    def initialize(self):
        self.elementlist = []  # Elementlist of aquifer
        d0 = 1.0 / (self.c * self.T)
        d0[:-1] += 1.0 / (self.c[1:] * self.T[:-1])
        dp1 = -1.0 / (self.c[1:] * self.T[1:])
        dm1 = -1.0 / (self.c[1:] * self.T[:-1])
        A = np.diag(dm1, -1) + np.diag(d0, 0) + np.diag(dp1, 1)
        w, v = np.linalg.eig(A)
        # take the real part for the rare occasion that the eig 
        # function returns a complex answer with very small 
        # imaginary part
        w = w.real
        v = v.real
        # sort lab in decending order, hence w in ascending order
        index = np.argsort(abs(w))
        w = w[index]
        v = v[:, index]
        if self.ilap:
            self.lab = np.zeros(self.naq)
            self.lab[1:] = 1.0 / np.sqrt(w[1:])
            self.zeropluslab = (
                self.lab
            )  # to be deprecated when new lambda is fully implemented
            v[:, 0] = self.T / np.sum(self.T)  # first column is normalized T
        else:
            self.lab = 1.0 / np.sqrt(w)
        self.eigvec = v
        self.coef = np.linalg.solve(v, np.diag(np.ones(self.naq))).T

    def add_element(self, e):
        self.elementlist.append(e)
        if isinstance(e, ConstantStar):
            self.hstar = e.hstar

    def isinside(self, x, y):
        raise Exception("Must overload AquiferData.isinside()")

    def storeinput(self, frame):
        self.inputargs, _, _, self.inputvalues = inspect.getargvalues(frame)

    def findlayer(self, z):
        """
        Returns layer-number, layer-type and model-layer-number"""
        if z > self.z[0]:
            modellayer, ltype = -1, "above"
            layernumber = None
        elif z < self.z[-1]:
            modellayer, ltype = len(self.layernumber), "below"
            layernumber = None
        else:
            modellayer = np.argwhere((z <= self.z[:-1]) & (z >= self.z[1:]))[0, 0]
            layernumber = self.layernumber[modellayer]
            ltype = self.ltype[modellayer]
        return layernumber, ltype, modellayer


class Aquifer(AquiferData):
    def __init__(self, model, kaq, c, z, npor, ltype):
        AquiferData.__init__(self, model, kaq, c, z, npor, ltype)
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

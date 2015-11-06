import numpy as np
import inspect # Used for storing the input
from element import Element
from equation import HeadEquation, MscreenWellEquation
from scipy.special import k0, k1

class WellBase(Element):
    def __init__(self, model, xw=0, yw=0, Qw=100.0, rw=0.1, \
                 res=0.0, layers=0, name='WellBase', label=None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers,\
                         name=name, label=label)
        self.Nparam = len(self.pylayers) # Defined here and not in Element as other elements can have multiple parameters per layers
        self.xw = float(xw)
        self.yw = float(yw)
        self.Qw = np.atleast_1d(Qw)
        self.rw = float(rw)
        self.res = float(res)
        self.model.add_element(self)
    def __repr__(self):
        return self.name + ' at ' + str((self.xw, self.yw))
    def initialize(self):
        self.xc = np.array([self.xw + self.rw])
        self.yc = np.array([self.yw])
        self.Ncp = 1
        self.aq = self.model.aq.find_aquifer_data(self.xw, self.yw)
        self.aq.add_element(self)
        self.parameters = np.empty((self.Nparam, 1))
        self.parameters[:,0] = self.Qw
        self.resfac = self.res / (2 * np.pi * self.rw * self.aq.Haq[self.pylayers])
    def potinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nparam, aq.Naq))
        if aq == self.aq:
            pot = np.zeros(aq.Naq)
            r = np.sqrt( (x-self.xw)**2 + (y-self.yw)**2 )
            if r < self.rw: r = self.rw  # If at well, set to at radius
            if aq.ltype[0] == 'a':
                pot[0] = np.log(r / self.rw) / (2 * np.pi)
                pot[1:] = -k0(r / aq.lab) / (2 * np.pi)
            else:
                pot[:] = k0(r / aq/lab) / (2 * np.pi)
            rv[:] = self.aq.coef[self.pylayers] * pot
        return rv
    def disinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.Nparam, aq.Naq))
        if aq == self.aq:
            qx = np.zeros(aq.Naq)
            qy = np.zeros(aq.Naq)
            rsq = (x - self.xw)**2 + (y - self.yw)**2
            r = np.sqrt(rsq)
            xminxw = x - self.xw
            yminyw = y - self.yw
            if r < self.rw:
                r = self.rw
                rsq = self.rwsq
                xminxw = self.rw
                yminyw = 0.0
            if aq.ltype[0] == 'a':
                qx[0] = -1 / (2 * np.pi) * xminxw / rsq
                qy[0] = -1 / (2 * np.pi) * yminyw / rsq 
                kone = k1(r / aq.lab)
                qx[1:] = -kone * xminxw / (r * aq.lab) / (2 * np.pi)
                qy[1:] = -kone * yminyw / (r * aq.lab) / (2 * np.pi)
            else:
                kone = k1(r / aq.lab)
                qx[:] = -kone * xminxw / (r * aq.lab) / (2 * np.pi)
                qy[:] = -kone * yminyw / (r * aq.lab) / (2 * np.pi)
            rv[0] = self.aq.coef[self.pylayers] * qx
            rv[1] = self.aq.coef[self.pylayers] * qy   
        return rv
    def headinside(self):
        h = self.model.head(self.xw + self.rw, self.yw, layers=self.pylayers)
        return h - self.resfac * self.parameters[:,0]
    
class Well(WellBase, MscreenWellEquation):
    def __init__(self, model, xw=0, yw=0, Qw=100.0, rw=0.1, \
                 res=0.0, layers=0, label=None):
        self.storeinput(inspect.currentframe())
        WellBase.__init__(self, model, xw, yw, 0.0, rw, res,\
                          layers=layers, name='Well', label=label)
        self.Qc = float(Qw)
        if self.Nlayers == 1:
            self.Nunknowns = 0
        else:
            self.Nunknowns = self.Nparam
    def initialize(self):
        WellBase.initialize(self)
    def setparams(self, sol):
        self.parameters[:,0] = sol
    
class HeadWell(WellBase, HeadEquation):
    def __init__(self, model, xw=0, yw=0, hw=10.0, rw=0.1, \
                 res=0.0, layers=0, label=None):
        self.storeinput(inspect.currentframe())
        WellBase.__init__(self, model, xw, yw, 0.0, rw, res,\
                          layers=layers, name='HeadWell', label=label)
        self.hc = hw
        self.Nunknowns = self.Nparam
    def initialize(self):
        WellBase.initialize(self)
        self.pc = self.hc * self.aq.T[self.pylayers] # Needed in solving
    def setparams(self, sol):
        self.parameters[:,0] = sol
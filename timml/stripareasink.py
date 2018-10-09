import numpy as np
import inspect  # Used for storing the input
from .element import Element

__all__ = ['StripAreaSink']

class StripAreaSink(Element):
    def __init__(self, model, xleft=-1, xright=1, N=0.001, layer=0, name='StripAreaSink', label=None):
        Element.__init__(self, model, nparam=1, nunknowns=0, layers=layer, \
                         name=name, label=label)
        self.xleft = xleft
        self.xright = xright
        self.N = N
        self.model.add_element(self)

    def __repr__(self):
        return self.name + ' at ' + str((self.xc, self.yc))

    def initialize(self):
        self.xc = 0.5 * (self.xleft + self.xright)
        self.yc = 0
        self.L = self.xright - self.xleft
        self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        self.aq.add_element(self)
        self.parameters = np.array([[self.N]])
        if self.aq.ilap:
            self.lab = self.aq.lab[1:]
            self.A = -self.aq.coef[self.layers, 1:] * self.lab ** 2 / 2
            self.B = self.A * (np.exp(-self.L / self.lab) - 1)
            self.plabsq = self.aq.coef[self.layers, 1:] * self.lab ** 2
        else:
            print('StripAreaSink cannot be added to semi-confined system')

    def potinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        if aq == self.aq:
            if x < self.xleft:
                rv[0, 0] = -(self.xleft - self.xc) * (x - self.xc) + \
                           self.L ** 2 / 8
                rv[0, 1:] = self.B * np.exp((x - self.xleft) / self.lab)
            elif x > self.xright:
                rv[0, 0] = -(self.xright - self.xc) * (x - self.xc) + \
                           self.L ** 2 / 8
                rv[0, 1:] = self.B * np.exp(-(x - self.xright) / self.lab)
            else:
                rv[0, 0] = -0.5 * (x ** 2 - 2 * self.xc * x + self.xc ** 2)
                rv[0, 1:] = self.A * (np.exp(-(x - self.xleft) / self.lab) +
                                      np.exp((x - self.xright) / self.lab)) + \
                                      self.plabsq
        return rv
    
    def disvecinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq == self.aq:
            if x < self.xleft:
                rv[0, 0, 0] = self.xleft - self.xc
                rv[0, 0, 1:] = -self.B / self.lab * np.exp((x - self.xleft) / self.lab)
            elif x > self.xright:
                rv[0, 0, 0] = self.xright - self.xc
                rv[0, 0, 1:] = self.B / self.lab * np.exp(-(x - self.xright) / self.lab)
            else:
                rv[0, 0, 0] = x - self.xc
                rv[0, 0, 1:] = self.A / self.lab * (
                                   np.exp(-(x - self.xleft) / self.lab) -
                                   np.exp((x - self.xright) / self.lab)) 
        return rv
    
    def K1RI0r(self, rin):
        rv = np.zeros(len(self.lab))
        if self.islarge.any():
            index = (self.R - rin) / self.labbig < 10
            if index.any():
                r = rin / self.labbig[index]
                R = self.R / self.labbig[index]
                rv[self.islarge * index] = np.sqrt(1 / (4 * r * R)) * np.exp(r - R) * \
                                   (1 + 3 / (8 * R) - 15 / (128 * R ** 2) + 315 / (3072 * R ** 3)) * \
                                   (1 + 1 / (8 * r) +  9 / (128 * r ** 2) + 225 / (3072 * r ** 3))
        if ~self.islarge.any():
            index = (self.R - rin) / self.labsmall < 10
            if index.any():
                r = rin / self.labsmall[index]
                rv[~self.islarge * index] = self.k1Roverlab[index] * i0(r)
        return rv
    
    def K1RI1r(self, rin):
        rv = np.zeros(len(self.lab))
        if self.islarge.any():
            index = (self.R - rin) / self.labbig < 10
            if index.any():
                r = rin / self.labbig[index]
                R = self.R / self.labbig[index]
                rv[self.islarge * index] = np.sqrt(1 / (4 * r * R)) * np.exp(r - R) * \
                                   (1 + 3 / (8 * R) - 15 / (128 * R ** 2) + 315 / (3072 * R ** 3)) * \
                                   (1 - 3 / (8 * r) - 15 / (128 * r ** 2) - 315 / (3072 * r ** 3))
        if ~self.islarge.any():
            index = (self.R - rin) / self.labsmall < 10
            if index.any():
                r = rin / self.labsmall[index]
                rv[~self.islarge * index] = self.k1Roverlab[index] * i1(r)
        return rv
    
    def I1RK0r(self, rin):
        rv = np.zeros(len(self.lab))
        if self.islarge.any():
            index = (rin - self.R) / self.labbig < 10
            if index.any():
                r = rin / self.labbig[index]
                R = self.R / self.labbig[index]
                rv[self.islarge * index] = np.sqrt(1 / (4 * r * R)) * np.exp(R - r) * \
                                   (1 - 3 / (8 * R) - 15 / (128 * R ** 2) - 315 / (3072 * R ** 3)) * \
                                   (1 - 1 / (8 * r) +  9 / (128 * r ** 2) - 225 / (3072 * r ** 3))
        if ~self.islarge.any():
            index = (self.R - rin) / self.labsmall < 10
            if index.any():
                r = rin / self.labsmall[index]
                rv[~self.islarge * index] = self.i1Roverlab[index] * k0(r)
        return rv
    
    def I1RK1r(self, rin):
        rv = np.zeros(len(self.lab))
        if self.islarge.any():
            index = (rin - self.R) / self.labbig < 10
            if index.any():
                r = rin / self.labbig[index]
                R = self.R / self.labbig[index]
                rv[self.islarge * index] = np.sqrt(1 / (4 * r * R)) * np.exp(R - r) * \
                                   (1 - 3 / (8 * R) - 15 / (128 * R ** 2) - 315 / (3072 * R ** 3)) * \
                                   (1 + 3 / (8 * r) - 15 / (128 * r ** 2) + 315 / (3072 * r ** 3))
        if ~self.islarge.any():
            index = (self.R - rin) / self.labsmall < 10
            if index.any():
                r = rin / self.labsmall[index]
                rv[~self.islarge * index] = self.i1Roverlab[index] * k1(r)
        return rv
    
    def qztop(self, x, y):
        rv = 0.0
        if np.sqrt((x - self.xc) ** 2 + (y - self.yc) ** 2) <= self.R:
            rv = -self.parameters[0, 0]  # minus cause the parameter is the infiltration rate
        return rv

    def changetrace(self, xyzt1, xyzt2, aq, layer, ltype, modellayer, direction, hstepmax):
        changed = False
        terminate = False
        xyztnew = 0
        eps = 1e-8
        r1sq = (xyzt1[0] - self.xc) ** 2 + (xyzt1[1] - self.yc) ** 2
        r2sq = (xyzt2[0] - self.xc) ** 2 + (xyzt2[1] - self.yc) ** 2
        if (r1sq < self.Rsq and r2sq > self.Rsq ) or (r1sq > self.Rsq and r2sq < self.Rsq):
            changed = True
            x1, y1 = xyzt1[0:2]
            x2, y2 = xyzt2[0:2]
            a = (x2 - x1) ** 2 + (y2 - y1) ** 2
            b = 2 * ((x2 - x1) * (x1 - self.xc) + (y2 - y1) * (y1 - self.yc))
            c = self.xc ** 2 + self.yc ** 2 + x1 ** 2 + y1 ** 2 - 2 * (self.xc * x1 +self.yc * y1) - self.Rsq
            u1 = (-b - np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
            u2 = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
            if u1 > 0:
                u = u1 * (1.0 + eps) # Go just beyond circle
            else:
                u = u2 * (1.0 + eps) # Go just beyond circle
            xn = x1 + u * (x2 - x1)
            yn = y1 + u * (y2 - y1)
            zn = xyzt1[2] + u * (xyzt2[2] - xyzt1[2])
            xyztnew = xyzt1 + u * (xyzt2 - xyzt1)
        return changed, terminate, xyztnew

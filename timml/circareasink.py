import numpy as np
import inspect  # Used for storing the input
from .element import Element
from scipy.special import k0, k1, i0, i1

__all__ = ['CircAreaSink']

class CircAreaSink(Element):
    def __init__(self, model, xc=0, yc=0, R=1, N=0.001, layer=0, name='CircAreasink', label=None):
        Element.__init__(self, model, nparam=1, nunknowns=0, layers=layer, \
                         name=name, label=label)
        self.xc = xc
        self.yc = yc
        self.R = R
        self.N = N
        self.model.add_element(self)

    def __repr__(self):
        return self.name + ' at ' + str((self.xc, self.yc))

    def initialize(self):
        self.Rsq = self.R ** 2
        self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        self.aq.add_element(self)
        self.parameters = np.array([[self.N]])
        self.Rlarge = 500.0  # If R/lab > Rlarge, then we use asymptotic approximation to compute potential
        if self.aq.ilap:
            self.lab = self.aq.lab[1:]
            self.A = -self.aq.coef[self.layers, 1:] * self.R * self.lab
            self.B = self.aq.coef[self.layers, 1:] * self.R * self.lab
            self.C = self.aq.coef[self.layers, 1:] * self.lab ** 2
        else:
            self.lab = self.aq.lab
            self.A = -self.aq.coef[self.layers] * self.R * self.lab
            self.B = self.aq.coef[self.layers] * self.R * self.lab
            self.C = self.aq.coef[self.layers] * self.lab ** 2
        self.islarge = self.R / self.lab > self.Rlarge
        self.labsmall = self.lab[~self.islarge]
        self.labbig = self.lab[self.islarge]
        self.k1Roverlab = k1(self.R / self.labsmall)
        self.i1Roverlab = i1(self.R / self.labsmall)

    def potinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        if aq == self.aq:
            r = np.sqrt((x - self.xc) ** 2 + (y - self.yc) ** 2)
            if r <= self.R:
                if aq.ilap:
                    rv[0, 0] = 0.25 * (self.Rsq - r ** 2)
                    rv[0, 1:] = self.A * self.K1RI0r(r) + self.C
                else:
                    rv[0] = self.A * self.K1RI0r(r) + self.C
            else:
                if aq.ilap:
                    rv[0, 0] = -0.5 * self.Rsq * np.log(r / self.R)
                    rv[0, 1:] = self.B * self.I1RK0r(r)
                else:
                    rv[0] = self.B * self.I1RK0r(r)
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq == self.aq:
            r = np.sqrt((x - self.xc) ** 2 + (y - self.yc) ** 2)
            if r <= self.R:
                if r > 1e-12:  # otherwise zero
                    if aq.ilap:
                        rv[0, 0, 0] = (x - self.xc) / 2
                        rv[1, 0, 0] = (y - self.yc) / 2
                        K1RI1r = self.K1RI1r(r)
                        rv[0, 0, 1:] = -self.A * K1RI1r * (x - self.xc) / (r * self.lab)
                        rv[1, 0, 1:] = -self.A * K1RI1r * (y - self.yc) / (r * self.lab)
                    else:
                        K1RI1r = self.K1RI1r(r) / self.lab
                        rv[0, 0] = -self.A * K1RI1r * (x - self.xc) / r
                        rv[1, 0] = -self.A * K1RI1r * (y - self.yc) / r
            else:
                if aq.ilap:
                    rv[0, 0, 0] = 0.5 * self.Rsq * (x - self.xc) / r ** 2
                    rv[1, 0, 0] = 0.5 * self.Rsq * (y - self.yc) / r ** 2
                    I1RK1r = self.I1RK1r(r)
                    rv[0, 0, 1:] = self.B * I1RK1r * (x - self.xc) / (r * self.lab)
                    rv[1, 0, 1:] = self.B * I1RK1r * (y - self.yc) / (r * self.lab)
                else:
                    I1RK1r = self.I1RK1r(r)
                    rv[0, 0] = self.B * I1RK1r * (x - self.xc) / (r * self.lab)
                    rv[1, 0] = self.B * I1RK1r * (y - self.yc) / (r * self.lab)
        else:
            r = np.sqrt((x - self.xc) ** 2 + (y - self.yc) ** 2)
            if r <= self.R:
                raise 'CircAreaSink should add constant inside inhomogeneity'
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
        message = None
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
        return changed, terminate, xyztnew, message

from __future__ import division
import numpy as np
import inspect  # Used for storing the input
from element import Element
from scipy.special import k0, k1, i0, i1

class CircAreaSink(Element):
    def __init__(self, model, xc=0, yc=0, R=1, N=0.001, layer=0, name='CircAreasink', label=None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layer, \
                         name=name, label=label)
        self.xc = float(xc)
        self.yc = float(yc)
        self.R = float(R)
        self.N = float(N)
        self.model.add_element(self)

    def __repr__(self):
        return self.name + ' at ' + str((self.xc, self.yc))

    def initialize(self):
        self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        self.aq.add_element(self)
        self.parameters = np.array([[self.N]])
        self.Rlarge = 500.0  # If R/lab > Rlarge, then we use asymptotic approximation to compute potential
        if self.aq.ilap:
            self.lab = self.aq.lab[1:]
            self.A = -self.aq.coef[self.pylayers, 1:] * self.R * self.lab
            self.B =  self.aq.coef[self.pylayers, 1:] * self.R * self.lab
            self.C =  self.aq.coef[self.pylayers, 1:] * self.lab ** 2
        else:
            self.lab = self.aq.lab
            self.A = -self.aq.coef[self.pylayers] * self.R * self.lab
            self.B =  self.aq.coef[self.pylayers] * self.R * self.lab
            self.C =  self.aq.coef[self.pylayers] * self.lab ** 2
        self.islarge = self.R / self.lab > self.Rlarge
        self.labsmall = self.lab[~self.islarge]
        self.labbig = self.lab[self.islarge]
        self.k1Roverlab = k1(self.R / self.labsmall)
        self.i1Roverlab = i1(self.R / self.labsmall)

    def potinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nparam, aq.Naq))
        if aq == self.aq:
            r = np.sqrt((x - self.xc) ** 2 + (y - self.yc) ** 2)
            if r <= self.R:
                if aq.ilap:
                    rv[0, 0] = 0.25 * (self.R ** 2 - r ** 2)
                    rv[0, 1:] = self.A * self.K1RI0r(r) + self.C
                else:
                    rv[0] = self.A * self.K1RI0r(r) + self.C
            else:
                if aq.ilap:
                    rv[0, 0] = -0.5 * self.R ** 2 * np.log(r / self.R)
                    rv[0, 1:] = self.B * self.I1RK0r(r)
                else:
                    rv[0] = self.B * self.I1RK0r(r)
        return rv
    
    def disinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.Nparam, aq.Naq))
        if aq == self.aq:
            r = np.sqrt((x - self.xc) ** 2 + (y - self.yc) ** 2)
            if r <= self.R:
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
                    rv[0, 0, 0] = 0.5 * self.R ** 2 * (x - self.xc) / r ** 2
                    rv[1, 0, 0] = 0.5 * self.R ** 2 * (y - self.yc) / r ** 2
                    I1RK1r = self.K1RI1r(r)
                    rv[0, 0, 1:] = self.B * I1RK1r * (x - self.xc) / (r * self.lab)
                    rv[1, 0, 1:] = self.B * I1RK1r * (y - self.yc) / (r * self.lab)
                else:
                    I1RK1r = self.K1RI1r(r)
                    rv[0, 0] = self.B * I1RK1r * (x - self.xc) / (r * self.lab)
                    rv[1, 0] = self.B * I1RK1r * (y - self.yc) / (r * self.lab)
        else:
            r = np.sqrt((x - self.xc) ** 2 + (y - self.yc) ** 2)
            if r <= self.R:
                raise 'CircAreaSink should add constant inside inhomogeneity'
        return rv

    
    def dischargeInfluence(self,aq,x,y,z=0,t=0):
        rvx = zeros((1,aq.Naquifers),'d'); rvy = zeros((1,aq.Naquifers),'d')
        rsq = (x-self.xp)**2 + (y-self.yp)**2
        r = sqrt(rsq)
        if r < 1e-6:  # Else evaluation at center blows up
            r = 1e-6
            rsq = r**2
        if r <= self.Rp:
            rvx[0,0] = (x - self.xp) / 2.0
            rvy[0,0] = (y - self.yp) / 2.0
        else:
            rvx[0,0] = self.Rp**2 * (x-self.xp) / (2.0 * rsq)
            rvy[0,0] = self.Rp**2 * (y-self.yp) / (2.0 * rsq)
        if self.aquiferParent == aq:
            for i in range(self.aquiferParent.Naquifers - 1):
                rolab = r / self.aquiferParent.lab[i]; Rolab = self.Rp / self.aquiferParent.lab[i]
                if r <= self.Rp:
                    if self.isLargeCircle[i] == 0:
                        dis = self.K1Rolab[i] * scipy.special.i1(rolab) / self.aquiferParent.lab[i]
                        rvx[0,i+1] = dis * (x-self.xp) / r
                        rvy[0,i+1] = dis * (y-self.yp) / r
                    else:
                        if (Rolab - rolab) < 10.0: # zero after 10 lambda
                            K1RI1r = sqrt( 1./(4.0 * rolab * Rolab) ) * exp( rolab - Rolab ) * \
                                (1. + 3./(8*Rolab) - 15./(128*Rolab**2) + 315./(3072*Rolab**3) ) * \
                                (1. - 3./(8*rolab) - 15./(128*rolab**2) - 315/(3072*rolab**3) )
                            rvx[0,i+1] = K1RI1r * (x-self.xp) / ( r * self.aquiferParent.lab[i] )
                            rvy[0,i+1] = K1RI1r * (y-self.yp) / ( r * self.aquiferParent.lab[i] )
                else:
                    if self.isLargeCircle[i] == 0:                    
                        dis = self.I1Rolab[i] * scipy.special.k1(rolab) / self.aquiferParent.lab[i]
                        rvx[0,i+1] = dis * (x-self.xp) / r
                        rvy[0,i+1] = dis * (y-self.yp) / r
                    else:
                        if (rolab - Rolab) < 10.0: # zero after 10 lambda
                            I1RK1r = sqrt( 1./(4.0 * rolab * Rolab) ) * exp( Rolab - rolab ) * \
                                (1. - 3./(8*Rolab) - 15./(128*Rolab**2) - 315./(3072*Rolab**3) ) * \
                                (1. + 3./(8*rolab) - 15./(128*rolab**2) + 315./(3072*rolab**3) )
                            rvx[0,i+1] = I1RK1r * (x-self.xp) / ( r * self.aquiferParent.lab[i] )
                            rvy[0,i+1] = I1RK1r * (y-self.yp) / ( r * self.aquiferParent.lab[i] )
            # Only if aq == aquiferParent, otherwise not
            rvx = self.coef * rvx
            rvy = self.coef * rvy
        return [rvx,rvy]
    
    def K1RI0r(self, rin):
        rv = np.zeros(len(self.lab))
        if self.islarge.any():
            index = (self.R - rin) / self.labbig < 10
            if index.any():
                r = rin / self.labbig[index]
                R = self.R / self.labbig[index]
                rv[self.islarge][index] = np.sqrt(1 / (4 * r * R)) * np.exp(r - R) * \
                                   (1 + 3 / (8 * R) - 15 / (128 * R ** 2) + 315 / (3072 * R ** 3)) * \
                                   (1 + 1 / (8 * r) +  9 / (128 * r ** 2) + 225 / (3072 * r ** 3))
        if ~self.islarge.any():
            index = (self.R - rin) / self.labsmall < 10
            if index.any():
                r = rin / self.labsmall[index]
                rv[~self.islarge][index] = self.k1Roverlab[index] * i0(r)
        return rv
    
    def K1RI1r(self, rin):
        rv = np.zeros(len(self.lab))
        if self.islarge.any():
            index = (self.R - rin) / self.labbig < 10
            if index.any():
                r = rin / self.labbig[index]
                R = self.R / self.labbig[index]
                rv[self.islarge][index] = np.sqrt(1 / (4 * r * R)) * np.exp(r - R) * \
                                   (1 + 3 / (8 * R) - 15 / (128 * R ** 2) + 315 / (3072 * R ** 3)) * \
                                   (1 - 3 / (8 * r) - 15 / (128 * r ** 2) - 315 / (3072 * r ** 3))
        if ~self.islarge.any():
            index = (self.R - rin) / self.labsmall < 10
            if index.any():
                r = rin / self.labsmall[index]
                rv[~self.islarge][index] = self.k1Roverlab[index] * i1(r)
        return rv
    
    def I1RK0r(self, rin):
        rv = np.zeros(len(self.lab))
        if self.islarge.any():
            index = (rin - self.R) / self.labbig < 10
            if index.any():
                r = rin / self.labbig[index]
                R = self.R / self.labbig[index]
                rv[self.islarge][index] = np.sqrt(1 / (4 * r * R)) * np.exp(R - r) * \
                                   (1 - 3 / (8 * R) - 15 / (128 * R ** 2) - 315 / (3072 * R ** 3)) * \
                                   (1 - 1 / (8 * r) +  9 / (128 * r ** 2) - 225 / (3072 * r ** 3))
        if ~self.islarge.any():
            index = (self.R - rin) / self.labsmall < 10
            if index.any():
                r = rin / self.labsmall[index]
                rv[~self.islarge][index] = self.i1Roverlab[index] * k0(r)
        return rv
    
    def I1RK1r(self, rin):
        rv = np.zeros(len(self.lab))
        if self.islarge.any():
            index = (rin - self.R) / self.labbig < 10
            if index.any():
                r = rin / self.labbig[index]
                R = self.R / self.labbig[index]
                rv[self.islarge][index] = np.sqrt(1 / (4 * r * R)) * np.exp(R - r) * \
                                   (1 - 3 / (8 * R) - 15 / (128 * R ** 2) - 315 / (3072 * R ** 3)) * \
                                   (1 + 3 / (8 * r) - 15 / (128 * r ** 2) + 315 / (3072 * r ** 3))
        if ~self.islarge.any():
            index = (self.R - rin) / self.labsmall < 10
            if index.any():
                r = rin / self.labsmall[index]
                rv[~self.islarge][index] = self.i1Roverlab[index] * k1(r)
        return rv

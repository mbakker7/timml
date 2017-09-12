import numpy as np
from .linesink import LineSinkHoBase
from .equation import HeadDiffEquation2, DisvecDiffEquation2
from .controlpoints import controlpoints
# needed for testing
#from .equation import DisvecDiffEquation


class IntHeadDiffLineSink(LineSinkHoBase, HeadDiffEquation2):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 order=0, ndeg=3, layers=0, label=None, addtomodel=True,
                 aq=None, aqin=None, aqout=None):
        LineSinkHoBase.__init__(self, model, x1, y1, x2, y2, Qls=0, \
                                layers=list(range(model.aq.Naq)), order=order, \
                                name='IntHeadDiffLineSink', label=label, \
                                addtomodel=addtomodel, aq=aq)
        self.inhomelement = True
        self.ndeg = ndeg
        self.Xleg, self.wleg = np.polynomial.legendre.leggauss(self.ndeg)
        self.nunknowns = self.nparam
        self.aqin = aqin
        self.aqout = aqout

    def initialize(self):
        LineSinkHoBase.initialize(self)
        self.xcin, self.ycin = controlpoints(self.ncp - 1, self.z1, self.z2,
                                             eps=1e-6, include_ends=True)
        self.xcout, self.ycout = controlpoints(self.ncp - 1, self.z1, self.z2,
                                               eps=-1e-6, include_ends=True)
        if self.aqin is None:
            self.aqin = self.model.aq.find_aquifer_data(self.xcin[1],
                                                        self.ycin[1])
        if self.aqout is None:
            self.aqout = self.model.aq.find_aquifer_data(self.xcout[1],
                                                         self.ycout[1])

    def setparams(self, sol):
        self.parameters[:, 0] = sol
        
    def changetrace(self, xyzt1, xyzt2, layer, ltype):
        changed = False
        xyztnew = 0
        if (ltype == 'a'):
            eps = 1e-8
            za = xyzt1[0] + xyzt1[1] * 1j
            zb = xyzt2[0] + xyzt2[1] * 1j
            Za = (2 * za - (self.z1 + self.z2)) / (self.z2 - self.z1)
            Zb = (2 * zb - (self.z1 + self.z2)) / (self.z2 - self.z1)
            if Za.imag * Zb.imag < 0:
                Xa, Ya = Za.real, Za.imag
                Xb, Yb = Zb.real, Zb.imag
                X = Xa - Ya * (Xb - Xa) / (Yb - Ya)
                if abs(X) <= 1:
                    Z = X + eps * np.sign(Yb) * 1j  # steps to side of Yb
                    z = 0.5 * ((self.z2 - self.z1) * Z + self.z1 + self.z2)
                    xnew, ynew = z.real, z.imag
                    dold = abs(za - zb)
                    dnew = np.sqrt((xnew - xyzt1[0]) ** 2 + (ynew - xyzt1[1]) ** 2)
                    znew = xyzt1[2] + dnew / dold * (xyzt2[2] - xyzt1[2])
                    tnew = xyzt1[3] + dnew / dold * (xyzt2[3] - xyzt1[3])
                    xyztnew = np.array([xnew, ynew, znew, tnew])
                    return True, xyztnew
        return changed, xyztnew


class IntFluxDiffLineSink(LineSinkHoBase, DisvecDiffEquation2):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 order=0, ndeg=3, label=None, addtomodel=True, aq=None,
                 aqin=None, aqout=None):
        LineSinkHoBase.__init__(self, model, x1, y1, x2, y2, Qls=0, \
                                layers=list(range(model.aq.Naq)), order=order,
                                name='IntFluxDiffLineSink', label=label, \
                                addtomodel=addtomodel, aq=aq)
        self.inhomelement = True
        self.ndeg = ndeg
        self.Xleg, self.wleg = np.polynomial.legendre.leggauss(self.ndeg)
        self.nunknowns = self.nparam
        self.aqin = aqin
        self.aqout = aqout

    def initialize(self):
        LineSinkHoBase.initialize(self)
        self.xcin, self.ycin = controlpoints(self.ncp - 1, self.z1, self.z2,
                                             eps=1e-6, include_ends=True)
        self.xcout, self.ycout = controlpoints(self.ncp - 1, self.z1, self.z2,
                                               eps=-1e-6, include_ends=True)
        if self.aqin is None:
            self.aqin = self.model.aq.find_aquifer_data(self.xcin[0],
                                                        self.ycin[0])
        if self.aqout is None:
            self.aqout = self.model.aq.find_aquifer_data(self.xcout[0],
                                                         self.ycout[0])

    def setparams(self, sol):
        self.parameters[:, 0] = sol
        
    def changetrace(self, xyzt1, xyzt2, layer, ltype):
        changed = False
        xyztnew = 0
        if (ltype == 'a'):
            eps = 1e-8
            za = xyzt1[0] + xyzt1[1] * 1j
            zb = xyzt2[0] + xyzt2[1] * 1j
            Za = (2 * za - (self.z1 + self.z2)) / (self.z2 - self.z1)
            Zb = (2 * zb - (self.z1 + self.z2)) / (self.z2 - self.z1)
            if Za.imag * Zb.imag < 0:
                Xa, Ya = Za.real, Za.imag
                Xb, Yb = Zb.real, Zb.imag
                X = Xa - Ya * (Xb - Xa) / (Yb - Ya)
                if abs(X) <= 1:
                    Z = X + eps * np.sign(Yb) * 1j  # steps to side of Yb
                    z = 0.5 * ((self.z2 - self.z1) * Z + self.z1 + self.z2)
                    xnew, ynew = z.real, z.imag
                    dold = abs(za - zb)
                    dnew = np.sqrt((xnew - xyzt1[0]) ** 2 + (ynew - xyzt1[1]) ** 2)
                    znew = xyzt1[2] + dnew / dold * (xyzt2[2] - xyzt1[2])
                    tnew = xyzt1[3] + dnew / dold * (xyzt2[3] - xyzt1[3])
                    xyztnew = np.array([xnew, ynew, znew, tnew])
                    return True, xyztnew
        return changed, xyztnew

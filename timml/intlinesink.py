"""Integrated line-sink elements.

Provides higher-order solutions and integrated conditions along line-sink elements. 
Used by inhomogeneities.
"""
import numpy as np

from .controlpoints import controlpoints
from .equation import (
    DisvecDiffEquation2,
    HeadDiffEquation2,
    IntDisVecEquation,
    IntLeakyWallEquation,
)
from .linesink import LineSinkHoBase


class IntHeadDiffLineSink(LineSinkHoBase, HeadDiffEquation2):
    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        order=0,
        ndeg=3,
        layers=None,
        label=None,
        addtomodel=True,
        aq=None,
        aqin=None,
        aqout=None,
    ):
        if layers is None:
            layers = np.arange(model.aq.naq)
        LineSinkHoBase.__init__(
            self,
            model,
            x1,
            y1,
            x2,
            y2,
            Qls=0,
            layers=layers,
            order=order,
            name="IntHeadDiffLineSink",
            label=label,
            addtomodel=addtomodel,
            aq=aq,
        )
        self.inhomelement = True
        self.ndeg = ndeg
        self.Xleg, self.wleg = np.polynomial.legendre.leggauss(self.ndeg)
        self.nunknowns = self.nparam
        self.aqin = aqin
        self.aqout = aqout

    def initialize(self):
        LineSinkHoBase.initialize(self)
        self.xcin, self.ycin = controlpoints(
            self.ncp - 1, self.z1, self.z2, eps=1e-6, include_ends=True
        )
        self.xcout, self.ycout = controlpoints(
            self.ncp - 1, self.z1, self.z2, eps=-1e-6, include_ends=True
        )
        if self.aqin is None:
            self.aqin = self.model.aq.find_aquifer_data(self.xcin[1], self.ycin[1])
        if self.aqout is None:
            self.aqout = self.model.aq.find_aquifer_data(self.xcout[1], self.ycout[1])

    def setparams(self, sol):
        self.parameters[:, 0] = sol

    # def changetrace(self, xyzt1, xyzt2, layer, ltype):
    def changetrace(
        self, xyzt1, xyzt2, aq, layer, ltype, modellayer, direction, hstepmax
    ):
        changed = False
        terminate = False
        xyztnew = 0
        message = None
        if ltype == "a":
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
                    changed = True
        return changed, terminate, [xyztnew], message


class IntFluxDiffLineSink(LineSinkHoBase, DisvecDiffEquation2):
    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        layers=None,
        order=0,
        ndeg=3,
        label=None,
        addtomodel=True,
        aq=None,
        aqin=None,
        aqout=None,
    ):
        if layers is None:
            layers = list(range(model.aq.naq))
        LineSinkHoBase.__init__(
            self,
            model,
            x1,
            y1,
            x2,
            y2,
            Qls=0,
            layers=layers,
            order=order,
            name="IntFluxDiffLineSink",
            label=label,
            addtomodel=addtomodel,
            aq=aq,
        )
        self.inhomelement = True
        self.ndeg = ndeg
        self.Xleg, self.wleg = np.polynomial.legendre.leggauss(self.ndeg)
        self.nunknowns = self.nparam
        self.aqin = aqin
        self.aqout = aqout

    def initialize(self):
        LineSinkHoBase.initialize(self)
        self.xcin, self.ycin = controlpoints(
            self.ncp - 1, self.z1, self.z2, eps=1e-6, include_ends=True
        )
        self.xcout, self.ycout = controlpoints(
            self.ncp - 1, self.z1, self.z2, eps=-1e-6, include_ends=True
        )
        if self.aqin is None:
            self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], self.ycin[0])
        if self.aqout is None:
            self.aqout = self.model.aq.find_aquifer_data(self.xcout[0], self.ycout[0])

    def setparams(self, sol):
        self.parameters[:, 0] = sol

    # def changetrace(self, xyzt1, xyzt2, layer, ltype):
    def changetrace(
        self, xyzt1, xyzt2, aq, layer, ltype, modellayer, direction, hstepmax
    ):
        changed = False
        terminate = False
        xyztnew = 0
        message = None
        if ltype == "a":
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
                    changed = True
                    # return True, False, xyztnew
        return changed, terminate, [xyztnew], message


class IntFluxLineSink(LineSinkHoBase, IntDisVecEquation):
    """Element to set numerically integrated flux along linesink to 0.

    Used in BuildingPit element.
    """

    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        layers=None,
        order=0,
        ndeg=3,
        label=None,
        addtomodel=True,
        aq=None,
        aqin=None,
        aqout=None,
    ):
        if layers is None:
            layers = list(range(model.aq.naq))
        LineSinkHoBase.__init__(
            self,
            model,
            x1,
            y1,
            x2,
            y2,
            Qls=0,
            layers=layers,
            order=order,
            name="IntFluxDiffLineSink",
            label=label,
            addtomodel=addtomodel,
            aq=aq,
        )
        self.inhomelement = True
        self.ndeg = ndeg
        self.Xleg, self.wleg = np.polynomial.legendre.leggauss(self.ndeg)
        self.nunknowns = self.nparam
        self.aqin = aqin
        self.aqout = aqout

    def initialize(self):
        LineSinkHoBase.initialize(self)

        # recalculated with ncp - 1 points instead of ncp points
        self.xcin, self.ycin = controlpoints(
            self.ncp - 1, self.z1, self.z2, eps=1e-6, include_ends=True
        )
        self.xcout, self.ycout = controlpoints(
            self.ncp - 1, self.z1, self.z2, eps=-1e-6, include_ends=True
        )

        if self.aqin is None:
            self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], self.ycin[0])
        if self.aqout is None:
            self.aqout = self.model.aq.find_aquifer_data(self.xcout[0], self.ycout[0])
        # set control points x,y depending on which side of
        # the boundary the element is for
        if self.aq == self.aqin:
            self.xc = self.xcin
            self.yc = self.ycin
        if self.aq == self.aqout:
            self.xc = self.xcout
            self.yc = self.ycout

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class LeakyIntHeadDiffLineSink(LineSinkHoBase, IntLeakyWallEquation):
    """Element to set numerically integrated head along linesink.

    Set integrated head equal to:

        Qnormal = H * (headin - headout) / res

    Used in LeakyBuildingPit element.
    """

    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        res=np.inf,
        layers=None,
        order=0,
        ndeg=3,
        label=None,
        addtomodel=True,
        aq=None,
        aqin=None,
        aqout=None,
    ):
        if layers is None:
            layers = list(range(model.aq.naq))
        LineSinkHoBase.__init__(
            self,
            model,
            x1,
            y1,
            x2,
            y2,
            Qls=0,
            layers=layers,
            order=order,
            name="LeakyIntHeadDiffLineSink",
            label=label,
            addtomodel=addtomodel,
            aq=aq,
        )
        self.res = res
        self.inhomelement = True
        self.ndeg = ndeg
        self.Xleg, self.wleg = np.polynomial.legendre.leggauss(self.ndeg)
        self.nunknowns = self.nparam
        self.aqin = aqin
        self.aqout = aqout

    def initialize(self):
        LineSinkHoBase.initialize(self)

        # recalculated with ncp - 1 points instead of ncp points
        self.xcin, self.ycin = controlpoints(
            self.ncp - 1, self.z1, self.z2, eps=1e-6, include_ends=True
        )
        self.xcout, self.ycout = controlpoints(
            self.ncp - 1, self.z1, self.z2, eps=-1e-6, include_ends=True
        )

        if self.aqin is None:
            self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], self.ycin[0])
        if self.aqout is None:
            self.aqout = self.model.aq.find_aquifer_data(self.xcout[0], self.ycout[0])
        # set control points x,y depending on which side of
        # the boundary the element is for
        if self.aq == self.aqin:
            self.xc = self.xcin
            self.yc = self.ycin
        if self.aq == self.aqout:
            self.xc = self.xcout
            self.yc = self.ycout

        # set resistance factor
        self.resfac = np.atleast_2d(self.aq.Haq[self.layers] / self.res).T

    def setparams(self, sol):
        self.parameters[:, 0] = sol

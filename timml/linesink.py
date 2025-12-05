"""Line-sink elements.

Implements line-sinks and line-sink strings for head- or discharge-specified features.

Example::

    HeadLineSink(ml, x1=-10, y1=0, x2=10, y2=0, hls=5, layers=0)
"""

import inspect  # Used for storing the input
import warnings

import matplotlib.pyplot as plt
import numpy as np

from . import bessel
from .controlpoints import controlpoints, strengthinf_controlpoints
from .element import Element
from .equation import HeadEquation

__all__ = [
    "HeadLineSink",
    "HeadLineSinkString",
    "Ditch",
    "DitchString",
    "LineSinkDitch",  # Deprecated alias for Ditch
    "LineSinkDitchString",  # Deprecated alias for DitchString
    "CollectorWell",
    "RadialCollectorWell",
]


class LineSinkChangeTrace:
    def changetrace(
        self,
        xyzt1,
        xyzt2,
        aq,
        layer,
        ltype,
        modellayer,
        direction,
        hstepmax,
        verbose=False,
    ):
        changed = False
        terminate = False
        xyztnew = 0
        message = None
        if ltype == "a":
            if True:
                # if (layer == self.layers).any():# in layer where line-sink is screened
                # not needed anymore, I thin this is all taken care of with
                # checking Qn1 and Qn2
                if verbose:
                    print("hello changetrace")
                    print("xyz1:", xyzt1[:-1])
                    print("xyz2:", xyzt2[:-1])
                x1, y1, z1, t1 = xyzt1
                x2, y2, z2, t2 = xyzt2
                eps = 1e-8
                za = x1 + y1 * 1j
                zb = x2 + y2 * 1j
                Za = (2 * za - (self.z1 + self.z2)) / (self.z2 - self.z1)
                Zb = (2 * zb - (self.z1 + self.z2)) / (self.z2 - self.z1)
                if verbose:
                    print("Za", Za)
                    print("Zb", Zb)
                if Za.imag * Zb.imag < 0:
                    Xa, Ya = Za.real, Za.imag
                    Xb, Yb = Zb.real, Zb.imag
                    X = Xa - Ya * (Xb - Xa) / (Yb - Ya)
                    if verbose:
                        print("X", X)
                    if abs(X) <= 1:  # crosses line-sink
                        if verbose:
                            print("crosses line-sink")
                        Znew1 = X - eps * np.sign(Yb) * 1j  # steps to side of Ya
                        Znew2 = X + eps * np.sign(Yb) * 1j  # steps to side of Yb
                        znew1 = 0.5 * ((self.z2 - self.z1) * Znew1 + self.z1 + self.z2)
                        znew2 = 0.5 * ((self.z2 - self.z1) * Znew2 + self.z1 + self.z2)
                        xnew1, ynew1 = znew1.real, znew1.imag
                        xnew2, ynew2 = znew2.real, znew2.imag
                        if Ya < 0:
                            theta = self.theta_norm_out
                        else:
                            theta = self.theta_norm_out + np.pi
                        Qx1, Qy1 = self.model.disvec(xnew1, ynew1)[:, layer] * direction
                        Qn1 = Qx1 * np.cos(theta) + Qy1 * np.sin(theta)
                        Qx2, Qy2 = self.model.disvec(xnew2, ynew2)[:, layer] * direction
                        Qn2 = Qx2 * np.cos(theta) + Qy2 * np.sin(theta)
                        if verbose:
                            print("xnew1, ynew1:", xnew1, ynew1)
                            print("xnew2, ynew2:", xnew2, ynew2)
                            print("Qn1, Qn2", Qn1, Qn2)
                            print("Qn2 > Qn1:", Qn2 > Qn1)
                        if Qn1 < 0:
                            # trying to cross line-sink that infiltrates, stay on
                            # bottom, don't terminate
                            if verbose:
                                print("change 1")
                            xnew = xnew1
                            ynew = ynew1
                            dold = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                            dnew = np.sqrt((x1 - xnew) ** 2 + (y1 - ynew) ** 2)
                            znew = z1 + dnew / dold * (z2 - z1)
                            tnew = t1 + dnew / dold * (t2 - t1)
                            changed = True
                            xyztnew = [np.array([xnew, ynew, znew, tnew])]
                        elif Qn2 < 0:  # all water is taken out, terminate
                            if verbose:
                                print("change 2")
                            xnew = xnew2
                            ynew = ynew2
                            dold = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                            dnew = np.sqrt((x1 - xnew) ** 2 + (y1 - ynew) ** 2)
                            znew = z1 + dnew / dold * (z2 - z1)
                            tnew = t1 + dnew / dold * (t2 - t1)
                            changed = True
                            terminate = True
                            xyztnew = [np.array([xnew, ynew, znew, tnew])]
                        elif Qn2 > Qn1:  # line-sink infiltrates
                            if verbose:
                                print("change 3")
                            xnew = xnew2
                            ynew = ynew2
                            dold = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                            dnew = np.sqrt((x1 - xnew) ** 2 + (y1 - ynew) ** 2)
                            znew = z1 + dnew / dold * (
                                z2 - z1
                            )  # elevation just before jump
                            tnew = t1 + dnew / dold * (t2 - t1)
                            Qbelow = (znew - aq.z[modellayer + 1]) / aq.Haq[layer] * Qn1
                            znew2 = aq.z[modellayer + 1] + Qbelow / Qn2 * aq.Haq[layer]
                            changed = True
                            xyztnew = [
                                np.array([xnew, ynew, znew, tnew]),
                                np.array([xnew, ynew, znew2, tnew]),
                            ]
                        else:  # line-sink takes part of water out
                            if verbose:
                                print("change 4")
                            xnew = xnew2
                            ynew = ynew2
                            dold = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                            dnew = np.sqrt((x1 - xnew) ** 2 + (y1 - ynew) ** 2)
                            znew = z1 + dnew / dold * (
                                z2 - z1
                            )  # elevation just before jump
                            tnew = t1 + dnew / dold * (t2 - t1)
                            Qbelow = (znew - aq.z[modellayer + 1]) / aq.Haq[layer] * Qn1
                            if Qbelow > Qn2:  # taken out
                                terminate = True
                                xyztnew = [np.array([xnew, ynew, znew, tnew])]
                            else:
                                znew2 = (
                                    aq.z[modellayer + 1] + Qbelow / Qn2 * aq.Haq[layer]
                                )
                                xyztnew = [
                                    np.array([xnew, ynew, znew, tnew]),
                                    np.array([xnew, ynew, znew2, tnew]),
                                ]
                            changed = True
        if terminate:
            message = "reached element of type linesink"
            if self.label is not None:
                message += " ({lab})".format(lab=self.label)
        if verbose:
            print("xyztnew, changed", xyztnew, changed)
        return changed, terminate, xyztnew, message


class LineSinkBase(LineSinkChangeTrace, Element):
    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        Qls=100.0,
        res=0,
        wh=1,
        layers=0,
        name="LineSinkBase",
        label=None,
        addtomodel=True,
    ):
        Element.__init__(
            self, model, nparam=1, nunknowns=0, layers=layers, name=name, label=label
        )
        self.nparam = len(self.layers)
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.x2 = float(x2)
        self.y2 = float(y2)
        self.Qls = np.atleast_1d(Qls)
        self.res = float(res)
        self.wh = wh
        self.addtomodel = addtomodel
        if self.addtomodel:
            self.model.add_element(self)

    def __repr__(self):
        return (
            self.name
            + " from "
            + str((self.x1, self.y1))
            + " to "
            + str((self.x2, self.y2))
        )

    def initialize(self):
        self.xc = np.array([0.5 * (self.x1 + self.x2)])
        self.yc = np.array([0.5 * (self.y1 + self.y2)])
        self.ncp = 1
        self.z1 = self.x1 + 1j * self.y1
        self.z2 = self.x2 + 1j * self.y2
        self.L = np.abs(self.z1 - self.z2)
        self.theta_norm_out = (
            np.arctan2(self.y2 - self.y1, self.x2 - self.x1) + np.pi / 2
        )
        self.order = 0  # This is for uniform strength only
        self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        if self.addtomodel:
            self.aq.add_element(self)
        self.parameters = np.empty((self.nparam, 1))
        self.parameters[:, 0] = self.Qls / self.L
        if self.wh == "H":
            self.whfac = self.aq.Haq[self.layers]
        elif self.wh == "2H":
            self.whfac = 2.0 * self.aq.Haq[self.layers]
        elif np.isscalar(self.wh):
            self.whfac = self.wh * np.ones(self.nlayers)
        self.resfac = self.aq.Haq[self.layers] * self.res / self.whfac

    def potinf(self, x, y, aq=None):
        """Can be called with only one x,y value."""
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        if aq == self.aq:
            pot = np.zeros(aq.naq)
            pot[:] = bessel.bessel.potbeslsho(
                float(x), float(y), self.z1, self.z2, aq.lab, 0, aq.ilap, aq.naq
            )
            rv[:] = self.aq.coef[self.layers] * pot
        return rv

    def disvecinf(self, x, y, aq=None):
        """Can be called with only one x,y value."""
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq == self.aq:
            qxqy = np.zeros((2, aq.naq))
            qxqy[:, :] = bessel.bessel.disbeslsho(
                float(x), float(y), self.z1, self.z2, aq.lab, 0, aq.ilap, aq.naq
            )
            rv[0] = self.aq.coef[self.layers] * qxqy[0]
            rv[1] = self.aq.coef[self.layers] * qxqy[1]
        return rv

    def discharge(self):
        # returns the discharge in each layer
        Q = np.zeros(self.aq.naq)
        Q[self.layers] = self.parameters[:, 0] * self.L
        return Q

    def plot(self, layer=None):
        if (layer is None) or (layer in self.layers):
            plt.plot([self.x1, self.x2], [self.y1, self.y2], "k")


class HeadLineSinkZero(LineSinkBase, HeadEquation):
    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        hls=1.0,
        res=0,
        wh=1,
        layers=0,
        label=None,
        addtomodel=True,
    ):
        self.storeinput(inspect.currentframe())
        LineSinkBase.__init__(
            self,
            model,
            x1,
            y1,
            x2,
            y2,
            Qls=0,
            res=res,
            wh=wh,
            layers=layers,
            name="HeadLineSink",
            label=label,
            addtomodel=addtomodel,
        )
        self.hc = np.atleast_1d(float(hls))
        self.nunknowns = self.nparam

    def initialize(self):
        LineSinkBase.initialize(self)
        self.pc = self.hc * self.aq.T[self.layers]  # Needed in solving

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class LineSinkHoBase(LineSinkChangeTrace, Element):
    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        Qls=0.0,
        layers=0,
        order=0,
        name="LineSinkHoBase",
        label=None,
        addtomodel=True,
        aq=None,
        dely=0,
        # zcinout=None,
    ):
        Element.__init__(
            self, model, nparam=1, nunknowns=0, layers=layers, name=name, label=label
        )
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.x2 = float(x2)
        self.y2 = float(y2)
        self.Qls = np.atleast_1d(Qls)
        self.order = order
        self.nparam = self.nlayers * (self.order + 1)
        self.addtomodel = addtomodel
        if addtomodel:
            self.model.add_element(self)
        self.aq = aq
        self.dely = dely
        # self.zcinout = zcinout

    def __repr__(self):
        return (
            self.name
            + " from "
            + str((self.x1, self.y1))
            + " to "
            + str((self.x2, self.y2))
        )

    def initialize(self):
        self.ncp = self.order + 1
        self.z1 = self.x1 + 1j * self.y1
        self.z2 = self.x2 + 1j * self.y2
        self.L = np.abs(self.z1 - self.z2)
        # bchanged minus to plus
        self.theta_norm_out = (
            np.arctan2(self.y2 - self.y1, self.x2 - self.x1) + np.pi / 2.0
        )
        self.cosnorm = np.cos(self.theta_norm_out) * np.ones(self.ncp)
        self.sinnorm = np.sin(self.theta_norm_out) * np.ones(self.ncp)
        # array of ncp by nlayers * (order + 1)
        self.strengthinf = strengthinf_controlpoints(self.ncp, self.nlayers)
        #
        self.xc, self.yc = controlpoints(
            self.ncp, self.z1, self.z2, eps=0, dely=self.dely
        )
        self.xcin, self.ycin = controlpoints(
            self.ncp, self.z1, self.z2, eps=1e-6, dely=self.dely
        )
        self.xcout, self.ycout = controlpoints(
            self.ncp, self.z1, self.z2, eps=-1e-6, dely=self.dely
        )
        # if self.zcinout is not None:
        #     self.xcin, self.ycin = controlpoints(
        #         self.ncp, self.zcinout[0], self.zcinout[1], eps=0
        #     )
        #     self.xcout, self.ycout = controlpoints(
        #         self.ncp, self.zcinout[2], self.zcinout[3], eps=0
        #     )
        # else:
        #     self.xcin, self.ycin = controlpoints(self.ncp, self.z1, self.z2, eps=1e-6)
        #     self.xcout, self.ycout = controlpoints(
        #         self.ncp, self.z1, self.z2, eps=-1e-6
        #     )
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.xc[0], self.yc[0])
        if self.addtomodel:
            self.aq.add_element(self)
        self.parameters = np.empty((self.nparam, 1))
        # Not sure if that needs to be here
        self.parameters[:, 0] = self.Qls / self.L

    def potinf(self, x, y, aq=None):
        """Can be called with only one x,y value.

        Returns
        -------
        array
            (nparam, self.aq.naq) with order
            order 0, layer[0]
            order 0, layer[1]
            ...
            order 1, layer[0]
            order 1, layer[1]
            etc.
        """
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        if aq == self.aq:
            # clever way of using a reshaped rv here
            potrv = rv.reshape((self.order + 1, self.nlayers, aq.naq))
            pot = np.zeros((self.order + 1, aq.naq))
            pot[:, :] = bessel.bessel.potbeslsv(
                float(x),
                float(y),
                self.z1,
                self.z2,
                aq.lab,
                self.order,
                aq.ilap,
                aq.naq,
            )
            potrv[:] = self.aq.coef[self.layers] * pot[:, np.newaxis, :]
        return rv

    def disvecinf(self, x, y, aq=None):
        """Can be called with only one x,y value.

        Returns
        -------
        array
            (nparam, self.aq.naq) with order
            order 0, layer[0]
            order 0, layer[1]
            ...
            order 1, layer[0]
            order 1, layer[1]
            etc.
        """
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq == self.aq:
            qxqyrv = rv.reshape((2, self.order + 1, self.nlayers, aq.naq))
            qxqy = np.zeros((2 * (self.order + 1), aq.naq))
            qxqy[:, :] = bessel.bessel.disbeslsv(
                float(x),
                float(y),
                self.z1,
                self.z2,
                aq.lab,
                self.order,
                aq.ilap,
                aq.naq,
            )
            qxqyrv[0, :] = (
                self.aq.coef[self.layers] * qxqy[: self.order + 1, np.newaxis, :]
            )
            qxqyrv[1, :] = (
                self.aq.coef[self.layers] * qxqy[self.order + 1 :, np.newaxis, :]
            )
        return rv

    def plot(self, layer=None, **kwargs):
        if (layer is None) or (layer in self.layers):
            plt.plot([self.x1, self.x2], [self.y1, self.y2], "k", **kwargs)

    def dischargeinf(self):
        # returns the unit contribution to the discharge in each layer
        # array of length nunknowns
        Qdisinf = np.zeros((self.order + 1, self.nlayers))
        for n in range(self.order + 1):
            Qdisinf[n] = (1 ** (n + 1) - (-1) ** (n + 1)) / (n + 1)
        rv = self.L / 2 * Qdisinf.ravel()
        return rv

    def discharge(self):
        # returns the discharge in each layer
        rv = np.zeros(self.aq.naq)
        Qls = self.parameters[:, 0] * self.dischargeinf()
        rv[self.layers] = np.sum(Qls.reshape(self.order + 1, self.nlayers), 0)
        return rv

    def headinside(self, icp=0):
        hinside = self.model.head(self.xc[icp], self.yc[icp], self.layers) - np.sum(
            self.resfac[icp] * self.parameters[:, 0], 1
        )
        return hinside


class HeadLineSink(LineSinkHoBase, HeadEquation):
    """Head-specified line-sink which may optionally have a width and resistance.

    Parameters
    ----------
    model : Model object
        Model to which the element is added
    x1 : scalar
        x-coordinate of fist point of line-sink
    y1 : scalar
        y-coordinate of fist point of line-sink
    x2 : scalar
        x-coordinate of second point of line-sink
    y2 : scalar
        y-coordinate of second point of line-sink
    hls : scalar, array or list
        head along line-sink
        if scalar: head is the same everywhere along line-sink
        if list or array of length 2: head at beginning and end of line-sink
        if list or array with length order + 1: heads at control points
    res : scalar (default is 0)
        resistance of line-sink
    wh : scalar or str
        distance over which water enters line-sink
        * if 'H': (default) distance is equal to the thickness of the aquifer layer
          (when flow comes mainly from one side)
        * if '2H': distance is twice the thickness of the aquifer layer
          (when flow comes from both sides)
        * if scalar: the width of the stream that partially penetrates the aquifer layer
    order : int (default is 0)
        polynomial order or inflow along line-sink
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    dely : float (default  is 0)
        distance between control points and line
    label : str or None
        label of element

    See Also
    --------
    :class:`.HeadLineSinkString`

    """

    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        hls=1.0,
        res=0,
        wh=1,
        order=0,
        layers=0,
        dely=0,
        label=None,
        name="HeadLineSink",
        addtomodel=True,
    ):
        self.storeinput(inspect.currentframe())
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
            name=name,
            label=label,
            addtomodel=addtomodel,
            dely=dely,
        )
        self.hls = np.atleast_1d(hls)
        self.res = res
        self.wh = wh
        self.nunknowns = self.nparam

    def initialize(self):
        LineSinkHoBase.initialize(self)
        if self.wh == "H":
            self.whfac = self.aq.Haq[self.layers]
        elif self.wh == "2H":
            self.whfac = 2.0 * self.aq.Haq[self.layers]
        elif np.isscalar(self.wh):
            self.whfac = self.wh * np.ones(self.nlayers)
        self.resfac = (
            np.tile(self.res / self.whfac, self.ncp) * self.strengthinf
        )  # this is resfach !
        self.resfac.shape = (self.ncp, self.nlayers, self.nunknowns)
        if len(self.hls) == 1:
            self.hc = self.hls * np.ones(self.nparam)
        elif len(self.hls) == self.ncp:  # head specified at control points
            self.hc = np.repeat(self.hls, self.nlayers)
        elif len(self.hls) == 2:
            s = np.sqrt((self.xc - self.x1) ** 2 + (self.yc - self.y1) ** 2)
            self.hc = np.interp(s, [0, self.L], self.hls)
            self.hc = np.repeat(self.hc, self.nlayers)

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class Ditch(HeadLineSink):
    """Line-sink with specified total discharge, and uniform but unknown head.

    Parameters
    ----------
    model : Model object
        Model to which the element is added
    x1 : scalar
        x-coordinate of fist point of line-sink
    y1 : scalar
        y-coordinate of fist point of line-sink
    x2 : scalar
        x-coordinate of second point of line-sink
    y2 : scalar
        y-coordinate of second point of line-sink
    Qls : scalar
        total discharge of the line-sink
    res : scalar (default is 0)
        resistance of line-sink
    wh : scalar or str
        distance over which water enters line-sink
        * if 'H': (default) distance is equal to the thickness of the aquifer layer
          (when flow comes mainly from one side)
        * if '2H': distance is twice the thickness of the aquifer layer
          (when flow comes from both sides)
        * if scalar: the width of the stream that partially penetrates the aquifer layer
    order : int (default is 0)
        polynomial order or inflow along line-sink
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    dely : float (default  is 0)
        distance between control points and line
    label: str or None
        label of element

    See Also
    --------
    :class:`.DitchString`

    """

    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        Qls=1,
        res=0,
        wh=1,
        order=0,
        layers=0,
        dely=0,
        label=None,
        addtomodel=True,
    ):
        self.storeinput(inspect.currentframe())
        HeadLineSink.__init__(
            self,
            model,
            x1,
            y1,
            x2,
            y2,
            hls=0,
            res=res,
            wh=wh,
            order=order,
            layers=layers,
            dely=dely,
            label=label,
            name="HeadLineSinkDitch",
            addtomodel=addtomodel,
        )
        self.Qls = Qls

    def initialize(self):
        HeadLineSink.initialize(self)

    def equation(self):
        mat, rhs = HeadLineSink.equation(self)
        for i in range(1, self.nunknowns):
            mat[i] -= mat[0]
            rhs[i] -= rhs[0]
        # first equation is sum of discharges equals Qls
        mat[0] = 0
        ieq = 0
        for e in self.model.elementlist:
            if e.nunknowns > 0:
                if e == self:
                    mat[0, ieq : ieq + e.nunknowns] = self.dischargeinf()
                    break
                ieq += e.nunknowns
        rhs[0] = self.Qls
        return mat, rhs

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class LineSinkStringBase2(Element):
    """Alternative implementation that loops through line-sinks to build equation.

    Has the advantage that it is easier to have different line-sinks in different layers
    and/or aquifers.

    """

    def __init__(
        self,
        model,
        xy,
        closed=False,
        layers=0,
        order=0,
        dely=0,
        name="LineSinkStringBase",
        label=None,
        aq=None,
    ):
        Element.__init__(
            self, model, nparam=1, nunknowns=0, layers=layers, name=name, label=label
        )
        self.xy = np.atleast_2d(xy).astype("d")
        if closed:
            self.xy = np.vstack((self.xy, self.xy[0]))
        self.x, self.y = self.xy[:, 0], self.xy[:, 1]
        self.order = order  # same for all segments in string
        self.dely = dely  # same for all segments in string
        self.lslist = []
        if self.xy.shape[1] == 2:
            self.nls = len(self.xy) - 1
        elif self.xy.shape[1] == 4:
            self.nls = len(self.xy)
        if self.layers.ndim == 1:
            if len(self.layers) == self.nls:
                self.layers = self.layers[:, np.newaxis]
            else:  # entire string in these layers
                self.layers = self.layers * np.ones(
                    (self.nls, len(self.layers)), dtype="int"
                )
        self.nlayers = len(self.layers[0])

    def __repr__(self):
        return self.name + " with nodes " + str(self.xy)

    def initialize(self):
        for ls in self.lslist:
            ls.initialize()
        self.aq = []
        for ls in self.lslist:
            if ls.aq not in self.aq:
                self.aq.append(ls.aq)
        for aq in self.aq:
            aq.add_element(self)
        # Same order for all elements in string
        # self.ncp = sum(ls.ncp for ls in self.lslist)
        self.nparam = sum(ls.nparam for ls in self.lslist)
        self.nunknowns = self.nparam
        # where are self.xls and self.yls used? self.xls and self.yls removed
        self.parameters = np.zeros((self.nparam, 1))

    def potinf(self, x, y, aq=None):
        """Compute unit potential influence of element.

        Returns
        -------
        array
            linesink 0, order 0, layer[0]
                        order 0, layer[1]
                        ...
                        order 1, layer[0]
                        order 1, layer[1]
                        ...
            linesink 1, order 0, layer[0]
                        order 0, layer[1]
                        ...
                        order 1, layer[0]
                        order 1, layer[1]
                        ...
        """
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nls, self.lslist[0].nparam, aq.naq))
        if aq in self.aq:
            for i, ls in enumerate(self.lslist):
                rv[i] = ls.potinf(x, y, aq)
        rv.shape = (self.nparam, aq.naq)
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nls, self.lslist[0].nparam, aq.naq))
        if aq in self.aq:
            for i, ls in enumerate(self.lslist):
                rv[:, i] = ls.disvecinf(x, y, aq)
        rv.shape = (2, self.nparam, aq.naq)
        return rv

    def dischargeinf(self):
        rv = np.zeros((self.nls, self.lslist[0].nparam))
        for i, ls in enumerate(self.lslist):
            rv[i] = ls.dischargeinf()
        return rv.ravel()

    def discharge_per_linesink(self):
        """Discharge of the linesinks in each layer.

        Returns
        -------
        rv : np.array
            array of shape (nlay, nlinesinks)
        """
        Qls = self.parameters[:, 0] * self.dischargeinf()
        Qls.shape = (self.nls, self.nlayers, self.order + 1)
        Qls = Qls.sum(axis=2)
        rv = np.zeros((self.model.aq.naq, self.nls))
        for i, q in enumerate(Qls):
            rv[self.layers[i], i] += q
        return rv

    def discharge(self):
        """Discharge of the element in each layer."""
        rv = np.zeros(self.aq[0].naq)
        Qls = self.parameters[:, 0] * self.dischargeinf()
        Qls.shape = (self.nls, self.nlayers, self.order + 1)
        Qls = np.sum(Qls, 2)
        for i, q in enumerate(Qls):
            rv[self.layers[i]] += q
            # rv[self.layers] = np.sum(Qls.reshape(self.nls * (self.order + 1),
            # self.nlayers), 0)
        return rv

    def headinside(self):
        return self.lslist[0].headinside()[0]

    def changetrace(
        self, xyzt1, xyzt2, aq, layer, ltype, modellayer, direction, hstepmax
    ):
        changed = False
        terminate = False
        xyztnew = 0
        message = None
        for ls in self.lslist:
            changed, terminate, xyztnew, message = ls.changetrace(
                xyzt1, xyzt2, aq, layer, ltype, modellayer, direction, hstepmax
            )
            if changed or terminate:
                return changed, terminate, xyztnew, message
        return changed, terminate, xyztnew, message

    def plot(self, layer=None, **kwargs):
        for ls in self.lslist:
            ls.plot(layer=layer, **kwargs)


class HeadLineSinkString(LineSinkStringBase2):
    """String of head-specified line-sinks with optional width and resistance.

    Parameters
    ----------
    model : Model object
        Model to which the element is added
    xy : array or list
        list or array of (x,y) pairs of coordinates of end-points of
        line-sinks in string
    hls : scalar, array or list
        head along string
        if scalar: head is the same everywhere along the string
        if list or array of length 2: head at beginning and end of string
        if list or array with same length as xy: heads at nodes, which
        may contain nans, except for first and last point
    res : scalar (default is 0)
        resistance of line-sink
    wh : scalar or str
        distance over which water enters line-sink
        if 'H': (default) distance is equal to the thickness of the aquifer
        layer (when flow comes mainly from one side)
        if '2H': distance is twice the thickness of the aquifer layer (when
        flow comes from both sides)
        if scalar: the width of the stream that partially penetrates the
        aquifer layer
    order : int (default is 0)
        order of all line-sinks in string
    dely : float (default  is 0)
        distance between control points and line
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    label: str or None

    See Also
    --------
    :class:`.HeadLineSink`
    """

    def __init__(
        self,
        model,
        xy=None,
        hls=0,
        res=0,
        wh=1,
        order=0,
        layers=0,
        dely=0,
        label=None,
        name="HeadLineSinkString",
    ):
        if xy is None:
            xy = [(-1, 0), (1, 0)]
        self.storeinput(inspect.currentframe())
        LineSinkStringBase2.__init__(
            self,
            model,
            xy,
            closed=False,
            layers=layers,
            order=order,
            name=name,
            label=label,
            dely=dely,
            aq=None,
        )
        self.hls = np.atleast_1d(hls)
        self.res = res
        self.wh = wh
        self.model.add_element(self)
        # TO DO: TEST FOR DIFFERENT AQUIFERS AND LAYERS

    def initialize(self):
        if len(self.hls) == 1:  # one value
            self.hls = self.hls * np.ones(self.nls + 1)  # at all nodes
        elif len(self.hls) == 2:  # values at beginning and end
            L = np.sqrt(
                (self.x[1:] - self.x[:-1]) ** 2 + (self.y[1:] - self.y[:-1]) ** 2
            )
            s = np.hstack((0, np.cumsum(L)))
            self.hls = np.interp(s, [0, s[-1]], self.hls)
        elif len(self.hls) == len(self.x):  # nodes may contain nan values
            if np.isnan(self.hls).any():
                L = np.sqrt(
                    (self.x[1:] - self.x[:-1]) ** 2 + (self.y[1:] - self.y[:-1]) ** 2
                )
                s = np.hstack((0, np.cumsum(L)))
                self.hls = np.interp(
                    s, s[~np.isnan(self.hls)], self.hls[~np.isnan(self.hls)]
                )
        else:
            print("Error: hls entry not supported in HeadLineSinkString")
        self.lslist = []  # start with empty list
        for i in range(self.nls):
            if self.label is not None:
                lslabel = self.label + "_" + str(i)
            else:
                lslabel = self.label
            if self.xy.shape[1] == 2:
                x1, x2 = self.x[i : i + 2]
                y1, y2 = self.y[i : i + 2]
            elif self.xy.shape[1] == 4:
                x1, y1, x2, y2 = self.xy[i]
            self.lslist.append(
                HeadLineSink(
                    self.model,
                    x1=x1,
                    y1=y1,
                    x2=x2,
                    y2=y2,
                    hls=self.hls[i : i + 2],
                    res=self.res,
                    wh=self.wh,
                    layers=self.layers[i],
                    order=self.order,
                    dely=self.dely,
                    label=lslabel,
                    addtomodel=False,
                )
            )
        LineSinkStringBase2.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol
        # assign parameters to individual linesinks
        i = 0
        for ls in self.lslist:
            ls.parameters[:, 0] = sol[i : i + ls.nparam]
            i += ls.nparam

    def equation(self):
        mat = np.empty((self.nunknowns, self.model.neq))
        rhs = np.empty(self.nunknowns)
        ieq = 0
        for ls in self.lslist:
            matls, rhsls = ls.equation()
            neq = len(rhsls)
            mat[ieq : ieq + neq] = matls
            rhs[ieq : ieq + neq] = rhsls
            ieq += neq
        # include resistance by computing position of coefficients in matrix
        # and subtracting resistance terms
        iself = self.model.elementlist.index(self)
        jcol = np.sum(e.nunknowns for e in self.model.elementlist[:iself])
        irow = 0
        for ls in self.lslist:
            for icp in range(ls.ncp):
                mat[irow : irow + ls.nlayers, jcol : jcol + ls.nunknowns] -= ls.resfac[
                    icp
                ]
                irow += ls.nlayers
            jcol += ls.nunknowns
        return mat, rhs


class DitchString(HeadLineSinkString):
    """String of Ditches with specified discharge and uniform unknown head.

    Parameters
    ----------
    model : Model object
        Model to which the element is added
    xy : array or list
        list or array of (x,y) pairs of coordinates of end-points of
        line-sinks in string
    Qls : scalar
        total discharge of the string
    res : scalar (default is 0)
        resistance of line-sinks in string
    wh : float or str
        distance over which water enters line-sink.
           * if 'H': (default) distance is equal to the thickness of the aquifer
             (when flow comes mainly from one side)
           * if '2H': distance is twice the thickness of the aquifer layer
             (when flow comes from both sides)
           * if scalar: the width of the stream that partially penetrates the aquifer
    order : int (default is 0)
        polynomial order or inflow along each line-sink in string
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    dely : float (default  is 0)
        distance between control points and line
    label: str or None
        label of element

    See Also
    --------
    :class:`.Ditch`
    """

    def __init__(
        self,
        model,
        xy=None,
        Qls=1,
        res=0,
        wh=1,
        order=0,
        layers=0,
        dely=0,
        label=None,
    ):
        if xy is None:
            xy = [(-1, 0), (1, 0)]
        self.storeinput(inspect.currentframe())
        HeadLineSinkString.__init__(
            self,
            model,
            xy=xy,
            hls=0,
            res=res,
            wh=wh,
            order=order,
            dely=dely,
            layers=layers,
            label=label,
            name="DitchString",
        )
        self.Qls = Qls

    def initialize(self):
        HeadLineSinkString.initialize(self)

    def equation(self):
        mat, rhs = HeadLineSinkString.equation(self)
        for i in range(1, self.nunknowns):
            mat[i] -= mat[0]
            rhs[i] -= rhs[0]
        # first equation is sum of discharges equals Qls
        mat[0] = 0
        ieq = 0
        for e in self.model.elementlist:
            if e.nunknowns > 0:
                if e == self:
                    mat[0, ieq : ieq + self.nunknowns] = self.dischargeinf()
                    break
                ieq += e.nunknowns
        rhs[0] = self.Qls
        return mat, rhs


class LineSinkDitch(Ditch):
    """Deprecated alias for :class:`.Ditch`.
    
    .. deprecated:: 
        Use :class:`.Ditch` instead. This alias will be removed in a future version.
    """

    def __init__(self, *args, **kwargs):
        warnings.warn(
            "LineSinkDitch is deprecated. Use Ditch instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)


class LineSinkDitchString(DitchString):
    """Deprecated alias for :class:`.DitchString`.
    
    .. deprecated:: 
        Use :class:`.DitchString` instead. This alias will be removed in a future version.
    """

    def __init__(self, *args, **kwargs):
        warnings.warn(
            "LineSinkDitchString is deprecated. Use DitchString instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)


class CollectorWell(DitchString):
    """Collector well: collection of line sinks with a specified total discharge.

    Collection of (discontinuous) line sinks with specified total discharge
    and unknown but uniform head.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    xy : np.array
        array of shape (N, 4) with start and end coordinates of the line sinks
        on each row: [(x1, y1, x2, y2), ...]
    Qw : float
        total discharge of the collector well
    rw : float
        radius of the collector well arms
    res : float
        resistance of the well screen
    layers : int, array or list
        layer (int) or layers (list or array) where well is screened
    order : int
        order of the line sink elements
    label : string, optional
        label of the collector well

    Examples
    --------
    Create a collector well with two arms::

        ml = timml.Model3D(kaq=10, z=np.arange(20, -1, -2), kzoverkh=0.1)
        xy = [(1, 0, 10, 0), (0, 1, 0, 10)]
        w = timml.CollectorWell(ml, xy=xy, Qw=1000, layers=np.arange(5, 10))
        ml.solve()

    """

    def __init__(
        self,
        model,
        xy,
        Qw=100.0,
        rw=0.1,
        res=0.0,
        layers=0,
        order=0,
        label=None,
    ):
        super().__init__(
            model,
            xy,
            Qls=Qw,
            res=res,
            layers=layers,
            order=order,
            dely=rw,
            label=label,
        )
        self.name = "CollectorWell"


class RadialCollectorWell(CollectorWell):
    """Radial collector well.

    Collection of (discontinuous) line sinks in a radial pattern with specified
    total discharge and unknown but uniform head.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    x : float
        x-coordinate of the center of the collector well
    y : float
        y-coordinate of the center of the collector well
    L : float
        length of each arm
    narms : int
        number of arms
    rcaisson : float
        radius of the caisson
    rw : float
        radius of the arms
    nls : int
        number of line sinks per arm
    Qw : float
        total discharge of the collector well
    res : float
        resistance of the arms
    layers : int, array or list
        layer(s) in which the well is screened
    label : string, optional
        label of the collector well

    Examples
    --------
    Create a radial collector well with 5 arms::

        ml = timml.Model3D(kaq=10, z=np.arange(20, -1, -2), kzoverkh=0.1)
        w = timml.RadialCollectorWell(ml, x=0, y=0, narms=5, nls=10, angle=0,
                                      rcaisson=2.0, rw=0.1, Qw=1000, layers=5)
        ml.solve()

    """

    def __init__(
        self,
        model,
        x=0,
        y=0,
        narms=5,
        nls=10,
        L=10.0,
        angle=0,
        rcaisson=1.0,
        rw=0.1,
        Qw=100.0,
        res=0.0,
        layers=0,
        label=None,
    ):
        if np.isscalar(angle):
            angle = np.deg2rad(angle) + np.linspace(0, 2 * np.pi, narms + 1)[:-1]
        else:
            angle = np.deg2rad(angle)
        if np.isscalar(L):
            L = L * np.ones(narms)
        if np.isscalar(nls):
            nls = nls * np.ones(narms, dtype="int")
        if np.isscalar(layers):
            layers = layers * np.ones(narms, dtype="int")
        xy, layers = self.compute_xy(x, y, rcaisson, narms, nls, L, angle, layers)
        super().__init__(
            model,
            xy,
            Qw=Qw,
            rw=rw,
            res=res,
            layers=layers,
            label=label,
        )
        self.name = "RadialCollectorWell"
        self.iarm = np.cumsum(np.hstack((0, nls)))  # starting index of arms

    def compute_xy(self, x, y, rcaisson, narms, nls, L, angle, layer_arms):
        """Compute the x,y-coordinates array for the radial collector well.

        Parameters
        ----------
        x : float
            x-coordinate of the center of the collector well
        y : float
            y-coordinate of the center of the collector well
        narms : int
            number of arms
        nls : int or array
            number of line sinks per arm
        L : float or array
            length of each arm
        angle : float or array
            angle of first arm or of each arm
        rcaisson : float
            radius of the caisson

        Returns
        -------
        xy : np.array
            array of shape (N, 4) with start and end coordinates of the line sinks
            on each row: [(x1, y1, x2, y2), ...]
        """
        xy = np.empty((np.sum(nls), 4))
        layers = np.empty(np.sum(nls), dtype="int")
        for i in range(narms):
            x = rcaisson * np.cos(angle[i]) + np.linspace(0, L[i], nls[i] + 1) * np.cos(
                angle[i]
            )
            y = rcaisson * np.sin(angle[i]) + np.linspace(0, L[i], nls[i] + 1) * np.sin(
                angle[i]
            )
            i0 = np.sum(nls[:i])
            xy[i0 : i0 + nls[i], 0] = x[:-1]
            xy[i0 : i0 + nls[i], 1] = y[:-1]
            xy[i0 : i0 + nls[i], 2] = x[1:]
            xy[i0 : i0 + nls[i], 3] = y[1:]
            layers[i0 : i0 + nls[i]] = layer_arms[i]
        return xy, layers

    def discharge_per_arm(self):
        """Discharge of each arm.

        Returns
        -------
        array (length number of arms)
            Total discharge of each arm
        """
        dis = np.zeros(len(self.iarm) - 1)
        for i in range(len(dis)):
            for j in range(self.iarm[i], self.iarm[i + 1]):
                dis[i] += self.lslist[j].discharge().sum()
        return dis

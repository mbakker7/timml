import inspect  # Used for storing the input
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np

from timml import bessel
from timml.controlpoints import controlpoints, strengthinf_controlpoints
from timml.element import Element
from timml.equation import HeadEquation
from timml.util import refine_n_segments

__all__ = [
    "LineSinkBase",
    "HeadLineSinkZero",
    "HeadLineSink",
    "LineSinkDitch",
    "HeadLineSinkString",
    "LineSinkDitchString",
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
                #            if (layer == self.layers).any():  # in layer where line-sink is screened
                # not needed anymore, I thin this is all taken care of with checking Qn1 and Qn2
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
                        if (
                            Qn1 < 0
                        ):  # trying to cross line-sink that infiltrates, stay on bottom, don't terminate
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
        refine_level=1,
    ):
        _input = {k: v for k, v in locals().items() if k not in ["self", "model"]}
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
        self.refine_level = refine_level
        if self.addtomodel:
            self.model.add_element(self)
        self._input = _input

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

    def _refine(self, n=None):
        if n is None:
            n = self.refine_level

        # refine xy
        xy = np.array([(self.x1, self.y1), (self.x2, self.y2)])
        xyr, _ = refine_n_segments(xy, "line", n_segments=n)

        # get input arguments
        input_args = deepcopy(self._input)
        cls = input_args.pop("__class__", self.__class__)
        input_args["model"] = self.model
        input_args["refine_level"] = 1  # set to 1 to prevent further refinement
        input_args["addtomodel"] = False  # these are internally created elements

        # build new elements
        refined_elements = []
        Qls = input_args.pop("Qls")
        L = np.sqrt((xyr[1:, 0] - xyr[:-1, 0]) ** 2 + (xyr[1:, 1] - xyr[:-1, 1]) ** 2)
        for ils in range(n):
            (input_args["x1"], input_args["y1"]) = xyr[ils]
            (input_args["x2"], input_args["y2"]) = xyr[ils + 1]
            # distribute discharge evenly according to new segment lengths
            input_args["Qls"] = Qls * L[ils] / np.sum(L)
            refined_elements.append(cls(**input_args))
        return refined_elements


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
        zcinout=None,
        refine_level=1,
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
        self.zcinout = zcinout
        self.refine_level = refine_level

    def __repr__(self):
        return (
            self.name
            + " from "
            + str((self.x1, self.y1))
            + " to "
            + str((self.x2, self.y2))
        )

    def initialize(self, addtoaq=True):
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
        self.xc, self.yc = controlpoints(self.ncp, self.z1, self.z2, eps=0)
        if self.zcinout is not None:
            self.xcin, self.ycin = controlpoints(
                self.ncp, self.zcinout[0], self.zcinout[1], eps=0
            )
            self.xcout, self.ycout = controlpoints(
                self.ncp, self.zcinout[2], self.zcinout[3], eps=0
            )
        else:
            self.xcin, self.ycin = controlpoints(self.ncp, self.z1, self.z2, eps=1e-6)
            self.xcout, self.ycout = controlpoints(
                self.ncp, self.z1, self.z2, eps=-1e-6
            )
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.xc[0], self.yc[0])
        # also respect addtomodel here to prevent sub-elements (e.g. parts of
        # HeadLineSinkString) from being added to the aquifer elementlists
        if (addtoaq is None and self.addtomodel) or addtoaq:
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

    def plot(self, layer=None):
        if (layer is None) or (layer in self.layers):
            plt.plot([self.x1, self.x2], [self.y1, self.y2], "k")

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
        hinside = (
            self.model.head(self.xc[icp], self.yc[icp])[self.layers[0]]
            - np.sum(self.strengthinf[icp] * self.parameters[:, 0])
            * self.res
            / self.whfac
        )
        return hinside


class HeadLineSink(LineSinkHoBase, HeadEquation):
    """Create a head-specified line-sink which may optionally have a width and
    resistance.

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
        if 'H': (default) distance is equal to the thickness of the aquifer layer (when flow comes mainly from one side)
        if '2H': distance is twice the thickness of the aquifer layer (when flow comes from both sides)
        if scalar: the width of the stream that partially penetrates the aquifer layer
    order : int (default is 0)
        polynomial order or inflow along line-sink
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    label: str or None
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
        label=None,
        name="HeadLineSink",
        addtomodel=True,
        refine_level=1,
    ):
        _input = {k: v for k, v in locals().items() if k not in ["self", "model"]}
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
            refine_level=refine_level,
        )
        self.hls = np.atleast_1d(hls)
        self.res = res
        self.wh = wh
        self.nunknowns = self.nparam
        self._input = _input

    def initialize(self, addtoaq=True):
        LineSinkHoBase.initialize(self, addtoaq=addtoaq)
        if self.wh == "H":
            self.whfac = self.aq.Haq[self.layers]
        elif self.wh == "2H":
            self.whfac = 2.0 * self.aq.Haq[self.layers]
        elif np.isscalar(self.wh):
            self.whfac = self.wh * np.ones(self.nlayers)
        self.resfac = np.tile(self.res / self.whfac, self.ncp) * self.strengthinf
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

    def _refine(self, n=None):
        if n is None:
            n = self.refine_level
        # refine xy
        xy = np.array([(self.x1, self.y1), (self.x2, self.y2)])
        xyr, _ = refine_n_segments(xy, "line", n_segments=n)
        # get input arguments
        input_args = deepcopy(self._input)
        cls = input_args.pop("__class__", self.__class__)
        input_args["model"] = self.model
        input_args["refine_level"] = 1  # set to 1 to prevent further refinement
        input_args["addtomodel"] = False
        # build new elements
        refined_elements = []
        hls = np.atleast_1d(input_args.pop("hls"))
        s = np.sqrt((xyr[:, 0] - xyr[0, 0]) ** 2 + (xyr[:, 1] - xyr[0, 1]) ** 2)
        for ils in range(n):
            (input_args["x1"], input_args["y1"]) = xyr[ils]
            (input_args["x2"], input_args["y2"]) = xyr[ils + 1]
            # refine head-specification if possible
            if len(hls) == 1:
                input_args["hls"] = hls
            elif len(hls) == 2:
                input_args["hls"] = np.interp(s, [0, s[-1]], self.hls)[ils : ils + 2]
            elif len(hls) == self.order + 1:
                # this is never well-defined, so raise error
                raise NotImplementedError(
                    "Cannot refine HeadLineSink when hls is defined at control points."
                )
            refined_elements.append(cls(**input_args))
        return refined_elements


class LineSinkDitch(HeadLineSink):
    """Class to create a line-sink for which the total discharge is specified, and for
    which the head along the line-sink is uniform but unknown.

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
        if 'H': (default) distance is equal to the thickness of the aquifer layer (when flow comes mainly from one side)
        if '2H': distance is twice the thickness of the aquifer layer (when flow comes from both sides)
        if scalar: the width of the stream that partially penetrates the aquifer layer
    order : int (default is 0)
        polynomial order or inflow along line-sink
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    label: str or None
        label of element

    See Also
    --------

    :class:`.LineSinkDitchString`
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


class LineSinkStringBase(Element):
    """Original implementation.

    Used for boundaries of inhomogenieities.
    """

    def __init__(
        self,
        model,
        xy,
        closed=False,
        layers=0,
        order=0,
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
        self.order = order
        self.aq = aq
        self.lslist = []
        self.x, self.y = self.xy[:, 0], self.xy[:, 1]
        self.nls = len(self.x) - 1
        for i in range(self.nls):
            if label is not None:
                lslabel = label + "_" + str(i)
            else:
                lslabel = label
            self.lslist.append(
                LineSinkHoBase(
                    model,
                    x1=self.x[i],
                    y1=self.y[i],
                    x2=self.x[i + 1],
                    y2=self.y[i + 1],
                    Qls=0.0,
                    layers=layers,
                    order=order,
                    label=lslabel,
                    addtomodel=False,
                    aq=aq,
                )
            )

    def __repr__(self):
        return self.name + " with nodes " + str(self.xy)

    def initialize(self):
        for ls in self.lslist:
            ls.initialize()
        # Same order for all elements in string
        self.ncp = self.nls * self.lslist[0].ncp
        self.nparam = self.nls * self.lslist[0].nparam
        self.nunknowns = self.nparam
        self.xls = np.empty((self.nls, 2))
        self.yls = np.empty((self.nls, 2))
        for i, ls in enumerate(self.lslist):
            self.xls[i, :] = [ls.x1, ls.x2]
            self.yls[i, :] = [ls.y1, ls.y2]
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(
                self.lslist[0].xc, self.lslist[0].yc
            )
        self.parameters = np.zeros((self.nparam, 1))
        # As parameters are only stored for the element not the list,
        # we need to combine the following
        self.xc = np.array([ls.xc for ls in self.lslist]).flatten()
        self.yc = np.array([ls.yc for ls in self.lslist]).flatten()
        self.xcin = np.array([ls.xcin for ls in self.lslist]).flatten()
        self.ycin = np.array([ls.ycin for ls in self.lslist]).flatten()
        self.xcout = np.array([ls.xcout for ls in self.lslist]).flatten()
        self.ycout = np.array([ls.ycout for ls in self.lslist]).flatten()
        self.cosnorm = np.array([ls.cosnorm for ls in self.lslist]).flatten()
        self.sinnorm = np.array([ls.sinnorm for ls in self.lslist]).flatten()
        self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], self.ycin[0])
        self.aqout = self.model.aq.find_aquifer_data(self.xcout[0], self.ycout[0])

    def potinf(self, x, y, aq=None):
        """Compute the unit potential influence of the element.

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
        for i in range(self.nls):
            rv[i] = self.lslist[i].potinf(x, y, aq)
        rv.shape = (self.nparam, aq.naq)
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nls, self.lslist[0].nparam, aq.naq))
        for i in range(self.nls):
            rv[:, i] = self.lslist[i].disvecinf(x, y, aq)
        rv.shape = (2, self.nparam, aq.naq)
        return rv

    def changetrace(
        self, xyzt1, xyzt2, aq, layer, ltype, modellayer, direction, hstepmax
    ):
        changed = False
        terminate = False
        xyztnew = 0
        message = None
        for ls in self.lslist:
            changed, terminate, xyztnew, message = ls.changetrace(
                xyzt1, xyzt2, aq, layer, ltype, modellayer, direction
            )
            if changed or terminate:
                return changed, terminate, xyztnew, message
        return changed, terminate, xyztnew, message

    def plot(self, layer=None):
        if (layer is None) or (layer in self.layers):
            plt.plot(self.x, self.y, "k")


class HeadLineSinkStringOLd(LineSinkStringBase, HeadEquation):
    def __init__(
        self, model, xy=[(-1, 0), (1, 0)], hls=0.0, layers=0, order=0, label=None
    ):
        self.storeinput(inspect.currentframe())
        LineSinkStringBase.__init__(
            self,
            model,
            xy,
            closed=False,
            layers=layers,
            order=order,
            name="HeadLineSinkString",
            label=label,
            aq=None,
        )
        self.hls = np.atleast_1d(hls)
        self.model.add_element(self)

    def initialize(self):
        LineSinkStringBase.initialize(self)
        self.aq.add_element(self)
        # self.pc = np.array([ls.pc for ls in self.lslist]).flatten()
        if len(self.hls) == 1:
            self.pc = self.hls * self.aq.T[self.layers] * np.ones(self.nparam)
        elif len(self.hls) == self.nls:  # head specified at centers
            self.pc = (self.hls[:, np.newaxis] * self.aq.T[self.layers]).flatten()
        elif len(self.hls) == 2:
            L = np.array([ls.L for ls in self.lslist])
            Ltot = np.sum(L)
            xp = np.zeros(self.nls)
            xp[0] = 0.5 * L[0]
            for i in range(1, self.nls):
                xp[i] = xp[i - 1] + 0.5 * (L[i - 1] + L[i])
            self.hls = np.interp(xp, [0, Ltot], self.hls)
            self.pc = (self.hls[:, np.newaxis] * self.aq.T[self.layers]).flatten()
        else:
            print("Error: hls entry not supported")
        self.resfac = 0.0

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
        name="LineSinkStringBase",
        label=None,
        aq=None,
        refine_level=1,
    ):
        Element.__init__(
            self, model, nparam=1, nunknowns=0, layers=layers, name=name, label=label
        )
        self.xy = np.atleast_2d(xy).astype("d")
        if closed:
            self.xy = np.vstack((self.xy, self.xy[0]))
        self.order = order  # same for all segments in string
        self.lslist = []
        self.x, self.y = self.xy[:, 0], self.xy[:, 1]
        self.nls = len(self.x) - 1
        if self.layers.ndim == 1:
            if len(self.layers) == self.nls:
                self.layers = self.layers[:, np.newaxis]
            else:  # entire string in these layers
                self.layers = self.layers * np.ones(
                    (self.nls, len(self.layers)), dtype=int
                )
        self.nlayers = len(self.layers[0])
        self.refine_level = refine_level

        # set _x, _y and _layers so that _refine can potentially overwrite
        # these attributes without modifying original input
        self._xy = self.xy.copy()
        self._x = self.x.copy()
        self._y = self.y.copy()
        self._layers = self.layers.copy()

    def __repr__(self):
        return self.name + " with nodes " + str(self._xy)

    def initialize(self):
        for ls in self.lslist:
            ls.initialize(addtoaq=False)
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
            rv[self._layers[i], i] += q
        return rv

    def discharge(self):
        """Discharge of the element in each layer."""

        rv = np.zeros(self.aq[0].naq)
        Qls = self.parameters[:, 0] * self.dischargeinf()
        Qls.shape = (self.nls, self.nlayers, self.order + 1)
        Qls = np.sum(Qls, 2)
        for i, q in enumerate(Qls):
            rv[self._layers[i]] += q
            # rv[self.layers] = np.sum(Qls.reshape(self.nls * (self.order + 1), self.nlayers), 0)
        return rv

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

    def plot(self, layer=None):
        if (layer is None) or (layer in self._layers):
            plt.plot(self.x, self.y, "k")

    def _reset(self):
        self._xy = self.xy.copy()
        self._x = self.x.copy()
        self._y = self.y.copy()
        self._layers = self.layers.copy()

class HeadLineSinkString(LineSinkStringBase2):
    """Class to create a string of head-specified line-sinks which may optionally have a
    width and resistance.

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
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    label: str or None
    refine_level : int, optional
        partition each linesink reach into refine_level segments, default is 1, which
        means no refinement is applied.

    See Also
    --------

    :class:`.HeadLineSink`
    """

    def __init__(
        self,
        model,
        xy=[(-1, 0), (1, 0)],
        hls=0,
        res=0,
        wh=1,
        order=0,
        layers=0,
        label=None,
        name="HeadLineSinkString",
        refine_level=1,
    ):
        self.storeinput(inspect.currentframe())
        LineSinkStringBase.__init__(
            self,
            model,
            xy,
            closed=False,
            layers=layers,
            order=order,
            name=name,
            label=label,
            aq=None,
            refine_level=refine_level,
        )
        # TODO: TEST FOR DIFFERENT AQUIFERS AND LAYERS
        self.res = res
        self.wh = wh
        self.model.add_element(self)
        self.hls = np.atleast_1d(hls)

        # set hls for computation, copy so that _refine
        # can potentially overwrite these attributes without modifying original input
        self._hls = self.hls.copy()

    def initialize(self):
        if len(self._hls) == 1:  # one value
            self._hls = self._hls * np.ones(self.nls + 1)  # at all nodes
        elif len(self._hls) == 2:  # values at beginning and end
            L = np.sqrt(
                (self._x[1:] - self._x[:-1]) ** 2 + (self._y[1:] - self._y[:-1]) ** 2
            )
            s = np.hstack((0, np.cumsum(L)))
            self._hls = np.interp(s, [0, s[-1]], self._hls)
        elif len(self._hls) == len(self.x):  # nodes may contain nan values
            if np.isnan(self._hls).any():
                L = np.sqrt(
                    (self._x[1:] - self._x[:-1]) ** 2
                    + (self._y[1:] - self._y[:-1]) ** 2
                )
                s = np.hstack((0, np.cumsum(L)))
                self._hls = np.interp(
                    s, s[~np.isnan(self._hls)], self._hls[~np.isnan(self._hls)]
                )
        else:
            print("Error: hls entry not supported in HeadLineSinkString")

        self.lslist = []  # start with empty list
        for i in range(self.nls):
            if self.label is not None:
                lslabel = self.label + "_" + str(i)
            else:
                lslabel = self.label
            self.lslist.append(
                HeadLineSink(
                    self.model,
                    x1=self._x[i],
                    y1=self._y[i],
                    x2=self._x[i + 1],
                    y2=self._y[i + 1],
                    hls=self._hls[i : i + 2],
                    res=self.res,
                    wh=self.wh,
                    layers=self._layers[i],
                    order=self.order,
                    label=lslabel,
                    addtomodel=False,
                )
            )
        LineSinkStringBase.initialize(self)

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
        # fix to include resistance
        # this is not pretty but works
        # not sure how to change the design to make this nicer
        # I guess the additional matrix can be pre-computed and stored
        jcol = 0
        for e in self.model.elementlist:
            if e == self:
                break
            elif e.nunknowns > 0:
                jcol += e.nunknowns
        irow = 0
        for ls in self.lslist:
            for icp in range(ls.ncp):
                mat[irow : irow + ls.nlayers, jcol : jcol + ls.nunknowns] -= ls.resfac[
                    icp
                ]
                irow += ls.nlayers
            jcol += ls.nunknowns
        return mat, rhs

    def _refine(self, n=None):
        if n is None:
            n = self.refine_level
        xyr, reindexer = refine_n_segments(self.xy, "line", n_segments=n)

        # update attributes
        self._xy = xyr
        self._x = xyr[:, 0]
        self._y = xyr[:, 1]

        # refining logic for hls
        if len(self.hls) == 1:
            self._hls = self.hls.copy()  # initialize will set hls for each segment
        elif len(self.hls == 2):
            self._hls = self.hls.copy()  # initialize will set hls for each segment
        elif len(self.hls) == len(self.x):
            self._hls = self.hls[reindexer]  # reindex user-specified hls
            mask = ~(np.diff(reindexer, prepend=[-1]).astype(bool))
            self._hls[mask] = np.nan

        self._layers = self.layers[reindexer]
        self.nlayers = len(self._layers[0])
        self.nls = len(self._x) - 1
        return [self]

    def _reset(self):
        super()._reset()
        self._hls = self.hls.copy()


class LineSinkDitchString(HeadLineSinkString):
    """Class to create a string of LineSinkDitch elements for which the total discharge
    of the string is specified, and for which the head along the entire string is
    uniform but unknown.

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
    wh : scalar or str
        distance over which water enters line-sink
        if 'H': (default) distance is equal to the thickness of the aquifer layer (when flow comes mainly from one side)
        if '2H': distance is twice the thickness of the aquifer layer (when flow comes from both sides)
        if scalar: the width of the stream that partially penetrates the aquifer layer
    order : int (default is 0)
        polynomial order or inflow along each line-sink in string
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    label: str or None
        label of element

    See Also
    --------

    :class:`.LineSinkDitch`
    """

    def __init__(
        self,
        model,
        xy=[(-1, 0), (1, 0)],
        Qls=1,
        res=0,
        wh=1,
        order=0,
        layers=0,
        label=None,
        refine_level=1,
    ):
        self.storeinput(inspect.currentframe())
        HeadLineSinkString.__init__(
            self,
            model,
            xy=xy,
            hls=0,
            res=res,
            wh=wh,
            order=order,
            layers=layers,
            label=label,
            name="LineSinkDitchString",
            refine_level=refine_level,
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

    def setparams(self, sol):
        self.parameters[:, 0] = sol

    def _refine(self, n=None):
        if n is None:
            n = self.refine_level
        xyr, _ = refine_n_segments(self.xy, "line", n_segments=n)
        # update attributes
        self._x = xyr[:, 0]
        self._y = xyr[:, 1]
        return [self]


class LineSinkContainer(Element):
    """
    Container class for bunch of line-sinks
    Required attributes:
    lslist: list of line-sinks
    nls: total number of line-sinks
    """

    def __init__(
        self,
        model,
        layers=0,
        order=0,
        name="LineSinkStringContainer",
        label=None,
        aq=None,
    ):
        Element.__init__(
            self,
            model,
            nparam=1,
            nunknowns=0,
            layers=layers,
            name=name,
            label=label,
        )
        self.order = order

    #         self.xy = np.atleast_2d(xy).astype('d')
    #         if closed: self.xy = np.vstack((self.xy, self.xy[0]))
    #         self.order = order  # same for all segments in string
    #         self.lslist = []
    #         self.x, self.y = self.xy[:, 0], self.xy[:, 1]
    #         self.nls = len(self.x) - 1
    #         if self.layers.ndim == 1:
    #             if len(self.layers) == self.nls:
    #                 self.layers = self.layers[:, np.newaxis]
    #             else:  # entire string in these layers
    #                 self.layers = self.layers * \
    #                     np.ones((self.nls, len(self.layers)), dtype='int')
    #         self.nlayers = len(self.layers[0])

    def __repr__(self):
        return self.name

    def initialize(self):
        self.layers = []
        for ls in self.lslist:
            ls.initialize()
            if ls.layers[0] not in self.layers:  # works if ls in 1 layer
                self.layers.append(ls.layers[0])
        self.aq = []
        for ls in self.lslist:
            if ls.aq not in self.aq:
                self.aq.append(ls.aq)
        for aq in self.aq:
            aq.add_element(self)
        self.nparam = sum(ls.nparam for ls in self.lslist)
        self.nunknowns = self.nparam
        self.parameters = np.zeros((self.nparam, 1))

    def potinf(self, x, y, aq=None):
        """Compute the unit potential influence of the element.

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

    def discharge(self):
        """Discharge of the element in each layer."""
        # TBD: it is unclear what this function does. shouldn't it return
        # values only for the layers where the element is screened?
        rv = np.zeros(self.aq[0].naq)
        Qls = self.parameters[:, 0] * self.dischargeinf()
        Qls.shape = (self.nls, self.nlayers, self.order + 1)
        Qls = np.sum(Qls, 2)
        for i, q in enumerate(Qls):
            rv[self.layers[i]] += q
            # rv[self.layers] = np.sum(Qls.reshape(self.nls * (self.order + 1), self.nlayers), 0)
        return rv

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

    def plot(self, layer=None):
        if (layer is None) or (layer in self.layers):
            for i in range(len(self.xls)):
                plt.plot(self.xls[i], self.yls[i], "k")


class HeadLineSinkContainer(LineSinkContainer):
    """Class to create a container of head-specified line-sinks which may optionally
    have a width and resistance.

    Parameters
    ----------

    model : Model object
        Model to which the element is added
    xydict : dictionary
        dictionary with lists or arrays of (x,y) pairs of coordinates of
        end-points of line-sinks that are stringed together
    hls : scalar
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
    laydict : scalar, list or array
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
        xydict={0: [(-1, 0), (1, 0)]},
        hls=0,
        res=0,
        wh=1,
        order=0,
        laydict={0: 0},
        label=None,
        name="HeadLineSinkContainer",
    ):
        self.storeinput(inspect.currentframe())
        LineSinkContainer.__init__(
            self, model, layers=0, order=order, name=name, label=label, aq=None
        )
        self.xydict = xydict
        self.hls = hls
        self.laydict = laydict
        self.res = res
        self.wh = wh
        self.model.add_element(self)
        # TODO: TEST FOR DIFFERENT AQUIFERS AND LAYERS

    def initialize(self):
        self.lslist = []
        self.xls = []  # used for layout
        self.yls = []
        for key in self.xydict.keys():
            xy = np.atleast_1d(self.xydict[key])
            layers = np.atleast_1d(self.laydict[key])
            if len(layers) == 1:
                layers = layers * np.ones(len(xy) - 1, dtype="int")
            self.xls.append(xy[:, 0])
            self.yls.append(xy[:, 1])
            for i in range(len(xy) - 1):
                if self.label is not None:
                    lslabel = self.label + "_" + str(i)
                else:
                    lslabel = self.label
                x1, y1 = xy[i]
                x2, y2 = xy[i + 1]
                ls = HeadLineSink(
                    self.model,
                    x1=x1,
                    y1=y1,
                    x2=x2,
                    y2=y2,
                    hls=self.hls,
                    res=self.res,
                    wh=self.wh,
                    layers=layers[i],
                    order=self.order,
                    label=lslabel,
                    addtomodel=False,
                )
                self.lslist.append(ls)
        self.nls = len(self.lslist)
        LineSinkContainer.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol

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
        # fix to include resistance
        # this is not pretty but works
        # not sure how to change the design to make this nicer
        # I guess the additional matrix can be pre-computed and stored
        jcol = 0
        for e in self.model.elementlist:
            if e == self:
                break
            elif e.nunknowns > 0:
                jcol += e.nunknowns
        irow = 0
        for ls in self.lslist:
            for icp in range(ls.ncp):
                mat[irow : irow + ls.nlayers, jcol : jcol + ls.nunknowns] -= ls.resfac[
                    icp
                ]
                irow += ls.nlayers
            jcol += ls.nunknowns
        return mat, rhs

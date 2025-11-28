"""1D line-doublet elements.

Doublet features for modeling barriers in a cross-section.

Example::
    
    ImpLineDoublet1D(ml, xld=0, layers=0)

"""
import inspect  # Used for storing the input

import matplotlib.pyplot as plt
import numpy as np

from .element import Element
from .equation import DisvecEquation, LeakyWallEquation

class LineDoublet1D(Element):
    def __init__(
        self,
        model,
        xld,
        delp=1,
        layers=0,
        name="LineDoublet1DBase",
        label=None,
        addtomodel=True,
        res=0,
        aq=None,
    ):
        Element.__init__(
            self, model, nparam=1, nunknowns=0, layers=layers, name=name, label=label
        )
        self.xld = float(xld)
        self.delp = np.atleast_1d(delp)
        self.res = float(res)
        self.aq = aq
        self.addtomodel = addtomodel
        if self.addtomodel:
            self.model.add_element(self)
        self.nparam = self.nlayers

    def __repr__(self):
        return self.name + " at " + str(self.xld) + " in layers: " + str(self.layers)

    def initialize(self):
        self.xc = np.array([self.xld])
        self.yc = np.zeros(1)
        self.ncp = 1
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.xc[0], self.yc[0])
        if self.addtomodel:
            self.aq.add_element(self)
        self.parameters = np.empty((self.nparam, 1))
        self.parameters[:, 0] = self.delp
        self.theta_norm_out = np.zeros(1)
        self.cosnorm = np.cos(self.theta_norm_out) * np.ones(self.ncp)
        self.sinnorm = np.sin(self.theta_norm_out) * np.ones(self.ncp)
        self.resfac = self.aq.Haq[self.layers] / self.res

    def potinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, 0)
        rv = np.zeros((self.nparam, aq.naq))
        if aq == self.aq:
            pot = np.zeros(aq.naq)
            if aq.ilap:
                if x - self.xld < 0.0:
                    # pot[0] = -0.5 * (x - self.xld - 1)  # so that pot = 0.5 at x=xld
                    pot[0] = -0.5  # so that pot = 0.5 at x=xld
                    pot[1:] = -0.5 * np.exp((x - self.xld) / aq.lab[1:])
                elif x - self.xld >= 0.0:
                    # pot[0] = 0.5 * (x - self.xld + 1)
                    pot[0] = 0.5
                    pot[1:] = 0.5 * np.exp(-(x - self.xld) / aq.lab[1:])
            else:
                if x - self.xld < 0.0:
                    pot[:] = -0.5 * np.exp((x - self.xld) / aq.lab)
                elif x - self.xld >= 0.0:
                    pot[:] = 0.5 * np.exp(-(x - self.xld) / aq.lab)
            rv[:] = self.aq.coef[self.layers] * pot
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, 0)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq == self.aq:
            qx = np.zeros(aq.naq)
            if aq.ilap:
                if x - self.xld < 0.0:
                    qx[0] = 0.0
                    qx[1:] = 0.5 / aq.lab[1:] * np.exp((x - self.xld) / aq.lab[1:])
                elif x - self.xld >= 0.0:
                    qx[0] = 0.0
                    qx[1:] = 0.5 / aq.lab[1:] * np.exp(-(x - self.xld) / aq.lab[1:])
            else:
                if x - self.xld < 0.0:
                    qx[:] = 0.5 / aq.lab * np.exp((x - self.xld) / aq.lab)
                elif x - self.xld >= 0.0:
                    qx[:] = 0.5 / aq.lab * np.exp(-(x - self.xld) / aq.lab)
            rv[0] = self.aq.coef[self.layers] * qx
        return rv

    def plot(self, ax=None):
        if ax is None:
            _, ax = plt.subplots()
        aq = self.model.aq.find_aquifer_data(self.xld, 0.0)
        ax.plot(
            [self.xld, self.xld],
            [aq.zaqtop[self.layers[0]], aq.zaqbot[self.layers[-1]]],
            "k-",
        )


class ImpLineDoublet1D(LineDoublet1D, DisvecEquation):
    """Create 1D impermeable wall."""

    def __init__(self, model, xld=0, layers=0, label=None):
        self.storeinput(inspect.currentframe())
        LineDoublet1D.__init__(
            self,
            model,
            xld,
            delp=0,
            layers=layers,
            name="ImpLineDoublet1D",
            label=label,
            addtomodel=True,
            res=np.inf,
            aq=None,
        )
        self.nunknowns = self.nparam

    def initialize(self):
        LineDoublet1D.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class LeakyLineDoublet1D(LineDoublet1D, LeakyWallEquation):
    """Create an infinitely long leaky or impermeable wall.

    Parameters
    ----------
    model : Model object
        Model to which the element is added
    xld : scalar
        x-location of line-doublet
    hls : scalar
        head in line-sink
    res : scalar (default is 0)
        resistance of leaky wall. use np.inf to create impermeable wall
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    label: str or None
        label of element
    """

    tiny = 1e-6

    def __init__(self, model, xld=0, res=np.inf, layers=0, label=None):
        self.storeinput(inspect.currentframe())
        LineDoublet1D.__init__(
            self,
            model,
            xld,
            delp=0,
            layers=layers,
            name="LeakyLineDoublet1D",
            label=label,
            addtomodel=True,
            res=res,
            aq=None,
        )
        self.nunknowns = self.nparam

    def initialize(self):
        LineDoublet1D.initialize(self)
        self.xcin = self.xc - self.tiny
        self.xcout = self.xc + self.tiny
        self.ycin = self.yc
        self.ycout = self.yc

    def setparams(self, sol):
        self.parameters[:, 0] = sol

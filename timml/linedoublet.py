import inspect  # Used for storing the input

import matplotlib.pyplot as plt
import numpy as np

from . import bessel
from .controlpoints import controlpoints
from .element import Element
from .equation import DisvecEquation, LeakyWallEquation

__all__ = [
    "ImpLineDoublet",
    "ImpLineDoubletString",
    "LeakyLineDoublet",
    "LeakyLineDoubletString",
]


class LineDoubletHoBase(Element):
    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        delp=0.0,
        res=0.0,
        layers=0,
        order=0,
        name="LineDoubletHoBase",
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
        self.delp = np.atleast_1d(delp).astype("d")
        self.res = np.atleast_1d(res).astype("d")
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
        self.thetaNormOut = (
            np.arctan2(self.y2 - self.y1, self.x2 - self.x1) - np.pi / 2.0
        )
        self.cosnorm = np.cos(self.thetaNormOut) * np.ones(self.ncp)
        self.sinnorm = np.sin(self.thetaNormOut) * np.ones(self.ncp)
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
        self.resfac = self.aq.Haq[self.layers] / self.res
        # also respect addtomodel here to prevent sub-elements (e.g. parts of
        # LineDoubletString) from being added to the aquifer elementlists
        if (addtoaq is None and self.addtomodel) or addtoaq:
            self.aq.add_element(self)
        self.parameters = np.empty((self.nparam, 1))
        # Not sure if this needs to be here
        self.parameters[:, 0] = self.delp

    def potinf(self, x, y, aq=None):
        """Can be called with only one x,y value
        Returns array(nparam, self.aq.naq) with order
        order 0, layer[0]
        order 0, layer[1]
        ...
        order 1, layer[0]
        order 1, layer[1]
        etc
        """
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        if aq == self.aq:
            potrv = rv.reshape(
                (self.order + 1, self.nlayers, aq.naq)
            )  # clever way of using a reshaped rv here
            pot = np.zeros((self.order + 1, aq.naq))
            pot[:, :] = bessel.bessel.potbesldv(
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
        """Can be called with only one x,y value
        Returns array(nparam, self.aq.naq) with order
        order 0, layer[0]
        order 0, layer[1]
        ...
        order 1, layer[0]
        order 1, layer[1]
        etc
        """
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq == self.aq:
            qxqyrv = rv.reshape(
                (2, self.order + 1, self.nlayers, aq.naq)
            )  # clever way of using a reshaped rv here
            qxqy = np.zeros((2 * (self.order + 1), aq.naq))
            qxqy[:, :] = bessel.bessel.disbesldv(
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

    def _refine(self, n=None):
        if n is None:
            n = self.refine_level
        # refine xy
        xy = np.array([(self.x1, self.y1), (self.x2, self.y2)])
        xyr, _ = refine_n_segments(xy, "line", n_segments=n)
        # get input args
        input_args = deepcopy(self._input)
        cls = input_args.pop("__class__", self.__class__)
        input_args["model"] = self.model
        input_args["refine_level"] = 1  # set to 1 to prevent further refinement
        input_args["addtomodel"] = False
        # build new elements
        refined_elements = []
        for ils in range(n):
            (input_args["x1"], input_args["y1"]) = xyr[ils]
            (input_args["x2"], input_args["y2"]) = xyr[ils + 1]
            refined_elements.append(cls(**input_args))
        return refined_elements


class ImpLineDoublet(LineDoubletHoBase, DisvecEquation):
    """
    Create a segment of an impermeable wall, which is
    simulated with a line-doublet

    Parameters
    ----------

    model : Model object
        Model to which the element is added
    x1 : scalar
        x-coordinate of fist point of line-doublet
    y1 : scalar
        y-coordinate of fist point of line-doublet
    x2 : scalar
        x-coordinate of second point of line-doublet
    y2 : scalar
        y-coordinate of second point of line-doublet
    order : int (default is 0)
        polynomial order of potential jump along line-doublet
        (head jump if transmissivity is equal on each side of wall)
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    label: str or None
        label of element

    See Also
    --------

    :class:`.ImpLineDoubletString`

    """

    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        order=0,
        layers=0,
        label=None,
        addtomodel=True,
        refine_level=1,
    ):
        _input = {k: v for k, v in locals().items() if k not in ["self", "model"]}
        self.storeinput(inspect.currentframe())
        LineDoubletHoBase.__init__(
            self,
            model,
            x1,
            y1,
            x2,
            y2,
            delp=0,
            res=np.inf,
            layers=layers,
            order=order,
            name="ImpLineDoublet",
            label=label,
            addtomodel=addtomodel,
            refine_level=refine_level,
        )
        self._input = _input
        self.nunknowns = self.nparam

    def initialize(self):
        LineDoubletHoBase.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class LeakyLineDoublet(LineDoubletHoBase, LeakyWallEquation):
    """
    Create a segment of a leaky wall, which is
    simulated with a line-doublet. The specific discharge through
    the wall is equal to the head difference across the wall
    divided by the resistance of the wall.

    Parameters
    ----------

    model : Model object
        Model to which the element is added
    x1 : scalar
        x-coordinate of fist point of line-doublet
    y1 : scalar
        y-coordinate of fist point of line-doublet
    x2 : scalar
        x-coordinate of second point of line-doublet
    y2 : scalar
        y-coordinate of second point of line-doublet
    res : scalar
        resistance of leaky wall
    order : int (default is 0)
        polynomial order of potential jump along line-doublet
        (head jump if transmissivity is equal on each side of wall)
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    label: str or None
        label of element

    See Also
    --------

    :class:`.LeakyLineDoubletString`

    """

    def __init__(
        self,
        model,
        x1=-1,
        y1=0,
        x2=1,
        y2=0,
        res=0,
        order=0,
        layers=0,
        label=None,
        addtomodel=True,
        refine_level=1,
    ):
        _input = {k: v for k, v in locals().items() if k not in ["self", "model"]}
        self.storeinput(inspect.currentframe())
        LineDoubletHoBase.__init__(
            self,
            model,
            x1,
            y1,
            x2,
            y2,
            delp=0,
            res=res,
            layers=layers,
            order=order,
            name="ImpLineDoublet",
            label=label,
            addtomodel=addtomodel,
            refine_level=refine_level,
        )
        self._input = _input
        self.nunknowns = self.nparam

    def initialize(self):
        LineDoubletHoBase.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class LineDoubletStringBase(Element):
    def __init__(
        self,
        model,
        xy,
        closed=False,
        layers=0,
        order=0,
        res=0,
        name="LineDoubletStringBase",
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
        self.ldlist = []
        self.x, self.y = self.xy[:, 0], self.xy[:, 1]
        self.Nld = len(self.x) - 1
        for i in range(self.Nld):
            self.ldlist.append(
                LineDoubletHoBase(
                    model,
                    x1=self.x[i],
                    y1=self.y[i],
                    x2=self.x[i + 1],
                    y2=self.y[i + 1],
                    delp=0.0,
                    res=res,
                    layers=layers,
                    order=order,
                    label=label,
                    addtomodel=False,
                    aq=aq,
                )
            )

    def __repr__(self):
        return self.name + " with nodes " + str(self.xy)

    def initialize(self):
        for ld in self.ldlist:
            ld.initialize()
        self.ncp = (
            self.Nld * self.ldlist[0].ncp
        )  # Same order for all elements in string
        self.nparam = self.Nld * self.ldlist[0].nparam
        self.nunknowns = self.nparam
        self.xld = np.empty((self.Nld, 2))
        self.yld = np.empty((self.Nld, 2))
        for i, ld in enumerate(self.ldlist):
            self.xld[i, :] = [ld.x1, ld.x2]
            self.yld[i, :] = [ld.y1, ld.y2]
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(
                self.ldlist[0].xc[0], self.ldlist[0].yc[0]
            )
        self.parameters = np.zeros((self.nparam, 1))
        ## As parameters are only stored for the element not the list, we need to combine the following
        self.xc = np.array([ld.xc for ld in self.ldlist]).flatten()
        self.yc = np.array([ld.yc for ld in self.ldlist]).flatten()
        self.xcin = np.array([ld.xcin for ld in self.ldlist]).flatten()
        self.ycin = np.array([ld.ycin for ld in self.ldlist]).flatten()
        self.xcout = np.array([ld.xcout for ld in self.ldlist]).flatten()
        self.ycout = np.array([ld.ycout for ld in self.ldlist]).flatten()
        self.cosnorm = np.array([ld.cosnorm for ld in self.ldlist]).flatten()
        self.sinnorm = np.array([ld.sinnorm for ld in self.ldlist]).flatten()
        self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], self.ycin[0])
        self.aqout = self.model.aq.find_aquifer_data(self.xcout[0], self.ycout[0])
        self.resfac = self.ldlist[0].resfac

    def potinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nld, self.ldlist[0].nparam, aq.naq))
        for i in range(self.Nld):
            rv[i] = self.ldlist[i].potinf(x, y, aq)
        rv.shape = (self.nparam, aq.naq)
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.Nld, self.ldlist[0].nparam, aq.naq))
        for i in range(self.Nld):
            rv[:, i] = self.ldlist[i].disvecinf(x, y, aq)
        rv.shape = (2, self.nparam, aq.naq)
        return rv

    def plot(self, layer=None):
        if (layer is None) or (layer in self.layers):
            plt.plot(self.x, self.y, "k")


class ImpLineDoubletString(LineDoubletStringBase, DisvecEquation):
    """
    Create a string of impermeable wall segements consisting
    of line-doublets

    Parameters
    ----------

    model : Model object
        Model to which the element is added
    xy : array or list
        list or array of (x,y) pairs of coordinates of end-points of
        the segements in the string
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    order : int (default is 0)
        polynomial order of potential jump along line-doublet
        (head jump if transmissivity is equal on each side of wall)
    label: str or None
        label of element

    See Also
    --------

    :class:`.ImpLineDoublet`

    """

    def __init__(self, model, xy=[(-1, 0), (1, 0)], layers=0, order=0, label=None):
        self.storeinput(inspect.currentframe())
        LineDoubletStringBase.__init__(
            self,
            model,
            xy,
            closed=False,
            res=np.inf,
            layers=layers,
            order=order,
            name="ImpLineDoubletString",
            label=label,
            aq=None,
        )
        self.model.add_element(self)

    def initialize(self):
        LineDoubletStringBase.initialize(self)
        self.aq.add_element(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class LeakyLineDoubletString(LineDoubletStringBase, LeakyWallEquation):
    """
    Create a string of leaky wall segements consisting
    of line-doublets

    Parameters
    ----------

    model : Model object
        Model to which the element is added
    xy : array or list
        list or array of (x,y) pairs of coordinates of end-points of
        the segements in the string
    res : scalar
        resistance of leaky wall
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers
    order : int (default is 0)
        polynomial order of potential jump along line-doublet
        (head jump if transmissivity is equal on each side of wall)
    label: str or None
        label of element

    See Also
    --------

    :class:`.ImpLineDoublet`

    """

    def __init__(
        self, model, xy=[(-1, 0), (1, 0)], res=np.inf, layers=0, order=0, label=None
    ):
        self.storeinput(inspect.currentframe())
        LineDoubletStringBase.__init__(
            self,
            model,
            xy,
            closed=False,
            layers=layers,
            order=order,
            res=res,
            name="ImpLineDoubletString",
            label=label,
            aq=None,
        )
        self.model.add_element(self)

    def initialize(self):
        LineDoubletStringBase.initialize(self)
        self.aq.add_element(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol

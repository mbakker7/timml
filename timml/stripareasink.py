"""Cross-section area-sink element.

Models uniform infiltration over along a strip (cross-section).

Example::

    XsectionAreaSink(ml, xleft=-50, xright=50, N=0.001, layer=0)
"""

import warnings

import matplotlib.pyplot as plt
import numpy as np

from timml.element import Element


class XsectionAreaSinkInhom(Element):
    """Create a cross-section area-sink in combination with an inhomogeneity.

    Notes
    -----
    Created automatically using XsectionMaq or Xsection3D.
    Can only be created if top boundary is confined.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    xleft : float
        left boundary of inhomogeneity (may not be -np.inf)
    xright : float
        right boundary of inhomogeneity (may not be np.inf)
    """

    def __init__(
        self,
        model,
        xleft=-1,
        xright=1,
        N=0.001,
        layer=0,
        name="XsectionAreaSink",
        label=None,
    ):
        super().__init__(
            model, nparam=1, nunknowns=0, layers=layer, name=name, label=label
        )
        self.xleft = xleft
        self.xright = xright
        self.N = N
        self.model.add_element(self)

    def __repr__(self):
        return self.name + " between " + str((self.xleft, self.xright))

    def initialize(self):
        self.xc = 0.5 * (self.xleft + self.xright)
        self.yc = 0
        self.L = self.xright - self.xleft
        self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        self.aq.add_element(self)
        self.parameters = np.array([[self.N]])
        if self.aq.ilap:
            self.lab = self.aq.lab[1:]
            self.A = -self.aq.coef[self.layers, 1:] * self.lab**2 / 2
            self.B = self.A * (np.exp(-self.L / self.lab) - 1)
            self.plabsq = self.aq.coef[self.layers, 1:] * self.lab**2
        else:
            print("XsectionAreaSink cannot be added to semi-confined system")

    def potinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        if aq == self.aq:
            if (x > self.xleft) and (x < self.xright):
                rv[0, 0] = -0.5 * (x**2 - 2 * self.xc * x + self.xc**2)
                rv[0, 1:] = self.plabsq
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq == self.aq:
            if (x > self.xleft) and (x < self.xright):
                rv[0, 0, 0] = x - self.xc
        return rv

    def qztop(self, x, y, aq=None):
        rv = 0.0
        if (x > self.xleft) and (x < self.xright):
            rv = -self.parameters[
                0, 0
            ]  # minus cause the parameter is the infiltration rate
        return rv

    def plot(self, ax=None, n_arrows=10, **kwargs):
        if ax is None:
            _, ax = plt.subplots()
        Lz = self.aq.z[0] - self.aq.z[-1]
        Lx = self.xright - self.xleft

        for i in np.linspace(self.xleft, self.xright, n_arrows):
            xtail = i
            ztail = self.aq.z[0] + Lz / 20.0
            dx = 0
            dy = -0.9 * Lz / 20.0
            ax.arrow(
                xtail,
                ztail,
                dx,
                dy,
                width=kwargs.pop("width", Lx / 300.0),
                length_includes_head=kwargs.pop("length_includes_head", True),
                head_width=kwargs.pop("head_width", 4 * Lx / 300.0),
                head_length=kwargs.pop("head_length", 0.4 * Lz / 20.0),
                color=kwargs.pop("color", "k"),
                joinstyle=kwargs.pop("joinstyle", "miter"),
                capstyle=kwargs.pop("capstyle", "projecting"),
            )
        return ax


class XsectionAreaSink(Element):
    def __init__(
        self,
        model,
        xleft=-1,
        xright=1,
        N=0.001,
        layer=0,
        name="XsectionAreaSink",
        label=None,
    ):
        warnings.warn(
            "XsectionAreaSink is only for testing purposes. It is recommended to add "
            "infiltration through XsectionMaq or Xsection3D and specifing 'N'.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(
            model, nparam=1, nunknowns=0, layers=layer, name=name, label=label
        )
        self.xleft = xleft
        self.xright = xright
        self.N = N
        self.model.add_element(self)

    def __repr__(self):
        return self.name + " at " + str((self.xc, self.yc))

    def initialize(self):
        self.xc = 0.5 * (self.xleft + self.xright)
        self.yc = 0
        self.L = self.xright - self.xleft
        self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        self.aq.add_element(self)
        self.parameters = np.array([[self.N]])
        if self.aq.ilap:
            self.lab = self.aq.lab[1:]
            self.A = -self.aq.coef[self.layers, 1:] * self.lab**2 / 2
            self.B = self.A * (np.exp(-self.L / self.lab) - 1)
            self.plabsq = self.aq.coef[self.layers, 1:] * self.lab**2
        else:
            print("XsectionAreaSink cannot be added to semi-confined system")

    def potinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        if aq == self.aq:
            if x < self.xleft:
                rv[0, 0] = -(self.xleft - self.xc) * (x - self.xc) + self.L**2 / 8
                rv[0, 1:] = self.B * np.exp((x - self.xleft) / self.lab)
            elif x > self.xright:
                rv[0, 0] = -(self.xright - self.xc) * (x - self.xc) + self.L**2 / 8
                rv[0, 1:] = self.B * np.exp(-(x - self.xright) / self.lab)
            else:
                rv[0, 0] = -0.5 * (x**2 - 2 * self.xc * x + self.xc**2)
                rv[0, 1:] = (
                    self.A
                    * (
                        np.exp(-(x - self.xleft) / self.lab)
                        + np.exp((x - self.xright) / self.lab)
                    )
                    + self.plabsq
                )
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
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
                rv[0, 0, 1:] = (
                    self.A
                    / self.lab
                    * (
                        np.exp(-(x - self.xleft) / self.lab)
                        - np.exp((x - self.xright) / self.lab)
                    )
                )
        return rv

    def qztop(self, x, y, aq=None):
        rv = 0.0
        if (x > self.xleft) and (x < self.xright):
            rv = -self.parameters[
                0, 0
            ]  # minus cause the parameter is the infiltration rate
        return rv

    def changetrace(
        self, xyzt1, xyzt2, aq, layer, ltype, modellayer, direction, hstepmax
    ):
        changed = False
        terminate = False
        xyztnew = 0
        message = None
        eps = 1e-4
        x1 = xyzt1[0]
        x2 = xyzt2[0]
        if (x1 < self.xleft and x2 >= self.xleft) or (
            x1 > self.xleft and x2 <= self.xleft
        ):
            dx1 = x1 - self.xleft
            dx2 = x2 - self.xleft
            # step just beyond boundary
            xyztnew = xyzt1 + (dx1 + np.sign(dx1) * eps) / (dx1 - dx2) * (xyzt2 - xyzt1)
            changed = True
            message = "exited or entered xsection area sink on the left side"
        elif (x1 < self.xright and x2 >= self.xright) or (
            x1 > self.xright and x2 <= self.xright
        ):  # type: ignore
            dx1 = x1 - self.xright
            dx2 = x2 - self.xright
            # step just beyond boundary
            xyztnew = xyzt1 + (dx1 + np.sign(dx1) * eps) / (dx1 - dx2) * (xyzt2 - xyzt1)
            changed = True
            message = "exited or entered area sink on the right side"
        return changed, terminate, xyztnew, message

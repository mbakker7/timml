import inspect  # Used for storing the input

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import k0, k1

from .element import Element
from .equation import MscreenWellEquation, PotentialEquation
from .trace import timtracelines

__all__ = ["WellBase", "Well", "HeadWell"]


class WellBase(Element):
    def __init__(
        self,
        model,
        xw=0,
        yw=0,
        Qw=100.0,
        rw=0.1,
        res=0.0,
        layers=0,
        name="WellBase",
        label=None,
    ):
        Element.__init__(
            self, model, nparam=1, nunknowns=0, layers=layers, name=name, label=label
        )
        # Defined here and not in Element as other elements can have multiple
        # parameters per layers
        self.nparam = len(self.layers)
        self.xw = float(xw)
        self.yw = float(yw)
        self.Qw = np.atleast_1d(Qw)
        self.rw = float(rw)
        self.res = float(res)
        self.model.add_element(self)

    def __repr__(self):
        return self.name + " at " + str((self.xw, self.yw))

    def initialize(self):
        self.xc = np.array([self.xw + self.rw])
        self.yc = np.array([self.yw])
        self.ncp = 1
        self.aq = self.model.aq.find_aquifer_data(self.xw, self.yw)
        self.aq.add_element(self)
        self.parameters = np.empty((self.nparam, 1))
        self.parameters[:, 0] = self.Qw
        self.resfac = self.res / (2 * np.pi * self.rw * self.aq.Haq[self.layers])

    def potinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        if aq == self.aq:
            pot = np.zeros(aq.naq)
            r = np.sqrt((x - self.xw) ** 2 + (y - self.yw) ** 2)
            if r < self.rw:
                r = self.rw  # If at well, set to at radius
            if aq.ilap:
                pot[0] = np.log(r / self.rw) / (2 * np.pi)
                pot[1:] = -k0(r / aq.lab[1:]) / (2 * np.pi)
            else:
                pot[:] = -k0(r / aq.lab) / (2 * np.pi)
            rv[:] = self.aq.coef[self.layers] * pot
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq == self.aq:
            qx = np.zeros(aq.naq)
            qy = np.zeros(aq.naq)
            rsq = (x - self.xw) ** 2 + (y - self.yw) ** 2
            r = np.sqrt(rsq)
            xminxw = x - self.xw
            yminyw = y - self.yw
            if r < self.rw:
                r = self.rw
                rsq = r ** 2
                xminxw = self.rw
                yminyw = 0.0
            if aq.ilap:
                qx[0] = -1 / (2 * np.pi) * xminxw / rsq
                qy[0] = -1 / (2 * np.pi) * yminyw / rsq
                kone = k1(r / aq.lab[1:])
                qx[1:] = -kone * xminxw / (r * aq.lab[1:]) / (2 * np.pi)
                qy[1:] = -kone * yminyw / (r * aq.lab[1:]) / (2 * np.pi)
            else:
                kone = k1(r / aq.lab)
                qx[:] = -kone * xminxw / (r * aq.lab) / (2 * np.pi)
                qy[:] = -kone * yminyw / (r * aq.lab) / (2 * np.pi)
            rv[0] = self.aq.coef[self.layers] * qx
            rv[1] = self.aq.coef[self.layers] * qy
        return rv

    def headinside(self):
        """The head inside the well

        Returns
        -------
        array (length number of screens)
            Head inside the well for each screen

        """

        h = self.model.head(self.xw + self.rw, self.yw, layers=self.layers)
        return h - self.resfac * self.parameters[:, 0]

    def discharge(self):
        """The discharge in each layer

        Returns
        -------
        array (length number of layers)
            Discharge in each screen with zeros for layers that are not
            screened

        """

        Q = np.zeros(self.aq.naq)
        Q[self.layers] = self.parameters[:, 0]
        return Q

    def changetrace(
        self, xyzt1, xyzt2, aq, layer, ltype, modellayer, direction, hstepmax
    ):
        changed = False
        terminate = False
        xyztnew = 0
        message = None
        if np.sqrt((xyzt2[0] - self.xw) ** 2 + (xyzt2[1] - self.yw) ** 2) < (
            hstepmax + self.rw
        ):
            if ltype == "a":
                if (layer == self.layers).any():  # in layer where well is screened
                    if (self.discharge()[layer] > 0 and direction > 0) or (
                        self.discharge()[layer] < 0 and direction < 0
                    ):
                        vx, vy, vz = self.model.velocity(*xyzt1[:-1])
                        tstep = np.sqrt(
                            (xyzt1[0] - self.xw) ** 2 + (xyzt1[1] - self.yw) ** 2
                        ) / np.sqrt(vx ** 2 + vy ** 2)
                        xnew = self.xw
                        ynew = self.yw
                        znew = xyzt1[2] + tstep * vz * direction
                        tnew = xyzt1[3] + tstep
                        xyztnew = np.array([xnew, ynew, znew, tnew])
                        changed = True
                        terminate = True
        if terminate:
            message = "reached element of type well"
            if self.label:
                message += " ({lab})".format(lab=self.label)
        return changed, terminate, [xyztnew], message

    def capzone(
        self,
        nt=10,
        zstart=None,
        hstepmax=10,
        vstepfrac=0.2,
        tmax=None,
        nstepmax=100,
        silent=".",
        *,
        metadata=False
    ):
        """Compute a capture zone

        Parameters
        ----------
        nt : int
            number of path lines
        zstart : scalar
            starting elevation of the path lines
        hstepmax : scalar
            maximum step in horizontal space
        vstepfrac : float
            maximum fraction of aquifer layer thickness during one step
        tmax : scalar
            maximum time
        nstepmax : scalar(int)
            maximum number of steps
        silent : boolean or string
            True (no messages), False (all messages), or '.'
            (print dot for each path line)

        Returns
        -------
        xyzt : list of arrays of x, y, z, and t values

        """
        xstart, ystart, zstart = self.capzonestart(nt, zstart)
        xyzt = timtracelines(
            self.model,
            xstart,
            ystart,
            zstart,
            -np.abs(hstepmax),
            vstepfrac=vstepfrac,
            tmax=tmax,
            nstepmax=nstepmax,
            silent=silent,
            metadata=metadata,
        )
        return xyzt

    def capzonestart(self, nt, zstart):
        eps = 1e-1
        angle = np.arange(eps, 2 * np.pi, 2 * np.pi / nt)
        xstart = self.xw + (1 + eps) * self.rw * np.cos(angle)
        ystart = self.yw + (1 + eps) * self.rw * np.sin(angle)
        if zstart is None:
            zstart = self.aq.zaqbot[self.layers[0]] + 0.5 * self.aq.Haq[self.layers[0]]
        zstart = zstart * np.ones(nt)
        return xstart, ystart, zstart

    def plot(self):
        plt.plot(self.xw, self.yw, "k.")

    def plotcapzone(
        self,
        nt=10,
        zstart=None,
        hstepmax=20,
        vstepfrac=0.2,
        tmax=365,
        nstepmax=100,
        silent=".",
        color=None,
        orientation="hor",
        win=[-1e30, 1e30, -1e30, 1e30],
        newfig=False,
        figsize=None,
        *,
        return_traces=False,
        metadata=False
    ):
        """Plot a capture zone

        Parameters
        ----------
        nt : int
            number of path lines
        zstart : scalar
            starting elevation of the path lines
        hstepmax : scalar
            maximum step in horizontal space
        vstepfrac : float
            maximum fraction of aquifer layer thickness during one step
        tmax : scalar
            maximum time
        nstepmax : scalar(int)
            maximum number of steps
        silent : boolean or string
            True (no messages), False (all messages), or '.'
            (print dot for each path line)
        color : color
        orientation : string
            'hor' for horizontal, 'ver' for vertical, or 'both' for both
        win : array_like (length 4)
            [xmin, xmax, ymin, ymax]
        newfig : boolean (default False)
            boolean indicating if new figure should be created
        figsize : tuple of integers, optional, default: None
            width, height in inches.

        """
        if not return_traces:
            metadata = True  # suppress future warning from timtraceline
        xstart, ystart, zstart = self.capzonestart(nt, zstart)
        traces = self.model.tracelines(
            xstart,
            ystart,
            zstart,
            hstepmax=-abs(hstepmax),
            vstepfrac=vstepfrac,
            tmax=tmax,
            nstepmax=nstepmax,
            silent=silent,
            color=color,
            orientation=orientation,
            win=win,
            newfig=newfig,
            figsize=figsize,
            return_traces=return_traces,
            metadata=metadata,
        )
        if return_traces:
            return traces


class Well(WellBase, MscreenWellEquation):
    """
    Well Class to create a well with a specified discharge. The well
    may be screened in multiple layers. The resistance of the screen may
    be specified. The head is computed such that the discharge :math:`Q_i`
    in layer :math:`i` is computed as

    .. math::
        Q_i = 2\pi r_w(h_i - h_w)/c

    where :math:`c` is the resistance of the well screen and :math:`h_w` is
    the head inside the well. The total discharge is distributed over the
    screens such that :math:`h_w` is the same in each screened layer.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    xw : float
        x-coordinate of the well
    yw : float
        y-coordinate of the well
    Qw : float
        total discharge of the well
    rw : float
        radius of the well
    res : float
        resistance of the well screen
    layers : int, array or list
        layer (int) or layers (list or array) where well is screened
    label : string or None (default: None)
        label of the well

    Examples
    --------
    >>> ml = Model3D(kaq=10, z=np.arange(20, -1, -2), kzoverkh=0.1)
    >>> Well(ml, 100, 200, 1000, layers=[0, 1, 2, 3])

    """

    def __init__(self, model, xw=0, yw=0, Qw=100.0, rw=0.1, \
                 res=0.0, layers=0, label=None):
        self.storeinput(inspect.currentframe())
        WellBase.__init__(self, model, xw, yw, Qw, rw, res, \
                          layers=layers, name="Well", label=label)
        self.Qc = float(Qw)
        if self.nlayers == 1:
            self.nunknowns = 0
        else:
            self.nunknowns = self.nparam

    def initialize(self):
        WellBase.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class HeadWell(WellBase, PotentialEquation):
    """
    HeadWell Class to create a well with a specified head inside the well.
    The well may be screened in multiple layers. The resistance of the screen
    may be specified. The head is computed such that the discharge :math:`Q_i`
    in layer :math:`i` is computed as

    .. math::
        Q_i = 2\pi r_w(h_i - h_w)/c

    where :math:`c` is the resistance of the well screen and :math:`h_w` is
    the head inside the well.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    xw : float
        x-coordinate of the well
    yw : float
        y-coordinate of the well
    hw : float
        head inside the well
    rw : float
        radius of the well
    res : float
        resistance of the well screen
    layers : int, array or list
        layer (int) or layers (list or array) where well is screened
    label : string (default: None)
        label of the well

    """

    def __init__(self, model, xw=0, yw=0, hw=10, rw=0.1, res=0, layers=0, label=None):
        self.storeinput(inspect.currentframe())
        WellBase.__init__(
            self,
            model,
            xw,
            yw,
            0.0,
            rw,
            res,
            layers=layers,
            name="HeadWell",
            label=label,
        )
        self.hc = hw
        self.nunknowns = self.nparam

    def initialize(self):
        WellBase.initialize(self)
        self.pc = self.hc * self.aq.T[self.layers]  # Needed in solving

    def setparams(self, sol):
        self.parameters[:, 0] = sol

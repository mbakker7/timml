import inspect  # Used for storing the input

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import k0, k1

from .element import Element
from .equation import HeadEquation, MscreenWellNoflowEquation
from .linesink import LineSinkDitchString
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
        xc=None,
        yc=None,
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
        self.xc = xc
        self.yc = yc
        self.model.add_element(self)

    def __repr__(self):
        return self.name + " at " + str((self.xw, self.yw))

    def initialize(self):
        if self.xc is None:
            self.xc = np.array([self.xw + self.rw])
        else:
            self.xc = np.atleast_1d(self.xc)
        if self.yc is None:
            self.yc = np.array([self.yw])
        else:
            self.yc = np.atleast_1d(self.yc)
        self.ncp = 1
        self.aq = self.model.aq.find_aquifer_data(self.xw, self.yw)
        self.aq.add_element(self)
        self.parameters = np.empty((self.nparam, 1))
        self.parameters[:, 0] = self.Qw
        self.resfac = self.res / (2 * np.pi * self.rw * self.aq.Haq[self.layers])
        self.resfac = self.resfac * np.identity(self.nlayers)
        self.resfac.shape = (
            self.ncp,
            self.nlayers,
            self.nlayers,  # changed to nlayers from nunknowns
        )  # required shape for HeadEquation

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
                rsq = r**2
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
        """The head inside the well.

        Returns
        -------
        array (length number of screens)
            Head inside the well for each screen
        """
        icp = 0  # there is only one control point
        hinside = self.model.head(self.xc[icp], self.yc[icp], self.layers) - np.sum(
            self.resfac[icp] * self.parameters[:, 0], 1
        )
        return hinside

    def discharge(self):
        """The discharge in each layer.

        Returns
        -------
        array (length number of layers)
            Discharge in each screen with zeros for layers that are not screened
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
                        ) / np.sqrt(vx**2 + vy**2)
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
        metadata=False,
    ):
        """Compute a capture zone.

        Parameters
        ----------
        nt : int
            number of path lines
        zstart : scalar or None
            starting elevation of the path lines, middle of aquifer if None
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

    def plot(self, layer=None):
        if (layer is None) or np.isin(layer, self.layers).any():
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
        win=None,
        newfig=False,
        figsize=None,
        *,
        return_traces=False,
        metadata=False,
    ):
        """Plot a capture zone.

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
        if win is None:
            win = [-1e30, 1e30, -1e30, 1e30]
        if not return_traces:
            metadata = True  # suppress future warning from timtraceline
        xstart, ystart, zstart = self.capzonestart(nt, zstart)
        traces = self.model.plots.tracelines(
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


class Well(WellBase):
    r"""Well Class to create a well with a specified discharge.

    Notes
    -----
    The well may be screened in multiple layers.
    The resistance of the screen may be specified.
    The head is computed such that the discharge :math:`Q_i` in layer :math:`i`
    is computed as.

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
    xc : float
        x-location of control point (default None, which puts it at xw)
    yc : float
        y-location of control point (default None, which puts it at yw + rw)

    Examples
    --------
    >>> ml = Model3D(kaq=10, z=np.arange(20, -1, -2), kzoverkh=0.1)
    >>> Well(ml, 100, 200, 1000, layers=[0, 1, 2, 3])
    """

    def __init__(
        self,
        model,
        xw=0,
        yw=0,
        Qw=100.0,
        rw=0.1,
        res=0.0,
        layers=0,
        label=None,
        xc=None,
        yc=None,
    ):
        self.storeinput(inspect.currentframe())
        WellBase.__init__(
            self,
            model,
            xw,
            yw,
            Qw,
            rw,
            res,
            layers=layers,
            name="Well",
            label=label,
            xc=xc,
            yc=yc,
        )
        self.hc = np.zeros(
            self.nlayers
        )  # needed as heads are same in all screened layers
        self.Qc = float(Qw)
        if self.nlayers == 1:
            self.nunknowns = 0
        else:
            self.nunknowns = self.nparam

    def initialize(self):
        WellBase.initialize(self)

    def equation(self):
        mat, rhs = HeadEquation.equation(self)
        for i in range(1, self.nunknowns):
            mat[i] -= mat[0]
            rhs[i] -= rhs[0]
        # first equation is sum of discharges equals Qw
        mat[0] = 0
        ieq = 0
        for e in self.model.elementlist:
            if e.nunknowns > 0:
                if e == self:
                    mat[0, ieq : ieq + e.nunknowns] = 1.0
                    break
                ieq += e.nunknowns
        rhs[0] = self.Qw
        return mat, rhs

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class HeadWell(WellBase, HeadEquation):
    r"""HeadWell Class to create a well with a specified head inside the well.

    Notes
    -----
    The well may be screened in multiple layers.
    The resistance of the screen may be specified.
    The head is computed such that the discharge :math:`Q_i` in layer :math:`i` is
    computed as:

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
    xc : float
        x-location of control point (default None, which puts it at xw)
    yc : float
        y-location of control point (default None, which puts it at yw + rw)
    """

    def __init__(
        self,
        model,
        xw=0,
        yw=0,
        hw=10,
        rw=0.1,
        res=0,
        layers=0,
        label=None,
        xc=None,
        yc=None,
    ):
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
            xc=xc,
            yc=yc,
        )
        self.nunknowns = self.nparam
        self.hc = hw * np.ones(self.nunknowns)

    def initialize(self):
        WellBase.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class LargeDiameterWell(WellBase, MscreenWellNoflowEquation):
    """Experimental class for radial flow to large diameter well.

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
    xc : float
        x-location of control point (default None, which puts it at xw)
    yc : float
        y-location of control point (default None, which puts it at yw + rw)

    Examples
    --------
    >>> ml = Model3D(kaq=10, z=np.arange(20, -1, -2), kzoverkh=0.1)
    >>> Well(ml, 100, 200, 1000, layers=[0, 1, 2, 3])
    """

    def __init__(
        self,
        model,
        xw=0,
        yw=0,
        Qw=100.0,
        rw=0.1,
        res=0.0,
        layers=0,
        label=None,
        xc=None,
        yc=None,
    ):
        self.storeinput(inspect.currentframe())
        WellBase.__init__(
            self,
            model,
            xw,
            yw,
            Qw,
            rw,
            res,
            layers=np.arange(model.aq.nlayers),
            name="Well",
            label=label,
            xc=xc,
            yc=yc,
        )
        self.Qc = float(Qw)
        self.screened = layers  # layers where well is screened
        self.nscreened = len(self.screened)
        if self.nlayers == 1:
            self.nunknowns = 0
        else:
            self.nunknowns = self.nparam

    def initialize(self):
        WellBase.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol

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
                pot[1:] = -k0(r / aq.lab[1:]) / (2 * np.pi) / k0(self.rw / aq.lab[1:])
            else:
                pot[:] = -k0(r / aq.lab) / (2 * np.pi) / k0(self.rw / aq.lab)
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
                rsq = r**2
                xminxw = self.rw
                yminyw = 0.0
            if aq.ilap:
                qx[0] = -1 / (2 * np.pi) * xminxw / rsq
                qy[0] = -1 / (2 * np.pi) * yminyw / rsq
                kone = k1(r / aq.lab[1:]) / k0(self.rw / aq.lab[1:])
                qx[1:] = -kone * xminxw / (r * aq.lab[1:]) / (2 * np.pi)
                qy[1:] = -kone * yminyw / (r * aq.lab[1:]) / (2 * np.pi)
            else:
                kone = k1(r / aq.lab) / k0(self.rw / aq.lab)
                qx[:] = -kone * xminxw / (r * aq.lab) / (2 * np.pi)
                qy[:] = -kone * yminyw / (r * aq.lab) / (2 * np.pi)
            rv[0] = self.aq.coef[self.layers] * qx
            rv[1] = self.aq.coef[self.layers] * qy
        return rv


class CollectorWell(LineSinkDitchString):
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
    >>> ml = timml.Model3D(kaq=10, z=np.arange(20, -1, -2), kzoverkh=0.1)
    >>> xy = [(1, 0, 10, 0), (0, 1, 0, 10)]
    >>> w = timml.CollectorWell(ml, xy=xy, Qw=1000, layers=np.arange(5, 10))
    >>> ml.solve()
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
    >>> ml = timml.Model3D(kaq=10, z=np.arange(20, -1, -2), kzoverkh=0.1)
    >>> w = timml.RadialCollectorWell(ml, x=0, y=0, narms=5, nls=10, angle=0,
    ... rcaisson=2.0, rw=0.1, Qw=1000, layers=5)
    >>> ml.solve()
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
        self.nls = nls
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

    # def compute_xyold(self, x, y, rcaisson, L, narms, nls):
    #     """Compute the x,y-coordinates array for the radial collector well.

    #     Parameters
    #     ----------
    #     x : float
    #         x-coordinate of the center of the collector well
    #     y : float
    #         y-coordinate of the center of the collector well
    #     rcaisson : float
    #         radius of the caisson
    #     L : float
    #         length of each arm
    #     narms : int
    #         number of arms
    #     nls : int
    #         number of line sinks per arm

    #     Returns
    #     -------
    #     xy : np.array
    #         array of shape (N, 4) with start and end coordinates of the line sinks
    #         on each row: [(x1, y1, x2, y2), ...]
    #     """
    #     xy = np.empty((narms * nls, 4))
    #     for i, theta in enumerate(np.arange(0, 2 * np.pi, 2 * np.pi / 5)):
    #         x = rcaisson * np.cos(theta) + np.linspace(0, L, nls + 1) * np.cos(theta)
    #         y = rcaisson * np.sin(theta) + np.linspace(0, L, nls + 1) * np.sin(theta)
    #         xy[i * nls : (i + 1) * nls, 0] = x[:-1]
    #         xy[i * nls : (i + 1) * nls, 1] = y[:-1]
    #         xy[i * nls : (i + 1) * nls, 2] = x[1:]
    #         xy[i * nls : (i + 1) * nls, 3] = y[1:]
    #     return xy

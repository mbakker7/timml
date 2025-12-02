import inspect  # Used for storing the input

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import k0, k1

from .element import Element
from .equation import HeadEquation, MscreenWellNoflowEquation
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
        addtomodel=True,
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
        if addtomodel:
            self.model.add_element(self)
        self.addtomodel = addtomodel

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
        if self.addtomodel:
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
        addtomodel=True,
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
            addtomodel=addtomodel,
        )
        self.nunknowns = self.nparam
        self.hc = hw * np.ones(self.nunknowns)

    def initialize(self):
        WellBase.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class TargetHeadWell(WellBase):
    r"""TargetHeadWell is a well with a specified head at (layer, x, y).

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
    rw : float
        radius of the well
    res : float
        resistance of the well screen
    layers : int, array or list
        layer (int) or layers (list or array) where well is screened
    hcp : float
        specified head at control point
    xcp : float
        x-coordinate where head is specified
    ycp : float
        y-coordinate where head is specified
    lcp : int
        layer in which head is specified, default is the first layer, layer 0
    label : string (default: None)
        label of the well
    """

    def __init__(
        self,
        model,
        xw=0,
        yw=0,
        rw=0.1,
        res=0,
        layers=0,
        hcp=10,
        xcp=10,
        ycp=10,
        lcp=0,
        label=None,
        addtomodel=True,
    ):
        self.storeinput(inspect.currentframe())
        super().__init__(
            model,
            xw,
            yw,
            0.0,
            rw,
            res,
            layers=layers,
            name="TargetHeadWell",
            label=label,
            addtomodel=addtomodel,
        )
        self.nunknowns = self.nparam
        self.hcp = hcp
        self.hc = hcp * np.ones(self.nunknowns)  # needed for HeadEquation
        self.xcp = xcp
        self.ycp = ycp
        self.lcp = np.atleast_1d(lcp)  # layer of control point for specified head

    def initialize(self):
        WellBase.initialize(self)

    def equation(self):
        mat, rhs = HeadEquation.equation(self)
        for i in range(1, self.nunknowns):
            mat[i] -= mat[0]
            rhs[i] -= rhs[0]
        # first equation is head at control point equals hc
        mat[0] = 0.0
        rhs[0] = self.hcp
        aq = self.model.aq.find_aquifer_data(self.xcp, self.ycp)
        ieq = 0
        for e in self.model.elementlist:
            if e.nunknowns > 0:
                mat[0:1, ieq : ieq + e.nunknowns] += (
                    e.potinflayers(self.xcp, self.ycp, self.lcp) / aq.Tcol[self.lcp]
                )
                ieq += e.nunknowns
            else:
                rhs[0] -= (
                    e.potentiallayers(self.xcp, self.ycp, self.lcp) / aq.T[self.lcp]
                )
        return mat, rhs

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


class WellStringBase(Element):
    """Base class for multiple connected wells.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    xy : array_like
        array of x and y coordinates of the well nodes
    layers : int, array or list
        layer (int) or layers (list or array) where well is screened
    name : string
        name of the element
    label : string or None (default: None)
        label of the well string
    aq : AquiferData object or None
        aquifer data object or None (default is None)
    """

    def __init__(self, model, xy, layers=0, name="WellStringBase", label=None, aq=None):
        super().__init__(model, nparam=1, nunknowns=0, layers=0, name=name, label=label)
        self.xy = np.array(xy)
        self.xw = self.xy[:, 0]
        self.yw = self.xy[:, 1]
        self.nw = len(self.xw)
        # convert layers to layers per well
        if isinstance(layers, (int, np.integer)):
            layers = layers * np.ones((self.nw, 1), dtype="int")
        elif isinstance(layers, (list, tuple)):
            # first try to coerce into an array
            try:
                nlayers = max(len(layers[i]) for i in range(self.nw))
                layers = np.array(layers) * np.ones((self.nw, nlayers), dtype="int")
            except ValueError:  # list/tuple of different size tuples
                pass
                # raise ValueError(
                #     "layers must be int, list of tuples of equal shape, or array of "
                #     "shape (nwells, nlayers)"
                # ) from e
            except TypeError:  # list/tuple of int
                layers = np.atleast_1d(layers) * np.ones((self.nw, 1), dtype="int")

        elif isinstance(layers, np.ndarray):
            if layers.ndim == 1:
                layers = layers * np.ones((self.nw, layers.shape[0]), dtype="int")
            else:
                assert layers.shape[0] == self.nw, (
                    "layers array must be shape (nwells, nlayers)"
                )
        self.layers = layers
        if isinstance(self.layers, np.ndarray):
            self.nlayers = max(self.layers[i : i + 1].shape[1] for i in range(self.nw))
        else:
            self.nlayers = max(len(self.layers[i]) for i in range(self.nw))
        self.wlist = []

    def __repr__(self):
        return self.name + " with nodes " + str(self.xy)

    def initialize(self):
        for w in self.wlist:
            w.initialize()

        self.aq = []
        for w in self.wlist:
            if w.aq not in self.aq:
                self.aq.append(w.aq)
        for aq in self.aq:
            aq.add_element(self)

        self.nparam = sum(ls.nparam for ls in self.wlist)
        self.nunknowns = self.nparam
        self.parameters = np.zeros((self.nparam, 1))

    def potinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        if aq in self.aq:
            j = 0
            for w in self.wlist:
                rv[j : j + w.nparam] = w.potinf(x, y, aq)
                j += w.nparam
        # rv.shape = (self.nparam, aq.naq)
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq in self.aq:
            j = 0
            for w in self.wlist:
                rv[:, j : j + w.nparam] = w.disvecinf(x, y, aq)
                j += w.nparam
        return rv

    def discharge(self):
        """Discharge of the element in each layer.

        Returns
        -------
        array (length naq)
            Discharge in each layer
        """
        return self.discharge_per_well().sum(axis=1)

    def discharge_per_well(self):
        """Discharge per well in each layer.

        Returns
        -------
        array (nlay, nwells)
            Discharge per well in each layer
        """
        Q = np.zeros((self.model.aq.naq, self.nw))
        j = 0
        for i, w in enumerate(self.wlist):
            Q[w.layers, i : i + 1] = self.parameters[j : j + w.nlayers]
            j += w.nlayers
        return Q

    def plot(self, layer=None):
        for iw, w in enumerate(self.wlist):
            if (layer is None) or (layer in self.layers[iw]):
                plt.plot(w.xw, w.yw, "k.")

    def equation(self):
        mat = np.zeros((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)
        ieq = 0
        for w in self.wlist:
            imat, irhs = w.equation()  # HeadWell equation
            mat[ieq : ieq + w.nunknowns] = imat
            rhs[ieq : ieq + w.nunknowns] = irhs
            ieq = ieq + w.nunknowns

        # include resistance by finding position of element in coefficient matrix
        # and subtracting resfac (resistance factor).
        iself = self.model.elementlist.index(self)
        jcol = np.sum(e.nunknowns for e in self.model.elementlist[:iself])
        irow = 0
        for w in self.wlist:
            mat[irow : irow + w.nlayers, jcol : jcol + w.nunknowns] -= w.resfac[
                0
            ]  # only one control point
            irow += w.nlayers
            jcol += w.nunknowns

        return mat, rhs

    def headinside(self):
        """The head inside the wells.

        Returns
        -------
        head: float
            head inside the well
        """
        return self.wlist[0].headinside()[0]


class WellString(WellStringBase):
    """
    WellString is a string of wells for which the total discharge is specified.

    The head inside the wells is equal but unknown, and the total discharge of the
    wells is specified.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    xy : list of tuples or np.ndarray
        list of (x, y) tuples or 2d array of x and y coordinates
        of the well nodes
    Qw : float
        total discharge of the well string
    rw : float
        radius of the wells
    res : float
        resistance of the well screens
    layers : int, array or list
        layer (int) or layers (list or array) where well is screened
    label : string or None (default: None)
        label of the well string
    """

    def __init__(
        self,
        model,
        xy,
        Qw=100,
        rw=0.1,
        res=0.0,
        layers=0,
        label=None,
    ):
        super().__init__(model, xy, layers=layers, name="WellString", label=label)
        self.Qw = float(Qw)
        self.rw = rw
        self.res = res
        self.model.add_element(self)

    def initialize(self):
        self.wlist = []
        for i in range(self.nw):
            self.wlist.append(
                HeadWell(
                    self.model,
                    xw=self.xw[i],
                    yw=self.yw[i],
                    hw=0.0,
                    rw=self.rw,
                    res=self.res,
                    layers=self.layers[i],
                    addtomodel=False,
                )
            )
        WellStringBase.initialize(self)

    def equation(self):
        mat, rhs = WellStringBase.equation(self)
        # print(rhs)
        for i in range(1, self.nunknowns):
            mat[i] -= mat[0]
            rhs[i] -= rhs[0]

        # first equation is sum of discharges equals Qw
        mat[0] = 0
        ieq = 0
        for e in self.model.elementlist:
            if e.nunknowns > 0:
                if e == self:
                    mat[0, ieq : ieq + self.nunknowns] = 1.0
                    break
                ieq += e.nunknowns
        rhs[0] = self.Qw
        return mat, rhs

    def setparams(self, sol):
        self.parameters[:, 0] = sol
        # assign parameters to individual wells
        i = 0
        for w in self.wlist:
            w.parameters[:, 0] = sol[i : i + w.nparam]
            i += w.nparam


class HeadWellString(WellStringBase):
    """
    HeadWellString is a string of wells for which the head is specified in the wells.

    This element is identical to a series of HeadWell elements.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    xy : list of tuples or np.ndarray
        list of (x, y) tuples or 2d array of x and y coordinates of the well nodes
    hw : float
        head in the well(s)
    rw : float
        radius of the wells
    res : float
        resistance of the well screens
    layers : int, array or list
        layer (int) or layers (list or array) where well is screened
    label : string or None (default: None)
        label of the well string
    """

    def __init__(
        self,
        model,
        xy,
        hw=10,
        rw=0.1,
        res=0.0,
        layers=0,
        label=None,
    ):
        super().__init__(model, xy, layers=layers, name="HeadWellString", label=label)

        self.hw = float(hw)
        self.rw = rw
        self.res = res
        self.model.add_element(self)

    def initialize(self):
        self.wlist = []
        for i in range(self.nw):
            self.wlist.append(
                HeadWell(
                    self.model,
                    xw=self.xw[i],
                    yw=self.yw[i],
                    hw=self.hw,
                    rw=self.rw,
                    res=self.res,
                    layers=self.layers[i],
                    addtomodel=False,
                )
            )
        super().initialize()

    def setparams(self, sol):
        self.parameters[:, 0] = sol
        # assign parameters to individual wells
        i = 0
        for w in self.wlist:
            w.parameters[:, 0] = sol[i : i + w.nparam]
            i += w.nparam


class TargetHeadWellString(WellStringBase):
    """
    A string of wells for which the head is specified at a point.

    The head in the wells is equal but unknown and the head at the control point
    (lc, xc, yc) must equal the specified head.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    xy : list of tuples or np.ndarray
        list of (x, y) tuples or 2d array of x and y coordinates of the wells
    rw : float
        radius of the wells
    res : float
        resistance of the well screens
    layers : int, array or list
        layer (int) or layers (list or array) where well is screened
    hcp : float
        head at the control point (lcp, xcp, ycp)
    xcp : float
        x-coordinate where head is specified
    ycp : float
        y-coordinate where head is specified
    lcp : int
        layer in which head is specified, default is the first layer, layer 0
    label : string or None (default: None)
        label of the well string
    """

    def __init__(
        self,
        model,
        xy,
        rw=0.1,
        res=0.0,
        layers=0,
        hcp=10,
        xcp=None,
        ycp=None,
        lcp=0,
        label=None,
    ):
        super().__init__(
            model, xy, layers=layers, name="TargetHeadWellString", label=label
        )

        self.rw = rw
        self.res = res
        self.hcp = float(hcp)
        self.xcp = xcp
        self.ycp = ycp
        self.lcp = np.atleast_1d(lcp)
        self.model.add_element(self)

    def initialize(self):
        self.wlist = []
        for i in range(self.nw):
            self.wlist.append(
                HeadWell(
                    self.model,
                    xw=self.xw[i],
                    yw=self.yw[i],
                    hw=0.0,
                    rw=self.rw,
                    res=self.res,
                    layers=self.layers[i],
                    addtomodel=False,
                )
            )
        super().initialize()

    def equation(self):
        mat, rhs = super().equation()
        for i in range(1, self.nunknowns):
            mat[i] -= mat[0]
            rhs[i] -= rhs[0]
        # first equation is head at control point equals hw
        mat[0] = 0
        rhs[0] = self.hcp
        aq = self.model.aq.find_aquifer_data(self.xcp, self.ycp)
        ieq = 0
        for e in self.model.elementlist:
            if e.nunknowns > 0:
                if e == self:
                    for w in self.wlist:
                        mat[0:1, ieq : ieq + w.nunknowns] += (
                            w.potinflayers(self.xcp, self.ycp, self.lcp)
                            / aq.Tcol[self.lcp]
                        )
                        ieq += w.nunknowns
                else:
                    mat[0:1, ieq : ieq + e.nunknowns] += (
                        e.potinflayers(self.xcp, self.ycp, self.lcp) / aq.Tcol[self.lcp]
                    )
                    ieq += e.nunknowns
            else:
                rhs[0] -= (
                    e.potentiallayers(self.xcp, self.ycp, self.lcp) / aq.T[self.lcp]
                )
        return mat, rhs

    def setparams(self, sol):
        self.parameters[:, 0] = sol
        # assign parameters to individual wells
        i = 0
        for w in self.wlist:
            w.parameters[:, 0] = sol[i : i + w.nparam]
            i += w.nparam

"""Model classes."""

import inspect  # Used for storing the input
import multiprocessing as mp
import warnings

import numpy as np
from scipy.integrate import quad_vec

from .aquifer import Aquifer, SimpleAquifer
from .aquifer_parameters import param_3d, param_maq
from .constant import ConstantStar
from .plots import PlotTim

__all__ = ["Model", "ModelMaq", "Model3D", "ModelXsection"]


class Model:
    """Create a model consisting of an arbitrary sequence of aquifers and leaky layers.

    Notes
    -----
    Use ModelMaq for regular sequence of aquifers and leaky layers.
    Use Model3D for multi-layer model of a single aquifer.

    Parameters
    ----------
    kaq : array
        hydraulic conductivity of each aquifer from the top down
    z : array
        elevation tops and bottoms of all layers
        layers may have zero thickness
    c : array
        resistance between two consecutive aquifer layers
        if ltype[0]='a': length is number of aquifers - 1
        if ltype[0]='l': length is number of aquifers
    npor : array
        porosity of all layers from the top down
    ltype : array of characters
        array indicating for each layer whether it is
        'a' aquifer layer
        'l' leaky layer
    """

    def __init__(self, kaq, c, z, npor, ltype):
        # All input variables are numpy arrays
        # That should be checked outside this function
        self.elementlist = []
        self.elementdict = {}  # only elements that have a label
        self.aq = Aquifer(self, kaq, c, z, npor, ltype)
        self.modelname = "ml"  # Used for writing out input

        self.plots = PlotTim(self)

    def initialize(self):
        # remove inhomogeneity elements (they are added again)
        self.elementlist = [e for e in self.elementlist if not e.inhomelement]
        self.aq.initialize()
        for e in self.elementlist:
            e.initialize()

    def add_element(self, e):
        self.elementlist.append(e)
        if e.label is not None:
            self.elementdict[e.label] = e

    def remove_element(self, e):
        """Remove element `e` from model."""
        if e.label is not None:
            self.elementdict.pop(e.label)
        self.elementlist.remove(e)

    def storeinput(self, frame):
        self.inputargs, _, _, self.inputvalues = inspect.getargvalues(frame)

    def potential(self, x, y, aq=None):
        if aq is None:
            aq = self.aq.find_aquifer_data(x, y)
        pot = np.zeros(aq.naq)
        for e in aq.elementlist:
            pot += e.potential(x, y, aq)
        rv = np.sum(pot * aq.eigvec, 1)
        if aq.ltype[0] == "l":
            # potential for head above leaky layer
            rv += aq.constantstar.potstar
        return rv

    def disvec(self, x, y, aq=None):
        """Discharge vector at `x`, `y`.

        Returns
        -------
        qxqy : array size (2, naq)
            first row is Qx in each aquifer layer, second row is Qy
        """
        if aq is None:
            aq = self.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, aq.naq))
        for e in aq.elementlist:
            rv += e.disvec(x, y, aq)
        rv = np.sum(rv[:, np.newaxis, :] * aq.eigvec, 2)
        return rv

    def normflux(self, x, y, theta):
        """Flux at point x, y in direction of angle theta.

        Parameters
        ----------
        x : float
        y : float
        theta : float
            angle for which to calculate flux,
            defined relative to positive x-axis.

        Returns
        -------
        flux
            flux in direction theta
        """
        qxqy = self.disvec(x, y)
        cosnorm = np.cos(theta)
        sinnorm = np.sin(theta)
        return cosnorm * qxqy[0] + sinnorm * qxqy[1]

    def _normflux_integrand(self, s, theta_norm, x1, y1):
        x = s * np.cos(theta_norm - np.pi / 2) + x1
        y = s * np.sin(theta_norm - np.pi / 2) + y1
        return self.normflux(x, y, theta_norm)

    def intnormflux_segment(self, x1, y1, x2, y2, method="legendre", ndeg=10):
        """Integrated normal (perpendicular) flux over specified line segment.

        Flux to the left is positive when going from (x1, y1) to (x2, y2).

        Parameters
        ----------
        x1 : float
        y1 : float
        x2 : float
        y2 : float
        method : str, optional
            integration method, either "quad" (numerical integration using scipy)
            or "legendre" (approximate integral using Gauss-Legendre quadrature),
            by default "legendre".
        ndeg : int, optional
            degree for legendre polynomial, by default 10,
            only used when method="legendre"

        Returns
        -------
        Qn : np.array
            integrated normal flux along specified line
        """
        z1 = x1 + y1 * 1j
        z2 = x2 + y2 * 1j
        normvec = (z2 - z1) / np.abs(z2 - z1) * np.exp(-np.pi * 1j / 2)
        theta_norm = np.angle(normvec) - np.pi
        L = np.abs(z2 - z1)
        if method == "quad":
            return quad_vec(
                self._normflux_integrand,
                0,
                L,
                args=(theta_norm, x1, y1),
            )[0]
        if method == "legendre":
            Xleg, wleg = np.polynomial.legendre.leggauss(ndeg)
            z = 0.5 * Xleg * (z2 - z1) + 0.5 * (z1 + z2)
            x = z.real
            y = z.imag
            qn = 0.0
            for i in range(ndeg):
                qn += wleg[i] * self.normflux(x=x[i], y=y[i], theta=theta_norm)
            return L * qn / 2.0

    def intnormflux(self, xy, method="legendre", ndeg=10):
        """Integrated normal (perpendicular) flux over polyline.

        Computes the flux per segment and per aquifer. Flux to the left is positive
        when going from (x1, y1) to (x2, y2).

        Parameters
        ----------
        xy : list [(x0, y0), (x1, y1),... , (xn, yn)] or 2D array
            if 2D-array, x in first column and y in second column
        method : str, optional
            integration method, either "quad" (numerical integration using scipy)
            or "legendre" (approximate integral using Gauss-Legendre quadrature),
            by default "legendre".
        ndeg : int, optional
            degree for legendre polynomial, by default 10,
            only used when method="legendre"

        Returns
        -------
        Qn : np.array of shape (naq, nsegments)
            integrated normal flux along specified polyline

        Examples
        --------
        Total flow across polyline can be obtained using:

        >>> np.sum(Qn)

        Total flow across segments summed over aquifers using

        >>> np.sum(Qn, axis=0)
        """
        xy = np.array(xy)  # convert to array
        if np.all(xy[-1] == xy[0]):
            Nsides = len(xy) - 1
        else:
            Nsides = len(xy)
        Qn = np.zeros((self.aq.naq, Nsides))
        for i in range(len(xy) - 1):
            x0, y0 = xy[i]
            x1, y1 = xy[i + 1]
            Qn[:, i] += self.intnormflux_segment(
                x0, y0, x1, y1, method=method, ndeg=ndeg
            )
        return Qn

    def qztop(self, x, y, aq=None):
        if aq is None:
            aq = self.aq.find_aquifer_data(x, y)
        rv = 0.0
        if aq.ltype[0] == "a":  # otherwise recharge cannot be added
            for e in aq.elementlist:
                rv += e.qztop(x, y, aq)
        return rv

    def head(self, x, y, layers=None, aq=None):
        """Head at `x`, `y`.

        Returns
        -------
        h : array length `naq` or `len(layers)`
            head in all `layers` (if not `None`),
            or all layers of aquifer (otherwise)
        """
        if aq is None:
            aq = self.aq.find_aquifer_data(x, y)
        rv = self.potential(x, y, aq) / aq.T
        if layers is None:
            return rv
        else:
            return rv[layers]

    def headgrid(self, xg, yg, layers=None, printrow=False):
        """Grid of heads.

        Parameters
        ----------
        xg : array
            x values of grid
        yg : array
            y values of grid
        layers : integer, list or array, optional
            layers for which grid is returned
        printrow : boolean, optional
            prints dot to screen for each row of grid if set to `True`

        Returns
        -------
        h : array size `nlayers, ny, nx`

        See Also
        --------
        :func:`~timml.model.Model.headgrid2`
        """
        nx, ny = len(xg), len(yg)
        if layers is None:
            Nlayers = self.aq.find_aquifer_data(xg[0], yg[0]).naq
        else:
            Nlayers = len(np.atleast_1d(layers))
        h = np.empty((Nlayers, ny, nx))
        for j in range(ny):
            if printrow:
                print(".", end="", flush=True)
            for i in range(nx):
                h[:, j, i] = self.head(xg[i], yg[j], layers)
        if printrow:
            print("", flush=True)
        return h

    def headgrid2(self, x1, x2, nx, y1, y2, ny, layers=None, printrow=False):
        """Grid of heads.

        Parameters
        ----------
        x1, x2, nx :
            x values are generated as linspace(x1, x2, nx)
        y1, y2, ny :
            y values are generated as linspace(y1, y2, ny)
        layers : integer, list or array, optional
            layers for which grid is returned
        printrow : boolean, optional
            prints dot to screen for each row of grid if set to `True`

        Returns
        -------
        h : array size `nlayers, ny, nx`

        See Also
        --------
        :func:`~timml.model.Model.headgrid`
        """
        xg, yg = np.linspace(x1, x2, nx), np.linspace(y1, y2, ny)
        return self.headgrid(xg, yg, layers=layers, printrow=printrow)

    def headalongline(self, x, y, layers=None):
        """Head along line or curve.

        Parameters
        ----------
        x : array
            x values of line
        y : array
            y values of line
        layers : integer, list or array, optional
            layers for which grid is returned

        Returns
        -------
        h : array size `nlayers, nx`
        """
        xg, yg = np.atleast_1d(x), np.atleast_1d(y)
        if layers is None:
            Nlayers = self.aq.find_aquifer_data(xg[0], yg[0]).naq
        else:
            Nlayers = len(np.atleast_1d(layers))
        nx = len(xg)
        if len(yg) == 1:
            yg = yg * np.ones(nx)
        h = np.zeros((Nlayers, nx))
        for i in range(nx):
            h[:, i] = self.head(xg[i], yg[i], layers)
        return h

    def disvecalongline(self, x, y, layers=None):
        """Compute discharge vector along line.

        Notes
        -----
        Assumes same number of layers for each x and y.
        Layers may be None or list of layers for which head is computed.

        Returns
        -------
        Qx
            [Nlayers,len(x)]
        Qy
            [Nlayers,len(x)]
        """
        xg, yg = np.atleast_1d(x), np.atleast_1d(y)
        if layers is None:
            nlayers = self.aq.find_aquifer_data(xg[0], yg[0]).naq
        else:
            nlayers = len(np.atleast_1d(layers))
        nx = len(xg)
        if len(yg) == 1:
            yg = yg * np.ones(nx)
        Qx = np.zeros((nlayers, nx))
        Qy = np.zeros((nlayers, nx))
        for i in range(nx):
            Qx[:, i], Qy[:, 1] = self.disvec(xg[i], yg[i], layers)
        return Qx, Qy

    #    def disvec_direction(self, s, x1, y1, cdirection):
    #        pass
    #
    #    def discharge_across_line(self, x1, y1, x2,  y2, layers=None):
    #        if layers is None:
    #            nlayers = self.aq.find_aquifer_data(x1, y1).naq
    #        else:
    #            nlayers = len(np.atleast_1d(layers))
    #        z1 = x1 + y1 * 1j
    #        z2 = x2 + y2 * 1j
    #        normvec = (z2 - z1) / np.abs(z2 - z1) * np.exp(-np.pi * 1j / 2)
    #        disvec = self.disvec(xg[i], yg[i], layers)

    def velocity(self, x, y, z):
        return self.velocomp(x, y, z)

    def velocomp(self, x, y, z, aq=None, layer_ltype=None):
        if aq is None:
            aq = self.aq.find_aquifer_data(x, y)
        if (z > aq.z[0]) or z < (aq.z[-1]):
            raise ValueError("z value not inside aquifer")
        if layer_ltype is None:
            layer, ltype, _ = aq.findlayer(z)
        else:
            layer, ltype = layer_ltype
        h = self.head(x, y, aq=aq)
        # qz between aquifer layers
        qzlayer = np.zeros(aq.naq + 1)
        qzlayer[1:-1] = (h[1:] - h[:-1]) / aq.c[1:]
        if aq.ltype[0] == "l":
            qzlayer[0] = (h[0] - aq.hstar) / aq.c[0]
        if ltype == "l":
            vz = qzlayer[layer] / aq.nporll[layer]
            vx = 0
            vy = 0
        else:
            qzbot = qzlayer[layer + 1]
            qztop = qzlayer[layer]
            if layer == 0:
                qztop += self.qztop(x, y)
            vz = (
                qzbot + (z - aq.zaqbot[layer]) / aq.Haq[layer] * (qztop - qzbot)
            ) / aq.nporaq[layer]
            qx, qy = self.disvec(x, y, aq=aq)
            vx = qx[layer] / (aq.Haq[layer] * aq.nporaq[layer])
            vy = qy[layer] / (aq.Haq[layer] * aq.nporaq[layer])
        return np.array([vx, vy, vz])

    def solve(self, printmat=0, sendback=0, silent=False):
        """Compute solution."""
        # Initialize elements
        self.initialize()
        # Compute number of equations
        self.neq = np.sum([e.nunknowns for e in self.elementlist])
        if self.neq == 0:
            return
        if silent is False:
            print(
                "Number of elements, Number of equations:",
                len(self.elementlist),
                ",",
                self.neq,
            )
        if self.neq == 0:
            if silent is False:
                print("No unknowns. Solution complete")
            return
        mat = np.empty((self.neq, self.neq))
        rhs = np.empty(self.neq)
        ieq = 0
        for e in self.elementlist:
            if e.nunknowns > 0:
                (
                    mat[ieq : ieq + e.nunknowns, :],
                    rhs[ieq : ieq + e.nunknowns],
                ) = e.equation()
                ieq += e.nunknowns
            if silent is False:
                print(".", end="", flush=True)
        if printmat:
            return mat, rhs
        sol = np.linalg.solve(mat, rhs)
        icount = 0
        for e in self.elementlist:
            if e.nunknowns > 0:
                e.setparams(sol[icount : icount + e.nunknowns])
                icount += e.nunknowns
        if silent is False:
            print()  # needed cause the dots are printed
            print("solution complete")
        elif (silent == "dot") or (silent == "."):
            print(".", end="", flush=True)
        if sendback:
            return sol
        return

    def solve_mp(self, nproc=4, printmat=0, sendback=0, silent=False):
        """Compute solution, multiprocessing implementation.

        Notes
        -----
        Estimated speedup approximately by factor of number of physical cores (virtual
        cores do not improve calculation time).
        """
        # Initialize elements
        self.initialize()
        # Compute number of equations
        self.neq = np.sum([e.nunknowns for e in self.elementlist])
        if self.neq == 0:
            return
        if silent is False:
            print(
                "Number of elements, Number of equations:",
                len(self.elementlist),
                ",",
                self.neq,
            )
        if self.neq == 0:
            if silent is False:
                print("No unknowns. Solution complete")
            return
        mat = np.empty((self.neq, self.neq))
        rhs = np.empty(self.neq)

        # start multiprocessing
        if nproc is None:
            nproc = (
                mp.cpu_count() - 1
            )  # make no. of processes equal to 1 less than no. of cores
        elif nproc > mp.cpu_count():
            print(
                "Given 'nproc' larger than no. of cores on machine. "
                f"Setting 'nproc' to {mp.cpu_count()}."
            )
            nproc = mp.cpu_count()

        pool = mp.Pool(processes=nproc)
        results = []
        for e in self.elementlist:
            if e.nunknowns > 0:
                results.append(pool.apply_async(e.equation))
            if silent is False:
                print(".", end="", flush=True)

        pool.close()
        pool.join()

        mat = np.empty((self.neq, self.neq))
        rhs = np.zeros(self.neq)

        ieq = 0

        for p in results:
            imat, irhs = p.get()
            mat[ieq : ieq + imat.shape[0], :] = imat
            rhs[ieq : ieq + irhs.shape[0]] = irhs
            ieq += imat.shape[0]

        # end multiprocessing

        if printmat:
            return mat, rhs
        sol = np.linalg.solve(mat, rhs)
        icount = 0
        for e in self.elementlist:
            if e.nunknowns > 0:
                e.setparams(sol[icount : icount + e.nunknowns])
                icount += e.nunknowns
        if silent is False:
            print()  # needed cause the dots are printed
            print("solution complete")
        elif (silent == "dot") or (silent == "."):
            print(".", end="", flush=True)
        if sendback:
            return sol
        return

    def write(self):
        rv = self.modelname + " = " + self.name + "(\n"
        for key in self.inputargs[1:]:  # The first argument (self) is ignored
            if isinstance(self.inputvalues[key], np.ndarray):
                rv += (
                    key
                    + " = "
                    + np.array2string(self.inputvalues[key], separator=",")
                    + ",\n"
                )
            elif isinstance(self.inputvalues[key], str):
                rv += key + " = '" + self.inputvalues[key] + "',\n"
            else:
                rv += key + " = " + str(self.inputvalues[key]) + ",\n"
        rv += ")\n"
        return rv

    def writemodel(self, fname):
        self.initialize()  # So that the model can be written without solving first
        f = open(fname, "w")
        f.write("from timml import *\n")
        f.write(self.write())
        for e in self.elementlist:
            f.write(e.write())
        f.close()

    def plot(self, *args, **kwargs):
        warnings.warn(
            "The 'ml.plot' method is deprecated. Use 'ml.plots.topview' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.plots.topview(**kwargs)

    def contour(self, *args, **kwargs):
        warnings.warn(
            "The 'ml.contour' method is deprecated. Use 'ml.plots.contour' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.plots.contour(*args, **kwargs)

    def vcontour(self, *args, **kwargs):
        warnings.warn(
            "The 'ml.vcontour' method is deprecated. Use 'ml.plots.vcontour' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.plots.vcontour(*args, **kwargs)

    def tracelines(self, *args, **kwargs):
        warnings.warn(
            "The 'ml.tracelines' method is deprecated. "
            "Use 'ml.plots.tracelines' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.plots.tracelines(*args, **kwargs)

    def vcontoursf1D(self, *args, **kwargs):
        warnings.warn(
            "The 'ml.vcontoursf1D' method is deprecated. "
            "Use 'ml.plots.vcontoursf1D' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.plots.vcontoursf1D(*args, **kwargs)


class ModelMaq(Model):
    """Create a model by specifying a mult-aquifer sequence of aquifer-leaky layer.

    Parameters
    ----------
    kaq : float, array or list
        Hydraulic conductivity of each aquifer from the top down.
        If float, hydraulic conductivity is the same in all aquifers.
    z : array or list
        Elevation of tops and bottoms of the aquifers from the top down.
        Leaky layers may have zero thickness.
           * if topboundary='conf': length is 2 * number of aquifers
           * if topboundary='semi': length is 2 * number of aquifers + 1
             as top of leaky layer on top of systems needs to be specified
    c : float, array or list
        Resistance of leaky layers from the top down.
           * if float, resistance is the same for all leaky layers
           * if topboundary='conf': length is number of aquifers - 1
           * if topboundary='semi': length is number of aquifers
    npor : float, array or list
        Porosity of all aquifers and leaky layers from the top down.
           * if float, porosity is the same for all layers
           * if topboundary='conf': length is 2 * number of aquifers - 1
           * if topboundary='semi': length is 2 * number of aquifers
    topboundary : string, 'conf' or 'semi' (default is 'conf')
        Indicates whether the topboundary is confined ('conf') or
        semi-confined ('semi').
    hstar : float or None (default is None)
        Head value above semi-confining top, only read if topboundary='semi'.

    Examples
    --------
    >>> ml = ModelMaq(kaq=[10, 20], z=[20, 12, 10, 0], c=1000)
    """

    def __init__(self, kaq=1, z=None, c=None, npor=0.3, topboundary="conf", hstar=None):
        if c is None:
            c = []
        if z is None:
            z = [1, 0]
        self.storeinput(inspect.currentframe())
        kaq, c, npor, ltype = param_maq(kaq, z, c, npor, topboundary)
        super().__init__(kaq, c, z, npor, ltype)
        self.name = "ModelMaq"
        if self.aq.ltype[0] == "l":
            ConstantStar(self, hstar, aq=self.aq)


class Model3D(Model):
    """Create a multi-layer model object consisting of stacked aquifer layers.

    The resistance between the layers is computed from the vertical hydraulic
    conductivity of the layers.

    Parameters
    ----------
    kaq : float, array or list
        hydraulic conductivity of each layer from the top down
        if float, hydraulic conductivity is the same in all aquifers
    z : array or list
        elevation of top of system followed by bottoms of all layers
        from the top down
        bottom of layer is automatically equal to top of layer below it
        length is number of aquifer layers + 1
    kzoverkh : float
        vertical anisotropy ratio vertical k divided by horizontal k
        if float, value is the same for all layers
        length is number of layers
    npor : float, array or list
        porosity of all aquifer layers
        from the top down
        if float, porosity is the same for all layers
        if topboundary='conf': length is number of layers
        if topboundary='semi': length is number of layers + 1
    topboundary : string, 'conf' or 'semi' (default is 'conf')
        indicating whether the top is confined ('conf') or
        semi-confined ('semi')
    topres : float
        resistance of top semi-confining layer (read if topboundary='semi')
    topthick: float
        thickness of top semi-confining layer (read if topboundary='semi')
    hstar : float or None (default is None)
        head value above semi-confining top (read if topboundary='semi')

    Notes
    -----
    For semi-confined aquifers, set top equal to 'semi' and provide:
       - topres: resistance of top
       - tophick: thickness of top
       - hstar: head above top

    Examples
    --------
    >>> ml = Model3D(kaq=10, z=np.arange(20, -1, -2), kzoverkh=0.1)
    """

    def __init__(
        self,
        kaq=1,
        z=None,
        kzoverkh=1,
        npor=0.3,
        topboundary="conf",
        topres=0,
        topthick=0,
        hstar=0,
    ):
        """Model3D.

        for semi-confined aquifers, set top equal to 'semi' and provide
        topres: resistance of top
        tophick: thickness of top
        hstar: head above top.
        """
        if z is None:
            z = [1, 0]
        self.storeinput(inspect.currentframe())
        kaq, c, npor, ltype = param_3d(kaq, z, kzoverkh, npor, topboundary, topres)
        if topboundary == "semi":
            z = np.hstack((z[0] + topthick, z))
        super().__init__(kaq, c, z, npor, ltype)
        self.name = "Model3D"
        if self.aq.ltype[0] == "l":
            ConstantStar(self, hstar, aq=self.aq)


class ModelXsection(Model):
    def __init__(self, naq=1):
        self.elementlist = []
        self.elementdict = {}  # only elements that have a label
        self.aq = SimpleAquifer(naq)
        self.modelname = "ml"  # Used for writing out input

        self.plots = PlotTim(self)
        self.name = "ModelXsection"

    def check_inhoms(self):
        """Check if number of aquifers in inhoms matches number of aquifers in model."""
        naqs = {}
        for inhom in self.aq.inhomlist:
            naqs[inhom.name] = inhom.naq
        check = np.array(list(naqs.values())) == self.aq.naq
        if not check.all():
            raise ValueError(
                f"Number of aquifers does not match {self.aq.naq}:\n{naqs}"
            )
        # # shared boundary check
        # # NOTE: does not deal with nested inhoms
        # xcoords = np.concatenate(
        #     [(inhom.x1, inhom.x2) for inhom in self.aq.inhomdict.values()]
        # )
        # xcoords.sort()
        # if not np.all(np.diff(xcoords[1:-1])[::2] < 1e-10):
        #     raise ValueError("Not all inhomogeneities have shared boundaries.")

    def initialize(self):
        self.check_inhoms()
        super().initialize()

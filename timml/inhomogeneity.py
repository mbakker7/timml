import inspect  # Used for storing the input
from copy import deepcopy
from warnings import warn

import numpy as np

from timml.aquifer import AquiferData
from timml.aquifer_parameters import param_3d, param_maq
from timml.constant import ConstantInside, ConstantStar
from timml.element import Element
from timml.intlinesink import (
    IntFluxDiffLineSink,
    IntFluxLineSink,
    IntHeadDiffLineSink,
    LeakyIntHeadDiffLineSink,
)
from timml.util import compute_z1z2, refine_n_segments


class PolygonInhom(AquiferData):
    tiny = 1e-8

    def __init__(
        self,
        model,
        xy,
        kaq,
        c,
        z,
        npor,
        ltype,
        hstar,
        N,
        order,
        ndeg,
        refine_level=1,
        addtomodel=True,
    ):
        self._input = {k: v for k, v in locals().items() if k not in ["self", "model"]}
        # All input variables except model should be numpy arrays
        # That should be checked outside this function):
        AquiferData.__init__(self, model, kaq, c, z, npor, ltype)
        self.order = order
        self.ndeg = ndeg
        self.hstar = hstar
        self.N = N
        self.addtomodel = addtomodel
        if self.addtomodel:
            self.inhom_number = self.model.aq.add_inhom(self)
        self.xy = xy
        self.z1, self.z2 = compute_z1z2(self.xy)
        self.Nsides = len(self.z1)
        Zin = 1e-6j
        Zout = -1e-6j
        self.zcin = 0.5 * (self.z2 - self.z1) * Zin + 0.5 * (
            self.z1 + self.z2
        )  # point at center on inside
        self.zcout = 0.5 * (self.z2 - self.z1) * Zout + 0.5 * (
            self.z1 + self.z2
        )  # point at center on outside
        self.zcenter = np.mean(self.z1)
        self.xcenter, self.ycenter = self.zcenter.real, self.zcenter.imag
        self.x = np.hstack((self.z1.real, self.z2[-1].real))
        self.y = np.hstack((self.z1.imag, self.z2[-1].imag))
        self.xmin = min(self.x)
        self.xmax = max(self.x)
        self.ymin = min(self.y)
        self.ymax = max(self.y)
        self.refine_level = refine_level

    def __repr__(self):
        return "PolygonInhom: " + str(list(zip(self.x, self.y)))

    def isinside(self, x, y):
        rv = 0
        if (
            (x >= self.xmin)
            and (x <= self.xmax)
            and (y >= self.ymin)
            and (y <= self.ymax)
        ):
            z = complex(x, y)
            bigZ = (2.0 * z - (self.z1 + self.z2)) / (self.z2 - self.z1)
            bigZmin1 = bigZ - 1.0
            bigZplus1 = bigZ + 1.0
            minAbsBigZmin1 = min(abs(bigZmin1))
            minAbsBigZplus1 = min(abs(bigZplus1))
            if minAbsBigZmin1 < self.tiny or minAbsBigZplus1 < self.tiny:
                rv = 1
                return rv
            angles = np.log(bigZmin1 / bigZplus1).imag
            angle = np.sum(angles)
            if angle > np.pi:
                rv = 1
        return rv

    def create_elements(self):
        aqin = self.model.aq.find_aquifer_data(self.zcin[0].real, self.zcin[0].imag)
        inhom_elements = []
        for i in range(self.Nsides):
            aqout = self.model.aq.find_aquifer_data(
                self.zcout[i].real, self.zcout[i].imag
            )
            if (aqout == self.model.aq) or (aqout.inhom_number > self.inhom_number):
                ls_in = IntHeadDiffLineSink(
                    self.model,
                    x1=self.x[i],
                    y1=self.y[i],
                    x2=self.x[i + 1],
                    y2=self.y[i + 1],
                    order=self.order,
                    ndeg=self.ndeg,
                    label=None,
                    addtomodel=True,
                    aq=aqin,
                    aqin=aqin,
                    aqout=aqout,
                )
                ls_out = IntFluxDiffLineSink(
                    self.model,
                    x1=self.x[i],
                    y1=self.y[i],
                    x2=self.x[i + 1],
                    y2=self.y[i + 1],
                    order=self.order,
                    ndeg=self.ndeg,
                    label=None,
                    addtomodel=True,
                    aq=aqout,
                    aqin=aqin,
                    aqout=aqout,
                )
                inhom_elements += [ls_in, ls_out]
        if aqin.ltype[0] == "a":  # add constant on inside
            c = ConstantInside(self.model, self.zcin.real, self.zcin.imag)
            c.inhomelement = True
            if self.N is not None:
                a = AreaSinkInhom(self.model, self.N, self.xcenter, aq=aqin)
                a.inhomelement = True
        if aqin.ltype[0] == "l":
            assert self.hstar is not None, "Error: hstar needs to be set"
            c = ConstantStar(self.model, self.hstar, aq=aqin)
            c.inhomelement = True
        inhom_elements += [c]
        return inhom_elements

    def _refine(self, n=None):
        if n is None:
            n = self.refine_level
        xyr, _ = refine_n_segments(self.xy, "polygon", n_segments=n)
        input_args = deepcopy(self._input)
        cls = input_args.pop("__class__", self.__class__)
        input_args["model"] = self.model
        # overwrite some input args for refined element
        input_args["xy"] = xyr
        input_args["refine_level"] = 1  # set to 1 to prevent further refinement
        input_args["addtomodel"] = False
        return cls(**input_args)


class PolygonInhomMaq(PolygonInhom):
    """Create a polygonal inhomogeneity.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    xy : array or list
        list or array of (x,y) pairs of coordinates of corners of the
        inhomogeneity. polygonal boundary is automatically closed (so first
        point is not repeated)
    kaq : float, array or list
        hydraulic conductivity of each aquifer from the top down
        if float, hydraulic conductivity is the same in all aquifers
    z : array or list
        elevation tops and bottoms of the aquifers from the top down
        leaky layers may have zero thickness
        if topboundary='conf': length is 2 * number of aquifers
        if topboundary='semi': length is 2 * number of aquifers + 1 as top
        of leaky layer on top of systems needs to be specified
    c : float, array or list
        resistance of leaky layers from the top down
        if float, resistance is the same for all leaky layers
        if topboundary='conf': length is number of aquifers - 1
        if topboundary='semi': length is number of aquifers
    npor : float, array or list
        porosity of all aquifers and leaky layers from the top down
        if float, porosity is the same for all layers
        if topboundary='conf': length is 2 * number of aquifers - 1
        if topboundary='semi': length is 2 * number of aquifers
    topboundary : string, 'conf' or 'semi' (default is 'conf')
        indicating whether the top is confined ('conf') or
        semi-confined ('semi')
    hstar : float or None (default is None)
        head value above semi-confining top, only read if topboundary='semi'
    N : float or None (default is None)
        infiltration rate (L/T) inside inhomogeneity. Only possible if
        topboundary='conf'
    order : int
        polynomial order of flux along each segment
    ndeg : int
        number of points used between two segments to numerically
        integrate normal discharge
    """

    tiny = 1e-8

    def __init__(
        self,
        model,
        xy,
        kaq=1,
        z=[1, 0],
        c=[],
        npor=0.3,
        topboundary="conf",
        hstar=None,
        N=None,
        order=3,
        ndeg=3,
        refine_level=1,
        addtomodel=True,
    ):
        _input = {k: v for k, v in locals().items() if k not in ["self", "model"]}
        if N is not None:
            assert (
                topboundary[:4] == "conf"
            ), "Error: infiltration can only be added if topboundary='conf'"
        self.storeinput(inspect.currentframe())
        (
            kaq,
            c,
            npor,
            ltype,
        ) = param_maq(kaq, z, c, npor, topboundary)
        PolygonInhom.__init__(
            self,
            model,
            xy,
            kaq,
            c,
            z,
            npor,
            ltype,
            hstar,
            N,
            order,
            ndeg,
            refine_level=refine_level,
            addtomodel=addtomodel,
        )
        self._input = _input


class PolygonInhom3D(PolygonInhom):
    """Model3D Class to create a multi-layer model object consisting of many aquifer
    layers. The resistance between the layers is computed from the vertical hydraulic
    conductivity of the layers.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    xy : array or list
        list or array of (x,y) pairs of coordinates of corners of the
        inhomogeneity. polygonal boundary is automatically closed (so first
        point is not repeated)
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
    N : float or None (default is None)
        infiltration rate (L/T) inside inhomogeneity. Only possible if
        topboundary='conf'
    order : int
        polynomial order of flux along each segment
    ndeg : int
        number of points used between two segments to numerically
        integrate normal discharge
    """

    def __init__(
        self,
        model,
        xy,
        kaq=1,
        z=[1, 0],
        kzoverkh=1,
        npor=0.3,
        topboundary="conf",
        topres=0,
        topthick=0,
        hstar=0,
        N=None,
        order=3,
        ndeg=3,
        refine_level=1,
        addtomodel=True,
    ):
        _input = {k: v for k, v in locals().items() if k not in ["self", "model"]}
        if N is not None:
            assert (
                topboundary[:4] == "conf"
            ), "Error: infiltration can only be added if topboundary='conf'"
        self.storeinput(inspect.currentframe())
        kaq, c, npor, ltype = param_3d(kaq, z, kzoverkh, npor, topboundary, topres)
        if topboundary == "semi":
            z = np.hstack((z[0] + topthick, z))
        PolygonInhom.__init__(
            self,
            model,
            xy,
            kaq,
            c,
            z,
            npor,
            ltype,
            hstar,
            N,
            order,
            ndeg,
            refine_level=refine_level,
            addtomodel=addtomodel,
        )
        self._input = _input


class BuildingPit(AquiferData):
    tiny = 1e-8

    def __init__(
        self,
        model,
        xy,
        kaq,
        c,
        z,
        ltype,
        hstar=None,
        npor=0.3,
        order=3,
        ndeg=3,
        layers=[0],
        refine_level=1,
        addtomodel=True,
    ):
        """Element to simulate a building pit with an impermeable wall.

        Parameters
        ----------
        model : Model object
            model to which the element is added
        xy : array or list
            list or array of (x,y) pairs of coordinates of corners of the
            inhomogeneity
            polygonal boundary is automatically closed (so first point
            is not repeated)
        kaq : float, array or list
            hydraulic conductivity of each aquifer from the top down
            if float, hydraulic conductivity is the same in all aquifers
        c : float, array or list
            resistance of leaky layers from the top down
            if float, resistance is the same for all leaky layers
            if topboundary='conf': length is number of aquifers - 1
            if topboundary='semi': length is number of aquifers
        z : array or list
            elevation tops and bottoms of the aquifers from the top down
            leaky layers may have zero thickness
            if topboundary='conf': length is 2 * number of aquifers
            if topboundary='semi': length is 2 * number of aquifers + 1 as top
            of leaky layer on top of systems needs to be specified
        ltype : list of string
            indicating whether layer is an aquifer ('a') or a leaky layer ('l').
        hstar : float or None (default is None)
            head value above semi-confining top, only read if topboundary='semi'
        npor : float, array or list
            porosity of all aquifers and leaky layers from the top down
            if float, porosity is the same for all layers
            if topboundary='conf': length is 2 * number of aquifers - 1
            if topboundary='semi': length is 2 * number of aquifers
        order : int
            polynomial order of flux along each segment
        ndeg : int
            number of points used between two segments to numerically
            integrate normal discharge
        layers : list or np.array
            layers in which impermeable wall is present.
        refine_level : int, optional
            refine element by partitioning each side into refine_level segments, default
            is 1, which means no refinement is applied.
        """
        self._input = {k: v for k, v in locals().items() if k not in ["self", "model"]}
        AquiferData.__init__(self, model, kaq, c, z, npor, ltype)
        self.order = order
        self.ndeg = ndeg
        self.xy = xy
        self.layers = np.atleast_1d(layers)  # layers with impermeable wall
        self.nonimplayers = list(
            set(range(self.model.aq.naq)) - set(self.layers)
        )  # layers without wall
        self.hstar = hstar
        self.refine_level = refine_level
        self.addtomodel = addtomodel
        if self.addtomodel:
            self.inhom_number = self.model.aq.add_inhom(self)

        # compute derived params
        self.compute_derived_params()

    def compute_derived_params(self):
        self.z1, self.z2 = compute_z1z2(self.xy)
        self.Nsides = len(self.z1)
        Zin = 1e-6j
        Zout = -1e-6j
        self.zcin = 0.5 * (self.z2 - self.z1) * Zin + 0.5 * (
            self.z1 + self.z2
        )  # point at center on inside
        self.zcout = 0.5 * (self.z2 - self.z1) * Zout + 0.5 * (
            self.z1 + self.z2
        )  # point at center on outside
        self.x = np.hstack((self.z1.real, self.z2[-1].real))
        self.y = np.hstack((self.z1.imag, self.z2[-1].imag))
        self.xmin = min(self.x)
        self.xmax = max(self.x)
        self.ymin = min(self.y)
        self.ymax = max(self.y)

    def __repr__(self):
        return (
            f"{self.__class__.__name__}: layers "
            + str(list(self.layers))
            + ", "
            + str(list(zip(self.x, self.y)))
        )

    def isinside(self, x, y):
        rv = 0
        if (
            (x >= self.xmin)
            and (x <= self.xmax)
            and (y >= self.ymin)
            and (y <= self.ymax)
        ):
            z = complex(x, y)
            bigZ = (2.0 * z - (self.z1 + self.z2)) / (self.z2 - self.z1)
            bigZmin1 = bigZ - 1.0
            bigZplus1 = bigZ + 1.0
            minAbsBigZmin1 = min(abs(bigZmin1))
            minAbsBigZplus1 = min(abs(bigZplus1))
            if minAbsBigZmin1 < self.tiny or minAbsBigZplus1 < self.tiny:
                rv = 1
                return rv
            angles = np.log(bigZmin1 / bigZplus1).imag
            angle = np.sum(angles)
            if angle > np.pi:
                rv = 1
        return rv

    def create_elements(self):
        aqin = self.model.aq.find_aquifer_data(self.zcin[0].real, self.zcin[0].imag)
        inhom_elements = []
        for i in range(self.Nsides):
            aqout = self.model.aq.find_aquifer_data(
                self.zcout[i].real, self.zcout[i].imag
            )
            if (aqout == self.model.aq) or (aqout.inhom_number > self.inhom_number):
                # Conditions for layers with impermeable walls
                ils_in = IntFluxLineSink(
                    self.model,
                    x1=self.x[i],
                    y1=self.y[i],
                    x2=self.x[i + 1],
                    y2=self.y[i + 1],
                    layers=self.layers,
                    order=self.order,
                    ndeg=self.ndeg,
                    label=None,
                    addtomodel=True,
                    aq=aqin,
                    aqin=aqin,
                    aqout=aqout,
                )
                ils_out = IntFluxLineSink(
                    self.model,
                    x1=self.x[i],
                    y1=self.y[i],
                    x2=self.x[i + 1],
                    y2=self.y[i + 1],
                    layers=self.layers,
                    order=self.order,
                    ndeg=self.ndeg,
                    label=None,
                    addtomodel=True,
                    aq=aqout,
                    aqin=aqin,
                    aqout=aqout,
                )
                inhom_elements += [ils_in, ils_out]

                if len(self.nonimplayers) > 0:
                    # use these conditions for layers without impermeable or leaky walls
                    ihdls_in = IntHeadDiffLineSink(
                        self.model,
                        x1=self.x[i],
                        y1=self.y[i],
                        x2=self.x[i + 1],
                        y2=self.y[i + 1],
                        layers=self.nonimplayers,
                        order=self.order,
                        ndeg=self.ndeg,
                        label=None,
                        addtomodel=True,
                        aq=aqin,
                        aqin=aqin,
                        aqout=aqout,
                    )
                    ihdls_out = IntFluxDiffLineSink(
                        self.model,
                        x1=self.x[i],
                        y1=self.y[i],
                        x2=self.x[i + 1],
                        y2=self.y[i + 1],
                        layers=self.nonimplayers,
                        order=self.order,
                        ndeg=self.ndeg,
                        label=None,
                        addtomodel=True,
                        aq=aqout,
                        aqin=aqin,
                        aqout=aqout,
                    )
                    inhom_elements += [ihdls_in, ihdls_out]

        if aqin.ltype[0] == "a":  # add constant on inside
            c = ConstantInside(self.model, self.zcin.real, self.zcin.imag)
            c.inhomelement = True
        if aqin.ltype[0] == "l":
            assert self.hstar is not None, "Error: hstar needs to be set"
            c = ConstantStar(self.model, self.hstar, aq=aqin)
            c.inhomelement = True
        inhom_elements += [c]

        return inhom_elements

    def _refine(self, n=None):
        if n is None:
            n = self.refine_level
        xyr, _ = refine_n_segments(self.xy, "polygon", n_segments=n)
        input_args = deepcopy(self._input)
        cls = input_args.pop("__class__")
        input_args["model"] = self.model
        # overwrite some input args for refined element
        input_args["xy"] = xyr
        input_args["refine_level"] = 1  # set to 1 to prevent further refinement
        input_args["addtomodel"] = False
        return cls(**input_args)


class BuildingPitMaq(BuildingPit):
    def __init__(
        self,
        model,
        xy,
        kaq=1.0,
        c=[],
        z=[1, 0],
        topboundary="conf",
        hstar=None,
        npor=0.3,
        order=3,
        ndeg=3,
        layers=[0],
        refine_level=1,
        addtomodel=True,
    ):
        """Element to simulate a building pit with an impermeable wall in ModelMaq.

        Parameters
        ----------
        model : Model object
            model to which the element is added
        xy : array or list
            list or array of (x,y) pairs of coordinates of corners of the
            inhomogeneity
            polygonal boundary is automatically closed (so first point
            is not repeated)
        kaq : float, array or list
            hydraulic conductivity of each aquifer from the top down
            if float, hydraulic conductivity is the same in all aquifers
        c : float, array or list
            resistance of leaky layers from the top down
            if float, resistance is the same for all leaky layers
            if topboundary='conf': length is number of aquifers - 1
            if topboundary='semi': length is number of aquifers
        z : array or list
            elevation tops and bottoms of the aquifers from the top down
            leaky layers may have zero thickness
            if topboundary='conf': length is 2 * number of aquifers
            if topboundary='semi': length is 2 * number of aquifers + 1 as top
            of leaky layer on top of systems needs to be specified
        topboundary : string, 'conf' or 'semi' (default is 'conf')
            indicating whether the top is confined ('conf') or
            semi-confined ('semi'). For a building pit, the 'conf'
            option is generally more applicable.
        hstar : float or None (default is None)
            head value above semi-confining top, only read if topboundary='semi'
        npor : float, array or list
            porosity of all aquifers and leaky layers from the top down
            if float, porosity is the same for all layers
            if topboundary='conf': length is 2 * number of aquifers - 1
            if topboundary='semi': length is 2 * number of aquifers
        order : int
            polynomial order of flux along each segment
        ndeg : int
            number of points used between two segments to numerically
            integrate normal discharge
        layers: list or np.array
            layers in which impermeable wall is present.
        refine_level : int, optional
            refine element by splitting up each side into refine_level segments, default
            is 1, which means no refinement is applied.
        """
        _input = {k: v for k, v in locals().items() if k not in ["self", "model"]}
        (kaq, c, npor, ltype) = param_maq(kaq, z, c, npor, topboundary)
        super().__init__(
            model=model,
            xy=xy,
            kaq=kaq,
            c=c,
            z=z,
            ltype=ltype,
            hstar=hstar,
            npor=npor,
            order=order,
            ndeg=ndeg,
            layers=layers,
            refine_level=refine_level,
            addtomodel=addtomodel,
        )
        self._input = _input


class BuildingPit3D(BuildingPit):
    def __init__(
        self,
        model,
        xy,
        kaq=1.0,
        kzoverkh=1.0,
        z=[1, 0],
        topboundary="conf",
        topres=0,
        topthick=0,
        hstar=None,
        npor=0.3,
        order=3,
        ndeg=3,
        layers=[0],
        refine_level=1,
        addtomodel=True,
    ):
        """Element to simulate a building pit with an impermeable wall in Model3D.

        Parameters
        ----------
        model : Model object
            model to which the element is added
        xy : array or list
            list or array of (x,y) pairs of coordinates of corners of the
            inhomogeneity
            polygonal boundary is automatically closed (so first point
            is not repeated)
        kaq : float, array or list
            hydraulic conductivity of each aquifer from the top down
            if float, hydraulic conductivity is the same in all aquifers
        kzoverkh : float, array or list
            vertical anisotropy ratio vertical k divided by horizontal k
            if float, value is the same for all layers, length is number of layers
        z : array or list
            elevation tops and bottoms of the aquifers from the top down
            leaky layers may have zero thickness
            if topboundary='conf': length is 2 * number of aquifers
            if topboundary='semi': length is 2 * number of aquifers + 1 as top
            of leaky layer on top of systems needs to be specified
        topboundary : string, 'conf' or 'semi' (default is 'conf')
            indicating whether the top is confined ('conf') or
            semi-confined ('semi'). For a building pit, the 'conf'
            option is generally more applicable.
        topres : float
            resistance of top leaky layer (read if topboundary="semi")
        topthick : float
            thickness of top semi-confining layer (read if topboundary='semi')
        hstar : float or None (default is None)
            head value above semi-confining top, only read if topboundary='semi'
        npor : float, array or list
            porosity of all aquifers and leaky layers from the top down
            if float, porosity is the same for all layers
            if topboundary='conf': length is 2 * number of aquifers - 1
            if topboundary='semi': length is 2 * number of aquifers
        order : int
            polynomial order of flux along each segment
        ndeg : int
            number of points used between two segments to numerically
            integrate normal discharge
        layers: list or np.array
            layers in which impermeable wall is present.
        refine_level : int, optional
            refine element by partitioning each side into refine_level segments, default
            is 1, which means no refinement is applied.
        """
        _input = {k: v for k, v in locals().items() if k not in ["self", "model"]}
        (kaq, c, npor, ltype) = param_3d(kaq, z, kzoverkh, npor, topboundary, topres)
        if topboundary == "semi":
            z = np.hstack((z[0] + topthick, z))
        super().__init__(
            model=model,
            xy=xy,
            kaq=kaq,
            c=c,
            z=z,
            ltype=ltype,
            hstar=hstar,
            npor=npor,
            order=order,
            ndeg=ndeg,
            layers=layers,
            refine_level=refine_level,
            addtomodel=addtomodel,
        )
        self._input = _input


class LeakyBuildingPit(BuildingPit):
    def __init__(
        self,
        model,
        xy,
        kaq,
        z,
        c,
        npor,
        ltype,
        hstar=None,
        order=3,
        ndeg=3,
        layers=[0],
        res=np.inf,
        refine_level=1,
        addtomodel=True,
    ):
        """Element to simulate a building pit with a leaky wall.

        Parameters
        ----------
        model : Model object
            model to which the element is added
        xy : array or list
            list or array of (x,y) pairs of coordinates of corners of the
            inhomogeneity
            polygonal boundary is automatically closed (so first point
            is not repeated)
        kaq : float, array or list
            hydraulic conductivity of each aquifer from the top down
            if float, hydraulic conductivity is the same in all aquifers
        c : float, array or list
            resistance of leaky layers from the top down
            if float, resistance is the same for all leaky layers
            if topboundary='conf': length is number of aquifers - 1
            if topboundary='semi': length is number of aquifers
        z : array or list
            elevation tops and bottoms of the aquifers from the top down
            leaky layers may have zero thickness
            if topboundary='conf': length is 2 * number of aquifers
            if topboundary='semi': length is 2 * number of aquifers + 1 as top
            of leaky layer on top of systems needs to be specified
        npor : float, array or list
            porosity of all aquifers and leaky layers from the top down
            if float, porosity is the same for all layers
            if topboundary='conf': length is 2 * number of aquifers - 1
            if topboundary='semi': length is 2 * number of aquifers
        ltype : list of string
            indicating whether layer is an aquifer ('a') or a leaky layer ('l').
        hstar : float or None (default is None)
            head value above semi-confining top, only read if topboundary='semi'
        order : int
            polynomial order of flux along each segment
        ndeg : int
            number of points used between two segments to numerically
            integrate normal discharge
        layers: list or np.array
            layers in which leaky wall is present.
        res : float or np.array, optional
            resistance of leaky wall, if passed as an array must be either
            shape (n_segments,) or (n_layers, n_segments). Default is np.inf,
            which simulates an impermeable wall.
        refine_level : int, optional
            refine element by partitioning each side into refine_level segments, default
            is 1, which means no refinement is applied.
        """
        super().__init__(
            model,
            xy,
            kaq=kaq,
            c=c,
            z=z,
            ltype=ltype,
            hstar=hstar,
            npor=npor,
            order=order,
            ndeg=ndeg,
            layers=layers,
            refine_level=refine_level,
            addtomodel=addtomodel,
        )
        self._input = {k: v for k, v in locals().items() if k not in ["self", "model"]}
        if isinstance(res, (int, float, np.integer)):
            # make 2D so indexing resistance works for all cases
            self.res = res * np.ones((1, self.Nsides))
        elif len(res) == self.Nsides and np.ndim(res) == 1:
            # make 2D so indexing resistance works for all cases
            self.res = np.atleast_2d(res)
        elif len(res) == len(self.layers):
            self.res = res
        else:
            raise ValueError(
                "Resistance `res` must be float or array of size "
                f"(n_segments={self.Nsides},) "
                f"or (n_layers={len(self.layers)}, n_segments={self.Nsides})."
            )
        if np.any(self.res < self.tiny):
            warn(
                f"Found resistances smaller than {self.tiny}, "
                f"these were replaced by {self.tiny}."
            )
            self.res[self.res < self.tiny] = self.tiny

    def create_elements(self):
        aqin = self.model.aq.find_aquifer_data(self.zcin[0].real, self.zcin[0].imag)
        inhom_elements = []
        for i in range(self.Nsides):
            aqout = self.model.aq.find_aquifer_data(
                self.zcout[i].real, self.zcout[i].imag
            )
            if (aqout == self.model.aq) or (aqout.inhom_number > self.inhom_number):
                # Conditions for layers with leaky walls
                ilhdls_in = LeakyIntHeadDiffLineSink(
                    self.model,
                    x1=self.x[i],
                    y1=self.y[i],
                    x2=self.x[i + 1],
                    y2=self.y[i + 1],
                    res=self.res[:, i],
                    layers=self.layers,
                    order=self.order,
                    ndeg=self.ndeg,
                    label=None,
                    addtomodel=False,
                    aq=aqin,
                    aqin=aqin,
                    aqout=aqout,
                )

                ilhdls_out = LeakyIntHeadDiffLineSink(
                    self.model,
                    x1=self.x[i],
                    y1=self.y[i],
                    x2=self.x[i + 1],
                    y2=self.y[i + 1],
                    res=self.res[:, i],
                    layers=self.layers,
                    order=self.order,
                    ndeg=self.ndeg,
                    label=None,
                    addtomodel=False,
                    aq=aqout,
                    aqin=aqin,
                    aqout=aqout,
                )
                inhom_elements += [ilhdls_in, ilhdls_out]

                if len(self.nonimplayers) > 0:
                    # use these conditions for layers without leaky walls
                    ihdls_in = IntHeadDiffLineSink(
                        self.model,
                        x1=self.x[i],
                        y1=self.y[i],
                        x2=self.x[i + 1],
                        y2=self.y[i + 1],
                        layers=self.nonimplayers,
                        order=self.order,
                        ndeg=self.ndeg,
                        label=None,
                        addtomodel=False,
                        aq=aqin,
                        aqin=aqin,
                        aqout=aqout,
                    )
                    ihdls_out = IntFluxDiffLineSink(
                        self.model,
                        x1=self.x[i],
                        y1=self.y[i],
                        x2=self.x[i + 1],
                        y2=self.y[i + 1],
                        layers=self.nonimplayers,
                        order=self.order,
                        ndeg=self.ndeg,
                        label=None,
                        addtomodel=False,
                        aq=aqout,
                        aqin=aqin,
                        aqout=aqout,
                    )
                    inhom_elements += [ihdls_in, ihdls_out]

        if aqin.ltype[0] == "a":  # add constant on inside
            c = ConstantInside(
                self.model, self.zcin.real, self.zcin.imag, addtomodel=False
            )
            c.inhomelement = True
        if aqin.ltype[0] == "l":
            assert self.hstar is not None, "Error: hstar needs to be set"
            c = ConstantStar(self.model, self.hstar, aq=aqin, addtomodel=False)
            c.inhomelement = True
        inhom_elements += [c]

        return inhom_elements

    def _refine(self, n=None):
        if n is None:
            n = self.refine_level
        xyr, reindexer = refine_n_segments(self.xy, "polygon", n_segments=n)
        input_args = deepcopy(self._input)
        cls = input_args.pop("__class__")
        input_args["model"] = self.model
        # overwrite some input args for refined element
        input_args["xy"] = xyr
        input_args["res"] = self.res[:, reindexer]
        input_args["refine_level"] = 1  # set to 1 to prevent further refinement
        input_args["addtomodel"] = False
        return cls(**input_args)


class LeakyBuildingPitMaq(LeakyBuildingPit):
    """Element to simulate a building pit with a leaky wall in ModelMaq.

    Parameters
    ----------
    model : Model object
        model to which the element is added
    xy : array or list
        list or array of (x,y) pairs of coordinates of corners of the
        inhomogeneity
        polygonal boundary is automatically closed (so first point
        is not repeated)
    kaq : float, array or list
        hydraulic conductivity of each aquifer from the top down
        if float, hydraulic conductivity is the same in all aquifers
    c : float, array or list
        resistance of leaky layers from the top down
        if float, resistance is the same for all leaky layers
        if topboundary='conf': length is number of aquifers - 1
        if topboundary='semi': length is number of aquifers
    z : array or list
        elevation tops and bottoms of the aquifers from the top down
        leaky layers may have zero thickness
        if topboundary='conf': length is 2 * number of aquifers
        if topboundary='semi': length is 2 * number of aquifers + 1 as top
        of leaky layer on top of systems needs to be specified
    topboundary : string, 'conf' or 'semi' (default is 'conf')
        indicating whether the top is confined ('conf') or
        semi-confined ('semi'). For a building pit, the 'conf'
        option is generally more applicable.
    hstar : float or None (default is None)
        head value above semi-confining top, only read if topboundary='semi'
    npor : float, array or list
        porosity of all aquifers and leaky layers from the top down
        if float, porosity is the same for all layers
        if topboundary='conf': length is 2 * number of aquifers - 1
        if topboundary='semi': length is 2 * number of aquifers
    order : int
        polynomial order of flux along each segment
    ndeg : int
        number of points used between two segments to numerically
        integrate normal discharge
    layers: list or np.array
        layers in which leaky wall is present.
    res : float or np.array, optional
        resistance of leaky wall, if passed as an array must be either
        shape (n_segments,) or (n_segments, n_layers). Default is np.inf,
        which simulates an impermeable wall.
    refine_level : int, optional
        refine element by partitioning each side into refine_level segments, default
        is 1, which means no refinement is applied.
    """

    def __init__(
        self,
        model,
        xy,
        kaq=1.0,
        z=[1, 0],
        c=[],
        npor=0.3,
        topboundary="conf",
        hstar=None,
        order=3,
        ndeg=3,
        layers=[0],
        res=np.inf,
        refine_level=1,
        addtomodel=True,
    ):
        _input = {k: v for k, v in locals().items() if k not in ["self"]}
        (kaq, c, npor, ltype) = param_maq(kaq, z, c, npor, topboundary)
        super().__init__(
            model=model,
            xy=xy,
            kaq=kaq,
            z=z,
            c=c,
            npor=npor,
            ltype=ltype,
            hstar=hstar,
            order=order,
            ndeg=ndeg,
            layers=layers,
            res=res,
            refine_level=refine_level,
            addtomodel=addtomodel,
        )
        self._input = _input


class LeakyBuildingPit3D(LeakyBuildingPit):
    def __init__(
        self,
        model,
        xy,
        kaq=1.0,
        z=[1, 0],
        kzoverkh=1.0,
        npor=0.3,
        topboundary="conf",
        topres=0.0,
        topthick=0.0,
        hstar=None,
        order=3,
        ndeg=3,
        layers=[0],
        res=np.inf,
        refine_level=1,
        addtomodel=True,
    ):
        """Element to simulate a building pit with a leaky wall in Model3D.

        Parameters
        ----------
        model : Model object
            model to which the element is added
        xy : array or list
            list or array of (x,y) pairs of coordinates of corners of the
            inhomogeneity
            polygonal boundary is automatically closed (so first point
            is not repeated)
        kaq : float, array or list
            hydraulic conductivity of each aquifer from the top down
            if float, hydraulic conductivity is the same in all aquifers
        kzoverkh : float, array or list
            vertical anisotropy ratio vertical k divided by horizontal k
            if float, value is the same for all layers, length is number of layers
        z : array or list
            elevation tops and bottoms of the aquifers from the top down
            leaky layers may have zero thickness
            if topboundary='conf': length is 2 * number of aquifers
            if topboundary='semi': length is 2 * number of aquifers + 1 as top
            of leaky layer on top of systems needs to be specified
        topboundary : string, 'conf' or 'semi' (default is 'conf')
            indicating whether the top is confined ('conf') or
            semi-confined ('semi'). For a building pit, the 'conf'
            option is generally more applicable.
        topres : float
            resistance of top leaky layer (read if topboundary="semi")
        topthick : float
            thickness of top semi-confining layer (read if topboundary='semi')
        hstar : float or None (default is None)
            head value above semi-confining top, only read if topboundary='semi'
        npor : float, array or list
            porosity of all aquifers and leaky layers from the top down
            if float, porosity is the same for all layers
            if topboundary='conf': length is 2 * number of aquifers - 1
            if topboundary='semi': length is 2 * number of aquifers
        order : int
            polynomial order of flux along each segment
        ndeg : int
            number of points used between two segments to numerically
            integrate normal discharge
        layers: list or np.array
            layers in which leaky wall is present.
        res : float or np.array, optional
            resistance of leaky wall, if passed as an array must be either
            shape (n_segments,) or (n_segments, n_layers). Default is np.inf,
            which simulates an impermeable wall.
        refine_level : int, optional
            refine element by partitioning each side into refine_level segments, default
            is 1, which means no refinement is applied.
        """
        _input = {k: v for k, v in locals().items() if k not in ["self"]}
        (kaq, c, npor, ltype) = param_3d(kaq, z, kzoverkh, npor, topboundary, topres)
        if topboundary == "semi":
            z = np.hstack((z[0] + topthick, z))
        super().__init__(
            model=model,
            xy=xy,
            kaq=kaq,
            z=z,
            c=c,
            npor=npor,
            ltype=ltype,
            hstar=hstar,
            order=order,
            ndeg=ndeg,
            layers=layers,
            res=res,
            refine_level=refine_level,
            addtomodel=addtomodel,
        )
        self._input = _input


class AreaSinkInhom(Element):
    def __init__(self, model, N, xc, name="AreaSinkInhom", label=None, aq=None):
        Element.__init__(
            self, model, nparam=1, nunknowns=0, layers=0, name=name, label=label
        )
        self.N = N
        self.xc = xc  # x-center of area-sink
        # self.nparam = 1  # Defined here and not in Element as other elements can have multiple parameters per layers
        # self.nunknowns = 0
        self.aq = aq
        self.model.add_element(self)

    def __repr__(self):
        return self.name

    def initialize(self):
        assert self.aq is not None, "Error: no aquifer passed"
        self.aq.add_element(self)
        self.plabsq = self.aq.coef[self.layers, 1:] * self.aq.lab[1:] ** 2
        self.parameters = np.atleast_2d(self.N)

    def potinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((1, aq.naq))
        if aq == self.aq:
            rv[0, 0] = -0.5 * (x - self.xc) ** 2
            rv[0, 1:] = self.plabsq
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq == self.aq:
            rv[0, 0, 0] = x - self.xc
        return rv

    def qztop(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = 0.0
        if aq == self.aq:
            rv = -self.parameters[0, 0]
        return rv

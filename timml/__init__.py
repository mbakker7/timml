"""Copyright (C), 2015, Mark Bakker.

Mark Bakker, Delft University of Technology
mark dot bakker at tudelft dot nl

TimML is a computer program for the simulation of steady-state multi-aquifer flow with
analytic elements and consists of a library of Python scripts and FORTRAN extensions.
"""

# ruff: noqa: F401
# --version number
__author__ = "Mark Bakker"
# Import all classes and functions
from timml import bessel
from timml.circareasink import CircAreaSink
from timml.constant import Constant, ConstantStar
from timml.inhomogeneity import (
    BuildingPit3D,
    BuildingPitMaq,
    LeakyBuildingPit3D,
    LeakyBuildingPitMaq,
    PolygonInhom3D,
    PolygonInhomMaq,
)
from timml.inhomogeneity1d import StripInhom3D, StripInhomMaq, Xsection3D, XsectionMaq
from timml.linedoublet import (
    ImpLineDoublet,
    ImpLineDoubletString,
    LeakyLineDoublet,
    LeakyLineDoubletString,
)
from timml.linedoublet1d import ImpLineDoublet1D, LeakyLineDoublet1D
from timml.linesink import (
    CollectorWell,
    HeadLineSink,
    # HeadLineSinkContainer,
    HeadLineSinkString,
    # HeadLineSinkZero,
    LineSinkBase,
    LineSinkDitch,
    LineSinkDitchString,
    RadialCollectorWell,
)
from timml.linesink1d import HeadLineSink1D, LineSink1D
from timml.model import Model, Model3D, ModelMaq, ModelXsection
from timml.stripareasink import XsectionAreaSink
from timml.trace import timtraceline, timtracelines
from timml.uflow import Uflow
from timml.version import __version__, show_versions
from timml.well import (
    HeadWell,
    HeadWellString,
    LargeDiameterWell,
    TargetHeadWell,
    TargetHeadWellString,
    Well,
    WellBase,
    WellString,
    WellStringBase,
)

# default bessel module is numba
bessel.set_bessel_method(method="numba")

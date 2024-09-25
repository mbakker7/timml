"""Copyright (C), 2015, Mark Bakker.

Mark Bakker, Delft University of Technology
mark dot bakker at tudelft dot nl

TimML is a computer program for the simulation of steady-state multiaquifer flow with
analytic elements and consists of a library of Python scripts and FORTRAN extensions.
"""
# ruff: noqa: F401
# from __future__ import division, print_function, absolute_import

# --version number
__name__ = "timml"
__author__ = "Mark Bakker"
# Import all classes and functions
from . import bessel
from .circareasink import CircAreaSink
from .constant import Constant, ConstantStar
from .inhomogeneity import (
    BuildingPit3D,
    BuildingPitMaq,
    LeakyBuildingPit3D,
    LeakyBuildingPitMaq,
    PolygonInhom3D,
    PolygonInhomMaq,
)
from timml.inhomogeneity1d import StripInhom3D, StripInhomMaq
from timml.linedoublet import (
    ImpLineDoublet,
    ImpLineDoubletString,
    LeakyLineDoublet,
    LeakyLineDoubletString,
)
from timml.linedoublet1d import ImpLineDoublet1D, LeakyLineDoublet1D
from timml.linesink import (
    HeadLineSink,
    HeadLineSinkContainer,
    HeadLineSinkString,
    HeadLineSinkZero,
    LineSinkBase,
    LineSinkDitch,
    LineSinkDitchString,
)
from timml.linesink1d import HeadLineSink1D, LineSink1D
from timml.model import Model, Model3D, ModelMaq
from timml.stripareasink import StripAreaSink
from timml.trace import timtraceline, timtracelines
from timml.uflow import Uflow
from timml.version import __version__
from timml.well import HeadWell, LargeDiameterWell, Well, WellBase

__all__ = [
    "CircAreaSink",
    "Constant",
    "ConstantStar",
    "BuildingPit3D",
    "BuildingPitMaq",
    "LeakyBuildingPit3D",
    "LeakyBuildingPitMaq",
    "PolygonInhom3D",
    "PolygonInhomMaq",
    "StripInhom3D",
    "StripInhomMaq",
    "ImpLineDoublet",
    "ImpLineDoubletString",
    "LeakyLineDoublet",
    "LeakyLineDoubletString",
    "ImpLineDoublet1D",
    "LeakyLineDoublet1D",
    "HeadLineSink",
    "HeadLineSinkContainer",
    "HeadLineSinkString",
    "HeadLineSinkZero",
    "LineSinkBase",
    "LineSinkDitch",
    "LineSinkDitchString",
    "HeadLineSink1D",
    "LineSink1D",
    "Model",
    "Model3D",
    "ModelMaq",
    "StripAreaSink",
    "timtraceline",
    "timtracelines",
    "Uflow",
    "__version__",
    "HeadWell",
    "LargeDiameterWell",
    "Well",
    "WellBase",
]

# default bessel module is numba
bessel.set_bessel_method(method="numba")

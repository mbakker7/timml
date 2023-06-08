'''
Copyright (C), 2015, Mark Bakker.
Mark Bakker, Delft University of Technology
mark dot bakker at tudelft dot nl

TimML is a computer program for the simulation of steady-state
multiaquifer flow with analytic elements and consists of a
library of Python scripts and FORTRAN extensions.
'''
#from __future__ import division, print_function, absolute_import

#--version number
__name__='timml'
__author__='Mark Bakker'
# Import all classes and functions
from .circareasink import CircAreaSink
from .constant import Constant, ConstantStar
from .inhomogeneity import (
    BuildingPit,
    LeakyBuildingPit,
    PolygonInhom3D,
    PolygonInhomMaq,
)
from .inhomogeneity1d import StripInhom3D, StripInhomMaq
from .linedoublet import (
    ImpLineDoublet,
    ImpLineDoubletString,
    LeakyLineDoublet,
    LeakyLineDoubletString,
)
from .linedoublet1d import ImpLineDoublet1D, LeakyLineDoublet1D
from .linesink import (
    HeadLineSink,
    HeadLineSinkContainer,
    HeadLineSinkString,
    HeadLineSinkZero,
    LineSinkBase,
    LineSinkDitch,
    LineSinkDitchString,
)
from .linesink1d import HeadLineSink1D, LineSink1D

from .model import Model, Model3D, ModelMaq
from .stripareasink import StripAreaSink
from .trace import timtraceline, timtracelines
from .uflow import Uflow
from .version import __version__
from .well import HeadWell, Well, WellBase

__all__ = [s for s in dir() if not s.startswith("_")]

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
from .version import __version__

# Import all classes and functions
from .model import ModelMaq, Model3D, Model
from .well import WellBase, Well, HeadWell
from .constant import Constant, ConstantStar
from .linesink import LineSinkBase, HeadLineSinkZero, HeadLineSink, \
                      LineSinkDitch, HeadLineSinkString, LineSinkDitchString, \
                      HeadLineSinkContainer
from .linedoublet import ImpLineDoublet, ImpLineDoubletString, \
                         LeakyLineDoublet, LeakyLineDoubletString
from .circareasink import CircAreaSink
from .inhomogeneity import PolygonInhomMaq, PolygonInhom3D, BuildingPit
from .inhomogeneity1d import StripInhomMaq, StripInhom3D
from .uflow import Uflow
from .trace import timtraceline, timtracelines
from .linesink1d import LineSink1D, HeadLineSink1D
from .linedoublet1d import ImpLineDoublet1D, LeakyLineDoublet1D
from .stripareasink import StripAreaSink

__all__ = [s for s in dir() if not s.startswith("_")]

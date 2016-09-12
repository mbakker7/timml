'''
Copyright (C), 2015, Mark Bakker.
Mark Bakker, Delft University of Technology
mark dot bakker at tudelft dot nl

TimML is a computer program for the simulation of steady-state
multiaquifer flow with analytic elements and consists of a
library of Python scripts and FORTRAN extensions.
'''

#--version number
__name__='timml'
__author__='Mark Bakker'
from version import __version__

# Import all classes and functions
from model import ModelMaq
from well import WellBase, Well, HeadWell
from constant import Constant
from linesink import HeadLineSink, HeadLineSinkHo, HeadLineSinkString
from inhomogeneity import PolygonInhomMaq
from uflow import Uflow
from util import timcontour




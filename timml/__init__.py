'''
Copyright (C), 2015, Mark Bakker.
Mark Bakker, Delft University of Technology
mark dot bakker at tudelft dot nl

TimML is a computer program for the simulation of steady-state
multiaquifer flow with analytic elements and consists of a
library of Python scripts and FORTRAN extensions.
'''
from __future__ import division, print_function, absolute_import

#--version number
__name__='timml'
__author__='Mark Bakker'
from .version import __version__

# Import all classes and functions
from .model import *
from .well import *
from .constant import *
from .linesink import *
from .linedoublet import *
from .circareasink import *
from .inhomogeneity import *
from .inhomogeneity1d import *
from .uflow import *
from .trace import *
from .linesink1d import *
from .linedoublet1d import *
from .stripareasink import *

__all__ = [s for s in dir() if not s.startswith("_")]

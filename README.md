[![Build Status](https://travis-ci.org/jentjr/timml.svg?branch=master)](https://travis-ci.org/jentjr/timml)
[![Build status](https://ci.appveyor.com/api/projects/status/h5y9fpjdb092kphg/branch/master?svg=true)](https://ci.appveyor.com/project/jentjr/timml/branch/master)

# TimML, A Multi-Layer, Analytic Element Model

## Introduction

TimML is a computer program for the modeling of steady-state multi-layer flow with analytic elements
and consists of a library of Python scripts and FORTRAN extensions.
TimML may be applied to an arbitrary number of aquifers and leaky layers.
The head, flow, and leakage between aquifers may be computed analytically at any point in the aquifer system.
The design of TimML is object-oriented and has been kept simple and flexible.
New analytic elements may be added to the code without making any changes in the existing part of the code.
TimML is coded in Python; use is made of FORTRAN extensions to improve performance.

## Installation

**Python versions:**

TimML requires **Python** > 3.5 and can be installed from PyPI.
The PyPI installation includes compiled versions of the FORTRAN extension
for both Windows and MacOS.


**Dependencies:**

TimML requires **NumPy** 1.12 (or higher) and **matplotlib** 2.0 (or higher). 

**For base Python distributions:**

To install TimML, open a command prompt and type:

    pip install timml

To update TimML type:

    pip install timml --upgrade

To uninstall TimML type:

    pip uninstall timml
    
## Documentation

* The manual is available from the docs directory or can be viewed [here](http://mbakker7.github.io/timml/docs/builddocs/html/index.html).
* Example Notebooks are available from the notebooks directory on github, of from [here](https://github.com/mbakker7/timml/tree/master/notebooks).

## TimML Version 5

TimML version 5 is a total rewrite and is not backwards compatible with previous TimML versions.
TimML version 5 is intended to be compatible with TTim.
TimML version 5 has many new features and elements, the code base is Python 3, and the object oriented design is much simpler.
TimML version 4 remains available through the timml4 branch.

## Release
TimML 5.0.1 - alpha release

## Citation

Some of the papers that you may want to cite when using TimML are

* Bakker, M., and O.D.L. Strack. 2003. Analytic Elements for Multiaquifer Flow. Journal of Hydrology, 271(1-4), 119-129.
* Bakker, M. 2006. An analytic element approach for modeling polygonal inhomogeneities in multi-aquifer systems. Advances in Water Resources, 29(10), 1546-1555.

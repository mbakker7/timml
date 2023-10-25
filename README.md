[![timml](https://github.com/mbakker7/timml/actions/workflows/ci.yml/badge.svg)](https://github.com/mbakker7/timml/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/mbakker7/timml/badge.svg?branch=master)](https://coveralls.io/github/mbakker7/timml?branch=master)
![PyPI](https://img.shields.io/pypi/v/timml?color=green)

# TimML, A Multi-Layer, Analytic Element Model

## Introduction

TimML is a computer program for the modeling of steady-state multi-layer flow with analytic elements
and consists of a library of Python scripts and FORTRAN extensions.
TimML may be applied to an arbitrary number of aquifers and leaky layers.
The head, flow, and leakage between aquifers may be computed analytically at any point in the aquifer system.
The design of TimML is object-oriented and has been kept simple and flexible.
New analytic elements may be added to the code without making any changes in the existing part of the code.
TimML is coded in Python and uses numba to speed up evaluation of the bessel line elements.

## Installation

**Python versions:**

TimML requires **Python** >= 3.7 and can be installed from PyPI.

**Dependencies:**

TimML requires **numpy** >=1.17, **scipy** >=1.5 and **matplotlib** >=3.1, **numba>=0.5**.

**Installation:**

To install TimML, open a command prompt and type:

    pip install timml

To update TimML type:

    pip install timml --upgrade

To uninstall TimML type:

    pip uninstall timml


## Documentation

* The documentation is hosted on [readthedocs](https://timml.readthedocs.io/).
* Example Notebooks are available from the notebooks directory on github, of from [here](https://github.com/mbakker7/timml/tree/master/notebooks).

## TimML Version 6

TimML version 6 has the same functionality as version 5, but doesn't depend on a fortran extension anymore, so installation is easy on all platforms.
TimML version 5 is a total rewrite and is not backwards compatible with previous TimML versions.
TimML version 5 is intended to be compatible with TTim.
TimML version 5 has many new features and elements, the code base is Python 3, and the object oriented design is much simpler.
TimML version 4 remains available through the timml4 branch.

## Release
TimML 6.0. First release that depends on numba and doesn't depend on fortran extension anymore. Code is now pure python (with numba for speed).

## Citation

Some of the papers that you may want to cite when using TimML are

* Bakker, M., and O.D.L. Strack. 2003. Analytic Elements for Multiaquifer Flow. Journal of Hydrology, 271(1-4), 119-129.
* Bakker, M. 2006. An analytic element approach for modeling polygonal inhomogeneities in multi-aquifer systems. Advances in Water Resources, 29(10), 1546-1555.

![image](/docs/_static/tim_logo.png)

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

TimML requires Python >= 3.8 and can be installed from PyPI.

**Dependencies:**

TimML requires:
* numpy 
* scipy
* matplotlib
* numba

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

## Latest release
TimML 0.6.5
* Improved documentation: new look, better organization, tutorials, how-to guides etc. Check it out [here](https://timml.readthedocs.io/)!
* New elements
    * Building pit elements for 3D (multi-layer single aquifer) models.
    * Large diamater wells (only for radial flow).
* Enhancements
    * Building pit leaky wall resistance can be set per layer and per side (for modeling leaks or gaps).
    * Return integrated normal flux per layer and per line segment.

## TimML Versions

* TimML version 0.6 has the same functionality as version 5, but doesn't depend on a fortran extension anymore, so installation is easy on all platforms.
* TimML version 0.5 
    * is a total rewrite and is not backwards compatible with previous TimML versions.
    * is intended to be compatible with TTim.
    * has many new features and elements, the code base is Python 3, and the object oriented design is much simpler.
* TimML version 0.4 remains available through the timml4 branch.

## Citation

Some of the papers that you may want to cite when using TimML are

* Bakker, M., and O.D.L. Strack. 2003. Analytic Elements for Multiaquifer Flow. Journal of Hydrology, 271(1-4), 119-129.
* Bakker, M. 2006. An analytic element approach for modeling polygonal inhomogeneities in multi-aquifer systems. Advances in Water Resources, 29(10), 1546-1555.

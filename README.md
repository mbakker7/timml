![image](/docs/_static/tim_logo_small.png)
[![timml](https://github.com/mbakker7/timml/actions/workflows/ci.yml/badge.svg)](https://github.com/mbakker7/timml/actions/workflows/ci.yml)
![PyPI](https://img.shields.io/pypi/v/timml?color=green)

> [!IMPORTANT]  
> TimML has moved! See [timflow](https://github.com/timflow-org/timflow). TimML is available through `timflow.steady`. `TimML` will no longer be developed. It will be maintained for the time being. All future development will focus on `timflow`.

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
* Example Notebooks are also available from the docs directory on github.

## Citation

Some of the papers that you may want to cite when using TimML are:

* Bakker, M., and O.D.L. Strack. 2003. Analytic Elements for Multiaquifer Flow. Journal of Hydrology, 271(1-4), 119-129.
* Bakker, M. 2006. An analytic element approach for modeling polygonal inhomogeneities in multi-aquifer systems. Advances in Water Resources, 29(10), 1546-1555.

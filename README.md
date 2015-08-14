# TimML, A Multi-Layer, Analytic Element Model

## Introduction

TimML is a computer program for the modeling of steady-state multi-layer flow with analytic elements
and consists of a library of Python scripts and FORTRAN extensions.
TimML may be applied to an arbitrary number of aquifers and leaky layers.
The head, flow, and leakage between aquifers may be computed analytically at any point in the aquifer system.
The design of TimML is object-oriented and has been kept simple and flexible.
New analytic elements may be added to the code without making any changes in the existing part of the code.
TimML is coded in Python; use is made of FORTRAN extensions to improve performance.

## TimML Changes

### Version 4.0.1
In version 4.0, the numbering of the layers in TimML has changed. Layers are numbered from the top down starting at number 0.
(in TimML versions prior to 4.0, layer numbers started at 1).

## Installation

**Python versions:**

TimML requires **Python** 2.7 and can be installed from PyPI.
The PyPI installation includes compiled versions of the FORTRAN extension
for both Windows and MaxOS.


**Dependencies:**

TimML requires **NumPy** 1.9 (or higher) and **matplotlib** 1.4 (or higher). T

**For base Python distributions:**

To install TimML, open a command prompt and type:

    pip install timml

To update TimML type:

    pip install timml --upgrade

To uninstall TimML type:

    pip uninstall timml

## Documentation

* The manual is available from the docs directory.
* Example Notebooks are available from the notebooks directory

# TimML, A Multi-Layer, Analytic Element Model

## Introduction

TimML is a computer program for the modeling of steady-state multi-layer flow with analytic elements
and consists of a library of Python scripts and FORTRAN extensions.
TimML may be applied to an arbitrary number of aquifers and leaky layers.
The head, flow, and leakage between aquifers may be computed analytically at any point in the aquifer system.
The design of TimML is object-oriented and has been kept simple and flexible.
New analytic elements may be added to the code without making any changes in the existing part of the code.
TimML is coded in Python; use is made of FORTRAN extensions to improve performance.

## TimML Version 5

TimML version 5 is a total rewrite and is not backwards compatible with previous TimML versions.
TimML version 5 is intended to be compatible with TTim.
TimML version 5 has many new features and elements, the code base is Python 3, and the object oriented design is much simpler.
TimML version 4 remains available through the timml4 branch.

## Release
The first release of TimML Version 5 is planned for the summer of 2017.

## Citation

Some of the papers that you may want to cite when using TimML are

* Bakker, M., and O.D.L. Strack. 2003. Analytic Elements for Multiaquifer Flow. Journal of Hydrology, 271(1-4), 119-129.
* Bakker, M. 2006. An analytic element approach for modeling polygonal inhomogeneities in multi-aquifer systems. Advances in Water Resources, 29(10), 1546-1555.

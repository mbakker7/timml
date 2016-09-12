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

### Version 4.0
In version 4.0, the numbering of the layers in TimML has changed. Layers are numbered from the top down starting at number 0.
(in TimML versions prior to 4.0, layer numbers started at 1).

## Installation

**Python versions:**

TimML requires **Python** 2.7 and can be installed from PyPI.
The PyPI installation includes compiled versions of the FORTRAN extension
for both Windows and MaxOS.


**Dependencies:**

TimML requires **NumPy** 1.9 (or higher) and **matplotlib** 1.4 (or higher). 

**For base Python distributions:**

To install TimML, open a command prompt and type:

    pip install timml

To update TimML type:

    pip install timml --upgrade

To uninstall TimML type:

    pip uninstall timml
    
**Testing installation:**

    ipython
    import timml.timmltest
    
An example model is imported and a contour plot is shown. When this is run from the regular Python prompt (not IPython), the
model is created and solved but the contour plot is probably not shown (depending on your default settings of matplotlib). 

## Documentation

* The manual is available from the docs directory or can be viewed [here](https://github.com/mbakker7/timml/blob/master/docs/timml.pdf).
Once you click on this link, you can click on the [Download] button to download the pdf file. 
* Example Notebooks are available from the git repository.

## Citation

Some of the papers that you may want to cite when using TimML are

* Bakker, M., and O.D.L. Strack. 2003. Analytic Elements for Multiaquifer Flow. Journal of Hydrology, 271(1-4), 119-129.
* Bakker, M. 2006. An analytic element approach for modeling polygonal inhomogeneities in multi-aquifer systems. Advances in Water Resources, 29(10), 1546-1555.
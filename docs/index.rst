.. timml documentation master file, created by
   sphinx-quickstart on Mon Jul 24 15:58:22 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

============
TimML
============
TimML is a computer program for the modeling of steady-state multi-layer flow with analytic elements
and consists of a library of Python scripts and FORTRAN extensions.
TimML may be applied to an arbitrary number of aquifers and leaky layers.
The head, flow, and leakage between aquifers may be computed analytically at any point in the aquifer system.
The design of TimML is object-oriented and has been kept simple and flexible.
New analytic elements may be added to the code without making any changes in the existing part of the code.
TimML is coded in Python; use is made of FORTRAN extensions to improve performance.

General Outline

Contents:

.. toctree::
    :maxdepth: 2
    :hidden:
    
    Elements <modules>

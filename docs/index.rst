.. timml documentation master file, created by
   sphinx-quickstart on Mon Jul 24 15:58:22 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

================================================
TimML - THIS DOCUMENTATION IS UNDER CONSTRUCTION
================================================
TimML is a computer program for the modeling of steady-state multi-layer flow with analytic elements
TimML may be applied to an arbitrary number of layers and arbitrary sequence of aquifers and leaky layers.
The Dupuit approximation is applied to aquifer layers, while flow in leaky layers is approximated as vertical.
The head, flow, and leakage between aquifer layers may be computed analytically at any point in the aquifer system.
The design of TimML is object-oriented and has been kept simple and flexible.
New analytic elements may be added to the code without making any changes in the existing part of the code.
TimML is coded in Python. Behind the scenes, use is made of FORTRAN extensions to improve performance.

Installation
------------
TimML is written for Python 3.

  installation instructions on pip
  
Main Approximations
-------------------

List of available elements
--------------------------
* Well

  * Discharge-specified well
  * Head-specified well
  * Multi-aquifer well. Well is screened in multiple layers and only total discharge is specified.
  
* Line-sink

  * Head-specified line-sink
  * String of head-specified line-sinks


  
.. toctree::
    :maxdepth: 3
    :hidden:
    
    Models <models/modelindex>
    Elements <aems>

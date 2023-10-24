Introduction
============

TimML is a Python pacakge for the modeling of steady-state multi-layer groundwater flow
with analytic elements.

TimML may be applied to an arbitrary number of layers and arbitrary sequence of
aquifers and leaky layers. The head, flow, and leakage between aquifer layers may be
computed analytically at any point in the aquifer system.

.. grid::

    .. grid-item-card:: User Guide
        :link: userguide/index
        :link-type: doc

        User guide on the basic concepts of TimML.

    .. grid-item-card:: Examples
        :link: examples/index
        :link-type: doc

        Examples of TimML usage.

    .. grid-item-card:: Code Reference
        :link: api/index
        :link-type: doc

        TimML code reference.


Quick Example
-------------

.. tab-set::

    .. tab-item:: Python

        In this example a well is modelled near a river in a single aquifer.

        .. code-block:: python

            # Import python packages
            import numpy as np
            import timml

            # Create model
            ml = timml.ModelMaq(kaq=10, z=[20, 0])
            
            # Add a river with a fixed water level
            yls = np.arange(-100.0, 101, 20)
            xls = 50.0 * np.ones_like(yls)
            river = timml.HeadLineSinkString(ml, xy=list(zip(xls, yls)), hls=0.0)
            
            # Add a well
            well = timml.Well(ml, 0.0, 0.0, rw=0.3, Qw=1000)
            
            # Solve model
            ml.solve()

            # Plot head contours at t=2 days
            ml.contour(win=[-30, 55, -30, 30], ngr=40, labels=True, decimals=1)
            

    .. tab-item:: Result

        .. figure:: _static/example_output.png
            :figwidth: 500px


Approximations
--------------

The Dupuit approximation is applied to aquifer layers, while flow in leaky layers is
approximated as vertical.



.. toctree::
   :maxdepth: 2
   :hidden:

    User Guide <userguide/index>
    Examples <examples/index>
    Code reference <api/index>
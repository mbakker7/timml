Introduction
============

TimML is a Python package for the modeling of steady-state multi-layer groundwater flow
with analytic elements.

TimML may be applied to an arbitrary number of layers and an arbitrary sequence of
aquifers and leaky layers. The head, flow, and leakage between aquifer layers may be
computed analytically at any point in the aquifer system. The Dupuit approximation is applied to flow in aquifer layers (i.e., the resistance to flow in the vertical direction is neglected), while flow in leaky layers is approximated as vertical.

.. grid::

    .. grid-item-card:: Tutorials
        :link: 00tutorials/index
        :link-type: doc

        Tutorials for getting started with TimML.

    .. grid-item-card:: Concepts
        :link: 02concepts/index
        :link-type: doc

        TimML basic concepts explained.

    .. grid-item-card:: How-to guides
        :link: 01howto/index
        :link-type: doc

        How-to guides for more advanced modeling with TimML.

.. grid::

    .. grid-item-card:: Examples
        :link: 03examples/index
        :link-type: doc

        TimML example notebooks.

    .. grid-item-card:: Cross-sections
        :link: 04xsections/index
        :link-type: doc

        Cross-sectional models explained.

    .. grid-item-card:: Code reference
        :link: 05api/index
        :link-type: doc

        TimML code reference.


Quick Example
-------------

.. tab-set::

    .. tab-item:: Python

        In this example a well is modelled near a river in a single aquifer.

        .. code-block:: python

            # import python packages
            import numpy as np
            import timml

            # create model
            ml = timml.ModelMaq(kaq=10, z=[20, 0]) # single layer model
            
            # add a river with a fixed water level
            yls = np.arange(-100, 101, 20) # 20 points, so 19 segments
            xls = 50 * np.ones_like(yls)
            river = timml.HeadLineSinkString(ml, xy=list(zip(xls, yls)), hls=0.0)
            
            # add a well
            well = timml.Well(ml, 0, 0, rw=0.3, Qw=1000)
            
            # solve model
            ml.solve()

            # plot head contours
            ml.contour(win=[-30, 55, -30, 30], ngr=40, labels=True, decimals=1)
            

    .. tab-item:: Result

        In this example a well is modelled near a river in a single aquifer.

        .. figure:: _static/example_output.png
            :figwidth: 500px


.. toctree::
   :maxdepth: 2
   :hidden:

    Tutorials <00tutorials/index>
    How-to guides <01howto/index>
    Concepts <02concepts/index>
    Examples <03examples/index>
    Cross-sections <04xsections/index>
    Code reference <05api/index>
    Cite <06about/index>
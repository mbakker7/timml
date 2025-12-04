Cross-sectional Modeling
========================

Cross-sectional (1D) models in TimML are used to simulate groundwater flow in problems
with symmetry along one axis. This section demonstrates how to build and use
cross-sectional models with strip inhomogeneities, infinitely long line-sinks, and
infinitely long line-doublets.

Elements
--------

Multi-layer cross-sectional (1D) models may consist of strip inhomogeneities,
infinitely long line-sinks, and infinitely long line-doublets. The cross-section is
parallel to the x-axis. Strip inhomogeneities may extend to infinity. The following
elements are available:

1. :class:`~timml.inhomogeneity1d.XsectionMaq` is a strip aquifer consisting of a
   regular sequence of aquifer - leaky layer - aquifer - leaky layer, aquifer, etc., like
   a ModelMaq.

2. :class:`~timml.inhomogeneity1d.Xsection3D` is a strip aquifer consisting of a
   stack of aquifer layers.

3. :class:`~timml.linesink1d.LineSink1D` is an infinitely long line-sink for which the
   discharge per unit length is specified

4. :class:`~timml.linesink1d.HeadLineSink1D` is an infinitely long line-sink for which
   the head and (optionally) a resistance is specified

5. :class:`~timml.linedoublet1d.LeakyLineDoublet1D` is an infinitely long leaky wall.
   An impermeable wall is created by specifying an infinitely large resistance

6. :class:`~timml.stripareasink.XsectionAreaSink` is a strip area-sink, which can only be
   added when the topboundary is confined (topboundary='conf').

Examples
--------

- :doc:`cross_section_models`

.. toctree::
    :maxdepth: 3
    :hidden:

    cross_section_models

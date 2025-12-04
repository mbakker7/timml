Cross-section Elements
======================

Cross-section (1D) models are used to simulate groundwater flow in problems with symmetry
along one axis. The cross-section is parallel to the x-axis. The following elements are
available for cross-section models:

Inhomogeneities
---------------

1. :class:`~timml.inhomogeneity1d.XsectionMaq` is a strip aquifer consisting of a regular
   sequence of aquifer - leaky layer - aquifer - leaky layer, aquifer, etc., like a
   ModelMaq.

2. :class:`~timml.inhomogeneity1d.Xsection3D` is a strip aquifer consisting of a stack of
   aquifer layers.


Line-sinks
----------

1. :class:`~timml.linesink1d.LineSink1D` is an infinitely long line-sink for which the
   discharge per unit length is specified.

2. :class:`~timml.linesink1d.HeadLineSink1D` is an infinitely long line-sink for which
   the head and (optionally) a resistance is specified.

Line-doublets
-------------

1. :class:`~timml.linedoublet1d.ImpLineDoublet1D` is an infinitely long impermeable wall.

2. :class:`~timml.linedoublet1d.LeakyLineDoublet1D` is an infinitely long leaky wall with
   specified resistance.

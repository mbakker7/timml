Xsection Models
---------------

Multi-layer cross-sectional models may consist of strip inhomogeneities, infinitely long line-sinks, and infinitely long line-doublets.
The cross-section is parallel to the x-axis. Strip inhomogeneities may extend to infinity. The following elements are available:

1. :class:`~timml.inhomogeneity1d.StripInhomMaq` is a strip aquifer consisting of a regular sequence of aquifer - leaky layer - aquifer - leaky layer, aquifer, etc., like a ModelMaq.

2. :class:`~timml.inhomogeneity1d.StripInhom3D` is a strip aquifer consisting of a stack of aquifer layers.

3. :class:`~timml.linesink1d.LineSink1D` is an infinitely long line-sink for which the discharge per unit length is specified

4. :class:`~timml.linesink1d.HeadLineSink1D` is an infinitely long line-sink for which the head and (optionally) a resistance is specified

5. :class:`~timml.linedoublet1d.LeakyLineDoublet1D` is an infinitely long leaky wall. An impermeable wall is created by specifying an infinitely large resistance

    
.. toctree::
    :maxdepth: 1
    :hidden:
    
    StripInhomMaq <stripinhommaq>
    StripInhom3D <stripinhom3d>
    LineSink1D <linesink1d>
    HeadLineSink1D <headlinesink1d>
    LeakyLineDoublet1D <leakylinedoublet1d>
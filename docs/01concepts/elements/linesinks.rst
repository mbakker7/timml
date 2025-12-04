Line-sinks
==========

Line-sinks are lines along which water is taken out of or put into the aquifer. The
inflow/outflow along the line-sink varies as a polynomial (the order of the polynomial
may be specified).

There are two types of line-sinks: line-sinks for which the head is specified along the
line-sink, and line-sinks for which the total discharge is specified and the head along
the line-sink is unknown but uniform.

Both types of line-sinks may have an entry resistance defined by the resistance
:math:`c` (dimension: time). The inflow :math:`\sigma_i` into the line-sink in layer
:math:`i` is a function of the head :math:`h_i` just outside the line-sink in layer
:math:`i` and the head :math:`h_{ls}` inside the line-sink:

    .. math::
        \sigma_i = w(h_i - h_{ls})/c
        
This equation is applied along each control point of the line-sink. :math:`w` is the
distance over which water infiltrates into the line-sink. The distance may be the width
of the stream, for example, in case of a partially penetrating stream. In case the
stream penetrates the aquifer layer fully, the distance may equal the thickness of the
aquifer layer (if water enters primarily from one side), or twice the aquifer thickness
(if water enters from both sides).

The two types of line-sinks may be entered one by one or as a string. For a string of
head-specified line-sinks, the head needs only be specified at the beginning and end of
the string and the head is interpolated between these values. Alternatively, the head
may be specified at the beginning and end of some or all line-sinks in the string For a
string of line-sink ditch elements, the head is uniform along the entire string while
the total discharge of all the line-sinks in the string is equal to the specified
value.

1. :class:`~timml.linesink.HeadLineSink` is a line-sink for which the head is specified
   along the line-sink.

2. :class:`~timml.linesink.HeadLineSinkString` is a string of head-specified line-sinks

3. :class:`~timml.linesink.LineSinkDitch` is a line-sink for which the head is unknown
   and uniform and the total discharge is specified

4. :class:`~timml.linesink.LineSinkDitchString` is a string of line-sink ditch elements.

5. :class:`~timml.linesink.CollectorWell` is a special line-sink representing a collector
   well with an arbitrary number of connected arms.

6. :class:`~timml.linesink.RadialCollectorWell` is a special line-sink representing a
   collector well with radial arms extending from a central caisson.

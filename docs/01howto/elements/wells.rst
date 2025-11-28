Wells
=====

There are two types of wells: wells for which the total discharge is specified and
wells for which the head inside the well is specified. Both types of wells may have an
entry resistance (skin effect) defined by the resistance :math:`c` (dimension: time).
The discharge :math:`Q_i` in layer :math:`i` is a function of the head :math:`h_i` in
layer :math:`i` just outside the well and the head :math:`h_w` inside the well:

    .. math::
        Q_i = 2\pi r_w H_i(h_i - h_w)/c

1. :class:`~timml.well.Well` is a well for which the total discharge is specified. The
   total discharge is distributed across the layers in which the well is screened such
   that the head inside the well is the same in each screened layer.

2. :class:`~timml.well.HeadWell` is a well for which the head inside the well is
   specified. The discharge in each layer is computed such that the head in all screened
   layers is equal to the specified head.

Additionally, there are more complex well types that consist of one or more Wells or
LineSinks.

1. :class:`~timml.well.WellString` is a well type that consists of multiple connected
   wells (e.g. multiple wells connected to the same pump) for which the total discharge is
   given. The total discharge is distributed across the connected wells and the head inside
   the connected wells is the same.

2. :class:`~timml.well.HeadWellString` is similar to the
   :class:`~timml.well.WellString`, but instead of specifying the total discharge, the
   head is specified at a specific point :math:`(x_c, y_c)`. 

3. :class:`~timml.well.CollectorWell` is a well type that consists of multiple
   connected LineSinks for which the total discharge is given. The total discharge is
   distributed across the connected LineSinks and the head inside the connected LineSinks
   is the same.`

4. :class:`~timml.well.RadialCollectorWell` is similar to the
   :class:`~timml.well.CollectorWell`, but instead of specifying the start and end points
   for each collector arm, the number of arms is specified, and these are automatically
   distributed in a radial pattern.

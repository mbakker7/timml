Wells
-----

There are two types of wells: wells for which the total discharge is specified and wells for which the head inside the well is specified.
Both types of wells may have an entry resistance (skin effect) defined by the resistance :math:`c` (dimension: time). The discharge :math:`Q_i` in layer :math:`i` is a function of the head
:math:`h_i` in layer :math:`i` just outside the well and the head :math:`h_w` inside the well:

    .. math::
        Q_i = 2\pi r_w(h_i - h_w)/c

1. :class:`~timml.well.Well` is a well for which the total discharge is specified. The total discharge is distributed across the layers in which the
well is screened such that the head inside the well is the same in each screened layer. 

2. :class:`~timml.well.HeadWell` is a well for which the head inside the well is specified. The discharge in each layer is computed such that
the head in all screened layers is equal to the specified head. 
    
.. toctree::
    :maxdepth: 1
    :hidden:
    
    Well <well>
    HeadWell <headwell>
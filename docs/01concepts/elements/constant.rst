Constant
========

Every model needs at least one head-specified condition. In addition, the head may be
specified at one point in one layer, provided this is in an area where the aquifer is
confined. This point is historically called the *reference point*. The *reference*
point is useful to control the behavior in the far field, especially in models where
the far field is not modeled explicitly.

1. :class:`~timml.constant.Constant` is an element that specifies the head at a specific
   location.

2. :class:`~timml.constant.ConstantStar` is a constant element for semi-confined aquifers
   that fixes the head at a reference point. 
Constant
--------
Every model needs at least one head-specified condition. In addition, the head may be specified at one point in one layer,
provided this is in an area where the aquifer is confined. This point is historically called the *reference point*.
The *reference* point is useful to control the behavior in the far field, especially in models where the far field is
not modeled explicitly.

.. autoclass:: timml.constant.Constant
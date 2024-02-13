import inspect  # Used for storing the input

import numpy as np

from .element import Element

__all__ = ["Uflow"]


class Uflow(Element):
    """Add uniform flow to the model.

    Uniform flow may only be added to a model of which the background aquifer system
    is confined.

    Parameters
    ----------
    model : Model object
        model to which the uniform flow is added
    slope : float
        slope of the head (head drop divided by the distance) in the
        direction of flow
    angle : float
        direction of flow in degrees (0 degrees is straight East,
        counter clock-wise is positive)
    label : string or None (default: None)
        label of the element
    """

    def __init__(self, model, slope, angle, label=None):
        assert model.aq.ilap, "TimML Error: Uflow can only be added to model with background confined aquifer"
        self.storeinput(inspect.currentframe())
        Element.__init__(
            self, model, nparam=2, nunknowns=0, layers=0, name="Uflow", label=label
        )
        self.slope = slope
        self.angle = angle
        self.model.add_element(self)

    def __repr__(self):
        return self.name + " with slope and angle: " + str((self.slope, self.angle))

    def initialize(self):
        self.aq = self.model.aq
        self.aq.add_element(self)
        self.Qx = self.slope * np.cos(self.angle * np.pi / 180) * np.sum(self.aq.T)
        self.Qy = self.slope * np.sin(self.angle * np.pi / 180) * np.sum(self.aq.T)
        self.parameters = np.array([[self.Qx], [self.Qy]])

    def potinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, aq.naq))
        if aq == self.aq:
            rv[0, 0] = -x
            rv[1, 0] = -y
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, 2, aq.naq))
        if aq == self.aq:
            rv[0, 0, 0] = 1.0
            rv[1, 1, 0] = 1.0
        return rv

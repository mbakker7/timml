import numpy as np
import inspect # Used for storing the input
from .element import Element

class Uflow(Element):
    def __init__(self, model, slope, angle,\
                 name='Uflow', label=None):
        assert model.aq.ilap, 'TimML Error: Uflow can only be added to model with background confined aquifer'
        self.storeinput(inspect.currentframe())
        Element.__init__(self, model, Nparam=2, Nunknowns=0, layers=0,\
                         name=name, label=label)
        self.slope = slope
        self.angle = angle
        self.model.add_element(self)
    def __repr__(self):
        return self.name + ' with slope and angle: ' + str((self.slope, self.angle))
    def initialize(self):
        self.aq = self.model.aq
        self.aq.add_element(self)
        self.Qx = self.slope * np.cos(self.angle * np.pi / 180) * np.sum(self.aq.T)
        self.Qy = self.slope * np.sin(self.angle * np.pi / 180) * np.sum(self.aq.T)
        self.parameters = np.array([[self.Qx], [self.Qy]])
    def potinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, aq.Naq))
        if aq == self.aq:
            rv[0, 0] = -x
            rv[1, 0] = -y
        return rv
    def disinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, 2, aq.Naq))
        rv[0, 0, 0] = 1.0
        rv[1, 1, 0] = 1.0
        return rv
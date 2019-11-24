import numpy as np
import inspect  # Used for storing the input
from .element import Element
from .equation import PotentialEquation

__all__ = ['Constant', 'ConstantStar']

class ConstantBase(Element, PotentialEquation):
    def __init__(self, model, xr=0, yr=0, hr=0.0, layer=0, \
                 name='ConstantBase', label=None, aq=None):
        self.storeinput(inspect.currentframe())
        Element.__init__(self, model, nparam=1, nunknowns=1, layers=layer, \
                         name=name, label=label)
        self.nparam = 1  # Defined here and not in Element as other elements can have multiple parameters per layers
        self.nunknowns = 0
        self.xr = xr
        self.yr = yr
        self.hr = hr
        self.aq = aq
        self.model.add_element(self)

    def __repr__(self):
        return self.name + ' at ' + str(
            (self.xr, self.yr)) + ' with head  ' + str(self.hr)

    def initialize(self):
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.xr, self.yr)
        self.aq.add_element(self)
        self.ncp = 1
        self.xc = np.array([self.xr])
        self.yc = np.array([self.yr])
        self.pc = self.hr * self.aq.T[self.layers]
        self.parameters = np.atleast_2d(self.pc)

    def potinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((1, aq.naq))
        if aq == self.aq:
            rv[0, 0] = 1
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, 1, aq.naq))
        return rv


class Constant(ConstantBase, PotentialEquation):
    """
    Specify the head at one point in the model in one layer.
    The head may only be specified in an area of the model where
    the aquifer system is confined.
    
    Parameters
    ----------
    model : Model object
        model to which the element is added
    xr : float
        x-coordinate of the point where the head is specified
    yr : float
        y-coordinate of the point where the head is specified
    hr : float
        specified head
    rw : float
        radius of the well
    layer : int
        layer where the head is specified
    label : string or None (default: None)
        label of the element
    
    """
    
    def __init__(self, model, xr=0, yr=0, hr=0.0, layer=0, label=None):
        self.storeinput(inspect.currentframe())
        ConstantBase.__init__(self, model, xr=xr, yr=yr, hr=hr, layer=layer, \
                              name='Constant', label=label)
        self.nunknowns = 1

    def initialize(self):
        ConstantBase.initialize(self)
        assert self.aq.ilap, 'Constant element added to area that is ' \
                             'semi-confined'
        self.resfac = np.zeros(1)  # required for HeadEquation
        self.strengthinf = np.zeros(1)  # required for HeadEquation

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class ConstantInside(Element):
    # Sets constant at points xc, yc equal to the average of the potential of all elements at points xc, yc
    # Used for the inside of an inhomogeneity
    def __init__(self, model, xc=0, yc=0, label=None):
        self.storeinput(inspect.currentframe())
        Element.__init__(self, model, nparam=1, nunknowns=1,
                         layers=list(range(model.aq.naq)), \
                         name='ConstantInside', label=label)
        self.xc = np.atleast_1d(xc)
        self.yc = np.atleast_1d(yc)
        self.parameters = np.zeros((1, 1))
        self.model.add_element(self)

    def __repr__(self):
        return self.name

    def initialize(self):
        self.aq = self.model.aq.find_aquifer_data(self.xc[0], self.yc[0])
        self.aq.add_element(self)
        self.ncp = len(self.xc)

    def potinf(self, x, y, aq=None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((1, aq.naq))
        if aq == self.aq:
            rv[0, 0] = 1
        return rv

    def disvecinf(self, x, y, aq=None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, 1, aq.naq))
        return rv

    def equation(self):
        mat = np.zeros((1, self.model.neq))
        rhs = np.zeros(1)  # Needs to be initialized to zero
        for icp in range(self.ncp):
            ieq = 0
            for e in self.model.elementlist:
                if e. nunknowns > 0:
                    if e != self:
                        mat[0:, ieq:ieq + e. nunknowns] += \
                            e.potinflayers(self.xc[icp], self.yc[icp],
                                           self.layers).sum(0)
                        ieq += e. nunknowns
                        # else:
                        #    mat[0, ieq:ieq+e. nunknowns] += -1
                else:
                    rhs[0] -= \
                        e.potentiallayers(self.xc[icp], self.yc[icp],
                                          self.layers).sum(0)
        return mat, rhs

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class ConstantStar(Element):
    def __init__(self, model, hstar=0.0, label=None, aq=None):
        self.storeinput(inspect.currentframe())
        Element.__init__(self, model, nparam=1, nunknowns=0, layers=0, \
                         name='ConstantStar', label=label)
        assert hstar is not None, 'a value for hstar needs to be specified'
        self.hstar = hstar
        self.aq = aq
        self.model.add_element(self)

    def __repr__(self):
        return self.name + ' with head  ' + str(self.hstar)

    def initialize(self):
        self.aq.add_element(self)
        self.aq.constantstar = self
        self.parameters = np.zeros((1, 1))
        self.potstar = self.hstar * self.aq.T

    def potinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((1, aq.naq))
        return rv

    def potentiallayers(self, x, y, layers, aq=None):
        '''Returns array of size len(layers) only used in building equations
        Defined here as it is the particular solution inside a semi-confined aquifer
        and cannot be added by using eigen vectors'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        pot = np.zeros(len(layers))
        if aq == self.aq:
            pot[:] = self.potstar[layers]
        return pot

    def disvecinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, 1, aq.naq))
        return rv

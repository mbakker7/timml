import numpy as np
import inspect  # Used for storing the input
from .element import Element
from .equation import HeadEquation, MscreenWellEquation, HeadDiffEquation, DisvecDiffEquation

__all__ = ['LineSink1D', 'HeadLineSink1D']

class LineSink1DBase(Element):

    def __init__(self, model, xls, sigls=1, layers=0, name="LineSink1DBase", label=None, 
                 addtomodel=True, res=0, wh=1, aq=None):
        Element.__init__(self, model, nparam=1, nunknowns=0, layers=layers,
                         name=name, label=label)
        self.xls = float(xls)
        self.sigls = np.atleast_1d(sigls)
        self.res = float(res)
        self.wh = wh
        self.aq = aq
        self.addtomodel = addtomodel
        if self.addtomodel: 
            self.model.add_element(self)
        self.nparam = self.nlayers
        self.tiny = 1e-12
        
    def __repr__(self):
        return self.name + " at " + str(self.xls) + " in layers: " + str(self.layers)
    
    def initialize(self):
        self.xc = np.array([self.xls])
        self.yc = np.zeros(1)
        self.ncp = 1
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.xc[0], self.yc[0])

        if self.addtomodel:
            self.aq.add_element(self)
        self.parameters = np.empty((self.nparam, 1))
        self.parameters[:, 0] = self.sigls
        self.theta_norm_out = np.zeros(1)
        self.cosnorm = np.cos(self.theta_norm_out) * np.ones(self.ncp)
        self.sinnorm = np.sin(self.theta_norm_out) * np.ones(self.ncp)
        if self.wh == 'H':
            self.wh = self.aq.Haq[self.layers]
        elif self.wh == '2H':
            self.wh = 2.0 * self.aq.Haq[self.layers]
        elif np.isscalar(self.wh):
            self.wh = self.wh * np.ones(self.nlayers)
        self.resfac = self.aq.T[self.layers] * self.res / self.wh

    def potinf(self, x, y, aq=None):
        if aq is None: 
            aq = self.model.aq.find_aquifer_data(x, 0)
        rv = np.zeros((self.nparam, aq.naq))
        if aq == self.aq:
            pot = np.zeros(aq.naq)
            if aq.ilap:
                if x - self.xls < 0.:
                    pot[0] = -0.5 * (x - self.xls - 1)  # so that pot = 0.5 at x=xls
                    pot[1:] = -0.5 * aq.lab[1:] * np.exp((x - self.xls) / aq.lab[1:])
                elif x - self.xls >= 0.:
                    pot[0] = 0.5 * (x - self.xls + 1)
                    pot[1:] = -0.5 * aq.lab[1:] * np.exp(-(x - self.xls) / aq.lab[1:])
            else:
                if x - self.xls < 0.:
                    pot[:] = -0.5 * aq.lab * np.exp((x - self.xls) / aq.lab)
                elif x - self.xls >= 0.:
                    pot[:] = -0.5 * aq.lab * np.exp(-(x - self.xls) / aq.lab)
            rv[:] = self.aq.coef[self.layers] * pot
        return rv
    
    def disvecinf(self, x, y, aq=None):
        if aq is None: 
            aq = self.model.aq.find_aquifer_data(x, 0)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq == self.aq:
            qx = np.zeros(aq.naq)
            if aq.ilap:
                if x - self.xls < 0.:
                    qx[0] = 0.5
                    qx[1:] = 0.5 * np.exp((x - self.xls) / aq.lab[1:])
                elif x - self.xls >= 0.:
                    qx[0] = -0.5
                    qx[1:] = -0.5 * np.exp(-(x - self.xls) / aq.lab[1:])
            else:
                if x - self.xls < 0.:
                    qx[:] = 0.5 * np.exp((x - self.xls) / aq.lab)
                elif x - self.xls >= 0.:
                    qx[:] = -0.5 * np.exp(-(x - self.xls) / aq.lab)
            rv[0] = self.aq.coef[self.layers] * qx
        return rv
    
    def discharge(self):
        """Discharge per unit length"""
        Q = np.zeros(self.aq.naq)
        Q[self.layers] = self.parameters[:, 0]
        return Q
    
class LineSink1D(LineSink1DBase, MscreenWellEquation):
    """
    Create an infinitely long line-sink with a given discharge
    per unit length. In case the line-sink is screened in multiple
    layers, the discharge is distributed over the layers such that
    the head inside the line-sink is constant in all layers.
    
    Parameters
    ----------
    
    model : Model object
        Model to which the element is added
    xls : scalar
        x-location of line-sink
    sigls : scalar
        discharge per unit length
    res : scalar (default is 0)
        resistance of line-sink
    wh : scalar or str
        distance over which water enters line-sink
        if 'H': (default) distance is equal to the thickness of the aquifer layer (when flow comes mainly from one side)
        if '2H': distance is twice the thickness of the aquifer layer (when flow comes from both sides)
        if scalar: the width of the stream that partially penetrates the aquifer layer
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers 
    label: str or None
        label of element

    """
    
    def __init__(self, model, xls=0, sigls=1, \
                 layers=0, label=None):
        self.storeinput(inspect.currentframe())
        LineSink1DBase.__init__(self, model, xls, sigls=0, layers=layers, \
                               name="Linesink1D", label=label, \
                               addtomodel=True, res=0, wh=1, aq=None)
        self.Qc = float(sigls)
        if self.nlayers == 1:
            self.nunknowns = 0
        else:
            self.nunknowns = self.nparam

    def initialize(self):
        LineSink1DBase.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol
    
class HeadLineSink1D(LineSink1DBase, HeadEquation):
    """
    Create an infinitely long line-sink with a given head.
    
    Parameters
    ----------
    
    model : Model object
        Model to which the element is added
    xls : scalar
        x-location of line-sink
    hls : scalar
        head in line-sink
    res : scalar (default is 0)
        resistance of line-sink
    wh : scalar or str
        distance over which water enters line-sink
        if 'H': (default) distance is equal to the thickness of the aquifer layer (when flow comes mainly from one side)
        if '2H': distance is twice the thickness of the aquifer layer (when flow comes from both sides)
        if scalar: the width of the stream that partially penetrates the aquifer layer
    layers : scalar, list or array
        layer(s) in which element is placed
        if scalar: element is placed in this layer
        if list or array: element is placed in all these layers 
    label: str or None
        label of element

    """
 
    def __init__(self, model, xls=0, hls=1, \
                 res=0, wh=1, layers=0, label=None):
        self.storeinput(inspect.currentframe())
        LineSink1DBase.__init__(self, model, xls, sigls=0, layers=layers, \
                               name="HeadLinesink1D", label=label, \
                               addtomodel=True, res=res, wh=wh, aq=None)
        self.hc = np.atleast_1d(hls) * np.ones(self.nlayers)
        self.nunknowns = self.nparam

    def initialize(self):
        LineSink1DBase.initialize(self)
        self.pc = self.hc * self.aq.T[self.layers]  # Needed in solving

    def setparams(self, sol):
        self.parameters[:, 0] = sol
    
class HeadDiffLineSink1D(LineSink1DBase, HeadDiffEquation):
    """HeadDiffLineSink1D for left side (xcout)
    """
    
    def __init__(self, model, xls, label=None, 
                 aq=None, aqin=None, aqout=None):
        LineSink1DBase.__init__(self, model, xls, sigls=0,
                                layers=np.arange(model.aq.naq), label=label, 
                                name='HeadDiffLineSink1D', addtomodel=True, aq=aq)
        self.inhomelement = True
        self.nunknowns = self.nparam
        self.aqin = aqin
        self.aqout = aqout
    
    def initialize(self):
        LineSink1DBase.initialize(self)
        self.xcout = self.xc + self.tiny * abs(self.xc) + self.tiny
        self.xcin = self.xc - self.tiny * abs(self.xc) - self.tiny
        self.ycout = np.zeros(1)
        self.ycin = np.zeros(1)
        if self.aqout is None:
            self.aqout = self.model.aq.find_aquifer_data(self.xcout[0], 0)
        if self.aqin is None:
            self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], 0)
            
    def setparams(self, sol):
        self.parameters[:, 0] = sol
        
class FluxDiffLineSink1D(LineSink1DBase, DisvecDiffEquation):
    """HeadDiffLineSink1D for left side (xcout)
    """
    
    def __init__(self, model, xls, label=None, 
                 aq=None, aqin=None, aqout=None):
        LineSink1DBase.__init__(self, model, xls, sigls=0,
                                layers=np.arange(model.aq.naq), label=label, 
                                name='FluxDiffLineSink1D', addtomodel=True, aq=aq)
        self.inhomelement = True
        self.nunknowns = self.nparam
        self.aqin = aqin
        self.aqout = aqout
    
    def initialize(self):
        LineSink1DBase.initialize(self)
        self.xcout = self.xc - self.tiny * abs(self.xc) - self.tiny
        self.xcin = self.xc + self.tiny * abs(self.xc) + self.tiny
        self.ycout = np.zeros(1)
        self.ycin = np.zeros(1)
        if self.aqout is None:
            self.aqout = self.model.aq.find_aquifer_data(self.xcout[0], 0)
        if self.aqin is None:
            self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], 0)
            
    def setparams(self, sol):
        self.parameters[:, 0] = sol

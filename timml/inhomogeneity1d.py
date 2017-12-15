import numpy as np
import inspect  # user for storing the input
from .aquifer import AquiferData
from .aquifer_parameters import param_maq, param_3d
from .constant import ConstantStar
from .linesink1d import HeadDiffLineSink1D, FluxDiffLineSink1D

__all__ = ['StripInhomMaq', 'StripInhom3D']

class StripInhom(AquiferData):
    tiny = 1e-12
    
    def __init__(self, model, x1, x2, kaq, c, z, npor, ltype, hstar):
        AquiferData.__init__(self, model, kaq, c, z, npor, ltype)
        self.hstar = hstar
        self.x1 = x1
        self.x2 = x2
        self.inhom_number = self.model.aq.add_inhom(self)
        self.addlinesinks = True  # Set to False not to add line-sinks
        
    def __repr__(self):
        return "Inhom1D: " + str(list([self.x1, self.x2]))
        
    def isinside(self, x, y):
        return (x >= self.x1) and (x < self.x2)

    def create_elements(self):
        # HeadDiff on right side, FluxDiff on left side
        if self.x1 == -np.inf:
            xin = self.x2 - self.tiny * abs(self.x2) - self.tiny
            xoutright = self.x2 + self.tiny * abs(self.x2) + self.tiny
            aqin = self.model.aq.find_aquifer_data(xin, 0)
            aqoutright = self.model.aq.find_aquifer_data(xoutright, 0)
            if self.addlinesinks:
                HeadDiffLineSink1D(self.model, self.x2, label=None, 
                                   aq=aqin, aqin=aqin, aqout=aqoutright)
        elif self.x2 == np.inf:
            xin = self.x1 + self.tiny * abs(self.x1) + self.tiny
            xoutleft = self.x1 - self.tiny * abs(self.x1) - self.tiny
            aqin = self.model.aq.find_aquifer_data(xin, 0)
            aqoutleft = self.model.aq.find_aquifer_data(xoutleft, 0)
            if self.addlinesinks:   
                FluxDiffLineSink1D(self.model, self.x1, label=None, 
                                   aq=aqin, aqin=aqin, aqout=aqoutleft)
        else:
            xin = 0.5 * (self.x1 + self.x2)
            xoutleft = self.x1 - self.tiny * abs(self.x1) - self.tiny
            xoutright = self.x2 + self.tiny * abs(self.x2) + self.tiny
            aqin = self.model.aq.find_aquifer_data(xin, 0)
            aqleft = self.model.aq.find_aquifer_data(xoutleft, 0)
            aqright = self.model.aq.find_aquifer_data(xoutright, 0)
            if self.addlinesinks:
                HeadDiffLineSink1D(self.model, self.x2, label=None, 
                                   aq=aqin, aqin=aqin, aqout=aqright)
                FluxDiffLineSink1D(self.model, self.x1, label=None, 
                                   aq=aqin, aqin=aqin, aqout=aqleft)
        if aqin.ltype[0] == 'l':
            assert self.hstar is not None, "Error: hstar needs to be set"
            c = ConstantStar(self.model, self.hstar, aq=aqin)
            c.inhomelement = True
        

class StripInhomMaq(StripInhom):
    
    def __init__(self, model, x1, x2, kaq=1, z=[1, 0], c=[], npor=0.3, 
                 topboundary='conf', hstar=None):
        self.storeinput(inspect.currentframe())
        kaq, c, npor, ltype, = param_maq(kaq, z, c, npor, topboundary)
        StripInhom.__init__(self, model, x1, x2, kaq, c, z, npor, ltype, hstar)


class StripInhom3D(StripInhom):
    
    def __init__(self, model, x1, x2, kaq, z=[1, 0], kzoverkh=1, npor=0.3, 
                 topboundary='conf', hstar=None, topres=None, topthick=0.):
        self.storeinput(inspect.currentframe())
        kaq, c, npor, ltype, = param_3d(kaq, z, kzoverkh, npor, topboundary, topres)
        if topbndary== 'semi':
            z = np.hstack((z[0] + topthick, z))
        StripInhom.__init__(self, model, x1, x2, kaq, c, z, npor, ltype, hstar)

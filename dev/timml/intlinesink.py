import numpy as np
from linesink import LineSinkHoBase
from equation import HeadDiffEquation2, DisvecDiffEquation2
from controlpoints import controlpoints

class IntHeadDiffLineSink(LineSinkHoBase, HeadDiffEquation2):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 order=0, ndeg=3, layers=0, label=None, addtomodel=True, aq=None, aqin=None, aqout=None):
        LineSinkHoBase.__init__(self, model, x1, y1, x2, y2, Qls=0, \
                 layers=range(model.aq.Naq), order=order, addconstant=False, \
                 name='IntHeadDiffLineSink', label=label, \
                 addtomodel=addtomodel, aq=aq)
        self.ndeg = ndeg
        self.Nunknowns = self.Nparam
        self.aqin = aqin
        self.aqout = aqout
    def initialize(self):
        LineSinkHoBase.initialize(self)
        self.xcin, self.ycin = controlpoints(self.Ncp-1, self.z1, self.z2, eps=1e-6, include_ends=True)
        self.xcout, self.ycout = controlpoints(self.Ncp-1, self.z1, self.z2, eps=-1e-6, include_ends=True)
        if self.aqin is None:
            self.aqin = self.model.aq.find_aquifer_data(self.xcin[1], self.ycin[1])
        if self.aqout is None:
            self.aqout = self.model.aq.find_aquifer_data(self.xcout[1], self.ycout[1])
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class IntFluxDiffLineSink(LineSinkHoBase, DisvecDiffEquation2):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 order=0, ndeg=3, label=None, addtomodel=True, aq=None, aqin=None, aqout=None):
        LineSinkHoBase.__init__(self, model, x1, y1, x2, y2, Qls=0, \
                 layers=range(model.aq.Naq), order= order, name='IntFluxLineSink', label=label, \
                 addtomodel=addtomodel, aq=aq)
        self.ndeg = ndeg
        self.Nunknowns = self.Nparam
        self.aqin = aqin
        self.aqout = aqout
    def initialize(self):
        LineSinkHoBase.initialize(self)
        self.xcin, self.ycin = controlpoints(self.Ncp-1, self.z1, self.z2, eps=1e-6, include_ends=True)
        self.xcout, self.ycout = controlpoints(self.Ncp-1, self.z1, self.z2, eps=-1e-6, include_ends=True)
        if self.aqin is None:
            self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], self.ycin[0])
        if self.aqout is None:
            self.aqout = self.model.aq.find_aquifer_data(self.xcout[0], self.ycout[0])
    def setparams(self, sol):
        self.parameters[:,0] = sol
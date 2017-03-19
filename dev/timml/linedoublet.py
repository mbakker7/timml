import numpy as np
import inspect  # Used for storing the input
from element import Element
from besselaesnew import besselaesnew
besselaesnew.initialize()
from controlpoints import controlpoints
from equation import DisvecEquation, LeakyWallEquation

class LineDoubletHoBase(Element):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, delp=0.0, res=0.0, \
                 layers=0, order=0, name='LineDoubletHoBase', \
                 label=None, addtomodel=True, aq=None, zcinout=None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers, \
                         name=name, label=label)
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.x2 = float(x2)
        self.y2 = float(y2)
        self.delp = np.atleast_1d(delp).astype('d')
        self.res = np.atleast_1d(res).astype('d')
        self.order = order
        self.Nparam = self.Nlayers * (self.order + 1)
        self.addtomodel = addtomodel
        if addtomodel: self.model.add_element(self)
        self.aq = aq
        self.zcinout = zcinout

    def __repr__(self):
        return self.name + ' from ' + str((self.x1, self.y1)) + ' to ' + str(
            (self.x2, self.y2))

    def initialize(self):
        self.Ncp = self.order + 1
        self.z1 = self.x1 + 1j * self.y1
        self.z2 = self.x2 + 1j * self.y2
        self.L = np.abs(self.z1 - self.z2)
        self.thetaNormOut = np.arctan2(self.y2 - self.y1,
                                       self.x2 - self.x1) - np.pi / 2.0
        self.cosnorm = np.cos(self.thetaNormOut) * np.ones(self.Ncp)
        self.sinnorm = np.sin(self.thetaNormOut) * np.ones(self.Ncp)
        #
        self.xc, self.yc = controlpoints(self.Ncp, self.z1, self.z2, eps=0)
        if self.zcinout is not None:
            self.xcin, self.ycin = controlpoints(self.Ncp, self.zcinout[0],
                                                 self.zcinout[1], eps=0)
            self.xcout, self.ycout = controlpoints(self.Ncp, self.zcinout[2],
                                                   self.zcinout[3], eps=0)
        else:
            self.xcin, self.ycin = controlpoints(self.Ncp, self.z1, self.z2,
                                                 eps=1e-6)
            self.xcout, self.ycout = controlpoints(self.Ncp, self.z1, self.z2,
                                                   eps=-1e-6)
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        self.resfac = self.aq.T[self.pylayers] / self.res
        if self.addtomodel:
            self.aq.add_element(self)
        self.parameters = np.empty((self.Nparam, 1))
        # Not sure if this needs to be here
        self.parameters[:, 0] = self.delp

    def potinf(self, x, y, aq=None):
        '''Can be called with only one x,y value
        Returns array(Nparam, self.aq.Naq) with order
        order 0, layer[0]
        order 0, layer[1]
        ...
        order 1, layer[0]
        order 1, layer[1]
        etc
        '''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nparam, aq.Naq))
        if aq == self.aq:
            potrv = rv.reshape((self.order + 1, self.Nlayers,
                                aq.Naq))  # clever way of using a reshaped rv here
            pot = np.zeros((self.order + 1, aq.Naq))
            pot[:, :] = besselaesnew.potbesldv(x, y, self.z1, self.z2, aq.lab,
                                               self.order, aq.ilap)
            potrv[:] = self.aq.coef[self.pylayers] * pot[:, np.newaxis, :]
        return rv

    def disinf(self, x, y, aq=None):
        '''Can be called with only one x,y value
        Returns array(Nparam, self.aq.Naq) with order
        order 0, layer[0]
        order 0, layer[1]
        ...
        order 1, layer[0]
        order 1, layer[1]
        etc
        '''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.Nparam, aq.Naq))
        if aq == self.aq:
            qxqyrv = rv.reshape((2, self.order + 1, self.Nlayers,
                                 aq.Naq))  # clever way of using a reshaped rv here
            qxqy = np.zeros((2 * (self.order + 1), aq.Naq))
            qxqy[:, :] = besselaesnew.disbesldv(x, y, self.z1, self.z2, aq.lab,
                                                self.order, aq.ilap)
            qxqyrv[0, :] = self.aq.coef[self.pylayers] * qxqy[:self.order + 1,
                                                         np.newaxis, :]
            qxqyrv[1, :] = self.aq.coef[self.pylayers] * qxqy[self.order + 1:,
                                                         np.newaxis, :]
        return rv

class ImpLineDoublet(LineDoubletHoBase, DisvecEquation):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 order=0, layers=0, label=None, addtomodel=True):
        self.storeinput(inspect.currentframe())
        LineDoubletHoBase.__init__(self, model, x1, y1, x2, y2, delp=0, \
                                layers=layers, order=order,
                                name='ImpLineDoublet', label=label, \
                                addtomodel=addtomodel)
        self.Nunknowns = self.Nparam

    def initialize(self):
        LineDoubletHoBase.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol
        
class LeakyLineDoublet(LineDoubletHoBase, LeakyWallEquation):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, res=0,\
                 order=0, layers=0, label=None, addtomodel=True):
        self.storeinput(inspect.currentframe())
        LineDoubletHoBase.__init__(self, model, x1, y1, x2, y2, delp=0, \
                                res=res, layers=layers, order=order,
                                name='ImpLineDoublet', label=label, \
                                addtomodel=addtomodel)
        self.Nunknowns = self.Nparam

    def initialize(self):
        LineDoubletHoBase.initialize(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class LineDoubletStringBase(Element):
    def __init__(self, model, xy, closed=False, layers=0, order=0, res=0,
                 name='LineDoubletStringBase', label=None, aq=None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers, \
                         name=name, label=label)
        self.xy = np.atleast_2d(xy).astype('d')
        if closed: self.xy = np.vstack((self.xy, self.xy[0]))
        self.order = order
        self.aq = aq
        self.ldlist = []
        self.x, self.y = self.xy[:, 0], self.xy[:, 1]
        self.Nld = len(self.x) - 1
        for i in range(self.Nld):
            self.ldlist.append(
                LineDoubletHoBase(model, x1=self.x[i], y1=self.y[i], x2=self.x[i + 1],
                                  y2=self.y[i + 1], delp=0.0, res=res, layers=layers,
                                  order=order, label=label, addtomodel=False, aq=aq))

    def __repr__(self):
        return self.name + ' with nodes ' + str(self.xy)

    def initialize(self):
        for ld in self.ldlist:
            ld.initialize()
        self.Ncp = self.Nld * self.ldlist[
            0].Ncp  # Same order for all elements in string
        self.Nparam = self.Nld * self.ldlist[0].Nparam
        self.Nunknowns = self.Nparam
        self.xld = np.empty((self.Nld, 2))
        self.yld = np.empty((self.Nld, 2))
        for i, ld in enumerate(self.ldlist):
            self.xld[i, :] = [ld.x1, ld.x2]
            self.yld[i, :] = [ld.y1, ld.y2]
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.ldlist[0].xc,
                                                      self.ldlist[0].yc)
        self.parameters = np.zeros((self.Nparam, 1))
        ## As parameters are only stored for the element not the list, we need to combine the following
        self.xc = np.array([ld.xc for ld in self.ldlist]).flatten()
        self.yc = np.array([ld.yc for ld in self.ldlist]).flatten()
        self.xcin = np.array([ld.xcin for ld in self.ldlist]).flatten()
        self.ycin = np.array([ld.ycin for ld in self.ldlist]).flatten()
        self.xcout = np.array([ld.xcout for ld in self.ldlist]).flatten()
        self.ycout = np.array([ld.ycout for ld in self.ldlist]).flatten()
        self.cosnorm = np.array([ld.cosnorm for ld in self.ldlist]).flatten()
        self.sinnorm = np.array([ld.sinnorm for ld in self.ldlist]).flatten()
        self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], self.ycin[0])
        self.aqout = self.model.aq.find_aquifer_data(self.xcout[0],
                                                     self.ycout[0])
        self.resfac = self.ldlist[0].resfac

    def potinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nld, self.ldlist[0].Nparam, aq.Naq))
        for i in range(self.Nld):
            rv[i] = self.ldlist[i].potinf(x, y, aq)
        rv.shape = (self.Nparam, aq.Naq)
        return rv

    def disinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.Nld, self.ldlist[0].Nparam, aq.Naq))
        for i in range(self.Nld):
            rv[:, i] = self.ldlist[i].disinf(x, y, aq)
        rv.shape = (2, self.Nparam, aq.Naq)
        return rv


class ImpLineDoubletString(LineDoubletStringBase, DisvecEquation):
    def __init__(self, model, xy=[(-1, 0), (1, 0)], delp=0.0, \
                 layers=0, order=0, label=None):
        self.storeinput(inspect.currentframe())
        LineDoubletStringBase.__init__(self, model, xy, closed=False,
                                    layers=layers, order=order, \
                                    name='ImpLineDoubletString', label=label,
                                    aq=None)
        self.model.add_element(self)

    def initialize(self):
        LineDoubletStringBase.initialize(self)
        self.aq.add_element(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol
        
class LeakyLineDoubletString(LineDoubletStringBase, LeakyWallEquation):
    def __init__(self, model, xy=[(-1, 0), (1, 0)], res=np.inf,\
                 layers=0, order=0, label=None):
        self.storeinput(inspect.currentframe())
        LineDoubletStringBase.__init__(self, model, xy, closed=False,
                                    layers=layers, order=order, res=res,\
                                    name='ImpLineDoubletString', label=label,
                                    aq=None)
        self.model.add_element(self)

    def initialize(self):
        LineDoubletStringBase.initialize(self)
        self.aq.add_element(self)

    def setparams(self, sol):
        self.parameters[:, 0] = sol

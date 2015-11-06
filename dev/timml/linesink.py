import numpy as np
import inspect # Used for storing the input
from element import Element
from equation import HeadEquation
from besselaes import potbeslsho, potbesonlylsho, disbeslsho, disbesonlylsho
from controlpoints import controlpoints

class LineSinkBase(Element):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, Qls=100.0, \
                 res=0, wh=1, layers=0, name='LineSinkBase', label=None, \
                 addtomodel=True):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers,\
                         name=name, label=label)
        self.Nparam = len(self.pylayers)
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.x2 = float(x2)
        self.y2 = float(y2)
        self.Qls = np.atleast_1d(Qls)
        self.res = float(res)
        self.wh = wh
        self.addtomodel = addtomodel
        if self.addtomodel: self.model.add_element(self)
        #self.xa,self.ya,self.xb,self.yb,self.np = np.zeros(1),np.zeros(1),np.zeros(1),np.zeros(1),np.zeros(1,'i')  # needed to call bessel.circle_line_intersection
    def __repr__(self):
        return self.name + ' from ' + str((self.x1, self.y1)) +' to '+str((self.x2, self.y2))
    def initialize(self):
        self.xc = np.array([0.5 * (self.x1 + self.x2)])
        self.yc = np.array([0.5 * (self.y1 + self.y2)])
        self.Ncp = 1
        self.z1 = self.x1 + 1j * self.y1
        self.z2 = self.x2 + 1j * self.y2
        self.L = np.abs(self.z1 - self.z2)
        self.order = 0 # This is for univform strength only
        self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        if self.addtomodel: self.aq.add_element(self)
        self.parameters = np.empty((self.Nparam, 1))
        self.parameters[:,0] = self.Qls / self.L
        if self.wh == 'H':
            self.wh = self.aq.Haq[self.pylayers]
        elif self.wh == '2H':
            self.wh = 2.0 * self.aq.Haq[self.pylayers]
        elif np.isscalar(self.wh):
            self.wh = self.wh * np.ones(self.Nlayers)
        self.resfac = self.aq.T[self.pylayers] * self.res / self.wh
    def potinf(self, x, y, aq = None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nparam, aq.Naq))
        if aq == self.aq:
            pot = np.zeros(aq.Naq)
            if aq.ltype[0] == 'a':
                potbeslsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                           aq.Naq, aq.zeropluslab, 0, pot)  # Call FORTRAN extension
            else:
                potbesonlylsho(x, y, self.x1, selr.y1, self.x2, self.y2, \
                               aq.Naq, aq.lab, 0, pot)
            rv[:] = self.aq.coef[self.pylayers] * pot
        return rv
    def disinf(self, x, y, aq = None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.Nparam, aq.Naq))
        if aq == self.aq:
            qx = np.zeros(aq.Naq)
            qy = np.zeros(aq.Naq)
            if aq.ltype[0] == 'a':
                disbeslsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                           aq.Naq, aq.zeropluslab, 0, qx, qy)
            else:
                disbesonlylsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                               aq.Naq, aq.lab, 0, qx, qy)
            rv[0] = self.aq.coef[self.pylayers] * qx
            rv[1] = self.aq.coef[self.pylayers] * qy
        return rv
    
class HeadLineSink(LineSinkBase, HeadEquation):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, hls=1.0, \
                 res=0, wh=1, layers=0, label=None, addtomodel=True):
        self.storeinput(inspect.currentframe())
        LineSinkBase.__init__(self, model, x1, y1, x2, y2, Qls=0,\
                 res=res, wh=wh, layers=layers, name='HeadLineSink', label=label, \
                 addtomodel=addtomodel)
        self.hc = hls
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineSinkBase.initialize(self)
        self.pc = self.hc * self.aq.T[self.pylayers] # Needed in solving
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class LineSinkHoBase(Element):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 Qls= 0.0, layers=0, order=0, addconstant=False, name='LineSinkHoBase', \
                 label=None, addtomodel=True, aq=None, zcinout=None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers,\
                         name=name, label=label)
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.x2 = float(x2)
        self.y2 = float(y2)
        self.Qls = np.atleast_1d(Qls)
        self.order = order
        self.addconstant = addconstant
        self.Nparam = self.Nlayers * (self.order + 1)
        if self.addconstant: self.Nparam += 1
        self.addtomodel = addtomodel
        if addtomodel: self.model.add_element(self)
        self.aq = aq
        self.zcinout = zcinout
    def __repr__(self):
        return self.name + ' from ' + str((self.x1, self.y1)) +' to '+str((self.x2, self.y2))
    def initialize(self):
        self.Ncp = self.order + 1
        if self.addconstant: self.Ncp += 1
        self.z1 = self.x1 + 1j * self.y1
        self.z2 = self.x2 + 1j * self.y2
        self.L = np.abs(self.z1 - self.z2)
        self.thetaNormOut = np.arctan2(self.y2 - self.y1, self.x2 - self.x1) - np.pi/2.0
        self.cosnorm = np.cos(self.thetaNormOut) * np.ones(self.Ncp)
        self.sinnorm = np.sin(self.thetaNormOut) * np.ones(self.Ncp)
        #
        self.xc, self.yc = controlpoints(self.Ncp, self.z1, self.z2, eps=0)
        if self.zcinout is not None:
            self.xcin, self.ycin = controlpoints(self.Ncp, self.zcinout[0], self.zcinout[1], eps=0)
            self.xcout, self.ycout = controlpoints(self.Ncp, self.zcinout[2], self.zcinout[3], eps=0)
        else:
            self.xcin, self.ycin = controlpoints(self.Ncp, self.z1, self.z2, eps=1e-6)
            self.xcout, self.ycout = controlpoints(self.Ncp, self.z1, self.z2, eps=-1e-6)
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        if self.addtomodel:
            self.aq.add_element(self)
        self.parameters = np.empty((self.Nparam, 1))
        self.parameters[:,0] = self.Qls / self.L  # Not sure if that needs to be here
    def potinf(self, x, y, aq = None):
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
            if self.addconstant:
                potrv = rv[:-1,:].reshape((self.order+1, self.Nlayers, aq.Naq))
            else:
                potrv = rv.reshape((self.order+1, self.Nlayers, aq.Naq))
            pot = np.zeros(aq.Naq)
            if aq.ltype[0] == 'a':
                for i in range(self.order+1): # This should be done inside FORTRAN extension
                    potbeslsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                               aq.Naq, aq.zeropluslab, i, pot)  # Call FORTRAN extension
                    potrv[i] = self.aq.coef[self.pylayers] * pot
                if self.addconstant:
                    rv[-1,0] = 1.0
            else:
                for i in range(self.order+1):
                    potbesonlylsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                                   aq.Naq + 1, aq.lab, i, pot)  # Need to specify aq.Naq + 1 (bug in besselaes)
                    potrv[i] = self.aq.coef[self.pylayers] * pot
        return rv
    def disinf(self, x, y, aq = None):
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
            rv.shape = (2, self.order+1, self.Nlayers, aq.Naq)
            qx = np.zeros(aq.Naq)
            qy = np.zeros(aq.Naq)
            if aq.ltype[0] == 'a':
                for i in range(self.order+1): # This should be done inside FORTRAN extension
                    disbeslsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                               aq.Naq, aq.zeropluslab, i, qx, qy)  # Call FORTRAN extension
                    rv[0, i] = self.aq.coef[self.pylayers] * qx
                    rv[1, i] = self.aq.coef[self.pylayers] * qy
            else:
                for i in range(self.order+1):
                    disbesonlylsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                                   aq.Naq + 1, aq.lab, i, qx, qy)  # Need to specify aq.Naq + 1 (bug in besselaes)
                    rv[0, i] = self.aq.coef[self.pylayers] * qx
                    rv[1, i] = self.aq.coef[self.pylayers] * qy
            rv.shape = (2, self.Nparam, aq.Naq)
        if self.addconstant:  # This is ugly here. Maybe I can find a better way
            rvnew = np.zeros((2, self.Nparam, aq.Naq))
            rv = np.vstack((rv, pot))
        return rv
    def disinfold(self, x, y, aq = None):
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
        rvx = np.zeros((self.Nparam, aq.Naq))
        rvy = np.zeros((self.Nparam, aq.Naq))
        if aq == self.aq:
            rvx.shape = (self.order+1, self.Nlayers, aq.Naq)
            rvy.shape = (self.order+1, self.Nlayers, aq.Naq)
            qx = np.zeros(aq.Naq)
            qy = np.zeros(aq.Naq)
            if aq.ltype[0] == 'a':
                for i in range(self.order+1): # This should be done inside FORTRAN extension
                    disbeslsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                               aq.Naq, aq.zeropluslab, i, qx, qy)  # Call FORTRAN extension
                    rvx[i] = self.aq.coef[self.pylayers] * qx
                    rvy[i] = self.aq.coef[self.pylayers] * qy
            else:
                for i in range(self.order+1):
                    disbesonlylsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                                   aq.Naq, aq.lab, i, qx, qy)
                    rvx[i] = self.aq.coef[self.pylayers] * qx
                    rvy[i] = self.aq.coef[self.pylayers] * qy
            rvx.shape = (self.Nparam, aq.Naq)
            rvy.shape = (self.Nparam, aq.Naq)
        return rvx, rvy
    
class HeadLineSinkHo(LineSinkHoBase, HeadEquation):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 hls=1.0, order=0, layers=0, label=None, addtomodel=True):
        self.storeinput(inspect.currentframe())
        LineSinkHoBase.__init__(self, model, x1, y1, x2, y2, Qls=0, \
                 layers=layers, order=order, name='HeadLineSinkHo', label=label, \
                 addtomodel=addtomodel)
        self.hc = hls
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineSinkHoBase.initialize(self)
        self.pc = self.hc * self.aq.T[self.pylayers] # Needed in solving
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class LineSinkStringBase(Element):
    def __init__(self, model, xy, closed=False, layers=0, order=0, name='LineSinkStringBase', label=None, aq=None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers, \
                         name=name, label=label)
        self.xy = np.atleast_2d(xy).astype('d')
        if closed: self.xy = np.vstack((self.xy, self.xy[0]))
        self.order = order
        self.aq = aq
        self.lslist = []
        self.x, self.y = self.xy[:,0], self.xy[:,1]
        self.Nls = len(self.x) - 1
        for i in range(self.Nls):
            self.lslist.append(LineSinkHoBase(model, \
                x1=self.x[i], y1=self.y[i], x2=self.x[i+1], y2=self.y[i+1], \
                Qls=0.0, layers=layers, order=order, label=label, addtomodel=False, aq=aq))
    def __repr__(self):
        return self.name + ' with nodes ' + str(self.xy)
    def initialize(self):
        for ls in self.lslist:
            ls.initialize()
        self.Ncp = self.Nls * self.lslist[0].Ncp  # Same order for all elements in string
        self.Nparam = self.Nls * self.lslist[0].Nparam
        self.Nunknowns = self.Nparam
        self.xls = np.empty((self.Nls,2))
        self.yls = np.empty((self.Nls,2))
        for i,ls in enumerate(self.lslist):
            self.xls[i,:] = [ls.x1,ls.x2]
            self.yls[i,:] = [ls.y1,ls.y2]
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.lslist[0].xc, self.lslist[0].yc)
        self.parameters = np.zeros((self.Nparam,1))
        ## As parameters are only stored for the element not the list, we need to combine the following
        self.xc = np.array([ls.xc for ls in self.lslist]).flatten()
        self.yc = np.array([ls.yc for ls in self.lslist]).flatten()
        self.xcin = np.array([ls.xcin for ls in self.lslist]).flatten()
        self.ycin = np.array([ls.ycin for ls in self.lslist]).flatten()
        self.xcout = np.array([ls.xcout for ls in self.lslist]).flatten()
        self.ycout = np.array([ls.ycout for ls in self.lslist]).flatten()
        self.cosnorm = np.array([ls.cosnorm for ls in self.lslist]).flatten()
        self.sinnorm = np.array([ls.sinnorm for ls in self.lslist]).flatten()
        self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], self.ycin[0])
        self.aqout = self.model.aq.find_aquifer_data(self.xcout[0], self.ycout[0])
    def potinf(self,x,y,aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data( x, y )
        rv = np.zeros((self.Nls, self.lslist[0].Nparam, aq.Naq))
        for i in range(self.Nls):
            rv[i] = self.lslist[i].potinf(x, y, aq)
        rv.shape = (self.Nparam, aq.Naq)
        return rv
    def disinf(self,x,y,aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data( x, y )
        rv = np.zeros((2, self.Nls, self.lslist[0].Nparam, aq.Naq))
        for i in range(self.Nls):
            rv[:,i] = self.lslist[i].disinf(x, y, aq)
        rv.shape = (2, self.Nparam, aq.Naq)
        return rv

class HeadLineSinkString(LineSinkStringBase, HeadEquation):
    def __init__(self, model, xy=[(-1,0), (1,0)], hls=0.0, \
                 layers=0, order=0, label=None):
        self.storeinput(inspect.currentframe())
        LineSinkStringBase.__init__(self, model, xy, closed=False, layers=layers, order=order, \
                                    name='HeadLineSinkString', label=label, aq=None)
        self.hls = hls
        self.model.add_element(self)
    def initialize(self):
        LineSinkStringBase.initialize(self)
        self.aq.add_element(self)
        self.pc = self.hls * self.aq.T[self.pylayers]
    def setparams(self, sol):
        self.parameters[:,0] = sol
import numpy as np
import matplotlib.pyplot as plt
import inspect  # Used for storing the input
from .element import Element
from .equation import HeadEquation
#from .besselaes import potbeslsho, potbesonlylsho, disbeslsho, disbesonlylsho
from .besselaesnew import *
besselaesnew.initialize()
from .controlpoints import controlpoints

class LineSinkChangeTrace:
    def changetrace(self, xyzt1, xyzt2, aq, layer, ltype, modellayer, direction, verbose=False):
        changed = False
        terminate = False
        xyztnew = 0
        if (ltype == 'a'):
            if (layer == self.pylayers).any():  # in layer where line-sink is screened
                if verbose:
                    print('hello changetrace')
                    print('xyz1:', xyzt1[:-1])
                    print('xyz2:', xyzt2[:-1])
                x1, y1, z1, t1 = xyzt1
                x2, y2, z2, t2 = xyzt2
                eps = 1e-8
                za = x1 + y1 * 1j
                zb = x2 + y2 * 1j
                Za = (2 * za - (self.z1 + self.z2)) / (self.z2 - self.z1)
                Zb = (2 * zb - (self.z1 + self.z2)) / (self.z2 - self.z1)
                if Za.imag * Zb.imag < 0:
                    Xa, Ya = Za.real, Za.imag
                    Xb, Yb = Zb.real, Zb.imag
                    X = Xa - Ya * (Xb - Xa) / (Yb - Ya)
                    if verbose: print('X', X)
                    if abs(X) <= 1:  # crosses line-sink
                        if verbose: print('crosses line-sink')
                        Znew1 = X - eps * np.sign(Yb) * 1j  # steps to side of Ya
                        Znew2 = X + eps * np.sign(Yb) * 1j  # steps to side of Yb
                        znew1 = 0.5 * ((self.z2 - self.z1) * Znew1 + self.z1 + self.z2)
                        znew2 = 0.5 * ((self.z2 - self.z1) * Znew2 + self.z1 + self.z2)
                        xnew1, ynew1 = znew1.real, znew1.imag
                        xnew2, ynew2 = znew2.real, znew2.imag
                        if Ya < 0:
                            theta = self.theta_norm_out
                        else:
                            theta = self.theta_norm_out + np.pi
                        Qx1, Qy1 = self.model.disvec(xnew1, ynew1)[:, layer] * direction
                        Qn1 = Qx1 * np.cos(theta) + Qy1 * np.sin(theta)
                        Qx2, Qy2 = self.model.disvec(xnew2, ynew2)[:, layer] * direction
                        Qn2 = Qx2 * np.cos(theta) + Qy2 * np.sin(theta)
                        if verbose:
                            print('xnew1, ynew1:', xnew1, ynew1)
                            print('xnew2, ynew2:', xnew2, ynew2)
                            print('Qn1, Qn2', Qn1, Qn2)
                            print('Qn2 > Qn1:', Qn2 > Qn1)
                        if Qn1 < 0:  # trying to cross line-sink that infiltrates, stay on bottom, don't terminate
                            if verbose: print('change 1')
                            xnew = xnew1
                            ynew = ynew1
                            dold = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                            dnew = np.sqrt((x1 - xnew) ** 2 + (y1 - ynew) ** 2)
                            znew = z1 + dnew / dold * (z2 - z1)
                            tnew = t1 + dnew / dold * (t2 - t1)
                            changed = True
                            xyztnew = [np.array([xnew, ynew, znew, tnew])]
                        elif Qn2 < 0:  # all water is taken out, terminate
                            if verbose: print('change 2')
                            xnew = xnew2
                            ynew = ynew2
                            dold = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                            dnew = np.sqrt((x1 - xnew) ** 2 + (y1 - ynew) ** 2)
                            znew = z1 + dnew / dold * (z2 - z1)
                            tnew = t1 + dnew / dold * (t2 - t1)
                            changed = True
                            terminate = True
                            xyztnew = [np.array([xnew, ynew, znew, tnew])]
                        elif Qn2 > Qn1:  # line-sink infiltrates
                            if verbose: print('change 3')
                            xnew = xnew2
                            ynew = ynew2
                            dold = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                            dnew = np.sqrt((x1 - xnew) ** 2 + (y1 - ynew) ** 2)
                            znew = z1 + dnew / dold * (z2 - z1)  # elevation just before jump
                            tnew = t1 + dnew / dold * (t2 - t1)
                            Qbelow = (znew - aq.z[modellayer + 1]) / aq.Haq[layer] * Qn1
                            znew2 = aq.z[modellayer + 1] + Qbelow / Qn2 * aq.Haq[layer]
                            changed = True
                            xyztnew = [np.array([xnew, ynew, znew, tnew]), np.array([xnew, ynew, znew2, tnew])]
                        else:  # line-sink takes part of water out
                            if verbose: print('change 4')
                            xnew = xnew2
                            ynew = ynew2
                            dold = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                            dnew = np.sqrt((x1 - xnew) ** 2 + (y1 - ynew) ** 2)
                            znew = z1 + dnew / dold * (z2 - z1)  # elevation just before jump
                            tnew = t1 + dnew / dold * (t2 - t1)
                            Qbelow = (znew - aq.z[modellayer + 1]) / aq.Haq[layer] * Qn1
                            if Qbelow > Qn2:  # taken out
                                terminate = True
                                xyztnew = [np.array([xnew, ynew, znew, tnew])]
                            else:
                                znew2 = aq.z[modellayer + 1] + Qbelow / Qn2 * aq.Haq[layer]
                                xyztnew = [np.array([xnew, ynew, znew, tnew]), np.array([xnew, ynew, znew2, tnew])]
                            changed = True
        return changed, terminate, xyztnew
    

class LineSinkBase(LineSinkChangeTrace, Element):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, Qls=100.0, \
                 res=0, wh=1, layers=0, name='LineSinkBase', label=None, \
                 addtomodel=True):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers, \
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
        # self.xa,self.ya,self.xb,self.yb,self.np = np.zeros(1),np.zeros(1),np.zeros(1),np.zeros(1),np.zeros(1,'i')  # needed to call bessel.circle_line_intersection

    def __repr__(self):
        return self.name + ' from ' + str((self.x1, self.y1)) + ' to ' + str(
            (self.x2, self.y2))

    def initialize(self):
        self.xc = np.array([0.5 * (self.x1 + self.x2)])
        self.yc = np.array([0.5 * (self.y1 + self.y2)])
        self.Ncp = 1
        self.z1 = self.x1 + 1j * self.y1
        self.z2 = self.x2 + 1j * self.y2
        self.L = np.abs(self.z1 - self.z2)
        self.theta_norm_out = np.arctan2(self.y2 - self.y1,
                                         self.x2 - self.x1) + np.pi / 2
        self.order = 0  # This is for uniform strength only
        self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        if self.addtomodel: self.aq.add_element(self)
        self.parameters = np.empty((self.Nparam, 1))
        self.parameters[:, 0] = self.Qls / self.L
        if self.wh == 'H':
            self.wh = self.aq.Haq[self.pylayers]
        elif self.wh == '2H':
            self.wh = 2.0 * self.aq.Haq[self.pylayers]
        elif np.isscalar(self.wh):
            self.wh = self.wh * np.ones(self.Nlayers)
        self.resfac = self.aq.T[self.pylayers] * self.res / self.wh

    def potinf(self, x, y, aq=None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nparam, aq.Naq))
        if aq == self.aq:
            pot = np.zeros(aq.Naq)
            pot[:] = besselaesnew.potbeslsho(x, y, self.z1, self.z2, aq.lab, 0,
                                             aq.ilap)
            rv[:] = self.aq.coef[self.pylayers] * pot
        return rv

    def disvecinf(self, x, y, aq=None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.Nparam, aq.Naq))
        if aq == self.aq:
            qxqy = np.zeros((2, aq.Naq))
            qxqy[:, :] = besselaesnew.disbeslsho(x, y, self.z1, self.z2, aq.lab,
                                                 0, aq.ilap)
            rv[0] = self.aq.coef[self.pylayers] * qxqy[0]
            rv[1] = self.aq.coef[self.pylayers] * qxqy[1]
        return rv
    
    def discharge(self):
        # returns the discharge in each layer
        Q = np.zeros(self.aq.Naq)
        Q[self.pylayers] = self.parameters[:, 0] * self.L
        return Q
    
    def plot(self):
        plt.plot([self.x1, self.x2], [self.y1, self.y2], 'k')

class HeadLineSink(LineSinkBase, HeadEquation):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, hls=1.0, \
                 res=0, wh=1, layers=0, label=None, addtomodel=True):
        self.storeinput(inspect.currentframe())
        LineSinkBase.__init__(self, model, x1, y1, x2, y2, Qls=0, \
                              res=res, wh=wh, layers=layers,
                              name='HeadLineSink', label=label, \
                              addtomodel=addtomodel)
        self.hc = hls
        self.Nunknowns = self.Nparam

    def initialize(self):
        LineSinkBase.initialize(self)
        self.pc = self.hc * self.aq.T[self.pylayers]  # Needed in solving

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class LineSinkHoBase(LineSinkChangeTrace, Element):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 Qls=0.0, layers=0, order=0, name='LineSinkHoBase', \
                 label=None, addtomodel=True, aq=None, zcinout=None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers, \
                         name=name, label=label)
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.x2 = float(x2)
        self.y2 = float(y2)
        self.Qls = np.atleast_1d(Qls)
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
        self.theta_norm_out = np.arctan2(self.y2 - self.y1,
                                       self.x2 - self.x1) + np.pi / 2.0  # changed minus to plus
        self.cosnorm = np.cos(self.theta_norm_out) * np.ones(self.Ncp)
        self.sinnorm = np.sin(self.theta_norm_out) * np.ones(self.Ncp)
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
        if self.addtomodel:
            self.aq.add_element(self)
        self.parameters = np.empty((self.Nparam, 1))
        # Not sure if that needs to be here
        self.parameters[:, 0] = self.Qls / self.L  

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
            # clever way of using a reshaped rv here
            potrv = rv.reshape((self.order + 1, self.Nlayers, aq.Naq)) 
            pot = np.zeros((self.order + 1, aq.Naq))
            pot[:, :] = besselaesnew.potbeslsv(x, y, self.z1, self.z2, aq.lab,
                                               self.order, aq.ilap)
            potrv[:] = self.aq.coef[self.pylayers] * pot[:, np.newaxis, :]
        return rv

    def disvecinf(self, x, y, aq=None):
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
            qxqyrv = rv.reshape((2, self.order + 1, self.Nlayers, aq.Naq))
            qxqy = np.zeros((2 * (self.order + 1), aq.Naq))
            qxqy[:, :] = besselaesnew.disbeslsv(x, y, self.z1, self.z2, aq.lab,
                                                self.order, aq.ilap)
            qxqyrv[0, :] = self.aq.coef[self.pylayers] * qxqy[:self.order + 1,
                                                         np.newaxis, :]
            qxqyrv[1, :] = self.aq.coef[self.pylayers] * qxqy[self.order + 1:,
                                                         np.newaxis, :]
        return rv
   
    def plot(self):
        plt.plot([self.x1, self.x2], [self.y1, self.y2], 'k')        
        
    def discharge(self):
        # returns the discharge in each layer
        rv = np.zeros(self.aq.Naq)
        Qdisinf = np.zeros((self.order + 1, self.Nlayers))
        for n in range(self.order + 1):
            Qdisinf[n] =  (1 ** (n + 1) - (-1) ** (n + 1)) / (n + 1)
        Qls = self.parameters[:, 0] * self.L / 2 * Qdisinf.ravel()
        rv[self.pylayers] = np.sum(Qls.reshape(self.order + 1, self.Nlayers), 0)
        return rv

class HeadLineSinkHo(LineSinkHoBase, HeadEquation):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 hls=1.0, order=0, layers=0, label=None, addtomodel=True):
        self.storeinput(inspect.currentframe())
        LineSinkHoBase.__init__(self, model, x1, y1, x2, y2, Qls=0, \
                                layers=layers, order=order,
                                name='HeadLineSinkHo', label=label, \
                                addtomodel=addtomodel)
        self.hc = hls
        self.Nunknowns = self.Nparam

    def initialize(self):
        LineSinkHoBase.initialize(self)
        self.resfac = 0.0  # Needs to be implemented
        self.pc = self.hc * self.aq.T[self.pylayers]  # Needed in solving

    def setparams(self, sol):
        self.parameters[:, 0] = sol


class LineSinkStringBase(Element):
    def __init__(self, model, xy, closed=False, layers=0, order=0,
                 name='LineSinkStringBase', label=None, aq=None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers, \
                         name=name, label=label)
        self.xy = np.atleast_2d(xy).astype('d')
        if closed: self.xy = np.vstack((self.xy, self.xy[0]))
        self.order = order
        self.aq = aq
        self.lslist = []
        self.x, self.y = self.xy[:, 0], self.xy[:, 1]
        self.Nls = len(self.x) - 1
        for i in range(self.Nls):
            self.lslist.append(LineSinkHoBase(model, \
                                              x1=self.x[i], y1=self.y[i],
                                              x2=self.x[i + 1],
                                              y2=self.y[i + 1], \
                                              Qls=0.0, layers=layers,
                                              order=order, label=label,
                                              addtomodel=False, aq=aq))

    def __repr__(self):
        return self.name + ' with nodes ' + str(self.xy)

    def initialize(self):
        for ls in self.lslist:
            ls.initialize()
        # Same order for all elements in string
        self.Ncp = self.Nls * self.lslist[0].Ncp  
        self.Nparam = self.Nls * self.lslist[0].Nparam
        self.Nunknowns = self.Nparam
        self.xls = np.empty((self.Nls, 2))
        self.yls = np.empty((self.Nls, 2))
        for i, ls in enumerate(self.lslist):
            self.xls[i, :] = [ls.x1, ls.x2]
            self.yls[i, :] = [ls.y1, ls.y2]
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.lslist[0].xc,
                                                      self.lslist[0].yc)
        self.parameters = np.zeros((self.Nparam, 1))
        # As parameters are only stored for the element not the list,
        # we need to combine the following
        self.xc = np.array([ls.xc for ls in self.lslist]).flatten()
        self.yc = np.array([ls.yc for ls in self.lslist]).flatten()
        self.xcin = np.array([ls.xcin for ls in self.lslist]).flatten()
        self.ycin = np.array([ls.ycin for ls in self.lslist]).flatten()
        self.xcout = np.array([ls.xcout for ls in self.lslist]).flatten()
        self.ycout = np.array([ls.ycout for ls in self.lslist]).flatten()
        self.cosnorm = np.array([ls.cosnorm for ls in self.lslist]).flatten()
        self.sinnorm = np.array([ls.sinnorm for ls in self.lslist]).flatten()
        self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], self.ycin[0])
        self.aqout = self.model.aq.find_aquifer_data(self.xcout[0],
                                                     self.ycout[0])

    def potinf(self, x, y, aq=None):
        '''
        linesink 0, order 0, layer[0]
                    order 0, layer[1]
                    ...
                    order 1, layer[0]
                    order 1, layer[1]
                    ...
        linesink 1, order 0, layer[0]
                    order 0, layer[1]
                    ...
                    order 1, layer[0]
                    order 1, layer[1]
                    ...
        '''    
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nls, self.lslist[0].Nparam, aq.Naq))
        for i in range(self.Nls):
            rv[i] = self.lslist[i].potinf(x, y, aq)
        rv.shape = (self.Nparam, aq.Naq)
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.Nls, self.lslist[0].Nparam, aq.Naq))
        for i in range(self.Nls):
            rv[:, i] = self.lslist[i].disvecinf(x, y, aq)
        rv.shape = (2, self.Nparam, aq.Naq)
        return rv
    
    def changetrace(self, xyzt1, xyzt2, aq, layer, ltype, modellayer, direction):
        changed = False
        terminate = False
        xyztnew = 0
        for ls in self.lslist:
            changed, terminate, xyztnew = ls.changetrace(xyzt1, xyzt2, aq, layer, ltype, modellayer, direction)
            if changed or terminate:
                return changed, terminate, xyztnew
        return changed, terminate, xyztnew
    
    def plot(self):
        plt.plot(self.x, self.y, 'k')
        
class HeadLineSinkString(LineSinkStringBase, HeadEquation):
    def __init__(self, model, xy=[(-1, 0), (1, 0)], hls=0.0, \
                 layers=0, order=0, label=None):
        self.storeinput(inspect.currentframe())
        LineSinkStringBase.__init__(self, model, xy, closed=False,
                                    layers=layers, order=order, \
                                    name='HeadLineSinkString', label=label,
                                    aq=None)
        self.hls = np.atleast_1d(hls)
        self.model.add_element(self)

    def initialize(self):
        LineSinkStringBase.initialize(self)
        self.aq.add_element(self)
        # self.pc = np.array([ls.pc for ls in self.lslist]).flatten()
        if len(self.hls) == 1:
            self.pc = self.hls * self.aq.T[self.pylayers] * np.ones(self.Nparam)
        elif len(self.hls) == self.Nls:  # head specified at centers
            self.pc = (self.hls[:, np.newaxis] * self.aq.T[self.pylayers]).flatten()
        elif len(self.hls) == 2:
            L = np.array([ls.L for ls in self.lslist])
            Ltot = np.sum(L)
            xp = np.zeros(self.Nls)
            xp[0] = 0.5 * L[0]
            for i in range(1, self.Nls):
                xp[i] = xp[i - 1] + 0.5 * (L[i - 1] + L[i])
            self.hls = np.interp(xp, [0, Ltot], self.hls)
            self.pc = (self.hls[:, np.newaxis] * self.aq.T[self.pylayers]).flatten()
        else:
            print('Error: hls entry not supported')
        self.resfac = 0.0

    def setparams(self, sol):
        self.parameters[:, 0] = sol
        
class LineSinkStringBase2(Element):
    '''
    Alternative implementation that loops through line-sinks to build equation
    '''
    def __init__(self, model, xy, closed=False, layers=0, order=0,
                 name='LineSinkStringBase', label=None, aq=None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers, \
                         name=name, label=label)
        self.xy = np.atleast_2d(xy).astype('d')
        if closed: self.xy = np.vstack((self.xy, self.xy[0]))
        self.order = order
        self.aq = aq
        self.lslist = []
        self.x, self.y = self.xy[:, 0], self.xy[:, 1]
        self.Nls = len(self.x) - 1
        for i in range(self.Nls):
            self.lslist.append(LineSinkHoBase(model, \
                                              x1=self.x[i], y1=self.y[i],
                                              x2=self.x[i + 1],
                                              y2=self.y[i + 1], \
                                              Qls=0.0, layers=layers,
                                              order=order, label=label,
                                              addtomodel=False, aq=aq))

    def __repr__(self):
        return self.name + ' with nodes ' + str(self.xy)

    def initialize(self):
        for ls in self.lslist:
            ls.initialize()
        # Same order for all elements in string
        self.Ncp = self.Nls * self.lslist[0].Ncp  
        self.Nparam = self.Nls * self.lslist[0].Nparam
        self.Nunknowns = self.Nparam
        self.xls = np.empty((self.Nls, 2))
        self.yls = np.empty((self.Nls, 2))
        for i, ls in enumerate(self.lslist):
            self.xls[i, :] = [ls.x1, ls.x2]
            self.yls[i, :] = [ls.y1, ls.y2]
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.lslist[0].xc,
                                                      self.lslist[0].yc)
        self.parameters = np.zeros((self.Nparam, 1))

    def potinf(self, x, y, aq=None):
        '''
        linesink 0, order 0, layer[0]
                    order 0, layer[1]
                    ...
                    order 1, layer[0]
                    order 1, layer[1]
                    ...
        linesink 1, order 0, layer[0]
                    order 0, layer[1]
                    ...
                    order 1, layer[0]
                    order 1, layer[1]
                    ...
        '''    
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nls, self.lslist[0].Nparam, aq.Naq))
        for i in range(self.Nls):
            rv[i] = self.lslist[i].potinf(x, y, aq)
        rv.shape = (self.Nparam, aq.Naq)
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.Nls, self.lslist[0].Nparam, aq.Naq))
        for i in range(self.Nls):
            rv[:, i] = self.lslist[i].disvecinf(x, y, aq)
        rv.shape = (2, self.Nparam, aq.Naq)
        return rv
    
    def changetrace(self, xyzt1, xyzt2, aq, layer, ltype, modellayer, direction):
        changed = False
        terminate = False
        xyztnew = 0
        for ls in self.lslist:
            changed, terminate, xyztnew = ls.changetrace(xyzt1, xyzt2, aq, layer, ltype, modellayer, direction)
            if changed or terminate:
                return changed, terminate, xyztnew
        return changed, terminate, xyztnew
    
    def plot(self):
        plt.plot(self.x, self.y, 'k')
        
class HeadLineSinkString2(LineSinkStringBase2):
    '''
    Not yet modified
    '''
    def __init__(self, model, xy=[(-1, 0), (1, 0)], hls=0.0, \
                 layers=0, order=0, label=None):
        self.storeinput(inspect.currentframe())
        LineSinkStringBase2.__init__(self, model, xy, closed=False,
                                    layers=layers, order=order, \
                                    name='HeadLineSinkString', label=label,
                                    aq=None)
        self.hls = np.atleast_1d(hls)
        self.model.add_element(self)

    def initialize(self):
        LineSinkStringBase2.initialize(self)
        self.aq.add_element(self)
        # self.pc = np.array([ls.pc for ls in self.lslist]).flatten()
        if len(self.hls) == 1:
            self.pc = self.hls * self.aq.T[self.pylayers] * np.ones(self.Nparam)
        elif len(self.hls) == self.Nls:  # head specified at centers
            self.pc = (self.hls[:, np.newaxis] * self.aq.T[self.pylayers]).flatten()
        elif len(self.hls) == 2:
            L = np.array([ls.L for ls in self.lslist])
            Ltot = np.sum(L)
            xp = np.zeros(self.Nls)
            xp[0] = 0.5 * L[0]
            for i in range(1, self.Nls):
                xp[i] = xp[i - 1] + 0.5 * (L[i - 1] + L[i])
            self.hls = np.interp(xp, [0, Ltot], self.hls)
            self.pc = (self.hls[:, np.newaxis] * self.aq.T[self.pylayers]).flatten()
        else:
            print('Error: hls entry not supported')
        self.resfac = 0.0

    def setparams(self, sol):
        self.parameters[:, 0] = sol

'''
Copyright (C), 2015, Mark Bakker.
TTim is distributed under the MIT license
'''

import numpy as np
import inspect # Used for storing the input
from scipy.special import k0
from besselaes import potbeslsho, potbesonlylsho

class ModelBase:
    def __init__(self, kaq, Haq, c, z, npor, ltype):
        # All input variables should be numpy arrays
        # That should be checked outside this function
        self.elementlist = []
        self.elementdict = {}
        self.aq = Aquifer(self, kaq, Haq, c, z, npor, ltype)
        self.modelname = 'ml' # Used for writing out input
    def initialize(self):
        self.aq.initialize()
        for e in self.elementlist:
            e.initialize()
    def add_element(self, e):
        self.elementlist.append(e)
        if e.label is not None: self.elementdict[e.label] = e
    def remove_element(self, e):
        if e.label is not None: self.elementdict.pop(e.label)
        self.elementlist.remove(e)
    def storeinput(self, frame):
        self.inputargs, _, _, self.inputvalues = inspect.getargvalues(frame)
    def potential(self, x, y, aq = None):
        if aq is None: aq = self.aq.find_aquifer_data(x,y)
        pot = np.zeros(aq.Naq)
        for e in self.elementlist:
            pot += e.potential(x, y, aq)
        rv = np.sum(pot * aq.eigvec, 1)
        return rv
    def head(self, x, y, layers = None, aq = None):
        if aq is None: aq = self.aq.find_aquifer_data(x,y)
        rv = self.potential(x, y, aq) / aq.T
        if layers is None:
            return rv
        else:
            return rv[layers]
    def headgrid(self, x1, x2, nx, y1, y2, ny, layers=None, printrow=False):
        '''Returns h[Nlayers,Ny,Nx]. If layers is None, all layers are returned'''
        xg, yg = np.linspace(x1,x2,nx), np.linspace(y1,y2,ny)
        if layers is None:
            Nlayers = self.aq.find_aquifer_data(xg[0],yg[0]).Naq
        else:
            Nlayers = len(np.atleast_1d(layers))
        h = np.empty((Nlayers, ny, nx))
        for j in range(ny):
            if printrow: print str(j)+' '
            for i in range(nx):
                h[:,j,i] = self.head(xg[i], yg[j], layers)
        return h
    def headalongline(self, x, y, layers=None):
        '''Returns head[Nlayers,len(x)]
        Assumes same number of layers for each x and y
        layers may be None or list of layers for which head is computed'''
        xg, yg = np.atleast_1d(x), np.atleast_1d(y)
        if layers is None:
            Nlayers = self.aq.find_aquifer_data(xg[0],yg[0]).Naq
        else:
            Nlayers = len(np.atleast_1d(layers))
        nx = len(xg)
        if len(yg) == 1:
            yg = yg * np.ones(nx)
        h = np.zeros((Nlayers, nx))
        for i in range(nx):
            h[:,i] = self.head(xg[i], yg[i], layers)
        return h
    def solve(self, printmat = 0, sendback = 0, silent = False):
        '''Compute solution'''
        # Initialize elements
        self.initialize()
        # Compute number of equations
        self.Neq = np.sum([e.Nunknowns for e in self.elementlist])
        if silent is False:
            print 'self.Neq ',self.Neq
        if self.Neq == 0:
            if silent is False: print 'No unknowns. Solution complete'
            return
        mat = np.empty((self.Neq,self.Neq))
        rhs = np.empty(self.Neq)
        ieq = 0
        for e in self.elementlist:
            if e.Nunknowns > 0:
                mat[ieq:ieq+e.Nunknowns, :], rhs[ieq:ieq+e.Nunknowns] = e.equation()
                ieq += e.Nunknowns
        if printmat:
            return mat,rhs
        sol = np.linalg.solve(mat, rhs)
        icount = 0
        for e in self.elementlist:
            if e.Nunknowns > 0:
                e.setparams(sol[icount:icount + e.Nunknowns])
                icount += e.Nunknowns
                # e.run_after_solve()
        if silent is False:
            print 'solution complete'
        elif (silent == 'dot') or (silent == '.'):
            print '.',
            sys.stdout.flush()  # Can be replaced with print with flush in Python 3.3
        if sendback:
            return sol
        return
        
def param_maq(kaq, z, c, npor, top):
    # Computes the parameters for a ModelBase from input for a maq model
    kaq = np.atleast_1d(kaq).astype('d')
    z = np.atleast_1d(z).astype('d')
    c = np.atleast_1d(c).astype('d')
    npor = np.atleast_1d(npor).astype('d')
    if top == 'conf':
        Naq = len(z) / 2
        ltype = np.array( list( (Naq-1) * 'al' + 'a' ) )
    else: # leaky layer on top
        Naq = (len(z) - 1) / 2
        ltype = np.array( list( Naq * 'la' ) )
    if len(kaq) == 1: kaq = kaq * np.ones(Naq)
    assert len(kaq) == Naq, 'Error: length of kaq needs to be 1 or' + str(Naq)
    H = z[:-1] - z[1:]
    assert np.all(H >= 0), 'Error: Not all layers thicknesses are non-negative' + str(H) 
    if top == 'conf':
        if len(c) == 1: c = c * np.ones(Naq - 1)
        if len(npor) == 1: npor = npor * np.ones(2 * Naq - 1)
        assert len(c) == Naq-1, 'Error: Length of c needs to be 1 or' + str(Naq-1)
        assert len(npor) == 2 * Naq - 1, 'Error: Length of npor needs to be 1 or' + str(2*Naq-1)
        Haq = H[::2]
        c = np.hstack((1e100,c))
    else: # leaky layer on top
        if len(c) == 1: c = c * np.ones(Naq)
        if len(npor) == 1: npor = npor * np.ones(2 * Naq)
        assert len(c) == Naq, 'Error: Length of c needs to be 1 or' + str(Naq)
        assert len(npor) == 2 * Naq, 'Error: Length of npor needs to be 1 or' + str(2*Naq)
        Haq = H[1::2]
    return kaq, Haq, c, npor, ltype
        
class ModelMaq(ModelBase):
    def __init__(self, kaq = 1, z = [1,0], c = [], npor = 0.3, top = 'conf'):
        self.storeinput(inspect.currentframe())
        kaq, Haq, c, npor, ltype = param_maq(kaq, z, c, npor, top)
        ModelBase.__init__(self, kaq, Haq, c, z, npor, ltype)
        self.name = 'ModelMaq'
    
class AquiferData:
    def __init__(self, model, kaq, Haq, c, z, npor, ltype):
        # All input variables except model should be numpy arrays
        # That should be checked outside this function
        self.model = model
        # Needed for heads
        self.kaq = kaq
        self.Naq = len(kaq)
        self.Haq = Haq
        self.T = self.kaq * self.Haq
        self.Tcol = self.T[:,np.newaxis]
        self.c = c
        # Needed for tracing
        self.z = z
        self.npor = npor
        self.ltype = ltype
    def initialize(self):
        d0 = 1.0 / (self.c * self.T)
        d0[:-1] += 1.0 / (self.c[1:] * self.T[:-1])
        dp1 = -1.0 / (self.c[1:] * self.T[1:])
        dm1 = -1.0 / (self.c[1:] * self.T[:-1])
        A = np.diag(dm1,-1) + np.diag(d0,0) + np.diag(dp1,1)
        w,v = np.linalg.eig(A)
        # sort lab in ascending order, hence lab in descending order
        index = np.argsort( abs(w) )
        w = w[index]; v = v[:,index]
        if self.ltype[0] == 'a':
            self.lab = 1.0 / np.sqrt(w[1:])
            self.zeropluslab = np.zeros(self.Naq)
            self.zeropluslab[1:] = self.lab
            v[:,0] = self.T / np.sum(self.T) # first column is normalized T
        else:
            self.lab = 1.0 / np.sqrt(w)
        self.eigvec = v
        self.coef = np.linalg.solve(v, np.diag(np.ones(self.Naq))).T
        
class Aquifer(AquiferData):
    def __init__(self, model, kaq, Haq, c, z, npor, ltype):
        AquiferData.__init__(self, model, kaq, Haq, c, z, npor, ltype)
        self.inhomlist = []
        self.area = 1e300 # Needed to find smallest inhom
    def initialize(self):
        AquiferData.initialize(self)  # cause we are going to call initialize for inhoms
        for inhom in self.inhomlist:
            inhom.initialize()
    def find_aquifer_data(self, x, y):
        rv = self
        for inhom in self.inhomlist:
            if inhom.is_inside(x,y):
                if inhom.area < rv.area:
                    rv = inhom
        return rv
    
class HeadEquationNores:
    def equation(self):
        '''Mix-in class that returns matrix rows for head-specified conditions. (really written as constant potential element)
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            rhs[istart:istart+self.Nlayers] = self.pc
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    e.potinflayers(self.xc[icp], self.yc[icp], self.pylayers)
                    ieq += e.Nunknowns
                else:
                    rhs[istart:istart+self.Nlayers] -= \
                    e.potentiallayers(self.xc[icp], self.yc[icp], self.pylayers)  # Pretty cool that this works, really
        return mat, rhs
    
class Element:
    def __init__(self, model, Nparam, Nunknowns, layers, name, label):
        self.model = model
        self.aq = None # Set in the initialization function
        self.Nparam = Nparam
        self.Nunknowns = Nunknowns
        self.pylayers = np.atleast_1d(layers)
        self.Nlayers = len(self.pylayers)
        self.name = name
        self.label = label
        if self.label is not None:
            assert self.label not in self.model.elementDict.keys(),\
            "TTim error: label " + self.label + " already exists"
    def initialzie(self):
        # must be overloade
        pass
    def potinf(self, x, y, aq = None):
        '''Returns array of size (Nparam, Naq)'''
        raise Exception('Must overload Element.potinf()')
    def potential(self, x, y, aq = None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        return np.sum(self.parameters * self.potinf(x, y, aq), 0)
    def potinflayers(self, x, y, pylayers, aq = None):
        '''Returns array of size (len(pylayers),Nparam)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data( x, y )
        pot = self.potinf(x, y, aq)  # Nparam rows, Naq cols
        rv = np.sum(pot[:,np.newaxis,:] * aq.eigvec, 2).T  # Transopose as the first axes needs to be the number of layers
        return rv[pylayers,:]
    def potentiallayers(self, x, y, pylayers, aq = None):
        '''Returns array of size len(pylayers)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        pot = np.sum(self.potential(x, y, aq) * aq.eigvec, 1 )
        return pot[pylayers]
    def setparams(self, sol):
        raise Exception('Must overload Element.setparams()')
    def storeinput(self,frame):
        self.inputargs, _, _, self.inputvalues = inspect.getargvalues(frame)

    
class WellBase(Element):
    def __init__(self, model, xw = 0, yw = 0, Qw = 100.0, rw = 0.1, \
                 res = 0.0, layers = 0, name = 'WellBase', label = None):
        Element.__init__(self, model, Nparam = 1, Nunknowns = 0, layers = layers,\
                         name = name, label = label)
        self.Nparam = len(self.pylayers) # Defined here and not in Element as other elements can have multiple parameters per layers
        self.xw = float(xw)
        self.yw = float(yw)
        self.Qw = np.atleast_1d(Qw)
        self.rw = float(rw)
        self.res = float(res)
        self.model.add_element(self)
    def __repr__(self):
        return self.name + ' at ' + str((self.xw, self.yw))
    def initialize(self):
        self.xc = np.array([self.xw + self.rw])
        self.yc = np.array([self.yw])
        self.Ncp = 1
        self.aq = self.model.aq.find_aquifer_data(self.xw, self.yw)
        self.parameters = np.empty((self.Nparam, 1))
        self.parameters[:,0] = self.Qw  # Not sure if that needs to be here
    def potinf(self, x, y, aq = None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nparam, aq.Naq))
        if aq == self.aq:
            pot = np.zeros(aq.Naq)
            r = np.sqrt( (x-self.xw)**2 + (y-self.yw)**2 )
            if r < self.rw: r = self.rw  # If at well, set to at radius
            if aq.ltype[0] == 'a':
                pot[0] = np.log(r / self.rw) / (2 * np.pi)
                pot[1:] = -k0(r / aq.lab) / (2 * np.pi)
            else:
                pot[:] = k0(r / aq/lab) / (2 * np.pi)
            rv[:] = self.aq.coef[self.pylayers] * pot
        return rv
    
class HeadWell(WellBase, HeadEquationNores):
    def __init__(self, model, xw = 0, yw = 0, hw = 10.0, rw = 0.1, \
                 layers = 0, label = None):
        self.storeinput(inspect.currentframe())
        WellBase.__init__(self, model, xw, yw, 0.0, rw,\
                          layers = layers, name = 'HeadWell', label = label)
        self.hc = hw
        self.Nunknowns = self.Nparam
    def initialize(self):
        WellBase.initialize(self)
        self.pc = self.hc * self.aq.T[self.pylayers] # Needed in solving
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class Constant(Element, HeadEquationNores):
    def __init__(self, model, xr = 0, yr = 0, hr = 0.0, \
                 layer = 0, label = None):
        Element.__init__(self, model, Nparam = 1, Nunknowns = 1, layers = layer,\
                         name = 'Constant', label = label)
        self.Nunknowns = 1
        self.xr = xr
        self.yr = yr
        self.hr = hr
        self.model.add_element(self)
    def initialize(self):
        self.aq = self.model.aq.find_aquifer_data(self.xr, self.yr)
        self.Ncp = 1
        self.xc = np.array([self.xr])
        self.yc = np.array([self.yr])
        self.pc = self.hr * self.aq.T[self.pylayers]
        self.parameters = np.empty((1, 1))
    def potinf(self, x, y, aq = None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((1, aq.Naq))
        if aq == self.aq:
            rv[0,0] = 1
        return rv
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class LineSinkBase(Element):
    def __init__(self, model, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
                 Qls = 100.0, layers = 0, name = 'LineSinkBase', label = None, \
                 addtomodel = True):
        Element.__init__(self, model, Nparam = 1, Nunknowns = 0, layers = layers,\
                         name = name, label = label)
        self.Nparam = len(self.pylayers)
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.x2 = float(x2)
        self.y2 = float(y2)
        self.Qls = np.atleast_1d(Qls)
        if addtomodel: self.model.add_element(self)
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
        self.parameters = np.empty((self.Nparam, 1))
        self.parameters[:,0] = self.Qls / self.L  # Not sure if that needs to be here
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
    
class HeadLineSink(LineSinkBase, HeadEquationNores):
    def __init__(self, model, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
                 hls = 1.0, layers = 0, label = None, addtomodel = True):
        self.storeinput(inspect.currentframe())
        LineSinkBase.__init__(self, model, x1, y1, x2, y2, Qls = 0,\
                 layers = layers, name = 'HeadLineSink', label = label, \
                 addtomodel = addtomodel)
        self.hc = hls
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineSinkBase.initialize(self)
        self.pc = self.hc * self.aq.T[self.pylayers] # Needed in solving
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class LineSinkHoBase(Element):
    def __init__(self, model, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
                 Qls = 100.0, layers = 0, order = 0, \
                 name = 'LineSinkBase', label = None, addtomodel = True):
        Element.__init__(self, model, Nparam = 1, Nunknowns = 0, layers = layers,\
                         name = name, label = label)
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.x2 = float(x2)
        self.y2 = float(y2)
        self.Qls = np.atleast_1d(Qls)
        self.order = order
        self.Nparam = self.Nlayers * (self.order + 1)
        if addtomodel: self.model.add_element(self)
        #self.xa,self.ya,self.xb,self.yb,self.np = np.zeros(1),np.zeros(1),np.zeros(1),np.zeros(1),np.zeros(1,'i')  # needed to call bessel.circle_line_intersection
    def __repr__(self):
        return self.name + ' from ' + str((self.x1, self.y1)) +' to '+str((self.x2, self.y2))
    def initialize(self):
        self.Ncp = self.order + 1
        self.z1 = self.x1 + 1j * self.y1
        self.z2 = self.x2 + 1j * self.y2
        self.L = np.abs(self.z1 - self.z2)
        #
        #thetacp = np.arange(np.pi, 0, -np.pi/self.Ncp) - 0.5 * np.pi/self.Ncp
        # The following works MUCH better for a uniform head along the line
        thetacp = np.linspace(np.pi, 0, self.Ncp+2)[1:-1]
        Zcp = np.zeros( self.Ncp, 'D' )
        Zcp.real = np.cos(thetacp)
        Zcp.imag = 1e-6  # control point just on positive site (this is handy later on)
        zcp = Zcp * (self.z2 - self.z1) / 2.0 + 0.5 * (self.z1 + self.z2)
        self.xc = zcp.real; self.yc = zcp.imag
        #
        self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
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
            rv.shape = (self.order+1, self.Nlayers, aq.Naq)
            pot = np.zeros(aq.Naq)
            if aq.ltype[0] == 'a':
                for i in range(self.order+1): # This should be done inside FORTRAN extension
                    potbeslsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                               aq.Naq, aq.zeropluslab, i, pot)  # Call FORTRAN extension
                    rv[i] = self.aq.coef[self.pylayers] * pot
            else:
                for i in range(self.order+1):
                    potbesonlylsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                                   aq.Naq, aq.lab, 0, pot)
                    rv[i] = self.aq.coef[self.pylayers] * pot
            rv.shape = (self.Nparam, aq.Naq)
        return rv
    
class HeadLineSinkHo(LineSinkHoBase, HeadEquationNores):
    def __init__(self, model, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
                 hls = 1.0, order = 0, layers = 0, label = None, addtomodel = True):
        self.storeinput(inspect.currentframe())
        LineSinkHoBase.__init__(self, model, x1, y1, x2, y2, Qls = 0, \
                 layers = layers, order = order, name = 'HeadLineSinkHo', label = label, \
                 addtomodel = addtomodel)
        self.hc = hls
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineSinkHoBase.initialize(self)
        self.pc = self.hc * self.aq.T[self.pylayers] # Needed in solving
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class LineSinkStringBase(Element):
    def __init__(self, model, layers = 0, order = 0, name = 'LineSinkStringBase', label = None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers = layers, \
                         name=name, label=label)
        self.order = order
        self.lslist = []
    def __repr__(self):
        return self.name + ' with nodes ' + str(zip(self.x,self.y))
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
        self.aq = self.model.aq.find_aquifer_data(self.lslist[0].xc, self.lslist[0].yc)
        self.parameters = np.zeros((self.Nparam,1))
        ## As parameters are only stored for the element not the list, we need to combine the following
        xc = []
        yc = []
        for ls in self.lslist:
            xc.extend(list(ls.xc))
            yc.extend(list(ls.yc))
        self.xc = np.array(xc)
        self.yc = np.array(yc)

        #self.resfac = self.ldList[0].resfac  # same for all elements in the list
        #self.xcneg, self.ycneg = np.zeros(self.Ncp), np.zeros(self.Ncp)
        #self.cosout, self.sinout = np.zeros(self.Ncp), np.zeros(self.Ncp)
        #for i,ld in enumerate(self.ldList):
        #    self.xc[i*ld.Ncp:(i+1)*ld.Ncp], self.yc[i*ld.Ncp:(i+1)*ld.Ncp] = ld.xc, ld.yc
        #    self.xcneg[i*ld.Ncp:(i+1)*ld.Ncp], self.ycneg[i*ld.Ncp:(i+1)*ld.Ncp] = ld.xcneg, ld.ycneg
        #    self.cosout[i*ld.Ncp:(i+1)*ld.Ncp], self.sinout[i*ld.Ncp:(i+1)*ld.Ncp] = ld.cosout, ld.sinout
    def potinf(self,x,y,aq=None):
        '''Returns array (Nunknowns,Nperiods)'''
        if aq is None: aq = self.model.aq.findAquiferData( x, y )
        rv = np.zeros((self.Nls, self.lslist[0].Nparam, aq.Naq))
        for i in range(self.Nls):
            rv[i] = self.lslist[i].potinf(x, y, aq)
        rv.shape = (self.Nparam, aq.Naq)
        return rv

class HeadLineSinkString(LineSinkStringBase, HeadEquationNores):
    def __init__(self, model, xy=[(-1,0), (1,0)], hls = 0.0, \
                 layers = 0, order = 0, label = None):
        self.storeinput(inspect.currentframe())
        LineSinkStringBase.__init__(self, model, layers = layers, \
                                    name = 'HeadLineSinkString', label = None)
        xy = np.atleast_2d(xy).astype('d')
        self.x, self.y = xy[:,0], xy[:,1]
        self.Nls = len(self.x) - 1
        for i in range(self.Nls):
            self.lslist.append(HeadLineSinkHo(model, \
                x1 = self.x[i], y1 = self.y[i], x2 = self.x[i+1], y2 = self.y[i+1], \
                hls = hls, layers = layers, order = order, label = label, addtomodel = False))
        self.model.add_element(self)
    def initialize(self):
        LineSinkStringBase.initialize(self)
        self.pc = self.lslist[0].pc  # same for all control points
    def setparams(self, sol):
        self.parameters[:,0] = sol
   


ml = ModelMaq(kaq = [1,2,3], z = [5,4,3,2,1,0], c = 1000)
#w = WellBase(ml, xw = 10, yw = 20, Qw = 100, layers = 0)
#w2 = WellBase(ml, xw = 0, yw = 0, Qw = [100,200], rw = 0.1, layers = [0,1])
#w2 = HeadWell(ml, 0, 0, hw = [10.0, 12.0], rw = 0.1, layers = [0,1])
#ls = LineSinkBase(ml, x1 = -20, y1 = -10, x2 = 40, y2 = 20, \
#                 Qls = 100.0, layers = 0)
#ls = HeadLineSink(ml, x1 = -20, y1 = -10, x2 = 40, y2 = 20, \
#                 hls = 12.0, layers = [0,1])
#rf = Constant(ml, -10, 20, 20, layer = 0)
#ls = LineSinkHoBase(ml, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
#                 Qls = 100.0, layers = [0,1], order = 3)
#ls = HeadLineSinkHo(ml, x1 = -20, y1 = -10, x2 = 40, y2 = 20, \
#                    hls = 12.0, layers = [0,2], order = 4)
lsstring = HeadLineSinkString(ml, xy = [(-10,-10), (0,0), (10,0), (10,10)], hls = 7, layers = [0,1], order = 5)
#ml.initialize()
ml.solve()
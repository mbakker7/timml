'''
Copyright (C), 2015, Mark Bakker.
TTim is distributed under the MIT license
'''

import numpy as np
import inspect # Used for storing the input
from scipy.special import k0, k1
from besselaes import potbeslsho, potbesonlylsho, disbeslsho, disbesonlylsho, potbesldho, disbesldho, omegalapldho, omegalaplsho

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
    def disvec(self, x, y, pylayers = None, aq = None):
        if aq is None: aq = self.aq.find_aquifer_data(x,y)
        rv = np.zeros((2, aq.Naq))
        for e in self.elementlist:
            rv += e.disvec(x, y, aq)
        rv = np.sum(rv[:,np.newaxis,:] * self.aq.eigvec, 2)
        return rv
    def disvecnorm(self, x, y, cosnorm, sinnorm):
        qxqy = self.disvec(x, y)
        return qxqy[0] * cosnorm + qxqy[1] * sinnorm
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
        if self.Neq == 0: return
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
        #
        self.area = 1e200 # Smaller than default of ml.aq so that inhom is found
    def initialize(self):
        self.elementlist = []  # Elementlist of aquifer
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
    def isinside(self, x, y):
        raise Exception('Must overload AquiferData.isinside()')
    def storeinput(self,frame):
        self.inputargs, _, _, self.inputvalues = inspect.getargvalues(frame)
    def add_element(self, e):
        self.elementlist.append(e)

class Aquifer(AquiferData):
    def __init__(self, model, kaq, Haq, c, z, npor, ltype):
        AquiferData.__init__(self, model, kaq, Haq, c, z, npor, ltype)
        self.inhomlist = []
        self.area = 1e300 # Needed to find smallest inhom
    def initialize(self):
        AquiferData.initialize(self)  # cause we are going to call initialize for inhoms
        for inhom in self.inhomlist:
            inhom.initialize()
    def add_inhom(self, inhom):
        self.inhomlist.append(inhom)
    def find_aquifer_data(self, x, y):
        rv = self
        for inhom in self.inhomlist:
            if inhom.isinside(x, y):
                if inhom.area < rv.area:
                    rv = inhom
        return rv
    
class PolygonInhom(AquiferData):
    tiny = 1e-8
    def __init__(self, model, xy, kaq, Haq, c, z, npor, ltype):
        # All input variables except model should be numpy arrays
        # That should be checked outside this function):        
        AquiferData.__init__(self, model, kaq, Haq, c, z, npor, ltype)
        self.model.aq.add_inhom(self)
        self.z1, self.z2 = compute_z1z2(xy)
        self.x, self.y = zip(*xy)
        self.xmin = min(self.x)
        self.xmax = max(self.x)
        self.ymin = min(self.y)
        self.ymax = max(self.y)
    def __repr__(self):
        return 'PolygonInhom'
    def isinside(self, x, y):
        rv = 0
        if (x >= self.xmin) and (x <= self.xmax) and (y >= self.ymin) and (y <= self.ymax):
            z = complex(x,y)
            bigZ = (2.0*z - (self.z1 + self.z2))/ (self.z2 - self.z1)
            bigZmin1 = bigZ - 1.0
            bigZplus1 = bigZ + 1.0
            minAbsBigZmin1 = min(abs(bigZmin1))
            minAbsBigZplus1 = min(abs(bigZplus1))
            if minAbsBigZmin1 < self.tiny or minAbsBigZplus1 < self.tiny:
                rv = 1
                return rv
            angles = np.log(bigZmin1 / bigZplus1).imag
            angle = np.sum(angles)
            if angle > np.pi: rv = 1
        return rv
    
class PolygonInhomMaq(PolygonInhom):
    tiny = 1e-8
    def __init__(self, model, xy, kaq = 1, z = [1,0], c = [], npor = 0.3, top = 'conf'):
        self.storeinput(inspect.currentframe())
        kaq, Haq, c, npor, ltype = param_maq(kaq, z, c, npor, top)
        PolygonInhom.__init__(self, model, xy, kaq, Haq, c, z, npor, ltype)
    
def compute_z1z2(xy):
    # Returns z1 and z2 of polygon, in clockwise order
    x,y = zip(*xy)
    if x[0] == x[-1] and y[0] == y[-1]:  # In case last point is repeated
        x = x[:-1]; y = y[:-1]
    z1 = np.array(x) + np.array(y) * 1j
    index = range(1,len(z1)) + [0]
    z2 = z1[index]
    Z = 1e-6j
    z = Z * (z2[0] - z1[0]) / 2.0 + 0.5 * (z1[0] + z2[0])
    bigZ = ( 2.0*z - (z1 + z2) )/ (z2 - z1)
    bigZmin1 = bigZ - 1.0
    bigZplus1 = bigZ + 1.0
    angle = np.sum(np.log(bigZmin1 / bigZplus1).imag)
    if angle < np.pi: # reverse order
        z1 = z1[::-1]
        z2 = z2[::-1]
    return z1, z2
    
class HeadEquation:
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
    
class DisvecEquation:
    def equation(self):
        '''Mix-in class that returns matrix rows for head-specified conditions. (really written as constant potential element)
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    qx, qy = e.disinflayers(self.xc[icp], self.yc[icp], self.pylayers)
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
                    ieq += e.Nunknowns
                else:
                    qx, qy = e.disveclayers(self.xc[icp], self.yc[icp], self.pylayers)
                    rhs[istart:istart+self.Nlayers] -= qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
        return mat, rhs
    
class DisvecEquationOut:
    def equation(self):
        '''Mix-in class that returns matrix rows for head-specified conditions. (really written as constant potential element)
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    qx, qy = e.disinflayers(self.xcout[icp], self.ycout[icp], self.pylayers)
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
                    ieq += e.Nunknowns
                else:
                    qx, qy = e.disveclayers(self.xcout[icp], self.ycout[icp], self.pylayers)
                    rhs[istart:istart+self.Nlayers] -= qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
        return mat, rhs
    
class HeadDiffEquation:
    def equation(self):
        '''Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    qx, qy = e.disinflayers(self.xcout[icp], self.ycout[icp], self.pylayers)
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    e.potinflayers(self.xcin[icp], self.ycin[icp], self.pylayers, aq=self.aqin) / self.aqin.Tcol - \
                    e.potinflayers(self.xcout[icp], self.ycout[icp], self.pylayers, aq=self.aqout) / self.aqout.Tcol
                    ieq += e.Nunknowns
                else:
                    rhs[istart:istart+self.Nlayers] -= \
                    e.potentiallayers(self.xcin[icp], self.ycin[icp], self.pylayers, aq=self.aqin) / self.aqin.T - \
                    e.potentiallayers(self.xcout[icp], self.ycout[icp], self.pylayers, aq=self.aqout) / self.aqout.T
        return mat, rhs
    
class HeadDiffEquation2:
    def equation(self):
        '''Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    qx, qy = e.disinflayers(self.xcout[icp], self.ycout[icp], self.pylayers)
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    e.potinflayers(self.xcin[icp], self.ycin[icp], self.pylayers, aq=self.aqin) / self.aqin.Tcol - \
                    e.potinflayers(self.xcout[icp], self.ycout[icp], self.pylayers, aq=self.aqout) / self.aqout.Tcol
                    ieq += e.Nunknowns
                else:
                    rhs[istart:istart+self.Nlayers] -= \
                    e.potentiallayers(self.xcin[icp], self.ycin[icp], self.pylayers, aq=self.aqin) / self.aqin.T - \
                    e.potentiallayers(self.xcout[icp], self.ycout[icp], self.pylayers, aq=self.aqout) / self.aqout.T
        return mat, rhs
    
class DisvecDiffEquation:
    def equation(self):
        '''Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    qxin, qyin = e.disinflayers(self.xcin[icp], self.ycin[icp], self.pylayers, aq=self.aqin)
                    qxout, qyout = e.disinflayers(self.xcout[icp], self.ycout[icp], self.pylayers, aq=self.aqout)
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    (qxin - qxout) * self.cosnorm[icp] + (qyin - qyout) * self.sinnorm[icp]
                    ieq += e.Nunknowns
                else:
                    qxin, qyin = e.disveclayers(self.xcin[icp], self.ycin[icp], self.pylayers, aq=self.aqin)
                    qxout, qyout = e.disveclayers(self.xcout[icp], self.ycout[icp], self.pylayers, aq=self.aqout)
                    rhs[istart:istart+self.Nlayers] -= (qxin - qxout) * self.cosnorm[icp] + (qyin - qyout) * self.sinnorm[icp]
        return mat, rhs
    
class DisvecDiffEquation2:
    def equation(self):
        '''Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    fluxin = self.intflux(e.disinflayers, self.xcin[icp], self.ycin[icp], \
                                          self.xcin[icp+1], self.ycin[icp+1], self.pylayers, ndeg=8, aq=self.aqin)
                    fluxout = self.intflux(e.disinflayers, self.xcout[icp], self.ycout[icp], \
                                          self.xcout[icp+1], self.ycout[icp+1], self.pylayers, ndeg=8, aq=self.aqout)
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = fluxin - fluxout
                    ieq += e.Nunknowns
                else:
                    fluxin = self.intflux(e.disveclayers, self.xcin[icp], self.ycin[icp], \
                                          self.xcin[icp+1], self.ycin[icp+1], self.pylayers, ndeg=8, aq=self.aqin)
                    fluxout = self.intflux(e.disveclayers, self.xcout[icp], self.ycout[icp], \
                                          self.xcout[icp+1], self.ycout[icp+1], self.pylayers, ndeg=8, aq=self.aqout)                    
                    rhs[istart:istart+self.Nlayers] -= fluxin - fluxout
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
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        pot = self.potinf(x, y, aq)  # Nparam rows, Naq cols
        rv = np.sum(pot[:,np.newaxis,:] * aq.eigvec, 2).T  # Transopose as the first axes needs to be the number of layers
        return rv[pylayers,:]
    def potentiallayers(self, x, y, pylayers, aq = None):
        '''Returns array of size len(pylayers)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        pot = np.sum(self.potential(x, y, aq) * aq.eigvec, 1 )
        return pot[pylayers]
    def disinf(self, x, y, aq = None):
        '''Returns array of size (2, Nparam, Naq)'''
        raise Exception('Must overload Element.disinf()')
    def disvec(self, x, y, aq = None):
        '''Returns array of size (2, Nparam, Naq)'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        return np.sum(self.parameters * self.disinf(x, y, aq), 1)
    def disvecnorm(self, x, y, cosnorm, sinnorm):
        qxqy = self.disvec(x, y)
        return qxqy[0] * cosnorm + qxqy[1] * sinnorm
    def disinflayers(self, x, y, pylayers, aq = None):
        '''Returns two arrays of size (len(pylayers),Nparam)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        qxqy = self.disinf(x, y, aq)  # Nparam rows, Naq cols
        qx = np.sum(qxqy[0,:,np.newaxis,:] * aq.eigvec, 2).T  # Transpose as the first axes needs to be the number of layers
        qy = np.sum(qxqy[1,:,np.newaxis,:] * aq.eigvec, 2).T
        return np.array((qx[pylayers], qy[pylayers]))
    def disveclayers(self, x, y, pylayers, aq = None):
        '''Returns two arrays of size len(pylayers)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        qxqy = self.disvec(x, y, aq)
        rv = np.sum(qxqy[:,np.newaxis,:] * aq.eigvec, 2)
        return rv[:,pylayers]
    def intflux(self, func, x1, y1, x2, y2, pylayers, ndeg=8, aq=None):
        if aq is None: print 'error, aquifer needs to be given'
        thetaNormOut = np.arctan2(y2 - y1, x2 - x1) - np.pi/2.0
        cosnorm = np.cos(thetaNormOut)
        sinnorm = np.sin(thetaNormOut)
        z1 = x1 + 1j * y1
        z2 = x2 + 1j * y2
        Xleg, wleg = np.polynomial.legendre.leggauss(ndeg)
        z = 0.5 * Xleg * (z2 - z1) + 0.5 * (z1 + z2)
        x = z.real
        y = z.imag
        qtot = 0.0
        for i in range(ndeg):
            qxqy = func(x=x[i], y=y[i], pylayers=pylayers, aq=aq)
            qtot += wleg[i] * (qxqy[0] * cosnorm + qxqy[1] * sinnorm)
        return qtot
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
        self.aq.add_element(self)
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
    def disinf(self, x, y, aq = None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.Nparam, aq.Naq))
        if aq == self.aq:
            qx = np.zeros(aq.Naq)
            qy = np.zeros(aq.Naq)
            rsq = (x - self.xw)**2 + (y - self.yw)**2
            r = np.sqrt(rsq)
            xminxw = x - self.xw
            yminyw = y - self.yw
            if r < self.rw:
                r = self.rw
                rsq = self.rwsq
                xminxw = self.rw
                yminyw = 0.0
            if aq.ltype[0] == 'a':
                qx[0] = -1 / (2 * np.pi) * xminxw / rsq
                qy[0] = -1 / (2 * np.pi) * yminyw / rsq 
                kone = k1(r / aq.lab)
                qx[1:] = -kone * xminxw / (r * aq.lab) / (2 * np.pi)
                qy[1:] = -kone * yminyw / (r * aq.lab) / (2 * np.pi)
            else:
                kone = k1(r / aq.lab)
                qx[:] = -kone * xminxw / (r * aq.lab) / (2 * np.pi)
                qy[:] = -kone * yminyw / (r * aq.lab) / (2 * np.pi)
            rv[0] = self.aq.coef[self.pylayers] * qx
            rv[1] = self.aq.coef[self.pylayers] * qy   
        return rv
    def disinfold(self, x, y, aq = None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rvx = np.zeros((self.Nparam, aq.Naq))
        rvy = np.zeros((self.Nparam, aq.Naq))
        if aq == self.aq:
            qx = np.zeros(aq.Naq)
            qy = np.zeros(aq.Naq)
            rsq = (x - self.xw)**2 + (y - self.yw)**2
            r = np.sqrt(rsq)
            xminxw = x - self.xw
            yminyw = y - self.yw
            if r < self.rw:
                r = self.rw
                rsq = self.rwsq
                xminxw = self.rw
                yminyw = 0.0
            if aq.ltype[0] == 'a':
                qx[0] = -1 / (2 * np.pi) * xminxw / rsq
                qy[0] = -1 / (2 * np.pi) * yminyw / rsq 
                kone = k1(r / aq.lab)
                qx[1:] = -kone * xminxw / (r * aq.lab) / (2 * np.pi)
                qy[1:] = -kone * yminyw / (r * aq.lab) / (2 * np.pi)
            else:
                kone = k1(r / aq.lab)
                qx[:] = -kone * xminxw / (r * aq.lab) / (2 * np.pi)
                qy[:] = -kone * yminyw / (r * aq.lab) / (2 * np.pi)
            rvx[:] = self.aq.coef[self.pylayers] * qx
            rvy[:] = self.aq.coef[self.pylayers] * qy   
        return rvx, rvy
    
class HeadWell(WellBase, HeadEquation):
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
        
class Constant(Element, HeadEquation):
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
        self.aq.add_element(self)
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
    def disinf(self, x, y, aq = None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, 1, aq.Naq))
        return rv
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
#class Constant2(Element, HeadEquationSpecial):
#    def __init__(self, model, xr=0, yr=0, ls=None, hr=0.0, \
#                 layer=0, label=None):
#        Element.__init__(self, model, Nparam=1, Nunknowns=1, layers=layer,\
#                         name='Constant', label=label)
#        self.Nunknowns = 1
#        self.xr = xr
#        self.yr = yr
#        self.hr = hr
#        self.model.add_element(self)
#    def initialize(self):
#        self.aq = self.model.aq.find_aquifer_data(self.xr, self.yr)
#        self.aq.add_element(self)
#        self.Ncp = 1
#        self.xc = np.array([self.xr])
#        self.yc = np.array([self.yr])
#        self.pc = self.hr * self.aq.T[self.pylayers]
#        self.parameters = np.empty((1, 1))
#    def potinf(self, x, y, aq = None):
#        '''Can be called with only one x,y value'''
#        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
#        rv = np.zeros((1, aq.Naq))
#        if aq == self.aq:
#            rv[0,0] = 1
#        return rv
#    def disinf(self, x, y, aq = None):
#        '''Can be called with only one x,y value'''
#        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
#        rv = np.zeros((2, 1, aq.Naq))
#        return rv
#    def setparams(self, sol):
#        self.parameters[:,0] = sol
        
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
    def disinf(self, x, y, aq = None):
        '''Can be called with only one x,y value'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rvx = np.zeros((2, self.Nparam, aq.Naq))
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
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 Qls= 0.0, layers=0, order=0, name='LineSinkHoBase', \
                 label=None, addtomodel=True, aq=None, zcinout=None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers,\
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
        return self.name + ' from ' + str((self.x1, self.y1)) +' to '+str((self.x2, self.y2))
    def initialize(self):
        self.Ncp = self.order + 1
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
                                   aq.Naq + 1, aq.lab, i, pot)  # Need to specify aq.Naq + 1 (bug in besselaes)
                    rv[i] = self.aq.coef[self.pylayers] * pot
            rv.shape = (self.Nparam, aq.Naq)
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
        
class IntfluxLineSinkHo(LineSinkHoBase, DisvecDiffEquation2):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 order=0, label=None, addtomodel=True, aq=None):
        self.storeinput(inspect.currentframe())
        LineSinkHoBase.__init__(self, model, x1, y1, x2, y2, Qls=0, \
                 layers=range(model.aq.Naq), order= order, name='IntfluxLineSinkHo', label=label, \
                 addtomodel=addtomodel, aq=aq)
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineSinkHoBase.initialize(self)
        self.xcin, self.ycin = controlpoints(self.Ncp-1, self.z1, self.z2, eps=1e-6, include_ends=True)
        self.xcout, self.ycout = controlpoints(self.Ncp-1, self.z1, self.z2, eps=-1e-6, include_ends=True)
        self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], self.ycin[0])
        self.aqout = self.model.aq.find_aquifer_data(self.xcout[0], self.ycout[0])
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
        
class HeadDiffLineSink(LineSinkHoBase, HeadDiffEquation):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 order=0, layers=0, label=None, addtomodel=True, aq=None):
        self.storeinput(inspect.currentframe())
        LineSinkHoBase.__init__(self, model, x1, y1, x2, y2, Qls=0, \
                 layers=range(model.aq.Naq), order=order, name='HeadDiffLineSink', label=label, \
                 addtomodel=addtomodel, aq=aq)
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineSinkHoBase.initialize(self)
        self.aqin = self.model.aq.find_aquifer_data(self.xcin[0], self.ycin[0])
        self.aqout = self.model.aq.find_aquifer_data(self.xcout[0], self.ycout[0])
    def setparams(self, sol):
        self.parameters[:,0] = sol

class HeadDiffLineSinkString(LineSinkStringBase, HeadDiffEquation):
    def __init__(self, model, xy=[(-1,0), (1,0)], closed=True, \
                 order=0, label=None, aq=None):
        self.storeinput(inspect.currentframe())
        LineSinkStringBase.__init__(self, model, xy, closed=closed, layers=layers, order=order, \
                                    name='HeadDiffLineSinkString', label=label, aq=aq)
        self.model.add_element(self)
    def initialize(self):
        LineSinkStringBase.initialize(self)
        self.aq.add_element(self)
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class DisvecDiffLineSinkString(LineSinkStringBase, DisvecDiffEquation):
    def __init__(self, model, xy=[(-1,0), (1,0)], closed=True,\
                 layers=0, order=0, label=None, aq=None):
        self.storeinput(inspect.currentframe())
        LineSinkStringBase.__init__(self, model, xy, closed=closed, layers=layers, order=order, \
                                    name='DisvecDiffLineSinkString', label=label, aq=aq)
        self.model.add_element(self)
    def initialize(self):
        LineSinkStringBase.initialize(self)
        self.aq.add_element(self)
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
def PolygonalInhomogeneity(model, xy, order=0):
        x1, y1 = xy[0]
        x2, y2 = xy[1]
        z1 = x1 + y1 * 1j
        z2 = x2 + y2 * 1j
        Zin = 1e-6j
        Zout = -1e-6j
        zin = 0.5 * (z2 - z1) * Zin + 0.5 * (z1 + z2)
        zout = 0.5 * (z2 - z1) * Zout + 0.5 * (z1 + z2)
        aqin = model.aq.find_aquifer_data(zin.real, zin.imag)
        aqout = model.aq.find_aquifer_data(zout.real, zout.imag)
        # 
        #ls_in = HeadDiffLineSinkString(model, xy=xy, closed=True, layers=range(model.aq.Naq), \
        #                               order=order, label=None, aq=aqin)
        #ls_out = DisvecDiffLineSinkString(model, xy=xy, closed=True, layers=range(model.aq.Naq), \
        #                                  order=order, label=None, aq=aqout)
        lsin = []
        lsout = []
        xy.append(xy[0])
        x,y = zip(*xy)
        for i in range(len(xy)-1):
            ls = HeadDiffLineSink(model, x1=x[i], y1=y[i], x2=x[i+1], y2=y[i+1], \
                                  order=order, label=None, addtomodel=True, aq=aqin)
            lsin.append(ls)
            ls = IntfluxLineSinkHo(model, x1=x[i], y1=y[i], x2=x[i+1], y2=y[i+1], \
                 order=order, label=None, addtomodel=True, aq=aqout)
            lsout.append(ls)
        return lsin, lsout

def PolygonalInhomogeneity_OK(model, xy, order=0):
        # Uses HeadDiffLineSinkString for inside
        # Works well for lakes
        x1, y1 = xy[0]
        x2, y2 = xy[1]
        z1 = x1 + y1 * 1j
        z2 = x2 + y2 * 1j
        Zin = 1e-6j
        Zout = -1e-6j
        zin = 0.5 * (z2 - z1) * Zin + 0.5 * (z1 + z2)
        zout = 0.5 * (z2 - z1) * Zout + 0.5 * (z1 + z2)
        aqin = model.aq.find_aquifer_data(zin.real, zin.imag)
        aqout = model.aq.find_aquifer_data(zout.real, zout.imag)
        # 
        ls_in = HeadDiffLineSinkString(model, xy=xy, closed=True, layers=range(model.aq.Naq), \
                                       order=order, label=None, aq=aqin)
        #ls_out = DisvecDiffLineSinkString(model, xy=xy, closed=True, layers=range(model.aq.Naq), \
        #                                  order=order, label=None, aq=aqout)
        lsout_list = []
        xy.append(xy[0])
        x,y = zip(*xy)
        for i in range(len(xy)-1):
            ls = IntfluxLineSinkHo(model, x1=x[i], y1=y[i], x2=x[i+1], y2=y[i+1], \
                 order=order, label=None, addtomodel=True, aq=aqout)
            lsout_list.append(ls)
        return ls_in, lsout_list
        
class LineSinkStringIn(LineSinkStringBase, HeadDiffEquation):
    def __init__(self, model, xy=[(-1,0), (1,0)], order=0, label=None, aq=None):
        LineSinkStringBase.__init__(self, model, layers=range(aq.Naq), order=order, \
                                    name='LineSinkStringIn', label=label, aq=aq)
        xy = np.atleast_2d(xy).astype('d')
        self.x, self.y = xy[:,0], xy[:,1]
        self.Nls = len(self.x) - 1
        for i in range(self.Nls):
            self.lslist.append(LineSinkHoBase(model, \
                x1=self.x[i], y1=self.y[i], x2=self.x[i+1], y2=self.y[i+1], \
                Q=0.0, layers=range(aq.Naq), order=order, label=label, addtomodel=False, aq=aq))
        self.model.add_element(self)
    def initialize(self):
        LineSinkStringBase.initialize(self)
        self.aq.add_element(self)
        self.pc = self.lslist[0].pc  # same for all control points
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class LineDoubletHoBase(Element):
    def __init__(self, model, x1=-1, y1=0, x2=1, y2=0, \
                 delp = 0.0, layers=0, order=0, name='LineDoubletHoBase', \
                 label=None, addtomodel=True, aq=None, zcinout=None):
        Element.__init__(self, model, Nparam=1, Nunknowns=0, layers=layers,\
                         name=name, label=label)
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.x2 = float(x2)
        self.y2 = float(y2)
        self.delp = np.atleast_1d(delp)
        self.order = order
        self.Nparam = self.Nlayers * (self.order + 1)
        self.addtomodel = addtomodel
        if addtomodel: self.model.add_element(self)
        self.aq = aq
        self.zcinout = zcinout
    def __repr__(self):
        return self.name + ' from ' + str((self.x1, self.y1)) +' to '+str((self.x2, self.y2))
    def initialize(self):
        self.Ncp = self.order + 1
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
        if self.aq is None:
            self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        if self.addtomodel:
            self.aq.add_element(self)
        self.parameters = np.zeros((self.Nparam, 1))
        self.parameters[:,0] = self.delp  # Not sure that needs to be here
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
                    potbesldho(x, y, self.x1, self.y1, self.x2, self.y2, \
                               aq.Naq, aq.zeropluslab, i, pot)  # Call FORTRAN extension
                    rv[i] = self.aq.coef[self.pylayers] * pot
            else:
                raise Exception('LineDoubletHoBase not implemented for semi-confined aquifer')
                for i in range(self.order+1):
                    potbesonlylsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                                   aq.Naq, aq.lab, 0, pot)
                    rv[i] = self.aq.coef[self.pylayers] * pot
            rv.shape = (self.Nparam, aq.Naq)
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
                    disbesldho(x, y, self.x1, self.y1, self.x2, self.y2, \
                               aq.Naq, aq.zeropluslab, i, qx, qy)  # Call FORTRAN extension
                    rv[0, i] = self.aq.coef[self.pylayers] * qx
                    rv[1, i] = self.aq.coef[self.pylayers] * qy
            else:
                raise Exception('Must overload AquiferData.isinside()')
                for i in range(self.order+1):
                    disbesonlylsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                                   aq.Naq, aq.lab, 0, qx, qy)
                    rv[0, i] = self.aq.coef[self.pylayers] * qx
                    rv[1, i] = self.aq.coef[self.pylayers] * qy
            rv.shape = (2, self.Nparam, aq.Naq)
        return rv
    
class ImpWall(LineDoubletHoBase, DisvecEquation):
    def __init__(self, model, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
                 order = 0, layers = 0, label = None, addtomodel = True):
        self.storeinput(inspect.currentframe())
        LineDoubletHoBase.__init__(self, model, x1, y1, x2, y2, delp = 0, \
                 layers = layers, order = order, name = 'ImpWall', label = label, \
                 addtomodel = addtomodel)
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineDoubletHoBase.initialize(self)
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
def inside_outside_polygon(x, y, eps = 1e-6):
    '''Returns closed polygon just on inside and just on outside of closed polygon'''
    x = np.hstack((x,x[0]))
    y = np.hstack((y,y[0]))
    N = len(x) - 1
    znew = np.zeros((N,2), 'D')
    for i in range(N):
        z1 = x[i] + y[i] * 1j
        z2 = x[i+1] + y[i+1] * 1j
        Z1 = -1.0 + eps * 1j
        Z2 = 1.0 + eps * 1j
        znew[i,0] = Z1 * 0.5 * (z2 - z1) + 0.5 * (z1 + z2)
        znew[i,1] = Z2 * 0.5 * (z2 - z1) + 0.5 * (z1 + z2)
    znew = np.vstack((znew, znew[0]))
    zcorners = np.zeros((N,2), 'D')
    for i in range(N):
        zcorners[i] = intersect(znew[i,0], znew[i,1], znew[i+1,0], znew[i+1,1])
    return zcorners
        
def intersect(z1,z2,z3,z4):
    Za = (z3 - 0.5 * (z1 + z2)) / (0.5 * (z2 - z1))
    Zb = (z4 - 0.5 * (z1 + z2)) / (0.5 * (z2 - z1))
    Ya = Za.imag
    Yb = Zb.imag
    a = -Ya / (Yb - Ya)
    Z = Za + a * (Zb - Za)
    z = 0.5 * (z2 - z1) * Z + 0.5 * (z1 + z2)
    return z
    
def intersect_old(x, y):
    x1, x2, x3, x4 = x
    y1, y2, y3, y4 = y
    A = np.array([[(x2 - x1), -(x4 - x3)], \
                  [(y2 - y1), -(y4 - y3)]])
    print np.linalg.det(A)
    rhs = np.array([(x3-x1),(y3-y1)])
    a, b = np.linalg.solve(A,rhs)
    print x1 + a * (x2 - x1)
    print y1 + a * (y2 - y1)
    print x3 + b * (x4 - x3)
    print y3 + b * (y4 - y3)
    return a,b
    
def numder(f, x, y, d):
    qx = (f(x-d, y) - f(x+d, y)) / (2*d)
    qy = (f(x, y-d) - f(x, y+d)) / (2*d)
    return qx, qy
    
def controlpoints(Ncp, z1, z2, eps=0, include_ends=False):
    #thetacp = np.arange(np.pi, 0, -np.pi/self.Ncp) - 0.5 * np.pi/self.Ncp
    # The following works MUCH better for a uniform head along the line
    thetacp = np.linspace(np.pi, 0, Ncp+2)[1:-1]
    if include_ends:
        Zcp = np.zeros(Ncp+2, 'D')
        Zcp[0] = -1
        Zcp[-1] = 1
        Zcp[1:-1] = np.cos(thetacp)
    else:
        #thetacp = np.arange(np.pi, 0, -np.pi/Ncp) - 0.5 * np.pi/Ncp
        Zcp = np.zeros(Ncp, 'D')
        Zcp.real = np.cos(thetacp)
    Zcp.imag = eps  # control point just on positive site (this is handy later on)
    zcp = Zcp * (z2 - z1) / 2.0 + 0.5 * (z1 + z2)
    return zcp.real, zcp.imag

class ImpLineSink(LineSinkHoBase, DisvecEquationOut):
    def __init__(self, model, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
                 order = 0, layers = 0, label=None, addtomodel=True, zcinout=None):
        self.storeinput(inspect.currentframe())
        LineSinkHoBase.__init__(self, model, x1, y1, x2, y2, Qls = 0, \
                 layers = layers, order = order, name = 'ImpLineSink', label = label, \
                 addtomodel = addtomodel, zcinout=zcinout)
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineSinkHoBase.initialize(self)
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class WellPsi(WellBase):
    def __init__(self, model, xw = 0, yw = 0, Qw = 10.0, rw = 0.1, \
                 layers = 0, label = None):
        self.storeinput(inspect.currentframe())
        WellBase.__init__(self, model, xw, yw, Qw, rw,\
                          layers = layers, name = 'PsiWell', label = label)
    def initialize(self):
        WellBase.initialize(self)
        self.zw = self.xw + 1j  * self.yw
    def psiinf(self, x, y, aq=None):
        z = x + 1j * y
        rv = 1.0 / (2 * np.pi) * np.log(z - self.zw)
        return np.array([[rv.imag]])
    def psi(self, x, y, aq = None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        return np.sum(self.parameters * self.psiinf(x, y, aq), 0)
    def psiinflayers(self, x, y, pylayers, aq = None):
        '''Returns array of size (len(pylayers),Nparam)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        psi = self.psiinf(x, y, aq)  # Nparam rows, Naq cols
        rv = np.sum(psi[:,np.newaxis,:] * aq.eigvec, 2).T  # Transopose as the first axes needs to be the number of layers
        return rv[pylayers,:]
    def psilayers(self, x, y, pylayers, aq = None):
        '''Returns array of size len(pylayers)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        psi = np.sum(self.psi(x, y, aq) * aq.eigvec, 1 )
        return psi[pylayers]
    
class LineDoubletPsi(LineDoubletHoBase):
    def __init__(self, model, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
                 order = 0, layers = 0, label = None, addtomodel = True):
        self.storeinput(inspect.currentframe())
        LineDoubletHoBase.__init__(self, model, x1, y1, x2, y2, delp = 0, \
                 layers = layers, order = order, name = 'ImpWall', label = label, \
                 addtomodel = addtomodel)
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineDoubletHoBase.initialize(self)
    def psiinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nparam, aq.Naq))
        if aq == self.aq:
            rv.shape = (self.order+1, self.Nlayers, aq.Naq)
            pot = np.zeros(self.order+1)
            psi = np.zeros(self.order+1)
            if aq.ltype[0] == 'a':
                for i in range(self.order+1): # This should be done inside FORTRAN extension
                    omegalapldho(x, y, self.x1, self.y1, self.x2, self.y2, \
                               self.order, pot, psi)  # Call FORTRAN extension
                    rv = self.aq.coef[self.pylayers] * psi
            else:
                raise Exception('LineDoubletHoBase not implemented for semi-confined aquifer')
            rv.shape = (self.Nparam, aq.Naq)
        return rv
    def psi(self, x, y, aq = None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        return np.sum(self.parameters * self.psiinf(x, y, aq), 0)
    def psiinflayers(self, x, y, pylayers, aq = None):
        '''Returns array of size (len(pylayers),Nparam)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        psi = self.psiinf(x, y, aq)  # Nparam rows, Naq cols
        rv = np.sum(psi[:,np.newaxis,:] * aq.eigvec, 2).T  # Transopose as the first axes needs to be the number of layers
        return rv[pylayers,:]
    def psilayers(self, x, y, pylayers, aq = None):
        '''Returns array of size len(pylayers)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        psi = np.sum(self.psi(x, y, aq) * aq.eigvec, 1 )
        return psi[pylayers]
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class LineSinkPsi(LineSinkHoBase):
    def __init__(self, model, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
                 order = 0, layers = 0, label = None, addtomodel=True, zcinout=None):
        self.storeinput(inspect.currentframe())
        LineSinkHoBase.__init__(self, model, x1, y1, x2, y2, Qls = 0, \
                 layers = layers, order = order, name = 'ImpWall', label = label, \
                 addtomodel = addtomodel, zcinout=zcinout)
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineSinkHoBase.initialize(self)
    def psiinf(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.Nparam, aq.Naq))
        if aq == self.aq:
            rv.shape = (self.order+1, self.Nlayers, aq.Naq)
            pot = np.zeros(self.order+1)
            psi = np.zeros(self.order+1)
            if aq.ltype[0] == 'a':
                for i in range(self.order+1): # This should be done inside FORTRAN extension
                    omegalaplsho(x, y, self.x1, self.y1, self.x2, self.y2, \
                               self.order, pot, psi)  # Call FORTRAN extension
                    rv = self.aq.coef[self.pylayers] * psi
            else:
                raise Exception('LineSinkHoBase not implemented for semi-confined aquifer')
            rv.shape = (self.Nparam, aq.Naq)
        return rv
    def psi(self, x, y, aq = None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        return np.sum(self.parameters * self.psiinf(x, y, aq), 0)
    def psiinflayers(self, x, y, pylayers, aq = None):
        '''Returns array of size (len(pylayers),Nparam)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        psi = self.psiinf(x, y, aq)  # Nparam rows, Naq cols
        rv = np.sum(psi[:,np.newaxis,:] * aq.eigvec, 2).T  # Transopose as the first axes needs to be the number of layers
        return rv[pylayers,:]
    def psilayers(self, x, y, pylayers, aq = None):
        '''Returns array of size len(pylayers)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        psi = np.sum(self.psi(x, y, aq) * aq.eigvec, 1 )
        return psi[pylayers]
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class ImpLineSink(LineSinkPsi, DisvecEquationOut):
    def __init__(self, model, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
                 order = 0, layers = 0, label=None, addtomodel=True, zcinout=None):
        self.storeinput(inspect.currentframe())
        LineSinkPsi.__init__(self, model, x1, y1, x2, y2, \
                 layers = layers, order = order, label = label, \
                 addtomodel = addtomodel, zcinout=zcinout)
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineSinkPsi.initialize(self)
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
class PsiDiffEquation:
    def equation(self):
        '''Mix-in class that returns matrix rows for head-specified conditions. (really written as constant potential element)
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    e.psiinflayers(self.xcout[icp], self.ycout[icp], self.pylayers) - \
                    e.psiinflayers(self.xcout[icp+1], self.ycout[icp+1], self.pylayers)
                    ieq += e.Nunknowns
                else:
                    rhs[istart:istart+self.Nlayers] -= \
                    e.psilayers(self.xcout[icp], self.ycout[icp], self.pylayers) - \
                    e.psilayers(self.xcout[icp+1], self.ycout[icp+1], self.pylayers)
        return mat, rhs
        
class ImpLineSink2(LineSinkPsi, PsiDiffEquation):
    def __init__(self, model, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
                 order = 0, layers = 0, label=None, addtomodel=True, zcinout=None):
        self.storeinput(inspect.currentframe())
        LineSinkPsi.__init__(self, model, x1, y1, x2, y2, \
                 layers = layers, order = order, label = label, \
                 addtomodel = addtomodel, zcinout=zcinout)
        self.Nunknowns = self.Nparam
    def initialize(self):
        LineSinkPsi.initialize(self)
        self.xcout, self.ycout = controlpoints(self.Ncp-1, self.z1, self.z2, eps=-1e-6, include_ends=True)
    def setparams(self, sol):
        self.parameters[:,0] = sol
        
def psi(model, x, y, aq = None):
        # Only for testing
        if aq is None: aq = model.aq.find_aquifer_data(x,y)
        psi = np.zeros(aq.Naq)
        for e in model.elementlist:
            psi += e.psi(x, y, aq)
        return psi
    

#ml = ModelMaq(kaq = [1,2,3], z = [5,4,3,2,1,0], c = 1000)
#xy = [(-5,0), (6,0), (4,8), (-5,8)]
#PolygonInhomMaq(ml, xy=xy, kaq = [10,2,30], z = [5,4,3,2,1,0], c = 1000)
#lsin, lsout = PolygonalInhomogeneity(ml, xy=xy, order=5)
#pi = PolygonInhomMaq(ml, xy = [(0,0), (10,0), (5,5)], kaq = 10, z = [5,4,3,2,1,0], c = 5000, npor = 0.3, top = 'conf')
#w = WellBase(ml, xw = 0, yw = -10, Qw = 100, layers = 0)
#w2 = WellBase(ml, xw = 0, yw = 0, Qw = [100,200], rw = 0.1, layers = [0,1])
#w2 = HeadWell(ml, 0, -10, hw = [10.0, 12.0], rw = 0.1, layers = [0,1])
#ls = LineSinkBase(ml, x1 = -20, y1 = -10, x2 = 40, y2 = 20, \
#                 Qls = 100.0, layers = 0)
#ls = HeadLineSink(ml, x1 = -20, y1 = -10, x2 = 40, y2 = 20, \
#                  hls = 12.0, layers = [0,1])
#rf = Constant(ml, -10, 20, 20, layer = 0)
#ls1 = LineSinkHoBase(ml, x1 = -1, y1 = 0, x2 = 1, y2 = 0, \
#                 Qls = 100.0, layers = [0,1], order = 3)
#ls1 = HeadLineSinkHo(ml, x1 = -20, y1 = -10, x2 = 40, y2 = 20, \
#                    hls = 12.0, layers = [0], order = 4)
#ls2 = HeadLineSinkHo(ml, x1 = -20, y1 = -30, x2 = 40, y2 = -10, \
#                    hls = 14.0, layers = [1], order = 4)
#lsstring = HeadLineSinkString(ml, xy = [(-10,-10), (0,0), (10,0), (10,10)], hls = 7, layers = [0,1], order = 5)
#ld = ImpWall(ml, -5, 0, 5, 0, order=3, layers=[0])
#ml.initialize()
#ml.solve()

#ml = ModelMaq(kaq = 1, z = [10,0])
#w = WellPsi(ml, 0, -10, 100)
#eps = 1e-6
#zcinout = [-5+eps*1j, 5+eps*1j, -5-eps*1j, 5-eps*1j]
#lsp = ImpLineSink2(ml, -5, 0, 5, 0, order = 3, zcinout=zcinout)
##lsp = ImpLineSink(ml, -5, 0, 5, 0, order = 3, zcinout=zcinout)
#ml.solve()

#def psi(x,y,layer):
#    return lsp.psiinf(x,y)[layer]
#psivec = np.vectorize(psi)

def intflux(func, x1, y1, x2, y2, ndeg=8, aq=None):
    thetaNormOut = np.arctan2(y2 - y1, x2 - x1) - np.pi/2.0
    cosnorm = np.cos(thetaNormOut)
    sinnorm = np.sin(thetaNormOut)
    z1 = x1 + 1j * y1
    z2 = x2 + 1j * y2
    Xleg, wleg = np.polynomial.legendre.leggauss(ndeg)
    z = 0.5 * Xleg * (z2 - z1) + 0.5 * (z1 + z2)
    x = z.real
    y = z.imag
    qtot = 0.0
    for i in range(ndeg):
        qxqy = func(x=x[i], y=y[i], aq=aq)
        qtot += wleg[i] * (qxqy[0] * cosnorm + qxqy[1] * sinnorm)
    return qtot * np.sqrt((x2 - x1) **2 + (y2 - y1) **2) / 2.0

#ml = ModelMaq(kaq = [1], z = [1,0])
#xy = [(-5,0), (5,0), (5,8), (-5,8)]
#p = PolygonInhomMaq(ml, xy=xy, kaq = [0.1], z = [1,0])
#lsin, lsout = PolygonalInhomogeneity(ml, xy=xy, order=7)
###pi = PolygonInhomMaq(ml, xy = [(0,0), (10,0), (5,5)], kaq = 10, z = [5,4,3,2,1,0], c = 5000, npor = 0.3, top = 'conf')
#w = WellBase(ml, xw = 0, yw = -10, Qw = 100, layers = 0)
#rf = Constant(ml, xr = 0, yr = -100, hr=20)
##ls1 = IntfluxLineSinkHo(ml, -2, 0, 2, 0, order=4, aq=ml.aq)
#ml.solve()
#
#ml = ModelMaq(kaq = [1], z = [1,0])
#xy = [(-5,0), (5,0), (5,8), (-5,8)]
#p = PolygonInhomMaq(ml, xy=xy, kaq = [0.2], z = [1,0])
#lsin, lsout = PolygonalInhomogeneity(ml, xy=xy, order=6)
###pi = PolygonInhomMaq(ml, xy = [(0,0), (10,0), (5,5)], kaq = 10, z = [5,4,3,2,1,0], c = 5000, npor = 0.3, top = 'conf')
#w = WellBase(ml, xw = 0, yw = -10, Qw = 100, layers = 0)
#rf = Constant(ml, xr = 0, yr = -100, hr=20)
##ls1 = IntfluxLineSinkHo(ml, -2, 0, 2, 0, order=4, aq=ml.aq)
#ml.solve()

#ml = ModelMaq(kaq = [1,2], z = [10,5,4,0], c=20)
#xy = [(-5,0), (5,0), (5,8), (-5,8)]
#p = PolygonInhomMaq(ml, xy=xy, kaq = [0.2,8], z = [10,5,4,0], c=20)
#lsin, lsout = PolygonalInhomogeneity(ml, xy=xy, order=7)
###pi = PolygonInhomMaq(ml, xy = [(0,0), (10,0), (5,5)], kaq = 10, z = [5,4,3,2,1,0], c = 5000, npor = 0.3, top = 'conf')
#w = WellBase(ml, xw = 0, yw = -10, Qw = 100, layers = 0)
#rf = Constant(ml, xr = 0, yr = -100, hr=20)
##ls1 = IntfluxLineSinkHo(ml, -2, 0, 2, 0, order=4, aq=ml.aq)
#ml.solve()

ml = ModelMaq(kaq = [1], z = [1,0])
xy = [(-3,0), (5,0), (5,8), (-5,8)]
p = PolygonInhomMaq(ml, xy=xy, kaq = [0.1], z = [2,1,0], c = 100, top = 'semi')
lsin, lsout = PolygonalInhomogeneity(ml, xy=xy, order=6)
##pi = PolygonInhomMaq(ml, xy = [(0,0), (10,0), (5,5)], kaq = 10, z = [5,4,3,2,1,0], c = 5000, npor = 0.3, top = 'conf')
w = WellBase(ml, xw = 0, yw = -10, Qw = 100, layers = 0)
rf = Constant(ml, xr = 0, yr = -100, hr=2)
#ls1 = IntfluxLineSinkHo(ml, -2, 0, 2, 0, order=4, aq=ml.aq)
ml.solve()
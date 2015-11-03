import numpy as np
import inspect # Used for storing the input
from aquifer import Aquifer
from aquifer_parameters import param_maq

class ModelBase:
    def __init__(self, kaq, Haq, c, z, npor, ltype):
        # All input variables are numpy arrays
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
    def potential(self, x, y, aq=None):
        if aq is None: aq = self.aq.find_aquifer_data(x,y)
        pot = np.zeros(aq.Naq)
        for e in aq.elementlist:
            pot += e.potential(x, y, aq)
        rv = np.sum(pot * aq.eigvec, 1)
        return rv
    def disvec(self, x, y, aq=None):
        if aq is None: aq = self.aq.find_aquifer_data(x,y)
        rv = np.zeros((2, aq.Naq))
        for e in aq.elementlist:
            rv += e.disvec(x, y, aq)
        rv = np.sum(rv[:,np.newaxis,:] * aq.eigvec, 2)
        return rv
    def head(self, x, y, layers=None, aq=None):
        if aq is None: aq = self.aq.find_aquifer_data(x,y)
        rv = self.potential(x, y, aq) / aq.T
        if layers is None:
            return rv
        else:
            return rv[layers]
    def headgrid(self, xg, yg, layers=None, printrow=False):
        '''Returns h[Nlayers,ny,nx]. If layers is None, all layers are returned'''
        nx, ny = len(xg), len(yg)
        if layers is None:
            Nlayers = self.aq.find_aquifer_data(xg[0],yg[0]).Naq
        else:
            Nlayers = len(np.atleast_1d(layers))
        h = np.empty((Nlayers, ny, nx))
        for j in range(ny):
            if printrow:
                print str(j)+' ',
                sys.stdout.flush()  # Can be replaced with print with flush in Python 3.3
            for i in range(nx):
                h[:,j,i] = self.head(xg[i], yg[j], layers)
        if printrow:
            print
            sys.stdout.flush()  # Can be replaced with print with flush in Python 3.3
        return h
    def headgrid2(self, x1, x2, nx, y1, y2, ny, layers=None, printrow=False):
        '''Returns h[Nlayers,ny,nx]. If layers is None, all layers are returned'''
        xg, yg = np.linspace(x1,x2,nx), np.linspace(y1,y2,ny)
        return self.headgrid(xg, yg, layers=layers, printrow=printrow)
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
    def solve(self, printmat=0, sendback=0, silent=False):
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
        mat = np.empty((self.Neq, self.Neq))
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
        if silent is False:
            print 'solution complete'
        elif (silent == 'dot') or (silent == '.'):
            print '.',
            sys.stdout.flush()  # Can be replaced with print with flush in Python 3.3
        if sendback:
            return sol
        return
        
class ModelMaq(ModelBase):
    def __init__(self, kaq = 1, z = [1,0], c = [], npor = 0.3, top = 'conf'):
        self.storeinput(inspect.currentframe())
        kaq, Haq, c, npor, ltype = param_maq(kaq, z, c, npor, top)
        ModelBase.__init__(self, kaq, Haq, c, z, npor, ltype)
        self.name = 'ModelMaq'
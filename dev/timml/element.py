import numpy as np
import inspect # Used for storing the input

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
        self.inhomelement = False  # elements used as part of an inhom boundary are tagged
        if self.label is not None:
            assert self.label not in self.model.elementdict.keys(),\
            "TTim error: label " + self.label + " already exists"
    def initialize(self):
        # must be overloaded
        pass
    def potinf(self, x, y, aq=None):
        '''Returns array of size (Nparam, Naq)'''
        raise Exception('Must overload Element.potinf()')
    def potential(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        return np.sum(self.parameters * self.potinf(x, y, aq), 0)
    def potinflayers(self, x, y, pylayers, aq=None):
        '''Returns array of size (len(pylayers),Nparam)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        pot = self.potinf(x, y, aq)  # Nparam rows, Naq cols
        rv = np.sum(pot[:,np.newaxis,:] * aq.eigvec, 2).T  # Transopose as the first axes needs to be the number of layers
        return rv[pylayers,:]
    def potentiallayers(self, x, y, pylayers, aq=None):
        '''Returns array of size len(pylayers)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        pot = np.sum(self.potential(x, y, aq) * aq.eigvec, 1 )
        return pot[pylayers]
    def disinf(self, x, y, aq=None):
        '''Returns array of size (2, Nparam, Naq)'''
        raise Exception('Must overload Element.disinf()')
    def disvec(self, x, y, aq=None):
        '''Returns array of size (2, Nparam, Naq)'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        return np.sum(self.parameters * self.disinf(x, y, aq), 1)
    def disinflayers(self, x, y, pylayers, aq=None):
        '''Returns two arrays of size (len(pylayers),Nparam)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        qxqy = self.disinf(x, y, aq)  # Nparam rows, Naq cols
        qx = np.sum(qxqy[0,:,np.newaxis,:] * aq.eigvec, 2).T  # Transpose as the first axes needs to be the number of layers
        qy = np.sum(qxqy[1,:,np.newaxis,:] * aq.eigvec, 2).T
        return np.array((qx[pylayers], qy[pylayers]))
    def disveclayers(self, x, y, pylayers, aq=None):
        '''Returns two arrays of size len(pylayers)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        qxqy = self.disvec(x, y, aq)
        rv = np.sum(qxqy[:,np.newaxis,:] * aq.eigvec, 2)
        return rv[:,pylayers]
    def intpot(self, func, x1, y1, x2, y2, pylayers, aq=None):
        if aq is None: print 'error, aquifer needs to be given'
        z1 = x1 + 1j * y1
        z2 = x2 + 1j * y2
        z = 0.5 * self.Xleg * (z2 - z1) + 0.5 * (z1 + z2)
        x = z.real
        y = z.imag
        pot = 0.0
        for i in range(self.ndeg):
            pot += self.wleg[i] * func(x=x[i], y=y[i], pylayers=pylayers, aq=aq)
        return pot
    def intflux(self, func, x1, y1, x2, y2, pylayers, aq=None):
        if aq is None: print 'error, aquifer needs to be given'
        thetaNormOut = np.arctan2(y2 - y1, x2 - x1) - np.pi/2.0
        cosnorm = np.cos(thetaNormOut)
        sinnorm = np.sin(thetaNormOut)
        z1 = x1 + 1j * y1
        z2 = x2 + 1j * y2
        z = 0.5 * self.Xleg * (z2 - z1) + 0.5 * (z1 + z2)
        x = z.real
        y = z.imag
        qtot = 0.0
        for i in range(self.ndeg):
            qxqy = func(x=x[i], y=y[i], pylayers=pylayers, aq=aq)
            qtot += self.wleg[i] * (qxqy[0] * cosnorm + qxqy[1] * sinnorm)
        return qtot
    def headinside(self):
        print 'headinside not implemented for this element'
    def setparams(self, sol):
        raise Exception('Must overload Element.setparams()')
    def storeinput(self,frame):
        self.inputargs, _, _, self.inputvalues = inspect.getargvalues(frame)
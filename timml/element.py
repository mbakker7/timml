import numpy as np
import inspect # Used for storing the input


class Element:
    def __init__(self, model, nparam, nunknowns, layers, name, label):
        self.model = model
        self.aq = None # Set in the initialization function
        self.nparam = nparam
        self.nunknowns = nunknowns
        self.layers = np.atleast_1d(layers)
        self.nlayers = len(self.layers)
        self.ncp = 0
        self.name = name
        self.label = label
        self.inhomelement = False  # elements used as part of an inhom boundary are tagged
        if self.label is not None:
            assert self.label not in list(self.model.elementdict.keys()), \
                "timml error: label " + self.label + " already exists"

    def initialize(self):
        # must be overloaded
        pass

    def potinf(self, x, y, aq=None):
        '''Returns array of size (nparam, naq)'''
        raise Exception('Must overload Element.potinf()')

    def potential(self, x, y, aq=None):
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        return np.sum(self.parameters * self.potinf(x, y, aq), 0)

    def potinflayers(self, x, y, layers, aq=None):
        '''Returns array of size (len(layers),nparam)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        pot = self.potinf(x, y, aq)  # nparam rows, naq cols
        rv = np.sum(pot[:,np.newaxis,:] * aq.eigvec, 2).T  # Transopose as the first axes needs to be the number of layers
        return rv[layers,:]

    def potentiallayers(self, x, y, layers, aq=None):
        '''Returns array of size len(layers)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        pot = np.sum(self.potential(x, y, aq) * aq.eigvec, 1 )
        return pot[layers]

    def disvecinf(self, x, y, aq=None):
        '''Returns array of size (2, nparam, naq)'''
        raise Exception('Must overload Element.disvecinf()')

    def disvec(self, x, y, aq=None):
        '''Returns array of size (2, nparam, naq)'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        return np.sum(self.parameters * self.disvecinf(x, y, aq), 1)

    def disvecinflayers(self, x, y, layers, aq=None):
        '''Returns two arrays of size (len(layers),nparam)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        qxqy = self.disvecinf(x, y, aq)  # nparam rows, naq cols
        qx = np.sum(qxqy[0,:,np.newaxis,:] * aq.eigvec, 2).T  # Transpose as the first axes needs to be the number of layers
        qy = np.sum(qxqy[1,:,np.newaxis,:] * aq.eigvec, 2).T
        return np.array((qx[layers], qy[layers]))

    def disveclayers(self, x, y, layers, aq=None):
        '''Returns two arrays of size len(layers)
        only used in building equations'''
        if aq is None: aq = self.model.aq.find_aquifer_data(x, y)
        qxqy = self.disvec(x, y, aq)
        rv = np.sum(qxqy[:,np.newaxis,:] * aq.eigvec, 2)
        return rv[:,layers]

    def intpot(self, func, x1, y1, x2, y2, layers, aq=None):
        if aq is None: print('error, aquifer needs to be given')
        z1 = x1 + 1j * y1
        z2 = x2 + 1j * y2
        z = 0.5 * self.Xleg * (z2 - z1) + 0.5 * (z1 + z2)
        x = z.real
        y = z.imag
        pot = 0.0
        for i in range(self.ndeg):
            pot += self.wleg[i] * func(x=x[i], y=y[i], layers=layers, aq=aq)
        return pot

    def intflux(self, func, x1, y1, x2, y2, layers, aq=None):
        if aq is None: print('error, aquifer needs to be given')
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
            qxqy = func(x=x[i], y=y[i], layers=layers, aq=aq)
            qtot += self.wleg[i] * (qxqy[0] * cosnorm + qxqy[1] * sinnorm)
        return qtot

    def headinside(self):
        print('headinside not implemented for this element')

    def setparams(self, sol):
        raise Exception('Must overload Element.setparams()')

    def storeinput(self,frame):
        self.inputargs, _, _, self.inputvalues = inspect.getargvalues(frame)

    #def stoptrace(self, xyz, layer, ltype, step, direction):
    #    return False, 0

    def changetrace(self, xyzt1, xyzt2, aq, layer, ltype, modellayer, direction, hstepmax):
        changed = False
        terminate = False
        xyztnew = 0
        message = None
        return changed, terminate, xyztnew, message

    def qztop(self, x, y):
        # given flux at top of aquifer system (as for area-sinks)
        return 0

    def plot(self, layer):
        pass
    
    def write(self):
        rv = self.name + '(' + self.model.modelname + ',\n'
        for key in self.inputargs[2:]:  # The first two are ignored
            if isinstance(self.inputvalues[key], np.ndarray):
                rv += key + ' = ' + np.array2string(self.inputvalues[key],
                                                    separator=',') + ',\n'
            elif isinstance(self.inputvalues[key], str):                
                rv += key + " = '" + self.inputvalues[key] + "',\n"
            else:
                rv += key + ' = ' + str(self.inputvalues[key]) + ',\n'
        rv += ')\n'
        return rv

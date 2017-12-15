"""
Model classes

"""

import numpy as np
import sys
import inspect  # Used for storing the input
from .aquifer import Aquifer
from .aquifer_parameters import param_maq, param_3d
from .constant import ConstantStar
from .util import PlotTim

__all__ = ['Model', 'ModelMaq', 'Model3D']

class Model(PlotTim):
    """
    Model Class to create a model object consisting of an arbitrary
    sequence of aquifer layers and leaky layers.
    Use ModelMaq for regular sequence of aquifer and leaky layers.
    Use Model3D for multi-layer model of single aquifer
    
    Parameters
    ----------
    kaq : array
        hydraulic conductivity of each aquifer from the top down
    z : array
        elevation tops and bottoms of all layers
        layers may have zero thickness
    c : array
        resistance between two consecutive aquifers
        if ltype[0]='a': length is number of aquifers - 1
        if ltype[0]='l': length is number of aquifers
    npor : array
        porosity of all layers from the top down
    ltype : array of characters
        array indicating for each layer whether it is
        'a' aquifer layer
        'l' leaky layer
    
    """
    
    def __init__(self, kaq, c, z, npor, ltype):
        # All input variables are numpy arrays
        # That should be checked outside this function
        self.elementlist = []
        self.elementdict = {}  # only elements that have a label
        self.aq = Aquifer(self, kaq, c, z, npor, ltype)
        self.modelname = 'ml'  # Used for writing out input

    def initialize(self):
        # remove inhomogeneity elements (they are added again)
        self.elementlist = [e for e in self.elementlist if not e.inhomelement]
        self.aq.initialize()
        for e in self.elementlist:
            e.initialize()

    def add_element(self, e):
        self.elementlist.append(e)
        if e.label is not None: self.elementdict[e.label] = e

    def remove_element(self, e):
        """Remove element `e` from model
        """
        
        if e.label is not None: self.elementdict.pop(e.label)
        self.elementlist.remove(e)

    def storeinput(self, frame):
        self.inputargs, _, _, self.inputvalues = inspect.getargvalues(frame)

    def potential(self, x, y, aq=None):
        if aq is None: aq = self.aq.find_aquifer_data(x, y)
        pot = np.zeros(aq.naq)
        for e in aq.elementlist:
            pot += e.potential(x, y, aq)
        rv = np.sum(pot * aq.eigvec, 1)
        if aq.ltype[0] == 'l':
            # potential for head above leaky layer
            rv += aq.constantstar.potstar
        return rv

    def disvec(self, x, y, aq=None):
        """Discharge vector at `x`, `y`
        
        Returns
        -------
        
        qxqy : array size (2, naq)
            first row is Qx in each aquifer layer, second row is Qy
        """
        
        if aq is None: aq = self.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, aq.naq))
        for e in aq.elementlist:
            rv += e.disvec(x, y, aq)
        rv = np.sum(rv[:, np.newaxis, :] * aq.eigvec, 2)
        return rv
    
    def qztop(self, x, y, aq=None):
        if aq is None: aq = self.aq.find_aquifer_data(x, y)
        rv = 0.0
        if aq.ltype[0] == 'a':  # otherwise recharge cannot be added
            for e in aq.elementlist:
                rv += e.qztop(x, y)
        return rv

    def head(self, x, y, layers=None, aq=None):
        """Head at `x`, `y`
        
        Returns
        -------
        
        h : array length `naq` or `len(layers)`
            head in all `layers` (if not `None`), or all layers of aquifer (otherwise)
        """
        
        if aq is None: aq = self.aq.find_aquifer_data(x, y)
        rv = self.potential(x, y, aq) / aq.T
        if layers is None:
            return rv
        else:
            return rv[layers]

    def headgrid(self, xg, yg, layers=None, printrow=False):
        """Grid of heads
        
        Parameters
        ----------
        xg : array
            x values of grid
        yg : array
            y values of grid
        layers : integer, list or array, optional
            layers for which grid is returned
        printrow : boolean, optional
            prints dot to screen for each row of grid if set to `True`
        
        Returns
        -------
        h : array size `nlayers, ny, ny`
        
        See also
        --------
        
        :func:`~timml.model.Model.headgrid2`

        """
        
        nx, ny = len(xg), len(yg)
        if layers is None:
            Nlayers = self.aq.find_aquifer_data(xg[0], yg[0]).naq
        else:
            Nlayers = len(np.atleast_1d(layers))
        h = np.empty((Nlayers, ny, nx))
        for j in range(ny):
            if printrow:
                print(str(j) + ' ', flush=True)
            for i in range(nx):
                h[:, j, i] = self.head(xg[i], yg[j], layers)
        if printrow:
            print('', flush=True)
        return h

    def headgrid2(self, x1, x2, nx, y1, y2, ny, layers=None, printrow=False):
        """Grid of heads
        
        Parameters
        ----------
        x1, x2, nx : 
            x values of grid
        y1, y2, ny : 
            y values of grid
        layers : integer, list or array, optional
            layers for which grid is returned
        printrow : boolean, optional
            prints dot to screen for each row of grid if set to `True`
        
        Returns
        -------
        h : array size `nlayers, ny, ny`

        """
        
        xg, yg = np.linspace(x1, x2, nx), np.linspace(y1, y2, ny)
        return self.headgrid(xg, yg, layers=layers, printrow=printrow)

    def headalongline(self, x, y, layers=None):
        '''Returns head[Nlayers,len(x)]
        Assumes same number of layers for each x and y
        layers may be None or list of layers for which head is computed'''
        xg, yg = np.atleast_1d(x), np.atleast_1d(y)
        if layers is None:
            Nlayers = self.aq.find_aquifer_data(xg[0], yg[0]).naq
        else:
            Nlayers = len(np.atleast_1d(layers))
        nx = len(xg)
        if len(yg) == 1:
            yg = yg * np.ones(nx)
        h = np.zeros((Nlayers, nx))
        for i in range(nx):
            h[:, i] = self.head(xg[i], yg[i], layers)
        return h
    
    def disvecalongline(self, x, y, layers=None):
        '''Returns Qx[Nlayers,len(x)], Qy[Nlayers,len(x)]
        Assumes same number of layers for each x and y
        layers may be None or list of layers for which head is computed'''
        xg, yg = np.atleast_1d(x), np.atleast_1d(y)
        if layers is None:
            Nlayers = self.aq.find_aquifer_data(xg[0], yg[0]).naq
        else:
            Nlayers = len(np.atleast_1d(layers))
        nx = len(xg)
        if len(yg) == 1:
            yg = yg * np.ones(nx)
        Qx = np.zeros((Nlayers, nx))
        Qy = np.zeros((Nlayers, nx))
        for i in range(nx):
            Qx[:, i], Qy[:, 1] = self.disvec(xg[i], yg[i], layers)
        return Qx, Qy
    
    def velocity(self, x, y, z):
        return self.velocomp(x, y, z)
    
    def velocomp(self, x, y, z, aq=None, layer_ltype=None):
        if aq is None: aq = self.aq.find_aquifer_data(x, y)
        assert z <= aq.z[0] and z >= aq.z[-1], "z value not inside aquifer"
        if layer_ltype is None:
            layer, ltype, dummy = aq.findlayer(z)
        else:
            layer, ltype = layer_ltype
        h = self.head(x, y, aq=aq)
        # qz between aquifer layers
        qzlayer = np.zeros(aq.naq + 1)
        qzlayer[1:-1] = (h[1:] - h[:-1]) / aq.c[1:]
        if aq.ltype[0] == 'l':
            qzlayer[0] = (h[0] - aq.hstar) / aq.c[0]
        if ltype == 'l':
            vz = qzlayer[layer] / aq.nporll[layer]
            vx = 0
            vy = 0
        else:
            qzbot = qzlayer[layer + 1]
            qztop = qzlayer[layer]
            if layer == 0:
                qztop += self.qztop(x, y)   
            vz = (qzbot + (z - aq.zaqbot[layer]) / aq.Haq[layer] * \
                 (qztop - qzbot)) / aq.nporaq[layer]        
            qx, qy = self.disvec(x, y, aq=aq)
            vx = qx[layer] / (aq.Haq[layer] * aq.nporaq[layer])
            vy = qy[layer] / (aq.Haq[layer] * aq.nporaq[layer])
        return np.array([vx, vy, vz])
               
    def solve(self, printmat=0, sendback=0, silent=False):
        '''Compute solution'''
        # Initialize elements
        self.initialize()
        # Compute number of equations
        self.neq = np.sum([e.nunknowns for e in self.elementlist])
        if self.neq == 0: return
        if silent is False:
            print('Number of elements, Number of equations:', len(
                self.elementlist), ',', self.neq)
        if self.neq == 0:
            if silent is False: print('No unknowns. Solution complete')
            return
        mat = np.empty((self.neq, self.neq))
        rhs = np.empty(self.neq)
        ieq = 0
        for e in self.elementlist:
            if e.nunknowns > 0:
                mat[ieq:ieq + e.nunknowns, :], rhs[ieq:ieq + e.nunknowns] = \
                e.equation()
                ieq += e.nunknowns
            if silent is False:
                print('.', end='', flush=True)
        if printmat:
            return mat, rhs
        sol = np.linalg.solve(mat, rhs)
        icount = 0
        for e in self.elementlist:
            if e.nunknowns > 0:
                e.setparams(sol[icount:icount + e.nunknowns])
                icount += e.nunknowns
        if silent is False:
            print()  # needed cause the dots are printed
            print('solution complete')
        elif (silent == 'dot') or (silent == '.'):
            print('.', end='', flush=True)
        if sendback:
            return sol
        return            
        
class ModelMaq(Model):
    """
    Create a Model object by specifying a mult-aquifer sequence of
    aquifer-leakylayer-aquifer-leakylayer-aquifer etc
    
    Parameters
    ----------
    kaq : float, array or list
        hydraulic conductivity of each aquifer from the top down
        if float, hydraulic conductivity is the same in all aquifers
    z : array or list
        elevation tops and bottoms of the aquifers from the top down
        leaky layers may have zero thickness
        if topboundary='conf': length is 2 * number of aquifers
        if topboundary='semi': length is 2 * number of aquifers + 1 as top
        of leaky layer on top of systems needs to be specified
    c : float, array or list
        resistance of leaky layers from the top down
        if float, resistance is the same for all leaky layers
        if topboundary='conf': length is number of aquifers - 1
        if topboundary='semi': length is number of aquifers
    npor : float, array or list
        porosity of all aquifers and leaky layers from the top down
        if float, porosity is the same for all layers
        if topboundary='conf': length is 2 * number of aquifers - 1
        if topboundary='semi': length is 2 * number of aquifers
    topboundary : string, 'conf' or 'semi' (default is 'conf')
        indicating whether the topboundary is confined ('conf') or
        semi-confined ('semi')
    hstar : float or None (default is None)
        head value above semi-confining top, only read if topboundary='semi'

    Examples
    --------
    >>> ml = ModelMaq(kaq=[10, 20], z=[20, 12, 10, 0], c=1000)
    
    """
    
    def __init__(self, kaq=1, z=[1, 0], c=[], npor=0.3, topboundary='conf', hstar=None):
        self.storeinput(inspect.currentframe())
        kaq, c, npor, ltype = param_maq(kaq, z, c, npor, topboundary)
        Model.__init__(self, kaq, c, z, npor, ltype)
        self.name = 'ModelMaq'
        if self.aq.ltype[0] == 'l':
            ConstantStar(self, hstar, aq=self.aq)
            
class Model3D(Model):
    """
    Model3D Class to create a multi-layer model object consisting of
    many aquifer layers. The resistance between the layers is computed
    from the vertical hydraulic conductivity of the layers.
    
    Parameters
    ----------
    kaq : float, array or list
        hydraulic conductivity of each layer from the top down
        if float, hydraulic conductivity is the same in all aquifers
    z : array or list
        elevation of top of system followed by bottoms of all layers
        from the top down
        bottom of layer is automatically equal to top of layer below it
        if topboundary='conf': length is number of layers + 1
        if topboundary='semi': length is number of layers + 2 as top
        of leaky layer on top of systems needs to be specified
    kzoverkh : float
        vertical anisotropy ratio vertical k divided by horizontal k
        if float, value is the same for all layers
        length is number of layers
    npor : float, array or list
        porosity of all aquifer layers from the top down
        if float, porosity is the same for all layers
        if topboundary='conf': length is number of layers
        if topboundary='semi': length is number of layers + 1
    topboundary : string, 'conf' or 'semi' (default is 'conf')
        indicating whether the top is confined ('conf') or
        semi-confined ('semi')
    topres : float
        resistance of top semi-confining layer, only read if topboundary='semi'
    topthick: float
        thickness of top semi-confining layer, only read if topboundary='semi'
    hstar : float or None (default is None)
        head value above semi-confining top, only read if topboundary='semi'

    Examples
    --------
    >>> ml = Model3D(kaq=10, z=np.arange(20, -1, -2), kzoverkh=0.1)
    
    """
    
    def __init__(self, kaq=1, z=[1, 0], kzoverkh=1, npor=0.3, topboundary='conf', topres=0, topthick=0, hstar=0):
        '''Model3D
        for semi-confined aquifers, set top equal to 'semi' and provide
        topres: resistance of top
        tophick: thickness of top
        hstar: head above top'''
        self.storeinput(inspect.currentframe())
        kaq, c, npor, ltype = param_3d(kaq, z, kzoverkh, npor, topboundary, topres)
        if topboundary == 'semi':
            z = np.hstack((z[0] + topthick, z))
        Model.__init__(self, kaq, c, z, npor, ltype)
        self.name = 'Model3D'
        if self.aq.ltype[0] == 'l':
            ConstantStar(self, hstar, aq=self.aq)


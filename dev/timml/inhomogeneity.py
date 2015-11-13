import numpy as np
import inspect # Used for storing the input
from aquifer import AquiferData
from aquifer_parameters import param_maq
from constant import ConstantInside, ConstantStar
from intlinesink import IntHeadDiffLineSink, IntFluxDiffLineSink

class PolygonInhom(AquiferData):
    tiny = 1e-8
    def __init__(self, model, xy, kaq, Haq, c, z, npor, ltype, hstar, order, ndeg):
        # All input variables except model should be numpy arrays
        # That should be checked outside this function):        
        AquiferData.__init__(self, model, kaq, Haq, c, z, npor, ltype)
        self.order = order
        self.ndeg = ndeg
        self.hstar = hstar
        self.inhom_number = self.model.aq.add_inhom(self)
        self.z1, self.z2 = compute_z1z2(xy)
        self.Nsides = len(self.z1)
        Zin = 1e-6j
        Zout = -1e-6j
        self.zcin = 0.5 * (self.z2 - self.z1) * Zin + 0.5 * (self.z1 + self.z2)  # point at center on inside
        self.zcout = 0.5 * (self.z2 - self.z1) * Zout + 0.5 * (self.z1 + self.z2)  #point at center on outside
        self.x = np.hstack((self.z1.real, self.z2[-1].real))
        self.y = np.hstack((self.z1.imag, self.z2[-1].imag))
        self.xmin = min(self.x)
        self.xmax = max(self.x)
        self.ymin = min(self.y)
        self.ymax = max(self.y)
    def __repr__(self):
        return 'PolygonInhom: ' + str(zip(self.x,self.y))
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
    def create_elements(self):
        aqin = self.model.aq.find_aquifer_data(self.zcin[0].real, self.zcin[0].imag)
        for i in range(self.Nsides):
            aqout = self.model.aq.find_aquifer_data(self.zcout[i].real, self.zcout[i].imag)
            if (aqout == self.model.aq) or (aqout.inhom_number > self.inhom_number):
                ls = IntHeadDiffLineSink(self.model, x1=self.x[i], y1=self.y[i], x2=self.x[i+1], y2=self.y[i+1], \
                                    order=self.order, ndeg=self.ndeg, label=None, addtomodel=True, \
                                    aq=aqin, aqin=aqin, aqout=aqout)
                ls = IntFluxDiffLineSink(self.model, x1=self.x[i], y1=self.y[i], x2=self.x[i+1], y2=self.y[i+1], \
                                    order=self.order, ndeg=self.ndeg, label=None, addtomodel=True, \
                                    aq=aqout, aqin=aqin, aqout=aqout)
        if aqin.ltype[0] == 'a':  # add constant on inside
            c = ConstantInside(self.model, self.zcin.real, self.zcin.imag)
            c.inhomelement = True
        if aqin.ltype[0] == 'l':
            assert self.hstar is not None, 'Error: hstar needs to be set'
            c = ConstantStar(self.model, self.hstar, aq=aqin)
            c.inhomelement = True

        
class PolygonInhomMaq(PolygonInhom):
    tiny = 1e-8
    def __init__(self, model, xy, kaq=1, z=[1,0], c=[], npor=0.3, top='conf', hstar=None, order=3, ndeg=3):
        self.storeinput(inspect.currentframe())
        kaq, Haq, c, npor, ltype = param_maq(kaq, z, c, npor, top)
        PolygonInhom.__init__(self, model, xy, kaq, Haq, c, z, npor, ltype, hstar, order, ndeg)
        
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
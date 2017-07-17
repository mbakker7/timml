from numpy import *
import numpy.linalg as linalg
import copy

class AquiferData:
    '''Base class for aquifer data
    Attributes provided on input:
    - k: row vector of hydraulic conductivities; length Naquifers
    - zb: row vector of bottom elevations of aquifers; length Naquifers
    - zt: row vector of top elevations of aquifers; length Naquifers
    - c: row vector of resistances of leaky layers; length Naquifers - 1 if confined, Naquifers if semi confined
    - n: row vector of porosities of aquifers; length Naquifers
    - nll: row vector of porosities of leaky layers; length Naquifers - 1 if confined, Naquifers if semi confined
    Optional input DEPRECATED:
    - type: 'conf' if confined (default), 'semi' if semi-confined
    - hstar: head above resistance layer on top of aquifer system (only useful for semi-confined), 0 by default
    Attibutes computed:
    - H: row vector of aquifer thicknesses; length Naquifers
    - T: row vector of transmissivities; length Naquifers
    - Tcol: column vector of transmissivities; length Naquifers
    - Tcomp: comprehensive transmissivity
    - lab: row vector of leakage factors; length Naquifers - 1
    - eigvec: matrix with normalized transmissivity in first column and eigenvectors of system matrix in remaining columns;
              size: Naquifers by Naquifers
    '''
    huge = 1e30
    conf = 1; semi = 2
    def __init__(self,k=[1],zb=[0],zt=[1],c=[],n=[],nll=[],semi=False,type='conf',hstar=0.0):
        # Note that it is important that hstar = 0.0 default. So for confined aquifer we can use same headToPotential!
        self.k = array(k,'d'); self.zb = atleast_1d(zb); self.zt = array(zt,'d')
        self.Naquifers = len(zb)
        self.n = array(n,'d')
        if type == 'conf':
            self.type = self.conf
            self.c = array([self.huge] + list(c) + [self.huge],'d')
            self.nll = array([0.0]+list(nll),'d')
        elif type == 'semi':
            self.type = self.semi
            self.c = array( list(c) + [self.huge], 'd' )
            self.nll = array(nll,'d')
        self.hstar = float(hstar)
        self.fakesemi = False
        if semi: self.fakesemi = True
        self.setCoefs()
    def __repr__(self):
	return 'AquiferData N,k,zb,zt,c: ' + str((self.Naquifers,list(self.k),list(self.zb),list(self.zt),list(self.c)))
    def setCoefs(self):
        self.H = self.zt - self.zb; self.T = self.k * self.H; self.zcenter = 0.5 * (self.zt+self.zb)
        self.Tcomp = sum(self.T); self.Tcol = transpose(self.T[newaxis,:])
        self.HLeakyLayer = zeros(self.Naquifers,'d')
        if self.Naquifers > 0: self.HLeakyLayer[1:] = self.zb[:-1] - self.zt[1:]
        assert min(self.H)>0, "TimML Input error: At least one thickness is negative; Enter zb and zt from top down"
        (self.lab, self.eigvec) = self.systemMatrix()
        if self.Naquifers > 1 or self.type == self.semi :
            self.threelabmaxsq = ( 3.0 * max(self.lab) )**2
            self.labsorted = sort(self.lab)[::-1]
            self.threelabsorted = 3.0 * self.labsorted
            self.threelabsqsorted = ( 3.0 * self.labsorted )**2
            self.zeropluslab = hstack(( array([0.0]), self.lab ))
        else:
            self.threelabmaxsq = 0; self.labsorted = zeros(0,'d'); self.threelabsqsorted = zeros(0,'d')
            self.zeropluslab = array([0.0])
        if len(self.n) == 0: self.n = 0.3 * ones(self.Naquifers,'d')  # Sets default n
        if len(self.nll) != self.Naquifers: self.nll = 0.3 * ones(self.Naquifers,'d')  # Sets default nll

        # Compute the jump in potential between layers due to unit infiltration in top aquifer
        Tcumsum = cumsum(self.eigvec[:,0])  # Cumulative sum array
	NthroughLeakyLayer = 1.0 - Tcumsum[:-1]  # Length Naquifers-1
	self.constVecHead = zeros(self.Naquifers,'d')
	for i in range(1,self.Naquifers):
            self.constVecHead[i] = self.constVecHead[i-1] - NthroughLeakyLayer[i-1] * self.c[i]
        constVecPotJump = self.constVecHead * self.T
	self.constVecPot = linalg.solve(self.eigvec,constVecPotJump)

    def systemMatrix(self):
        A = zeros((self.Naquifers,self.Naquifers),'d')
        for n in range(1,self.Naquifers):     #Starting at 0
            A[n][n-1] = -1. / (self.c[n] * self.T[n-1])
        for n in range(0,self.Naquifers-1):
            A[n][n+1] = -1. / (self.c[n+1] * self.T[n+1])
        for n in range(0,self.Naquifers):
            A[n][n] = 1/ (self.c[n] * self.T[n]) + 1 / (self.c[n+1] * self.T[n])
        (W,U) = linalg.eig(A)
        if self.type == self.conf:
            imin = argmin(W)
            W = hstack(( W[:imin],W[imin+1:] )) # Remove minimum value (zero)
            U = U[:,range(0,imin)+range(imin+1,self.Naquifers)]
            U = hstack(( self.Tcol/self.Tcomp, U ))         # Make first row nomalized transmissivity
        W = 1. / sqrt(W)
	if W.dtype == 'complex128':
	    if abs(W.imag).max() < 1e-12:  # Quick fix for possibly small imaginary parts
		print 'Warning: tiny imaginary part of W neglected'
	        W = W.real
	if U.dtype == 'complex128':
            if abs(U.imag).max() < 1e-12:
		print 'Warning: tiny imaginary part of W neglected'
	        U = U.real
        return (W,U)
##    def systemMatrixOld(self):
##        A = zeros((self.Naquifers,self.Naquifers),'d')
##        for n in range(1,self.Naquifers):     #Starting at 0
##            A[n][n-1] = -1. / (self.c[n] * self.T[n-1])
##        for n in range(0,self.Naquifers-1):
##            A[n][n+1] = -1. / (self.c[n+1] * self.T[n+1])
##        for n in range(0,self.Naquifers):
##            A[n][n] = 1/ (self.c[n] * self.T[n]) + 1 / (self.c[n+1] * self.T[n])
##        (W,U) = MLab.eig(A)
##        W = list(W); U = list(U)                                 # Convert to lists
##        imin = argmin(W); void = W.pop(imin); void = U.pop(imin) # Remove minimum value
##        W = array(W); U = array(U)                               # Convert back to arrrays
##        W = 1. / sqrt(W)
##        U = array( list([self.T/self.Tcomp]) + list(U) )         # Make first row nomalized transmissivity
##        U = transpose(U)                                         # Take transpose so that 
##        return (W,U)
    def headToPotential(self,pylayer,head):
        if self.type == self.conf:
            rv = head * self.T[pylayer]
        elif self.type == self.semi:
            rv = self.T[pylayer] * (head - self.hstar)
        return rv
    def potentialToHead(self,pylayer,pot):
        if self.type == self.conf:
            rv = pot / self.T[pylayer]
        elif self.type == self.semi:
            rv = pot / self.T[pylayer] + self.hstar
        return rv
    def potentialVectorToHeadVector(self,potVec):
        if self.type == self.conf:
            rv = potVec / self.T
        elif self.type == self.semi:
            rv = potVec / self.T + self.hstar
        return rv
    def isInside(self,x,y):
        return 0
    def inWhichPyLayer(self,z):
        '''Returns -9999 if above top of system, +9999 if below bottom of system, negative for in leaky layer.
        leaky layer -n is on top of aquifer n'''
        if z > self.zt[0]:
            return -9999
        for i in range(self.Naquifers-1):
            if z >= self.zb[i]:
                return i
            if z > self.zt[i+1]:
                return -i-1
        if z >= self.zb[self.Naquifers-1]:
            return self.Naquifers - 1
        return +9999
    def layout(self):
        return [0]
    def crossBoundary(self,pyLayer,xyz1,xyz2,aqOther,step,idir):
        '''Returns 4 values:
        changed
        stop
        xyznew
        backToOld is set to true if the step size is reduced
        '''
        return [0,0,0,0]
    
class Aquifer(AquiferData):
    '''Aquifer class, derived from AquiferData. This is the background aquifer
    Attributes:
    - inhomList: list of inhomogeneity elements (for now only PolygonInhom)
    - all attributes of AquiferData
    '''
    def __init__(self,k=[1],zb=[0],zt=[1],c=[],n=[],nll=[],semi=False,type='conf',hstar=0.0):
        AquiferData.__init__(self,k,zb,zt,c,n,nll,semi,type,hstar)
        self.inhomList = []
    def addInhom(self,aq):
        self.inhomList.append(aq)
    def findAquiferData(self,x,y):
        rv = self
        for inhom in self.inhomList:
            if inhom.isInside(x,y):
                rv = inhom
                return rv
        return rv
  
class PolygonInhom(AquiferData):
    '''Class containing data of a polygonal inhomogeneity. Derived from AquiferData
    Attributes provided on input:
    - modelParent: model it belongs to
    - k: row vector of hydraulic conductivities; length Naquifers
    - zb: row vector of bottom elevations of aquifers; length Naquifers
    - zt: row vector of top elevations of aquifers; length Naquifers
    - c: row vector of resistances of leaky layers; length Naquifers - 1
    - xylist: list of tuples of xy pairs: [(x1,y1),(x2,y2),...]; length Ncorners; first point is NOT repeated
    Keyword arguments:
    - n: list with porosities of the aquifers (if not provided set to 0.3)
    - nll: list with porosities of the leaky layers (if not provided set to 0.3)
    - semi: default is False. When set to True, the water in the top aquifer is treated as a water table above 
      a semi-confining layer. In this case, the number of aquifers must be 1 more than the number of aquifers
      in the part of the model where semi is set to False. 
    Attributes computed:
    - Ncorners: number of corners of polygon
    - z1: list of complex coordinates of begin points of polygon segments; length Ncorners-1
    - z2: list of complex coordinates of end points of polygon segments; length Ncorners-1
    - xmin: minimum x of bounding box
    - xmax: maximum x of bounding box
    - ymin: minimum y of bounding box
    - xmax: maximum x of bounding box
    '''
    tiny = 1e-8
    def __init__(self,modelParent,k,zb,zt,c,xylist,n=[],nll=[],semi=False,type='conf',hstar=0.0):        
        AquiferData.__init__(self,k,zb,zt,c,n,nll,semi,type,hstar)
        self.modelParent = modelParent; modelParent.aq.addInhom(self)
        xylist = check_direction(xylist)
        self.xylist = xylist; self.Ncorners = len(self.xylist)
        xcorners = []; ycorners = []; self.z1 = []
        for xy in self.xylist:
            xcorners = xcorners + [xy[0]]; ycorners = ycorners + [xy[1]]
            self.z1 = self.z1 + [complex(xy[0],xy[1])]
        self.z2 = copy.copy(self.z1); zfirst = self.z2.pop(0); self.z2 = self.z2 + [zfirst]
        self.z1 = array(self.z1); self.z2 = array(self.z2)
        self.xmin = min(xcorners); self.xmax = max(xcorners)
        self.ymin = min(ycorners); self.ymax = max(ycorners)
        self.thetaNormOut = arctan2( self.z2.imag - self.z1.imag, self.z2.real - self.z1.real ) - pi/2.0
    def __repr__(self):
	return 'PolygonInhom N,k,zb,zt,c,xylist: ' + str((self.Naquifers,list(self.k),list(self.zb),list(self.zt),list(self.c),self.xylist))
    def isInside(self,x,y):
        rv = 0
        if x >= self.xmin and x <= self.xmax and y >= self.ymin and y <= self.ymax:
            z = complex(x,y)
            bigZ = ( 2.0*z - (self.z1 + self.z2) )/ (self.z2 - self.z1)
            bigZmin1 = bigZ - 1.0; bigZplus1 = bigZ + 1.0
            minAbsBigZmin1 = min(abs(bigZmin1)); minAbsBigZplus1 = min(abs(bigZplus1))
            if minAbsBigZmin1 < self.tiny or minAbsBigZplus1 < self.tiny:
                    rv = 1
                    return rv
            angles = log( bigZmin1 / bigZplus1 ).imag
            angle = sum(angles)
            if angle > pi: rv = 1
        return rv
    def movePoint(self,zin):
        # Move point if at corner point; this doesn't happy very often, so it is not that efficient
        bigZ = ( 2.0*zin - (self.z1 + self.z2) )/ (self.z2 - self.z1)
        bigZmin1 = bigZ - 1.0
        imin = argmin( abs(bigZmin1) )
        if abs(bigZmin1[imin]) < self.tiny:  # Then we are inside, as by isInside routine
            znew = self.z2[imin] + complex(4.0,-1.0) * self.tiny * (self.z1[imin] - self.z2[imin])
            return znew
        # In the off chance that this doesn't do it, we need to check the other side
        bigZplus1 = bigZ + 1.0
        imin = argmin( abs(bigZplus1) )
        if abs(bigZplus1[imin]) < self.tiny:  # Then we are inside, as by isInside routine
            if imin > 0:
                imin = imin - 1  # Element to the left
            else:
                imin = self.Ncorners-1
            znew = self.z2[imin] + complex(4.0,-1.0) * self.tiny * (self.z1[imin] - self.z2[imin])
            return znew
        # If that doesn't do anything then there is no problem and we return the input variable
        return zin
    def layout(self):
        return zip(*self.xylist)
    def crossBoundary(self,pyLayer,xyz1,xyz2,aqOther,step,idir):
        changed = 0; stop = 0; xyznew = 0.0; backToOld = 0
        [x1,y1,z1] = xyz1; [x2,y2,z2] = xyz2
        z1in = complex(x1,y1); z2in = complex(x2,y2)
        bigZ1 = ( 2.0*z1in - (self.z1 + self.z2) )/ (self.z2 - self.z1)
        bigZ2 = ( 2.0*z2in - (self.z1 + self.z2) )/ (self.z2 - self.z1)
        for i in range(len(bigZ1)):
            Z1 = bigZ1[i]; Z2 = bigZ2[i]            
        # If point 1 on one side and point 2 on other side, find out if segment intesects line-doublet
            if (Z1.imag > 0 and Z2.imag < 0) or (Z1.imag < 0 and Z2.imag > 0):
                Xintersec = Z1.real + (0.0-Z1.imag) / (Z2.imag-Z1.imag) * (Z2.real-Z1.real)
                if abs(Xintersec) <= 1:  # There is an intersection point

                    changed = 1  # Point will be changed

                    # Check if number of aquifers is same on both sides, taking into account fakesemi
                    Naq = self.Naquifers
                    if self.fakesemi: Naq = Naq - 1
                    NaqOther = aqOther.Naquifers
                    if aqOther.fakesemi: NaqOther = NaqOther - 1
                    if Naq == NaqOther:
                        theta = self.thetaNormOut[i]
                        if Z1.imag < 0:  # old point on outside, Znew1 just out, Znew2 just in
                            Znew1 = complex(Xintersec,-1e-6); Znew2 = complex(Xintersec,1e-6)
                            if idir > 0: theta = self.thetaNormOut[i] + pi  # Such that theta points in supposed direction of flow
                        else:
                            Znew1 = complex(Xintersec,1e-6); Znew2 = complex(Xintersec,-1e-6)
                            if idir < 0: theta = self.thetaNormOut[i] + pi
                        znew1 = ( (self.z2[i] - self.z1[i])*Znew1 + (self.z1[i] + self.z2[i]) ) / 2.0  # New complex coordinate
                        znew2 = ( (self.z2[i] - self.z1[i])*Znew2 + (self.z1[i] + self.z2[i]) ) / 2.0
                        xnew1 = znew1.real; ynew1 = znew1.imag
                        xnew2 = znew2.real; ynew2 = znew2.imag
                        horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
                        horstepnew = sqrt( (xyz1[0]-xnew1)**2 + (xyz1[1]-ynew1)**2 )
                        zvertnew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )  # Vertical coordinate at new1
                        # Check if Z1 just before crossing. If not, set new point equal to znew1
                        # This avoids problems with stepping into inhom and stepping into new layer at the same time
                        if abs(Z1.imag) > 2e-6:
                            xyznew = array([xnew1,ynew1,zvertnew])
                            backToOld = 1
                            return [changed, stop, xyznew, backToOld]
                        disnorm = self.modelParent.dischargeNormVector(xnew1,ynew1,theta)[pyLayer]  # Normal discharge at new1
                        if disnorm < 0:
                            # Normal discharge at boundary points the other way; can only be if step size too large
                            # Half step to boundary
                            xyznew = 0.5 * ( xyz1 + array([xnew1,ynew1,zvertnew]) )
                            backToOld = 1
                        else:
                            pyLayernew = pyLayer
                            if Z1.imag < 0:  # old point on outside, step from aqOther into self
                                if self.fakesemi and not aqOther.fakesemi: pyLayernew = pyLayernew + 1
                                if not self.fakesemi and aqOther.fakesemi: pyLayernew = pyLayernew - 1
                                frac = ( zvertnew - aqOther.zb[pyLayer] ) / aqOther.H[pyLayer]
                                znew2 = self.zb[pyLayernew] + frac * self.H[pyLayernew]
                            else: # old point on inside, step from self to aqOther
                                if self.fakesemi and not aqOther.fakesemi: pyLayernew = pyLayernew - 1
                                if not self.fakesemi and aqOther.fakesemi: pyLayernew = pyLayernew + 1
                                frac = ( zvertnew - self.zb[pyLayer] ) / self.H[pyLayer]
                                znew2 = aqOther.zb[pyLayernew] + frac * aqOther.H[pyLayernew]
                            xyznew = array([xnew2,ynew2,znew2])

                        return [changed, stop, xyznew, backToOld]  # Didn't adjust time

                    else:
                        theta = self.thetaNormOut[i]
                        if Z1.imag < 0:  # old point on outside, Znew1 just out, Znew2 just in
                            Znew1 = complex(Xintersec,-1e-6); Znew2 = complex(Xintersec,1e-6)
                            if idir > 0: theta = self.thetaNormOut[i] + pi  # Such that theta points in supposed direction of flow
                        else:
                            Znew1 = complex(Xintersec,1e-6); Znew2 = complex(Xintersec,-1e-6)
                            if idir < 0: theta = self.thetaNormOut[i] + pi
                        znew1 = ( (self.z2[i] - self.z1[i])*Znew1 + (self.z1[i] + self.z2[i]) ) / 2.0  # New complex coordinate
                        znew2 = ( (self.z2[i] - self.z1[i])*Znew2 + (self.z1[i] + self.z2[i]) ) / 2.0
                        xnew1 = znew1.real; ynew1 = znew1.imag
                        xnew2 = znew2.real; ynew2 = znew2.imag
                        horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
                        horstepnew = sqrt( (xyz1[0]-xnew1)**2 + (xyz1[1]-ynew1)**2 )
                        zvertnew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )  # Vertical coordinate at new1
                    # This part only needed for different number of aquifers left and right.
                    # Modification for samen number of aquifers but change in thickness needs to be here
                        disnorm = self.modelParent.dischargeNormVector(xnew1,ynew1,theta)[pyLayer]  # Normal discharge at xnew1,ynew1
                        if disnorm < 0:
                            # Normal discharge at boundary points the other way; can only be if step size too large
                            # Half step to boundary
                            xyznew = 0.5 * ( xyz1 + array([xnew1,ynew1,zvertnew]) )
                            backToOld = 1
                        else:
                            dis = self.modelParent.dischargeNormBelowZ(xnew1,ynew1,zvertnew,theta)
                            [xnew,ynew,znew] = self.modelParent.zForGivenNormalDischarge\
                                            (zvertnew,xnew1,ynew1,xnew2,ynew2,theta,dis)
                            xyznew = array([xnew,ynew,znew])

                        return [changed, stop, xyznew, backToOld]  # Didn't adjust time
        return [changed, stop, xyznew, backToOld]

def check_direction(xylist):
    x,y = zip(*xylist)
    if x[0] == x[-1] and y[0] == y[-1]:  # In case last point is repeated
        x = x[:-1]; y = y[:-1]
    z1 = array(x) + array(y) * 1j
    index = range(1,len(z1)) + [0]
    z2 = z1[index]
    Z = 1e-6j
    z = Z * (z2[0] - z1[0]) / 2.0 + 0.5 * (z1[0] + z2[0])
    bigZ = ( 2.0*z - (z1 + z2) )/ (z2 - z1)
    bigZmin1 = bigZ - 1.0; bigZplus1 = bigZ + 1.0
    angle = sum( log( bigZmin1 / bigZplus1 ).imag )
    if angle < pi: # reverse order
        xylist = zip(x[::-1],y[::-1])
    else:
        xylist = zip(x,y) # Is corrected for repeated point at end
    return xylist
        
class PolygonInhomComp(PolygonInhom):
    def __init__(self,modelParent,z=[1,0.5,0],kh=1.0,kzoverkh=1.0,xylist=[],n=[],nll=[]):
        z=atleast_1d(z,'d')
        zt = z[0:-1]
        zb = z[1:]
        H = zt - zb
        k = kh * ones(Naq)
        Hinter = 0.5 * ( H[0:-1] + H[1:] )
        c = Hinter / ( kzoverkh * 0.5*( k[0:-1] + k[1:] ) )
        PolygonInhom.__init__(self,modelParent,k,zb,zt,c,xylist,n,nll,type='conf',hstar=0.0)

class CircleInhomData(AquiferData):
    '''Class containing data of a circular inhomogeneity.
    Derived from AquiferData
    Attributes provided on input:
    - modelParent: model it belongs to
    - k: row vector of hydraulic conductivities; length Naquifers
    - zb: row vector of bottom elevations of aquifers; length Naquifers
    - zt: row vector of top elevations of aquifers; length Naquifers
    - c: row vector of resistances of leaky layers; length Naquifers - 1
    - xc: x-coordinate of center of circle
    - yc: y-coordinate of center of circle
    - R: radius of circle
    Attributes computed:
    - Rsq: square of radius
    '''
    def __init__(self,modelParent,k,zb,zt,c,xc,yc,R,n=[],nll=[],semi=False,type='conf',hstar=0.0):        
        AquiferData.__init__(self,k,zb,zt,c,n,nll,semi,type,hstar)
        self.modelParent = modelParent; modelParent.aq.addInhom(self)
        self.xc = float(xc); self.yc = float(yc); self.R = float(R); self.Rsq = R*R
    def __repr__(self):
	return 'CircleInhom N,k,zb,zt,c,xc,yc,R: ' + str((self.Naquifers,list(self.k),list(self.zb),list(self.zt),\
                                                          list(self.c),self.xc,self.yc,self.R))
    def isInside(self,x,y):
        return ( (x-self.xc)**2 + (y-self.yc)**2 ) < self.Rsq
    def layout(self):
        theta = arange(0,2*pi+1e-5,pi/50)
        return [ list( self.xc + self.R * cos(theta) ), list( self.yc + self.R * sin(theta) ) ]
    def crossBoundary(self,pyLayer,xyz1,xyz2,aqOther,step,idir):
        tol = 1e-8  # Tolerance used to step to edge of circle
        changed = 0; stop = 0; xyznew = 0.0; backToOld = 0
        [x1,y1,z1] = xyz1; [x2,y2,z2] = xyz2
        R1 = sqrt( (x1-self.xc)**2 + (y1-self.yc)**2 )
        R2 = sqrt( (x2-self.xc)**2 + (y2-self.yc)**2 )
        if ( R1 < self.R and R2 > self.R ) or ( R1 > self.R and R2 < self.R ):
            changed = 1  # Point will be changed
            # Check if number of aquifers is same on both sides, taking into account fakesemi
            # Not tested for anyting but Naq == NaqOther
            Naq = self.Naquifers
            if self.fakesemi: Naq = Naq - 1
            NaqOther = aqOther.Naquifers
            if aqOther.fakesemi: NaqOther = NaqOther - 1
            if Naq == NaqOther:
                a = (x2-x1)**2 + (y2-y1)**2
                b = 2.0 * ( (x2-x1) * (x1-self.xc) + (y2-y1) * (y1-self.yc) )
                if R1 < self.R:
                    if R1 < (1.0-tol)*self.R:  # first step to just inside circle
                        backToOld = 1
                        c = self.xc**2 + self.yc**2 + x1**2 + y1**2 - 2.0 * (self.xc*x1 +self.yc*y1) - (1.0-tol)*self.Rsq
                    else:
                        c = self.xc**2 + self.yc**2 + x1**2 + y1**2 - 2.0 * (self.xc*x1 +self.yc*y1) - (1.0+tol)*self.Rsq
                else:
                    if R1 > (1.0+tol)*self.R:  # first step to just outside circle
                        backToOld = 1
                        c = self.xc**2 + self.yc**2 + x1**2 + y1**2 - 2.0 * (self.xc*x1 +self.yc*y1) - (1.0+tol)*self.Rsq
                    else:
                        c = self.xc**2 + self.yc**2 + x1**2 + y1**2 - 2.0 * (self.xc*x1 +self.yc*y1) - (1.0-tol)*self.Rsq
                u1 = ( -b - sqrt(b**2 - 4.0*a*c) ) / (2.0*a)
                u2 = ( -b + sqrt(b**2 - 4.0*a*c) ) / (2.0*a)
                u = u1
                if u <= 0: u = u2
                xn = x1 + u * (x2-x1); yn = y1 + u * (y2-y1)
                zn = xyz1[2] + u * ( xyz2[2] - xyz1[2] )
                xyznew = array([xn,yn,zn],'d')
            else:
                print 'Naq not equal to Naqother. This case is not implemented'
        return [changed, stop, xyznew, backToOld]

class EllipseInhomData(AquiferData):
    '''Class containing data of elliptical inhomogeneity.
    Derived from AquiferData
    Attributes provided on input:
    - modelParent: model it belongs to
    - k: row vector of hydraulic conductivities; length Naquifers
    - zb: row vector of bottom elevations of aquifers; length Naquifers
    - zt: row vector of top elevations of aquifers; length Naquifers
    - c: row vector of resistances of leaky layers; length Naquifers - 1
    - xc: x-coordinate of center of ellipse
    - yc: y-coordinate of center of ellipse
    - along: length of long axis of ellipse
    - bshort: length of short axis of ellipse
    - angle: angle of long axis with horizontal
    Attributes computed:
    - cosal: cosine of alpha
    - sinal: sine of alpha
    - afoc: focal distance
    - etastar: eta value along ellipse
    '''
    def __init__(self,modelParent,k,zb,zt,c,xc,yc,along,bshort,alpha,n=[],nll=[],semi=False,type='conf',hstar=0.0):        
        AquiferData.__init__(self,k,zb,zt,c,n,nll,semi,type,hstar)
        self.modelParent = modelParent; modelParent.aq.addInhom(self)
        self.xc = float(xc); self.yc = float(yc); self.along = float(along); self.bshort = bshort;
        assert self.along > self.bshort, "TimML Input error: Long axis of ellipse must be larger than short axis"
        self.angle = alpha; self.cosal = cos(self.angle); self.sinal = sin(self.angle)
        self.afoc = sqrt( self.along**2 - self.bshort**2 )
        self.etastar = arccosh( self.along / self.afoc )
        self.zc = complex(self.xc,self.yc)
    def __repr__(self):
	return 'EllipseInhom N,k,zb,zt,c,xc,yc,a,b,angle: ' \
               + str((self.Naquifers,list(self.k),list(self.zb),list(self.zt),\
                        list(self.c),self.xc,self.yc,self.along,self.bshort,self.angle))
    def isInside(self,x,y):
        [eta,psi] = self.xytoetapsi(x,y)
        return eta < self.etastar
    def outwardNormal(self,x,y):
        [eta,psi] = self.xytoetapsi(x,y)
        alpha = arctan2( cosh(eta)*sin(psi), sinh(eta)*cos(psi) ) + self.angle
        return [cos(alpha), sin(alpha)]
    def outwardNormalAngle(self,x,y):
        [eta,psi] = self.xytoetapsi(x,y)
        alpha = arctan2( cosh(eta)*sin(psi), sinh(eta)*cos(psi) ) + self.angle
        return alpha
    def xytoetapsiold(self,xin,yin):
        xloc = xin - self.xc ; yloc = yin - self.yc
        x =  xloc * self.cosal + yloc * self.sinal
        y = -xloc * self.sinal + yloc * self.cosal
        r1 = sqrt( (x + self.afoc)**2 + y**2 )
        r2 = sqrt( (x - self.afoc)**2 + y**2 )
        eta = arccosh ( (r1+r2) / (2.0 * self.afoc) )
        arg = 2.0 * x / (r1+r2)
        if arg > 1.0: arg = 1.0  # This happens because of round-off error
        if arg < -1.0: arg = -1.0 
        psi = arccos( arg )
        if y < 0 :
            psi = 2.0 * pi - psi
        return [eta,psi]
    def xytoetapsi(self,xin,yin):
        z = complex(xin,yin)
        Z = (z - self.zc) * exp( -complex(0,1) * self.angle )
        tau = log( Z / self.afoc + sqrt( Z/self.afoc - 1 ) * sqrt( Z/self.afoc + 1 ) )
        return [tau.real,tau.imag]
    def etapsitoxy(self,eta,psi):
        xloc = self.afoc * cosh(eta) * cos(psi)
        yloc = self.afoc * sinh(eta) * sin(psi)
        x = xloc * self.cosal - yloc * self.sinal + self.xc
        y = xloc * self.sinal + yloc * self.cosal + self.yc
        return [x,y]
    def layout(self):
        theta = arange(0,2*pi+0.001,pi/50)
        return [ list( self.etapsitoxy(self.etastar,theta)[0] ), list( self.etapsitoxy(self.etastar,theta)[1] ) ]

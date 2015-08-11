
from numpy import *
import numpy.linalg as linalg
from mlelement import *
from mlholineelements import *
from besselaes import *

class PolyAreaSink(Element):
    tiny = 1e-6
    def __init__(self,modelParent,xylist,infil,label=None):
	Element.__init__(self,modelParent)
	self.xylist = xylist
	self.infil = float(infil)
        self.label = None
        self.type = 'polyareasink'
        self.setCoefs()
        self.modelParent.addElement(self)
    def __repr__(self):
	return 'PolyAreaSink xylist,infil: ' + str(self.xylist) + ', ' + str(self.infil)
    def setCoefs(self):
        self.parameters = array([[self.infil]])
        self.Nsides = len(self.xylist)
        # Store coordinates
        self.xp = zeros(self.Nsides+1,'d'); self.yp = zeros(self.Nsides+1,'d')
        self.z1 = zeros(self.Nsides,'D'); self.z2 = zeros(self.Nsides,'D')
        for i in range(self.Nsides):
            self.xp[i] = self.xylist[i][0]
            self.yp[i] = self.xylist[i][1]
            self.z1[i] = complex( self.xp[i], self.yp[i] )
        self.xp[-1] = self.xylist[0][0]
        self.yp[-1] = self.xylist[0][1]
        self.z2[:-1] = self.z1[1:]; self.z2[-1] = self.z1[0]
        self.lengths = abs(self.z2 - self.z1)
        self.xc = average(self.xp[:-1]); self.yc = average(self.yp[:-1])  # Center of area-sink
        self.xmin = min(self.xp)-2.0*self.tiny; self.xmax = max(self.xp)+2.0*self.tiny
        self.ymin = min(self.yp)-2.0*self.tiny; self.ymax = max(self.yp)+2.0*self.tiny
        self.xpmin = zeros(self.Nsides,'d'); self.ypmin = zeros(self.Nsides,'d')  # Used for quickly finding intesection
        self.xpmax = zeros(self.Nsides,'d'); self.ypmax = zeros(self.Nsides,'d') 
        for i in range(self.Nsides):
            self.xpmin[i] = min(self.xp[i],self.xp[i+1]); self.xpmax[i] = max(self.xp[i],self.xp[i+1])
            self.ypmin[i] = min(self.yp[i],self.yp[i+1]); self.ypmax[i] = max(self.yp[i],self.yp[i+1])
        # Set aquifer parent
	self.aquiferParent = self.modelParent.aq.findAquiferData(self.xc,self.yc)
	# Compute constant on inside
	Tcumsum = cumsum(self.aquiferParent.eigvec[:,0])
	NthroughLeakyLayer = 1.0 - Tcumsum[:-1]  # Length Naquifers-1
	self.constVecHead = zeros(self.aquiferParent.Naquifers,'d')
	for i in range(1,self.aquiferParent.Naquifers):
            self.constVecHead[i] = self.constVecHead[i-1] - NthroughLeakyLayer[i-1] * self.aquiferParent.c[i]
        self.constVecPotJump = self.constVecHead * self.aquiferParent.T
	self.constVecPot = linalg.solve(self.aquiferParent.eigvec,self.constVecPotJump)
	# Construct line-sinks and line-doublets
        self.lsList = []; self.ldListBes = []; self.ldListLap = []
        for i in range(self.Nsides):
            a = 0.5 * (self.xp[i+1] - self.xp[i])
            b = 0.5 * (self.xp[i] + self.xp[i+1])
            nout = -(self.z2[i] - self.z1[i]) / abs(self.z2[i]-self.z1[i]) * complex(0.0,1.0)
            noutx = nout.real
            sigma = - array([ b-self.xc, a ]) * noutx
            ls = LineSinkHo(self.modelParent,self.xp[i],self.yp[i],self.xp[i+1],self.yp[i+1],\
                            1,sigma,addToModel=0,Bessel=0)
            self.lsList.append(ls)
            if self.aquiferParent.Naquifers > 1:
                ld = LineDoubletHo(self.modelParent,self.xp[i],self.yp[i],self.xp[i+1],self.yp[i+1],\
                            0,-self.constVecPotJump[newaxis,:],addToModel=0,Bessel=1)
                self.ldListBes.append(ld)
            potjump = 0.5 * array([ (b-self.xc)**2, 2*a*(b-self.xc), a**2 ])
            ld = LineDoubletHo(self.modelParent,self.xp[i],self.yp[i],self.xp[i+1],self.yp[i+1],\
                            2,potjump,addToModel=0,Bessel=0)
            self.ldListLap.append(ld)
        self.potLapLsInf = zeros(2,'d')
        self.potLapLdInf = zeros(3,'d')
        self.potBesInf = zeros(self.aquiferParent.Naquifers,'d')
        self.disxLapLsInf = zeros(2,'d'); self.disyLapLsInf = zeros(2,'d')
        self.disxLapLdInf = zeros(3,'d'); self.disyLapLdInf = zeros(3,'d')
        self.disxBesInf = zeros(self.aquiferParent.Naquifers,'d'); self.disyBesInf = zeros(self.aquiferParent.Naquifers,'d')
    def potentialInfluence(self,aq,x,y):
        # Need to modify this for different aquifer types and different inhomogeneities!
        pot = zeros(aq.Naquifers,'d')
        [isInside,x,y] = self.isInside(x,y)
        for ls in self.lsList:
            potlaplsho(x,y,ls.x1,ls.y1,ls.x2,ls.y2,ls.order,self.potLapLsInf)
            pot[0] = pot[0] + self.potLapLsInf[0] * ls.coef[0,0] + self.potLapLsInf[1] * ls.coef[1,0]
        for ld in self.ldListBes:
            potbesldho(x,y,ld.x1,ld.y1,ld.x2,ld.y2,aq.Naquifers,aq.zeropluslab,0,self.potBesInf)
            pot = pot + ld.coef[0,:] * self.potBesInf
        for ld in self.ldListLap:
            potlapldho(x,y,ld.x1,ld.y1,ld.x2,ld.y2,ld.order,self.potLapLdInf)
            pot[0] = pot[0] + self.potLapLdInf[0] * ld.coef[0,0] + self.potLapLdInf[1] * ld.coef[1,0] + \
                     self.potLapLdInf[2] * ld.coef[2,0]
        if isInside:
            pot[0] = pot[0] - 0.5 * (x - self.xc)*(x - self.xc)
            if aq == self.aquiferParent:
                pot = pot + self.constVecPot
            else:
                pot[0] = pot[0] + self.constVecPot[0]
        return pot
    def potentialCollection(self,potsum,potadd,elementList,aq,x,y):
        for el in elementList:
            potadd[:] = 0.0
            [isInside,x,y] = el.isInside(x,y)
            for ls in el.lsList:
                potlaplsho(x,y,ls.x1,ls.y1,ls.x2,ls.y2,ls.order,el.potLapLsInf)
                potadd[0] = potadd[0] + el.potLapLsInf[0] * ls.coef[0,0] + el.potLapLsInf[1] * ls.coef[1,0]
            for ld in el.ldListBes:
                potbesldho(x,y,ld.x1,ld.y1,ld.x2,ld.y2,aq.Naquifers,aq.zeropluslab,0,el.potBesInf)
                potadd = potadd + ld.coef[0,:] * el.potBesInf
            for ld in el.ldListLap:
                potlapldho(x,y,ld.x1,ld.y1,ld.x2,ld.y2,ld.order,el.potLapLdInf)
                potadd[0] = potadd[0] + el.potLapLdInf[0] * ld.coef[0,0] + el.potLapLdInf[1] * ld.coef[1,0] + \
                         el.potLapLdInf[2] * ld.coef[2,0]
            if isInside:
                potadd[0] = potadd[0] - 0.5 * (x - el.xc)*(x - el.xc)
                if aq == el.aquiferParent:
                    potadd = potadd + el.constVecPot
                else:
                    potadd[0] = potadd[0] + el.constVecPot[0]
            potsum = potsum + el.parameters[0,0] * potadd
        return potsum
    def dischargeCollection(self,dissum,disadd,elementList,aq,x,y):
        for el in elementList:
            disadd[:,:] = 0.0
            [isInside,x,y] = el.isInside(x,y)
            for ls in el.lsList:
                dislaplsho(x,y,ls.x1,ls.y1,ls.x2,ls.y2,ls.order,el.disxLapLsInf,el.disyLapLsInf)
                disadd[0,0] = disadd[0,0] + el.disxLapLsInf[0] * ls.coef[0,0] + el.disxLapLsInf[1] * ls.coef[1,0]
                disadd[1,0] = disadd[1,0] + el.disyLapLsInf[0] * ls.coef[0,0] + el.disyLapLsInf[1] * ls.coef[1,0]
            for ld in el.ldListBes:
                disbesldho(x,y,ld.x1,ld.y1,ld.x2,ld.y2,aq.Naquifers,aq.zeropluslab,0,el.disxBesInf,el.disyBesInf)
                disadd[0,:] = disadd[0,:] + ld.coef[0,:] * el.disxBesInf
                disadd[1,:] = disadd[1,:] + ld.coef[0,:] * el.disyBesInf
            for ld in el.ldListLap:
                dislapldho(x,y,ld.x1,ld.y1,ld.x2,ld.y2,ld.order,el.disxLapLdInf,el.disyLapLdInf)
                disadd[0,0] = disadd[0,0] + el.disxLapLdInf[0] * ld.coef[0,0] + el.disxLapLdInf[1] * ld.coef[1,0] + \
                         el.disxLapLdInf[2] * ld.coef[2,0]
                disadd[1,0] = disadd[1,0] + el.disyLapLdInf[0] * ld.coef[0,0] + el.disyLapLdInf[1] * ld.coef[1,0] + \
                         el.disyLapLdInf[2] * ld.coef[2,0]
            if isInside:
                disadd[0,0] = disadd[0,0] + (x - el.xc)
            dissum = dissum + el.parameters[0,0] * disadd
        return dissum
    def potentialContribution(self,aq,x,y):
        '''Returns VECTOR of potentialContribution; doesn't multiply with eigenvectors'''
        return self.parameters[0,0] * self.potentialInfluence(aq,x,y)
    def dischargeInfluence(self,aq,x,y):
        # Need to modify this for different aquifer types and different inhomogeneities!
        disx = zeros(aq.Naquifers,'d'); disy = zeros(aq.Naquifers,'d')
        [isInside,x,y] = self.isInside(x,y)
        for ls in self.lsList:
            dislaplsho(x,y,ls.x1,ls.y1,ls.x2,ls.y2,ls.order,self.disxLapLsInf,self.disyLapLsInf)
            disx[0] = disx[0] + self.disxLapLsInf[0] * ls.coef[0,0] + self.disxLapLsInf[1] * ls.coef[1,0]
            disy[0] = disy[0] + self.disyLapLsInf[0] * ls.coef[0,0] + self.disyLapLsInf[1] * ls.coef[1,0]
        for ld in self.ldListBes:
            disbesldho(x,y,ld.x1,ld.y1,ld.x2,ld.y2,aq.Naquifers,aq.zeropluslab,0,self.disxBesInf,self.disyBesInf)
            disx = disx + ld.coef[0,:] * self.disxBesInf
            disy = disy + ld.coef[0,:] * self.disyBesInf
        for ld in self.ldListLap:
            dislapldho(x,y,ld.x1,ld.y1,ld.x2,ld.y2,ld.order,self.disxLapLdInf,self.disyLapLdInf)
            disx[0] = disx[0] + self.disxLapLdInf[0] * ld.coef[0,0] + self.disxLapLdInf[1] * ld.coef[1,0] + \
                     self.disxLapLdInf[2] * ld.coef[2,0]
            disy[0] = disy[0] + self.disyLapLdInf[0] * ld.coef[0,0] + self.disyLapLdInf[1] * ld.coef[1,0] + \
                     self.disyLapLdInf[2] * ld.coef[2,0]
        if isInside:
            disx[0] = disx[0] + (x - self.xc)
        return [disx,disy]
    def isInside(self,x,y):
        '''Checks if x,y is inside element and returns new points if at corner point'''
        isInside = 0
        if x >= self.xmin and x <= self.xmax and y >= self.ymin and y <= self.ymax:
            z = complex(x,y)
            zdis = abs(z-self.z1) / self.lengths  # Alas, we are comparing only to one side. two would be better
            if min(zdis) < self.tiny:
                for i in range(len(zdis)):
                    if zdis[i] < self.tiny:
                        z = z + self.tiny * self.lengths[i]
                        x = z.real; y=z.imag
                        break
            bigZ = ( 2.0*z - (array(self.z1)+array(self.z2)) )/ (array(self.z2)-array(self.z1))
            bigZmin1 = bigZ - 1.0; bigZplus1 = bigZ + 1.0
            angles = log( bigZmin1 / bigZplus1 ).imag
            angle = sum(angles)
            if angle > pi: isInside = 1
        return [isInside,x,y]
    def layout(self):
        rv = [ self.Nsides+1,[],[] ]
        rv[1] = list( self.xp )
        rv[2] = list( self.yp )
        return rv
    def check(self):
        print 'PolyAreaSink centered at '+str((self.xc,self.yc))+\
                  ' Infiltration '+str(self.parameters[0,0])+' has no unknown parameters'
        return None
    def qzTop(self,x,y):
        qz = 0.0
        if self.isInside(x,y)[0]:
            qz = -self.parameters[0,0]  # As infiltration is positive, but qz postive up
        return qz
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):
        changed = 0; stop = 0; xyznew = zeros(3,'d')
        [isin1,xdum,ydum] = self.isInside(xyz1[0],xyz1[1])
        [isin2,xdum,ydum] = self.isInside(xyz2[0],xyz2[1])
        if ( isin1 and not isin2 ) or ( not isin1 and isin2 ):
            [x1,y1,z1] = xyz1; [x2,y2,z2] = xyz2
            for i in range(self.Nsides):
                if x1 < self.xpmin[i] - step or x1 > self.xpmax[i] + step or\
                   y1 < self.ypmin[i] - step or y1 > self.ypmax[i] + step:
                   continue  # Not intersecting this one
                z1 = complex(x1,y1)
                z2 = complex(x2,y2)
                Z1 = (2.0*z1 - (self.z1[i]+self.z2[i]))/(self.z2[i]-self.z1[i])
                Z2 = (2.0*z2 - (self.z1[i]+self.z2[i]))/(self.z2[i]-self.z1[i])
                if abs(Z1.imag) < 1e-5:  # Point basically on edge, so no reason to change
                    break
                if abs(Z2.imag) < 1e-5:  # Point basically on edge, so no reason to change
                    break
                # If point 1 on one side and point 2 on other side, find out if segment intesects line-sink
                if (Z1.imag > 0 and Z2.imag < 0) or (Z1.imag < 0 and Z2.imag > 0):
                    Xintersec = Z1.real + (0.0-Z1.imag) / (Z2.imag-Z1.imag) * (Z2.real-Z1.real)
                    if abs(Xintersec) <= 1:
                        changed = 1
                        if isin1:  # First point is inside, so put new point outside
                            Znew = complex(Xintersec,-1e-6) # Put 0.000001*0.5*L outside
                        else:  # Stepping inside, so put new point inside
                            Znew = complex(Xintersec,1e-6)
                        znew = ( (self.z2[i]-self.z1[i])*Znew + (self.z1[i]+self.z2[i]) ) / 2.0
                        xnew = znew.real; ynew = znew.imag
                        horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
                        horstepnew = sqrt( (xyz1[0]-xnew)**2 + (xyz1[1]-ynew)**2 )
                        znew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )
                        xyznew[:] = [xnew,ynew,znew]
                        break  # We found the intersection point, so break loop
        return [changed, stop, xyznew]
    def distanceSquaredToElement(self,x,y):
        dissqlist = []
        for i in range(self.Nsides):
            Z = ( 2.0 * complex(x,y) - (self.z1[i] + self.z2[i]) ) / (self.z2[i] - self.z1[i])
            if abs(Z.real) <= 1.0:
                dissq = ( Z.imag * self.lengths[i] / 2.0 )**2
            elif Z.real < -1.0:
                dissq = ( x - self.xp[i] )**2 + ( y - self.yp[i] ) **2
            else:
                dissq = ( x - self.xp[i+1] )**2 + ( y - self.yp[i+1] ) **2
            dissqlist.append(dissq)
        return min(dissqlist)

'''
mllinesinkgenearl.py contains the LineSinkGen class.
This file is part of the TimML library and is distributed under
the GNU LPGL. See the TimML.py file for more details.
(c) Mark Bakker, 2002-2007
'''

from numpy import *
import numpy.linalg as linalg
from mlelement import *
from besselaes import *

class LineSinkGen(Element):
    def __init__(self,modelParent,x1,y1,x2,y2,sigma=[1.0],order=0,layers=[1],label=None,aquiferParentFixed=None,addToModel=1):
	Element.__init__(self,modelParent)
        self.x1 = float(x1); self.y1 = float(y1); self.x2 = float(x2); self.y2 = float(y2);
        self.sigma = array( sigma, 'd' )
        self.order = order
        self.layers = array(layers)
        assert len(self.sigma) == self.order + 1, 'TimML input error, sigma of LineSinkGen must be length order+1'
        self.aquiferParentFixed = aquiferParentFixed
        self.label = label
        self.type = 'linesinkgen'
        self.setCoefs()
        if addToModel:
            self.modelParent.addElement(self)
    def __repr__(self):
	return 'LineSinkGen z1,z2,sigma,order,layers: ' + str((self.x1,self.y1,self.x2,self.y2,self.sigma,self.order,list(self.layers)))
    def setCoefs(self):
        self.xy1 = (self.x1,self.y1); self.xy2 = (self.x2,self.y2)
        self.xmin = min(self.x1,self.x2); self.xmax = max(self.x1,self.x2)
        self.ymin = min(self.y1,self.y2); self.ymax = max(self.y1,self.y2)
        self.thetaNormOut = arctan2(self.y2-self.y1,self.x2-self.x1) - pi/2.0
        self.thetaNormIn = self.thetaNormOut + pi
        self.alpha = arctan2(self.y2-self.y1,self.x2-self.x1)
        self.cosalpha = cos(self.alpha); self.sinalpha = sin(self.alpha)
        self.z1 = complex(self.x1,self.y1); self.z2 = complex(self.x2,self.y2);
        self.xc = 0.5*(self.x1+self.x2); self.yc = 0.5*(self.y1+self.y2)
        self.L = abs(self.z2-self.z1); self.Lover4pi = self.L/(4.0*pi); self.Lover2 = self.L / 2.0
        self.Ndegree = self.order+1
        if self.aquiferParentFixed == None:
            self.aquiferParent = self.modelParent.aq.findAquiferData(self.xc,self.yc)  # Determined at xc,yc
        else:
            self.aquiferParent = self.aquiferParentFixed
      	self.pylayers = self.layers-1; self.NscreenedLayers = len(self.layers)
	if self.NscreenedLayers == 1:  # Screened in only one layer, no unknowns
            self.parameters = self.sigma[:,newaxis]
        else:
            self.parameters = zeros( (self.NscreenedLayers * self.Ndegree,1),'d')
        self.coef = ones((self.NscreenedLayers,self.aquiferParent.Naquifers),'d')
	if self.aquiferParent.Naquifers > 1:   # Multiple aquifers, must compute coefficients
            if self.aquiferParent.type == self.aquiferParent.conf:
                for i in range(self.NscreenedLayers):
                    pylayer = self.pylayers[i]
                    ramat = self.aquiferParent.eigvec[:,1:]  # All eigenvectors, but not first normalized transmissivity
                    ramat = vstack(( ramat[0:pylayer,:], ramat[pylayer+1:,:] ))  # Remove row pylayer
                    rb = - self.aquiferParent.eigvec[:,0]  # Store Tn vector in rb
                    rb = hstack(( rb[0:pylayer], rb[pylayer+1:] )) #  takes rowvector
                    self.coef[i,1:self.aquiferParent.Naquifers] = linalg.solve(ramat,rb)
            elif self.aquiferParent.type == self.aquiferParent.semi:
                print 'higher order line-sinks not yet implemented for semi-confined aquifer'
        self.Ncp = self.Ndegree
        thetacp = arange(pi,0,-pi/self.Ncp) - 0.5 * pi/self.Ncp
        Zcp = zeros( self.Ncp, 'D' )
        Zcp.real = cos(thetacp)
        Zcp.imag = 1e-6  # control point just on positive side (this is handy later on)
        zcp = Zcp * (self.z2 - self.z1) / 2.0 + 0.5 * (self.z1 + self.z2)
        self.xcp = zcp.real; self.ycp = zcp.imag
        self.paramxcoef = zeros( (self.Ndegree,self.aquiferParent.Naquifers), 'd' )
        for ior in range(self.Ndegree):
            self.paramxcoef[ior,:] = \
            sum( self.parameters[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,:] * self.coef, 0 )  # Parameters times coefficients
        self.potLapLsInf = zeros(self.Ndegree,'d'); self.potInf = zeros(self.aquiferParent.Naquifers,'d')
        self.disxLapLsInf = zeros(self.Ndegree,'d'); self.disyLapLsInf = zeros(self.Ndegree,'d');
        self.disxInf = zeros(self.aquiferParent.Naquifers,'d'); self.disyInf = zeros(self.aquiferParent.Naquifers,'d')
    def potentialInfluence(self,aq,x,y):
        '''Returns array of (NscreenedLayers*Order,aq.Naquifers)
        Sequence: first all order 0 linesinks, than all order 1, etc.'''
        rv = zeros((self.NscreenedLayers * self.Ndegree,aq.Naquifers),'d')
        # Order of if-s is slightly different than for wells, because Laplace part is not computed separately
        if self.aquiferParent.type == self.aquiferParent.conf:
            if self.aquiferParent == aq:  # Same confined aquifer
                lab = self.aquiferParent.zeropluslab
                for iorder in range(self.Ndegree):
                    potbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers,lab,iorder,self.potInf)
                    # Call FORTRAN extension
                    rv[iorder*self.NscreenedLayers:iorder*self.NscreenedLayers+self.NscreenedLayers,:] = \
                        self.coef * self.potInf
            else:  # Different confined aquifer, Laplace part only
                potlaplsho(x,y,self.x1,self.y1,self.x2,self.y2,self.order,self.potLapLsInf)
                for iorder in range(self.Ndegree):
                    rv[iorder*self.NscreenedLayers:iorder*self.NscreenedLayers+self.NscreenedLayers,0] = \
                        self.potLapLsInf[iorder]
        elif self.aquiferParent.type == self.aquiferParent.semi:
            print 'not yet implemented'
        return rv
    def dischargeInfluence(self,aq,x,y):
        '''Returns two arrays of (order,aq.Naquifers)'''
        rvx = zeros((self.NscreenedLayers * self.Ndegree,aq.Naquifers),'d')
        rvy = zeros((self.NscreenedLayers * self.Ndegree,aq.Naquifers),'d')
        if self.aquiferParent.type == self.aquiferParent.conf:
            if self.aquiferParent == aq:  # Same confined aquifer
                lab = self.aquiferParent.zeropluslab
                for ior in range(self.Ndegree):
                    disbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers,lab,ior,self.disxInf,self.disyInf)
                    # Call FORTRAN extension
                    rvx[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,:] = \
                        self.coef * self.disxInf
                    rvy[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,:] = \
                        self.coef * self.disyInf                
            else:  # Different confined aquifer, Laplace part only
                dislaplsho(x,y,self.x1,self.y1,self.x2,self.y2,self.order,self.disxLapLsInf,self.disyLapLsInf)
                for ior in range(self.Ndegree):
                    rvx[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,0] = \
                        self.disxLapLsInf[ior]
                    rvy[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,0] = \
                        self.disyLapLsInf[ior]
        elif self.aquiferParent.type == self.aquiferParent.semi:
            print 'not yet implemented'
        return [rvx,rvy]
    def potentialCollection(self,potsum,potadd,elementList,aq,x,y):
        for el in elementList:
            if el.aquiferParent.type == el.aquiferParent.conf:
                if el.aquiferParent == aq:  # Same confined aquifer
                    lab = el.aquiferParent.zeropluslab
                    for ior in range(el.Ndegree):
                        potbeslsho(x,y,el.x1,el.y1,el.x2,el.y2,el.aquiferParent.Naquifers,lab,ior,el.potInf)
                        potsum = potsum + el.paramxcoef[ior,:] * el.potInf
                elif aq.type == aq.conf:  # Different confined aquifer, Laplace part only
                    potlaplsho(x,y,el.x1,el.y1,el.x2,el.y2,el.order,el.potLapLsInf)  # Call FORTRAN extension
                    for ior in range(el.Ndegree):
                        potsum[0] = potsum[0] + el.paramxcoef[ior,0] * el.potLapLsInf[ior]
            elif el.aquiferParent.type == el.aquiferParent.semi:
                print 'not yet implemented'
        return potsum
    def dischargeCollection(self,dissum,disadd,elementList,aq,x,y):
        for el in elementList:
            if el.aquiferParent.type == el.aquiferParent.conf:
                if el.aquiferParent == aq:  # Same confined aquifer
                    lab = el.aquiferParent.zeropluslab
                    for ior in range(el.Ndegree):
                        disbeslsho(x,y,el.x1,el.y1,el.x2,el.y2,el.aquiferParent.Naquifers,lab,ior,disadd[0,:],disadd[1,:])
                        dissum = dissum + el.paramxcoef[ior,:] * disadd
                elif aq.type == aq.conf:  # Different confined aquifer, Laplace part only
                    dislaplsho(x,y,el.x1,el.y1,el.x2,el.y2,el.order,el.disxLapLsInf,el.disyLapLsInf)
                    for ior in range(el.Ndegree):
                        dissum[0,0] = dissum[0,0] + el.paramxcoef[ior,0] * el.disxLapLsInf[ior]
                        dissum[1,0] = dissum[1,0] + el.paramxcoef[ior,0] * el.disyLapLsInf[ior]
            elif el.aquiferParent.type == el.aquiferParent.semi:
                print 'not yet implemented'
        return dissum
    def layout(self):
        rv = [ 2,[],[] ]
        rv[1] = [ self.x1,self.x2 ]
        rv[2] = [ self.y1,self.y2 ]
        return rv
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):
        changed = 0; stop = 0; xyznew = 0.0
        if pyLayer == self.pylayers.any():  # In layer line-sink is screened. This is true if self.pylayers contains pyLayer
            [x1,y1,z1] = xyz1; [x2,y2,z2] = xyz2
            if self.NscreenedLayers == 1:  # Find intersection point
                z1 = complex(x1,y1)
                z2 = complex(x2,y2)
                Z1 = (2.0*z1 - (self.z1+self.z2))/(self.z2-self.z1)
                Z2 = (2.0*z2 - (self.z1+self.z2))/(self.z2-self.z1)
                # If point 1 on one side and point 2 on other side, find out if segment intesects line-sink
                if (Z1.imag > 0 and Z2.imag < 0) or (Z1.imag < 0 and Z2.imag > 0):
                    Xintersec = Z1.real + (0.0-Z1.imag) / (Z2.imag-Z1.imag) * (Z2.real-Z1.real)
                    if abs(Xintersec) <= 1:
                        changed = 1
                        if ( self.parameters[0,0] > 0 and idir > 0 ) or ( self.parameters[0,0] < 0 and idir < 0 ): # Takes water out, stop trace (should check for jump). Allows for zero to allow for line-sinks used elsewhere
                            Znew = complex(Xintersec,0.0)
                            znew = ( (self.z2-self.z1)*Znew + (self.z1+self.z2) ) / 2.0
                            xnew = znew.real; ynew = znew.imag
                            horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
                            horstepnew = sqrt( (xyz1[0]-xnew)**2 + (xyz1[1]-ynew)**2 )
                            znew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )
                            xyznew = array([xnew,ynew,znew]) 
                            stop = 1
                        else: # Line-sink puts water in (should jump streamline down)
                            if Z1.imag < 0:  # Is below line-sink, put above
                                Znew = complex(Xintersec,0.001) # Put 0.001*0.5*L off line-sink
                            else:  #Is above line-sink, put below
                                Znew = complex(Xintersec,0.001)
                            znew = ( (self.z2-self.z1)*Znew + (self.z1+self.z2) ) / 2.0
                            xnew = znew.real; ynew = znew.imag
                            xyznew = array([xnew,ynew,xyz1[2]])  # Should modify z
                            changed = 1
                elif abs(Z2.imag) < 2e-6 * step / self.L and abs(Z2.real) <= 1:  # Point 2 is on linesink
                    stop = 1
#            elif abs( sum(sum(self.parameters)) ) < 1e-10:  # Screened in multiple layers, but zero total discharge
            elif abs( sum(self.parameters) ) < 1e-10:  # Screened in multiple layers, but zero total discharge
                # Do nothing for now; when used for inhomogeneity this is taken care off by aquifer class
                pass
##                z1 = complex(x1,y1); z2 = complex(x2,y2)
##                Z1 = (2.0*z1 - (self.z1+self.z2))/(self.z2-self.z1)
##                Z2 = (2.0*z2 - (self.z1+self.z2))/(self.z2-self.z1)
##                # If point 1 on one side and point 2 on other side, find out if segment intesects line-sink
##                if (Z1.imag > 0 and Z2.imag < 0) or (Z1.imag < 0 and Z2.imag > 0):
##                    Xintersec = Z1.real + (0.0-Z1.imag) / (Z2.imag-Z1.imag) * (Z2.real-Z1.real)
##                    if abs(Xintersec) <= 1:  # There is an intersection point
##                        changed = 1
##                        theta = self.thetaNormOut
##                        if Z1.imag < 0:  # old point on outside, Znew1 just out, Znew2 just in
##                            Znew1 = complex(Xintersec,-1e-6); Znew2 = complex(Xintersec,1e-6)
##                            if idir > 0: theta = self.thetaNormOut + pi  # Such that theta points in supposed direction of flow
##                        else:
##                            Znew1 = complex(Xintersec,1e-6); Znew2 = complex(Xintersec,-1e-6)
##                            if idir < 0: theta = self.thetaNormOut + pi
##                        znew1 = ( (self.z2-self.z1)*Znew1 + (self.z1+self.z2) ) / 2.0  # New complex coordinate
##                        znew2 = ( (self.z2-self.z1)*Znew2 + (self.z1+self.z2) ) / 2.0
##                        xnew1 = znew1.real; ynew1 = znew1.imag
##                        xnew2 = znew2.real; ynew2 = znew2.imag
##                        horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
##                        horstepnew = sqrt( (xyz1[0]-xnew1)**2 + (xyz1[1]-ynew1)**2 )
##                        zvertnew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )  # Just to be confusing, this is z (vertical coordinate) 
##                        disnorm = self.modelParent.dischargeNormVector(xnew1,ynew1,theta)[pyLayer]  # Normal discharge at xnew1,ynew1
##                        if disnorm < 0:
##                            # Normal discharge at boundary points the other way; can only be if step size too large
##                            # Half step to boundary
##                            xyznew = 0.5 * ( xyz1 + array([xnew1,ynew1,zvertnew]) )
##                        else:
##                            dis = self.modelParent.dischargeNormBelowZ(xnew1,ynew1,zvertnew,theta)
##                            [xnew,ynew,znew] = self.modelParent.zForGivenNormalDischarge\
##                                            (zvertnew,xnew1,ynew1,xnew2,ynew2,theta,dis)
##                            xyznew = array([xnew,ynew,znew])
            else:
                stop = 1  # Seems to be only logical step. Stop traceline, there is nothing else we know
        return [changed, stop, xyznew]
    def nearElementOld(self,pyLayer,xyz1,xyz2,step):
        changed = 0; stop = 0; xyznew = zeros(3,'d')
        if pyLayer == self.pylayers.any():  # In layer line-sink is screened. This is true if self.pylayers contains pyLayer
            [x1,y1,z1] = xyz1; [x2,y2,z2] = xyz2
            if x1 < self.xmin - abs(step) or x1 > self.xmax + abs(step) or \
               y1 < self.ymin - abs(step) or y1 > self.ymax + abs(step):
               return [changed, stop, xyznew]  # Away from line-sink
            if self.NscreenedLayers == 1:  # Find intersection point
                z1 = complex(x1,y1)
                z2 = complex(x2,y2)
                Z1 = (2.0*z1 - (self.z1+self.z2))/(self.z2-self.z1)
                Z2 = (2.0*z2 - (self.z1+self.z2))/(self.z2-self.z1)
                # If point 1 on one side and point 2 on other side, find out if segment intesects line-sink
                if (Z1.imag > 0 and Z2.imag < 0) or (Z1.imag < 0 and Z2.imag > 0):
                    Xintersec = Z1.real + (0.0-Z1.imag) / (Z2.imag-Z1.imag) * (Z2.real-Z1.real)
                    if abs(Xintersec) <= 1:
                        changed = 1
                        if ( self.parameters[0,0] > 0 and step > 0 ) or ( self.parameters[0,0] < 0 and step < 0 ): # Takes water out, stop trace (should check for jump). Allows for zero to allow for line-sinks used elsewhere
                            Znew = complex(Xintersec,0.0)
                            znew = ( (self.z2-self.z1)*Znew + (self.z1+self.z2) ) / 2.0
                            xnew = znew.real; ynew = znew.imag
                            horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
                            horstepnew = sqrt( (xyz1[0]-xnew)**2 + (xyz1[1]-ynew)**2 )
                            znew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )
                            xyznew[:] = [xnew,ynew,znew] 
                            stop = 1
                            print 'reached line-sink'
                        else: # Line-sink puts water in (should jump streamline down)
                            if Z1.imag < 0:  # Is below line-sink, put above
                                Znew = complex(Xintersec,0.001) # Put 0.001*0.5*L off line-sink
                            else:  #Is above line-sink, put below
                                Znew = complex(Xintersec,0.001)
                            znew = ( (self.z2-self.z1)*Znew + (self.z1+self.z2) ) / 2.0
                            xnew = znew.real; ynew = znew.imag
                            xyznew[:] = [xnew,ynew,xyz1[2]]  # Should modify z
                            changed = 1
                elif abs(Z2.imag) < 2e-6 * abs(step) / self.L and abs(Z2.real) <= 1:  # Point 2 is on linesink
                    stop = 1
                    print 'reached line-sink'
#            elif abs( sum(sum(self.parameters)) ) < 1e-12:  # Screened in multiple layers, but zero total discharge
            elif abs( sum(self.parameters) ) < 1e-12:  # Screened in multiple layers, but zero total discharge
                # Do nothing for now. Is for now covered in inhomogeneity section.
                # Needs modification for cracks.
                pass
                
##                z1 = complex(x1,y1); z2 = complex(x2,y2)
##                Z1 = (2.0*z1 - (self.z1+self.z2))/(self.z2-self.z1)
##                Z2 = (2.0*z2 - (self.z1+self.z2))/(self.z2-self.z1)
##                # If point 1 on one side and point 2 on other side, find out if segment intesects line-sink
##                if (Z1.imag > 0 and Z2.imag < 0) or (Z1.imag < 0 and Z2.imag > 0):
##                    Xintersec = Z1.real + (0.0-Z1.imag) / (Z2.imag-Z1.imag) * (Z2.real-Z1.real)
##                    if abs(Xintersec) <= 1:  # There is an intersection point
##                        changed = 1
##                        theta = self.thetaNormOut
##                        if Z1.imag < 0:  # old point on outside, Znew1 just out, Znew2 just in
##                            Znew1 = complex(Xintersec,-1e-6); Znew2 = complex(Xintersec,1e-6)
##                            if step > 0: theta = self.thetaNormOut + pi  # Such that theta points in supposed direction of flow
##                        else:
##                            Znew1 = complex(Xintersec,1e-6); Znew2 = complex(Xintersec,-1e-6)
##                            if step < 0: theta = self.thetaNormOut + pi
##                        znew1 = ( (self.z2-self.z1)*Znew1 + (self.z1+self.z2) ) / 2.0  # New complex coordinate
##                        znew2 = ( (self.z2-self.z1)*Znew2 + (self.z1+self.z2) ) / 2.0
##                        xnew1 = znew1.real; ynew1 = znew1.imag
##                        xnew2 = znew2.real; ynew2 = znew2.imag
##                        horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
##                        horstepnew = sqrt( (xyz1[0]-xnew1)**2 + (xyz1[1]-ynew1)**2 )
##                        zvertnew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )  # Just to be confusing, this is z (vertical coordinate) 
##                        disnorm = self.modelParent.dischargeNormVector(xnew1,ynew1,theta)[pyLayer]  # Normal discharge at xnew1,ynew1
##                        if disnorm < 0:
##                            # Normal discharge at boundary points the other way; can only be if step size too large
##                            # Half step to boundary
##                            xyznew = 0.5 * ( xyz1 + array([xnew1,ynew1,zvertnew]) )
##                        else:
##                            dis = self.modelParent.dischargeNormBelowZ(xnew1,ynew1,zvertnew,theta)
##                            [xnew,ynew,znew] = self.modelParent.zForGivenNormalDischarge\
##                                            (zvertnew,xnew1,ynew1,xnew2,ynew2,theta,dis)
##                            xyznew[:] = [xnew,ynew,znew]
            else:
                stop = 1  # Seems to be only logical step. Stop traceline, there is nothing else we know
                print 'reached line-sink'
        return [changed, stop, xyznew]
    def distanceSquaredToElement(self,x,y):
        X = (x-self.xc) * self.cosalpha + (y-self.yc) * self.sinalpha
        if X < -self.Lover2:
            dissq = ( x - self.x1 )*( x - self.x1 ) + ( y - self.y1 )*( y - self.y1 )
        elif X > self.Lover2:
            dissq = ( x - self.x2 )*( x - self.x2 ) + ( y - self.y2 )*( y - self.y2 )
        else:
            Y = -(x-self.xc) * self.sinalpha + (y-self.yc) * self.cosalpha
            dissq = Y * Y
        return dissq
    def getMatrixRows(self,elementList):
        rows=[]
        if self.NscreenedLayers > 1:  # Otherwise no unknowns !
            for icp in range(self.Ncp):
                xc = self.xcp[icp]; yc = self.ycp[icp]
                rowlast = zeros((0),'d')
                lastpylayer = self.pylayers[self.NscreenedLayers-1]; Tlast = self.aquiferParent.T[lastpylayer]
                for e in elementList:
                    rowpart = e.getMatrixCoefficients(self.aquiferParent,lastpylayer,xc,yc,\
                                    lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                    rowlast = hstack(( rowlast, rowpart ))
                potlast = self.modelParent.potentialInLayer(self.aquiferParent,lastpylayer,xc,yc)
                for i in range(self.NscreenedLayers - 1):                
                    row = zeros(0,'d'); T = self.aquiferParent.T[self.pylayers[i]]
                    for e in elementList:
                        rowpart = e.getMatrixCoefficients(self.aquiferParent,self.pylayers[i],xc,yc,\
                                        lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                        row = hstack(( row, rowpart ))
                    row = Tlast * row - T * rowlast    # Multiply with T's and subtract last row
                    row = hstack(( row,  T * potlast - Tlast * \
                        self.modelParent.potentialInLayer(self.aquiferParent,self.pylayers[i],xc,yc) ))  # Add right hand side as last element
                    rows = rows + [row.tolist()]
                row = zeros(0,'d')
                # Last equation is sum of strength parameters for order icp = sigma[icp] (recall order=Ncp)
                rowsection = zeros(self.Ndegree*self.NscreenedLayers,'d')
                rowsection[icp*self.NscreenedLayers:icp*self.NscreenedLayers+self.NscreenedLayers] = 1.0
        # I added ,0 Does that work ?
                paramsum = sum(self.parameters[icp*self.NscreenedLayers:icp*self.NscreenedLayers+self.NscreenedLayers,0])
                for e in elementList:
                    if e == self:
                        row = hstack(( row, rowsection ))
                    else:
                        row = hstack(( row, e.getMatrixCoefficients(self.aqdum,self.ldum,self.xdum,self.ydum,\
                                                 lambda el,aqdum,ldum,xdum,ydum:el.zeroFunction(aqdum,ldum,xdum,ydum)) ))
                row = hstack(( row, self.sigma[icp] - paramsum ))
                rows = rows + [row.tolist()]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,xc,yc,func):
        rv = zeros(0,'d')
        if self.NscreenedLayers > 1:
            rv = func(self,aq,pylayer,xc,yc)
        return rv
    def takeParameters(self,xsol,icount):
        if self.NscreenedLayers > 1:
            for i in range(self.Ncp*self.NscreenedLayers):
                self.parameters[i,0] = self.parameters[i,0] + xsol[icount]
                icount = icount+1
            for ior in range(self.Ndegree):
                self.paramxcoef[ior,:] = \
                sum( self.parameters[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,:] * self.coef, 0 )  # Parameters times coefficients
        return icount
    def check(self):
        if self.NscreenedLayers == 1:
            print 'Linesink from '+str(self.xy1)+' to '+str(self.xy2)+' Layer: '+str(self.layers[0])+\
                  ' Strength: '+str(self.sigma)+'Order: '+str(self.order)+' has no unknown parameters'
        else:
            print 'Linesink from '+str(self.xy1)+' to '+str(self.xy2)
            for ior in range(self.Ndegree):
                sig = sum(self.parameters[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,0])
                print 'Specified sigma'+str(ior)+': '+str(self.sigma[ior])+' Computed sigma: '+str(sig)
            for icp in range(self.Ncp):
                print 'Control point '+str(icp)
                for i in range(self.NscreenedLayers):
                    print 'Layer '+str(self.layers[i])+' Head: '+\
                      str(self.modelParent.head(self.layers[i],self.xcp[icp],self.ycp[icp],self.aquiferParent))

class DoubleLineSinkGen(Element):
    def __init__(self,modelParent,x1,y1,x2,y2,order,aqin,aqout,label=None):
	Element.__init__(self,modelParent)
        self.x1 = float(x1); self.y1 = float(y1); self.x2 = float(x2); self.y2 = float(y2);
        assert order < 9, "TimML Input error: Maximum order is 8"
        self.order = order
        self.aqin = aqin; self.aqout = aqout
        if not aqin.fakesemi:
            assert self.aqin.Naquifers == self.aqout.Naquifers, 'TimML input error: DoubleLineSinkGen number of aquifers must be same in and out'
        else:
            assert self.aqin.Naquifers == self.aqout.Naquifers+1, 'TimML input error: DoubleLineSinkGen number of aquifers must be in = out+1'
        self.label = label
        self.type = 'doublelinesinkgen'
        self.setCoefs()
        self.modelParent.addElement(self)
    def __repr__(self):
	return 'DoubleLineSinkGen z1,z2,order: ' + str((self.x1,self.y1,self.x2,self.y2,self.order))
    def setCoefs(self):
        self.xy1 = (self.x1,self.y1); self.xy2 = (self.x2,self.y2)
        self.xmin = min(self.x1,self.x2); self.xmax = max(self.x1,self.x2)
        self.ymin = min(self.y1,self.y2); self.ymax = max(self.y1,self.y2)
        self.thetaNormOut = arctan2(self.y2-self.y1,self.x2-self.x1) - pi/2.0
        self.cosout = cos( self.thetaNormOut ); self.sinout = sin( self.thetaNormOut )
        self.thetaNormIn = self.thetaNormOut + pi
        self.alpha = arctan2(self.y2-self.y1,self.x2-self.x1)
        self.cosalpha = cos(self.alpha); self.sinalpha = sin(self.alpha)
        self.z1 = complex(self.x1,self.y1); self.z2 = complex(self.x2,self.y2);
        self.xc = 0.5*(self.x1+self.x2); self.yc = 0.5*(self.y1+self.y2)
        self.L = abs(self.z2-self.z1); self.Lover4pi = self.L/(4.0*pi); self.Lover2 = self.L / 2.0
        self.Ndegree = self.order+1
        self.aquiferParent = None # Actually means nothing here, so setting to None is save as it will throw error if used
        self.Naqin = self.aqin.Naquifers; self.Naqout = self.aqout.Naquifers
        self.paramin = zeros( ( self.Ndegree, self.Naqin ),'d'); self.Nparamin = self.Ndegree * (self.Naqin-1)  # as Laplace part is zero
        self.paramout = zeros( ( self.Ndegree, self.Naqout ),'d'); self.Nparamout = self.Ndegree * (self.Naqout-1)
        # I don't think there will be coefs
        self.Ncp = self.Ndegree
        self.Nparam = (self.Naqin - 1) * self.Ndegree + (self.Naqout - 1) * self.Ndegree
        # self.halfNun = self.Ncp * ( self.Naqin - 1 )  # All inside parameters are stored first
        thetacp = arange(pi,0,-pi/self.Ncp) - 0.5 * pi/self.Ncp
        Zcp = zeros( self.Ncp, 'D' )
        Zcp.real = cos(thetacp)
        Zcp.imag = 1e-6  # control point just on inside
        zcp = Zcp * (self.z2 - self.z1) / 2.0 + 0.5 * (self.z1 + self.z2)
        self.xcpin = zcp.real; self.ycpin = zcp.imag
        Zcp.imag = -1e-6  # control point just on outside
        zcp = Zcp * (self.z2 - self.z1) / 2.0 + 0.5 * (self.z1 + self.z2)
        self.xcpout = zcp.real; self.ycpout = zcp.imag
        # No paramxcoef either I guess
#        self.potInfin = zeros(self.Naqin,'d'); self.potInfout = zeros(self.Naqout,'d')
        self.potInfin = zeros(self.Naqin-1,'d'); self.potInfout = zeros(self.Naqout-1,'d')
        self.disxLapLsInf = zeros(self.Ndegree,'d'); self.disyLapLsInf = zeros(self.Ndegree,'d');
        self.disxInfin = zeros(self.Naqin-1,'d'); self.disyInfin = zeros(self.Naqin-1,'d')
        self.disxInfout = zeros(self.Naqout-1,'d'); self.disyInfout = zeros(self.Naqout-1,'d')
    def potentialInfluence(self,aq,x,y):
        '''Returns array of (2*NscreenedLayers*(Order+1),aq.Naquifers), with inside first
        Sequence: first all order 0 linesinks, than all order 1, etc.'''
        rv = zeros( (self.Ndegree, aq.Naquifers), 'd')
        # Order of if-s is slightly different than for wells, because Laplace part is not computed separately
        if self.aqin == aq:  # Inside aquifer
            for ior in range(self.Ndegree):
                potbesonlylsho(x,y,self.x1,self.y1,self.x2,self.y2,self.Naqin,self.aqin.lab,ior,self.potInfin)
                rv[ ior, 1: ] = self.potInfin[:]
        elif self.aqout == aq: # Outside aquifer
            for ior in range(self.Ndegree):
                potbesonlylsho(x,y,self.x1,self.y1,self.x2,self.y2,self.Naqout,self.aqout.lab,ior,self.potInfout)
                rv[ ior, 1: ] = self.potInfout[:]
        elif aq.type == aq.conf:  # Different confined aquifer, Laplace part only, which is zero
            dum = 0
        else:
            print 'not yet implemented, linesinkgen pot'
        return rv
    def dischargeInfluence(self,aq,x,y):
        '''Returns array of (2*NscreenedLayers*(Order+1),aq.Naquifers), with inside first
        Sequence: first all order 0 linesinks, than all order 1, etc.'''
        rvx = zeros( (self.Ndegree, aq.Naquifers), 'd')
        rvy = zeros( (self.Ndegree, aq.Naquifers), 'd')
        # Order of if-s is slightly different than for wells, because Laplace part is not computed separately
        if self.aqin == aq:  # Inside aquifer
            for ior in range(self.Ndegree):
                disbesonlylsho(x,y,self.x1,self.y1,self.x2,self.y2,self.Naqin,self.aqin.lab,ior,self.disxInfin,self.disyInfin)
                rvx[ ior, 1: ] = self.disxInfin[:]
                rvy[ ior, 1: ] = self.disyInfin[:]
        elif self.aqout == aq:
            for ior in range(self.Ndegree):
                disbesonlylsho(x,y,self.x1,self.y1,self.x2,self.y2,self.Naqout,self.aqout.lab,ior,self.disxInfout,self.disyInfout)
                rvx[ ior, 1: ] = self.disxInfout[:]
                rvy[ ior, 1: ] = self.disyInfout[:]
        elif aq.type == aq.conf:  # Different confined aquifer, Laplace part only, which is zero
            dum = 0
        else:
            print 'not yet implemented'
        return [rvx,rvy]
    def potentialContribution(self,aq,x,y):
        '''Returns array of potentialContribution. Needs to be overloaded cause there is inside and outside'''
        if aq == self.aqin:
            return sum( self.paramin * self.potentialInfluence(aq,x,y), 0 )
        elif aq == self.aqout:
            return sum( self.paramout * self.potentialInfluence(aq,x,y), 0 )
        else:  # Really all zeros, as Laplace part is zero
            return zeros(aq.Naquifers,'d')
    def dischargeContribution(self,aq,x,y):
        '''Returns matrix with two rowvectors of dischargeContributions Qx and Qy. Needs to be overloaded cause there is inside and outside'''
        disInf = self.dischargeInfluence(aq,x,y)
        if aq == self.aqin:
            return array([ sum( self.paramin * disInf[0], 0 ), \
                           sum( self.paramin * disInf[1], 0 ) ])
        elif aq == self.aqout:
            return array([ sum( self.paramout * disInf[0], 0 ), \
                           sum( self.paramout * disInf[1], 0 ) ])
        else:  # Really all zeros, as Laplace part is zero
            return zeros((2,aq.Naquifers),'d')
    def potentialInfluenceInLayer(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in pylayer as 1D array (1 value per parameter)'''
        rv = zeros( self.Nparamin + self.Nparamout, 'd' )
        if aq == self.aqin:
            potInf = self.potentialInfluence(aq,x,y)[:,1:]
            potInf = potInf * aq.eigvec[pylayer,1:]  # again excluding the Laplace part
            rv[:self.Nparamin] = potInf.ravel()
        elif aq == self.aqout:
            potInf = self.potentialInfluence(aq,x,y)[:,1:]
            potInf = potInf * aq.eigvec[pylayer,1:]  # again excluding the Laplace part
            rv[self.Nparamin:] = potInf.ravel()
        return rv
    def dischargeInfluenceInLayer(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in pylayer as 1D array (1 value per parameter)'''
        rvx,rvy = zeros( self.Nparamin + self.Nparamout, 'd' ), zeros( self.Nparamin + self.Nparamout, 'd' )
        if aq == self.aqin:
            disxInf,disyInf = self.dischargeInfluence(aq,x,y)
            disx = disxInf[:,1:] * aq.eigvec[pylayer,1:]  # again excluding the Laplace part
            disy = disyInf[:,1:] * aq.eigvec[pylayer,1:]
            rvx[:self.Nparamin] = disx.ravel(); rvy[:self.Nparamin] = disy.ravel()
        elif aq == self.aqout:
            disxInf,disyInf = self.dischargeInfluence(aq,x,y)
            disx = disxInf[:,1:] * aq.eigvec[pylayer,1:]  # again excluding the Laplace part
            disy = disyInf[:,1:] * aq.eigvec[pylayer,1:]
            rvx[self.Nparamin:] = disx.ravel(); rvy[self.Nparamin:] = disy.ravel()
        return [rvx,rvy]
    def potentialInfluenceAllLayers(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in pylayer as 1D array (1 value per parameter)'''
        rv = zeros( (aq.Naquifers, self.Nparamin + self.Nparamout), 'd' )
        if aq == self.aqin:
            potInf = self.potentialInfluence(aq,x,y)[:,1:]
            for i in range( aq.Naquifers ):
                potdum = potInf * aq.eigvec[i,1:]  # again excluding the Laplace part
                rv[i,:self.Nparamin] = potdum.ravel()
        elif aq == self.aqout:
            potInf = self.potentialInfluence(aq,x,y)[:,1:]
            for i in range( aq.Naquifers ):
                potdum = potInf * aq.eigvec[i,1:]  # again excluding the Laplace part
                rv[i,self.Nparamin:] = potdum.ravel()
        return rv
    def dischargeInfluenceAllLayers(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in pylayer as 1D array (1 value per parameter)'''
        rvx = zeros( (aq.Naquifers, self.Nparamin + self.Nparamout), 'd' )
        rvy = zeros( (aq.Naquifers, self.Nparamin + self.Nparamout), 'd' )
        if aq == self.aqin:
            disxInf,disyInf = self.dischargeInfluence(aq,x,y)
            for i in range( aq.Naquifers ):
                disx = disxInf[:,1:] * aq.eigvec[i,1:]  # again excluding the Laplace part
                disy = disyInf[:,1:] * aq.eigvec[i,1:]
                rvx[i,:self.Nparamin] = disx.ravel(); rvy[i,:self.Nparamin] = disy.ravel()
        elif aq == self.aqout:
            disxInf,disyInf = self.dischargeInfluence(aq,x,y)
            for i in range( aq.Naquifers ):
                disx = disxInf[:,1:] * aq.eigvec[i,1:]  # again excluding the Laplace part
                disy = disyInf[:,1:] * aq.eigvec[i,1:]
                rvx[i,self.Nparamin:] = disx.ravel(); rvy[i,self.Nparamin:] = disy.ravel()
        return [rvx,rvy]
    def potentialInfluenceSpecLayers(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in pylayer as 1D array (1 value per parameter)'''
        rv = zeros( (len(pylayer), self.Nparamin + self.Nparamout), 'd' )
        if aq == self.aqin:
            potInf = self.potentialInfluence(aq,x,y)[:,1:]
            for i in range(len(pylayer)):
                potdum = potInf * aq.eigvec[pylayer[i],1:]  # again excluding the Laplace part
                rv[i,:self.Nparamin] = potdum.ravel()
        elif aq == self.aqout:
            potInf = self.potentialInfluence(aq,x,y)[:,1:]
            for i in range(len(pylayer)):
                potdum = potInf * aq.eigvec[pylayer[i],1:]  # again excluding the Laplace part
                rv[i,self.Nparamin:] = potdum.ravel()
        return rv
    def dischargeInfluenceSpecLayers(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in pylayer as 1D array (1 value per parameter)'''
        rvx = zeros( (len(pylayer), self.Nparamin + self.Nparamout), 'd' )
        rvy = zeros( (len(pylayer), self.Nparamin + self.Nparamout), 'd' )
        if aq == self.aqin:
            disxInf,disyInf = self.dischargeInfluence(aq,x,y)
            for i in range( len(pylayer) ):
                disx = disxInf[:,1:] * aq.eigvec[pylayer[i],1:]  # again excluding the Laplace part
                disy = disyInf[:,1:] * aq.eigvec[pylayer[i],1:]
                rvx[i,:self.Nparamin] = disx.ravel(); rvy[i,:self.Nparamin] = disy.ravel()
        elif aq == self.aqout:
            disxInf,disyInf = self.dischargeInfluence(aq,x,y)
            for i in range( len(pylayer) ):
                disx = disxInf[:,1:] * aq.eigvec[pylayer[i],1:]  # again excluding the Laplace part
                disy = disyInf[:,1:] * aq.eigvec[pylayer[i],1:]
                rvx[i,self.Nparamin:] = disx.ravel(); rvy[i,self.Nparamin:] = disy.ravel()
        return [rvx,rvy]
    def layout(self):
        rv = [ 2,[],[] ]
        rv[1] = [ self.x1,self.x2 ]
        rv[2] = [ self.y1,self.y2 ]
        return rv
    def distanceSquaredToElement(self,x,y):
        X = (x-self.xc) * self.cosalpha + (y-self.yc) * self.sinalpha
        if X < -self.Lover2:
            dissq = ( x - self.x1 )*( x - self.x1 ) + ( y - self.y1 )*( y - self.y1 )
        elif X > self.Lover2:
            dissq = ( x - self.x2 )*( x - self.x2 ) + ( y - self.y2 )*( y - self.y2 )
        else:
            Y = -(x-self.xc) * self.sinalpha + (y-self.yc) * self.cosalpha
            dissq = Y * Y
        return dissq
    def getMatrixRows(self,elementList):
        rows=[]

        if self.Naqout > 1:        
            pylayerin = range(1,self.Naqin); pylayerout = range(1,self.Naqout)
            if self.aqin.fakesemi: pylayerin = range(2,self.Naqin)
            in1 = pylayerin[0]; in2 = pylayerin[-1]
            out1 = pylayerout[0]; out2 = pylayerout[-1]

            for icp in range(self.Ncp):
                xcpin = self.xcpin[icp]; ycpin = self.ycpin[icp]
                xcpout = self.xcpout[icp]; ycpout = self.ycpout[icp]

             #   pylayerout = pylayer; pylayerin = pylayer
             #   if self.aqin.fakesemi: pylayerin = pylayerin + 1

                rowin = (self.Naqout-1)*[[]]; rowout = (self.Naqout-1)*[[]]
                for e in elementList:  # headin = headout for pylayers 1 through Naq
                    rowpartin = e.getMatrixCoefficients(self.aqin,pylayerin,xcpin,ycpin,\
                                lambda el,aq,pylayer,x,y:el.potentialInfluenceSpecLayers(aq,pylayer,x,y))
                    rowpartout = e.getMatrixCoefficients(self.aqout,pylayerout,xcpout,ycpout,\
                                lambda el,aq,pylayer,x,y:el.potentialInfluenceSpecLayers(aq,pylayer,x,y))
                    if size(rowpartout) > 0: 
                        rowpartin = rowpartin.tolist()
                        rowpartout = rowpartout.tolist()
                        for pylayer in range(self.Naqout-1):
                            # Need to modify here for fake semi
                            rowin[pylayer] = rowin[pylayer] + rowpartin[pylayer]  # Note that append doesn't work here
                            rowout[pylayer] = rowout[pylayer] + rowpartout[pylayer]
                row = self.aqout.Tcol[out1:] * array(rowin) - self.aqin.Tcol[in1:] * array(rowout)
                rhs = \
                    self.aqin.T[in1:] * self.aqout.T[out1:] * ( self.aqout.hstar - self.aqin.hstar ) + \
                    self.aqin.T[in1:] * self.modelParent.potentialCollection(xcpout,ycpout,self.aqout)[out1:] - \
                    self.aqout.T[out1:] * self.modelParent.potentialCollection(xcpin,ycpin,self.aqin)[in1:]
                row = row.tolist()
                for ieq in range(0,self.Naqout-1):
                    row[ieq].append( rhs[ieq] )
                    rows.append( row[ieq] )


        pylayerin = range(1,self.Naqin); pylayerout = range(1,self.Naqout)
        if self.aqin.fakesemi: pylayerout = range(0,self.Naqout)
        in1 = pylayerin[0]; in2 = pylayerin[-1]
        out1 = pylayerout[0]; out2 = pylayerout[-1]

        for icp in range(self.Ncp):
            xcpin = self.xcpin[icp]; ycpin = self.ycpin[icp]
            xcpout = self.xcpout[icp]; ycpout = self.ycpout[icp]

##            for pylayer in range(1,self.Naqin):  # flowin = flowout is applied at one more layer in case of fakesemi
##                pylayerout = pylayer; pylayerin = pylayer
##                if self.aqin.fakesemi:  # so here it is different too
##                    pylayerout = pylayerout - 1
##                rowin = []; rowout = []
            rowin = (self.Naqin-1)*[[]]; rowout = (self.Naqin-1)*[[]]
            for e in elementList:  # flowin = flowout for pylayers 1 through Naq
                rowqxqyin = e.getMatrixCoefficients(self.aqin,pylayerin,xcpin,ycpin,\
                            lambda el,aq,pylayer,x,y:el.dischargeInfluenceSpecLayers(aq,pylayer,x,y))
                rowqxqyout = e.getMatrixCoefficients(self.aqout,pylayerout,xcpout,ycpout,\
                            lambda el,aq,pylayer,x,y:el.dischargeInfluenceSpecLayers(aq,pylayer,x,y))
                
                if size(rowqxqyin) > 0:
                    rowpartin = rowqxqyin[0] * self.cosout + rowqxqyin[1] * self.sinout
                    rowpartin = rowpartin.tolist()
                    for pylayer in range(self.Naqin-1):
                        # Need to modify here for fake semi
                        rowin[pylayer] = rowin[pylayer] + rowpartin[pylayer]
                if size(rowqxqyout) > 0:  # Does this ever happen ?????? I don't think we need to check...
                    rowpartout = rowqxqyout[0] * self.cosout + rowqxqyout[1] * self.sinout
                    rowpartout = rowpartout.tolist()
                    for pylayer in range(self.Naqin-1):
                        # Need to modify here for fake semi
                        rowout[pylayer] = rowout[pylayer] + rowpartout[pylayer]

            row = array(rowin) - array(rowout)
            row = row.tolist()
            rhs = self.modelParent.dischargeNormVector(xcpout,ycpout,self.thetaNormOut,self.aqout)[out1:] -\
                  self.modelParent.dischargeNormVector(xcpin,ycpin,self.thetaNormOut,self.aqin)[in1:]
            for ieq in range(0,self.Naqin-1):
                row[ieq].append( rhs[ieq] )
                rows.append( row[ieq] )
        return rows
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        for i in range( self.Ndegree ):
            self.paramin[i,1:] = self.paramin[i,1:] + xsol[ icount : icount + self.Naqin - 1 ]
            icount = icount + self.Naqin-1
        for i in range( self.Ndegree ):
            self.paramout[i,1:] = self.paramout[i,1:] + xsol[ icount : icount + self.Naqout - 1 ]
            icount = icount + self.Naqout-1
        return icount
    def check(self):
        print 'DoubleLineSinkGen from '+str(self.xy1)+' to '+str(self.xy2)
        for icp in range(self.Ncp):
            print 'Control point '+str(icp)
            for i in range(2,self.Naqout+1):
                #i_in = i
                #if self.aqin.fakesemi: i_in = i + 1
                print 'Layer '+str(i)+' Head in,out: '+\
                  str(self.modelParent.head(i,self.xcpin[icp],self.ycpin[icp],self.aqin)) + ' ' +\
                  str(self.modelParent.head(i,self.xcpout[icp],self.ycpout[icp],self.aqout))
            for i in range(0,self.Naqout):
                i_in = i
                if self.aqin.fakesemi: i_in = i + 1
                print 'Layer '+str(i+1)+' Qnorm in,out: '+\
                  str(self.modelParent.dischargeNormInLayer(self.aqin,i_in,self.xcpin[icp],self.ycpin[icp],self.thetaNormOut)) + ' ' +\
                  str(self.modelParent.dischargeNormInLayer(self.aqout,i,self.xcpout[icp],self.ycpout[icp],self.thetaNormOut))
        return None
    def zeroFunction(self,aqdum,ldum,xdum,ydum):
        '''Returns list of zeros of length number of parameters'''
        return list( zeros( self.Nparam,'d' ) )
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):  # For now, put new point just on other side of layer
        # Do nothing for now; moving inhomogeneity action to aquifer class
        changed = 0; stop = 0; xyznew = 0.0
##        [x1,y1,z1] = xyz1; [x2,y2,z2] = xyz2
##        z1 = complex(x1,y1); z2 = complex(x2,y2)
##        Z1 = (2.0*z1 - (self.z1+self.z2))/(self.z2-self.z1)
##        Z2 = (2.0*z2 - (self.z1+self.z2))/(self.z2-self.z1)
##        # If point 1 on one side and point 2 on other side, find out if segment intesects line-sink
##        if (Z1.imag > 0 and Z2.imag < 0) or (Z1.imag < 0 and Z2.imag > 0):
##            Xintersec = Z1.real + (0.0-Z1.imag) / (Z2.imag-Z1.imag) * (Z2.real-Z1.real)
##            if abs(Xintersec) <= 1:  # There is an intersection point
##                changed = 1
##                theta = self.thetaNormOut
##                if Z1.imag < 0:  # old point on outside, Znew1 just out, Znew2 just in
##                    Znew1 = complex(Xintersec,-1e-6); Znew2 = complex(Xintersec,1e-6)
##                    if idir > 0: theta = self.thetaNormOut + pi  # Such that theta points in supposed direction of flow
##                else:
##                    Znew1 = complex(Xintersec,1e-6); Znew2 = complex(Xintersec,-1e-6)
##                    if idir < 0: theta = self.thetaNormOut + pi
##                znew1 = ( (self.z2-self.z1)*Znew1 + (self.z1+self.z2) ) / 2.0  # New complex coordinate
##                znew2 = ( (self.z2-self.z1)*Znew2 + (self.z1+self.z2) ) / 2.0
##                xnew1 = znew1.real; ynew1 = znew1.imag
##                xnew2 = znew2.real; ynew2 = znew2.imag
##                horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
##                horstepnew = sqrt( (xyz1[0]-xnew1)**2 + (xyz1[1]-ynew1)**2 )
##                zvertnew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )  # Just to be confusing, this is z (vertical coordinate) 
##                # I think this is the new value
##                xyznew = array([xnew2,ynew2,zvertnew])
####                disnorm = self.modelParent.dischargeNormVector(xnew1,ynew1,theta)[pyLayer]  # Normal discharge at xnew1,ynew1
####                if disnorm < 0:
####                    # Normal discharge at boundary points the other way; can only be if step size too large
####                    # Half step to boundary
####                    xyznew = 0.5 * ( xyz1 + array([xnew1,ynew1,zvertnew]) )
####                else:
####                    dis = self.modelParent.dischargeNormBelowZ(xnew1,ynew1,zvertnew,theta)
####                    [xnew,ynew,znew] = self.modelParent.zForGivenNormalDischarge\
####                                    (zvertnew,xnew1,ynew1,xnew2,ynew2,theta,dis)
####                    xyznew = array([xnew,ynew,znew])
        return [changed, stop, xyznew]
    def layout(self):
        return [0]

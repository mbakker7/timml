'''
mlholinesink.py contains the LineSinkHoLap class.
It is a higher order Laplace line-sink. Only one layer can be specified
(where the boundary condition is applied), and it does not generate leakage.
This file is part of the TimML library and is distributed under
the GNU LPGL. See the TimML.py file for more details.
(c) Mark Bakker, 2002-2007
'''

from numpy import *
from mlelement import *
from besselaes import *

class LineSinkHoLap(Element):
    def __init__(self,modelParent,x1,y1,x2,y2,sigma=[1.0],order=0,layers=[1],label=None,aquiferParentFixed=None,addToModel=1):
	Element.__init__(self,modelParent)
        self.x1 = float(x1); self.y1 = float(y1); self.x2 = float(x2); self.y2 = float(y2);
        self.sigma = array( sigma, 'd' )
        self.order = order
        self.layers = array(layers)
        assert len(self.layers) == 1, 'TimML input error, LineSinkHoLap must be screened in one layer only'
        assert len(self.sigma) == self.order + 1, 'TimML input error, sigma of LineSinkHoLap must be length order+1'
        self.aquiferParentFixed = aquiferParentFixed
        self.label = label
        self.type = 'LineSinkHoLap'
        self.setCoefs()
        if addToModel:
            self.modelParent.addElement(self)
    def __repr__(self):
	return 'LineSink z1,z2,sigma,order,layers: ' + str((self.x1,self.y1,self.x2,self.y2,self.sigma,self.order,list(self.layers)))
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
        if self.aquiferParentFixed == None:
            self.aquiferParent = self.modelParent.aq.findAquiferData(self.xc,self.yc)  # Determined at xc,yc
        else:
            self.aquiferParent = self.aquiferParentFixed
        assert self.aquiferParent.type == self.aquiferParent.conf, 'TimML input error, LineSinkHoLap must be screened in confined aqufier system'
      	self.pylayers = self.layers-1; self.NscreenedLayers = len(self.layers)
        self.parameters = self.sigma[:,newaxis]
        self.orderpone = self.order + 1
        self.rangeorderpone = range(self.orderpone)
        self.Ncp = self.order + 1
        thetacp = arange(pi,0,-pi/self.Ncp) - 0.5 * pi/self.Ncp
        Zcp = zeros( self.Ncp, 'D' )
        Zcp.real = cos(thetacp)
        Zcp.imag = 1e-6  # control point just on positive site (this is handy later on)
        zcp = Zcp * (self.z2 - self.z1) / 2.0 + 0.5 * (self.z1 + self.z2)
        self.xcp = zcp.real; self.ycp = zcp.imag
        self.potInf = zeros((self.orderpone,self.aquiferParent.Naquifers),'d')
        self.potLap = zeros(self.orderpone,'d')
        self.disxInf = zeros((self.orderpone,self.aquiferParent.Naquifers),'d')
        self.disyInf = zeros((self.orderpone,self.aquiferParent.Naquifers),'d')
        self.disxLap = zeros(self.orderpone,'d'); self.disyLap = zeros(self.orderpone,'d')
    def potentialInfluence(self,aq,x,y):
        '''Returns array of (NscreenedLayers,aq.Naquifers)'''
        if self.aquiferParent.type == self.aquiferParent.conf:
            if self.aquiferParent == aq:  # Same confined aquifer
                potlaplsho(x,y,self.x1,self.y1,self.x2,self.y2,self.order,self.potLap)  # Call FORTRAN extension
                rv = self.potInf
                rv[:,0] = self.potLap[:]
            else:
                rv = zeros((self.orderpone,aq.Naquifers),'d')
                potlaplsho(x,y,self.x1,self.y1,self.x2,self.y2,self.order,self.potLap)  # Call FORTRAN extension
                rv[:,0] = self.potLap[:]
        elif self.aquiferParent.type == self.aquiferParent.semi:
            rv = zeros((self.orderpone,aq.Naquifers),'d')
        return rv
    def potentialCollection(self,potsum,potadd,elementList,aq,x,y):
        for el in elementList:
            if el.aquiferParent.type == el.aquiferParent.conf:
                # Always the same action, as long as it is confined !
                potlaplsho(x,y,el.x1,el.y1,el.x2,el.y2,el.order,el.potLap)  # Call FORTRAN extension
                for i in el.rangeorderpone:  # always faster if just 8 long (maximum order)
                    potsum[0] = potsum[0] + el.parameters[i,0] * el.potLap[i]
            # If semi-confined you don't do anything, as this is just Laplace
        return potsum
    def dischargeInfluence(self,aq,x,y):
        '''Returns two arrays of (order,aq.Naquifers)'''
        if self.aquiferParent.type == self.aquiferParent.conf:
            if self.aquiferParent == aq:  # Same confined aquifer
                dislaplsho(x,y,self.x1,self.y1,self.x2,self.y2,self.order,self.disxLap,self.disyLap)  # Call FORTRAN extension
                rvx = self.disxInf; rvy = self.disyInf
                rvx[:,0] = self.disxLap[:]; rvy[:,0] = self.disyLap[:]
            else:
                rvx = zeros((self.orderpone,aq.Naquifers),'d')
                rvy = zeros((self.orderpone,aq.Naquifers),'d')
                dislaplsho(x,y,self.x1,self.y1,self.x2,self.y2,self.order,self.disxLap,self.disyLap)  # Call FORTRAN extension
                rvx[:,0] = self.disxLap[:]; rvy[:,0] = self.disyLap[:]
        elif self.aquiferParent.type == self.aquiferParent.semi:
            rvx = zeros((self.orderpone,aq.Naquifers),'d')
            rvy = zeros((self.orderpone,aq.Naquifers),'d')
        return [rvx,rvy]
    def dischargeCollection(self,dissum,disadd,elementList,aq,x,y):
        for el in elementList:
            if el.aquiferParent.type == el.aquiferParent.conf:
                dislaplsho(x,y,el.x1,el.y1,el.x2,el.y2,el.order,el.disxLap,el.disyLap)  # Call FORTRAN extension
                for i in el.rangeorderpone:  # always faster if just 8 long (maximum order)
                    dissum[0,0] = dissum[0,0] + el.parameters[i,0] * el.disxLap[i]
                    dissum[1,0] = dissum[1,0] + el.parameters[i,0] * el.disyLap[i]
        return dissum
    def layout(self):
        rv = [ 2,[],[] ]
        rv[1] = [ self.x1,self.x2 ]
        rv[2] = [ self.y1,self.y2 ]
        return rv
    def check(self):
        print 'Linesink from '+str(self.xy1)+' to '+str(self.xy2)+' Layer: '+str(self.layers[0])+\
              ' Strength: '+str(self.parameters[0,0])+' has no unknown parameters'
        return None
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
            elif abs( sum(sum(self.parameters)) ) < 1e-10:  # Screened in multiple layers, but zero total discharge
                z1 = complex(x1,y1); z2 = complex(x2,y2)
                Z1 = (2.0*z1 - (self.z1+self.z2))/(self.z2-self.z1)
                Z2 = (2.0*z2 - (self.z1+self.z2))/(self.z2-self.z1)
                # If point 1 on one side and point 2 on other side, find out if segment intesects line-sink
                if (Z1.imag > 0 and Z2.imag < 0) or (Z1.imag < 0 and Z2.imag > 0):
                    Xintersec = Z1.real + (0.0-Z1.imag) / (Z2.imag-Z1.imag) * (Z2.real-Z1.real)
                    if abs(Xintersec) <= 1:  # There is an intersection point
                        changed = 1
                        theta = self.thetaNormOut
                        if Z1.imag < 0:  # old point on outside, Znew1 just out, Znew2 just in
                            Znew1 = complex(Xintersec,-1e-6); Znew2 = complex(Xintersec,1e-6)
                            if idir > 0: theta = self.thetaNormOut + pi  # Such that theta points in supposed direction of flow
                        else:
                            Znew1 = complex(Xintersec,1e-6); Znew2 = complex(Xintersec,-1e-6)
                            if idir < 0: theta = self.thetaNormOut + pi
                        znew1 = ( (self.z2-self.z1)*Znew1 + (self.z1+self.z2) ) / 2.0  # New complex coordinate
                        znew2 = ( (self.z2-self.z1)*Znew2 + (self.z1+self.z2) ) / 2.0
                        xnew1 = znew1.real; ynew1 = znew1.imag
                        xnew2 = znew2.real; ynew2 = znew2.imag
                        horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
                        horstepnew = sqrt( (xyz1[0]-xnew1)**2 + (xyz1[1]-ynew1)**2 )
                        zvertnew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )  # Just to be confusing, this is z (vertical coordinate) 
                        disnorm = self.modelParent.dischargeNormVector(xnew1,ynew1,theta)[pyLayer]  # Normal discharge at xnew1,ynew1
                        if disnorm < 0:
                            # Normal discharge at boundary points the other way; can only be if step size too large
                            # Half step to boundary
                            xyznew = 0.5 * ( xyz1 + array([xnew1,ynew1,zvertnew]) )
                        else:
                            dis = self.modelParent.dischargeNormBelowZ(xnew1,ynew1,zvertnew,theta)
                            [xnew,ynew,znew] = self.modelParent.zForGivenNormalDischarge\
                                            (zvertnew,xnew1,ynew1,xnew2,ynew2,theta,dis)
                            xyznew = array([xnew,ynew,znew])
            else:
                stop = 1  # Seems to be only logical step. Stop traceline, there is nothing else we know
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
    def distanceSquaredToElementOld(self,x,y):
        Z = ( 2.0 * complex(x,y) - (self.z1 + self.z2) ) / (self.z2 - self.z1)
        if abs(Z.real) <= 1.0:
            dissq = ( Z.imag * self.L / 2 )**2
        elif Z.real < -1.0:
            dissq = ( x - self.x1 )**2 + ( y - self.y1 ) **2
        else:
            dissq = ( x - self.x2 )**2 + ( y - self.y2 ) **2
        return dissq

class CompHeadBoundaryElement(LineSinkHoLap):
    '''Class for comprehensive head specified boundary element'''
    def __init__(self,modelParent,x1,y1,x2,y2,head,order,label=None):
        LineSinkHoLap.__init__(self,modelParent,x1,y1,x2,y2,zeros(order+1),order,layers=[1],addToModel=0)
        self.head = head
        self.label = label
        self.modelParent.addElement(self)
    def __repr__(self):
	return 'CompHeadBoundaryElement z1,z2,head,order: ' + \
               str((self.x1,self.y1,self.x2,self.y2,self.head,self.order))
    def getMatrixRows(self,elementList):
        rows=[]
        pylayer = self.pylayers[0]
        pot = self.aquiferParent.headToPotential(pylayer,self.head)
        for i in range(self.Ncp):
            row = []
            for e in elementList:
                rowpart = e.getMatrixCoefficients(self.aquiferParent,pylayer,self.xcp[i],self.ycp[i],\
                                                lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                row = row + rowpart.tolist()
            row = row + [ pot - self.modelParent.potentialInLayer(self.aquiferParent,pylayer,self.xcp[i],self.ycp[i]) ]
            rows = rows + [ row ]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        self.parameters[:,0] = self.parameters[:,0] + xsol[icount:icount+self.Ncp]
        icount = icount+self.Ncp
        return icount
    def check(self):
        print 'CompHeadBoundaryElement from '+str(self.xy1)+' to '+str(self.xy2)
        for i in range(self.Ncp):
            print 'Control point '+str((self.xcp[i],self.ycp[i]))+' Specified head: '+str(self.head)+\
                  ' Computed head: '+str(self.modelParent.head(self.layers[0],self.xcp[i],self.ycp[i]))+\
                  ' Strength: '+str(self.parameters[i,0])
        return None

class CompResBoundaryElement(LineSinkHoLap):
    '''Class for comprehensive resistance and head specified boundary element'''
    def __init__(self,modelParent,x1,y1,x2,y2,head,res,width,order,label=None):
        LineSinkHoLap.__init__(self,modelParent,x1,y1,x2,y2,zeros(order+1),order,layers=[1],addToModel=0)
        self.head = float(head); self.res = float(res); self.width = float(width)
        self.label = label
        self.modelParent.addElement(self)
    def __repr__(self):
	return 'CompResBoundaryElement z1,z2,head,res,width,order: ' + \
               str((self.x1,self.y1,self.x2,self.y2,self.head,self.res,self.width,self.order))
    def getMatrixRows(self,elementList):
        rows=[]
        pylayer = self.pylayers[0]
        pot = self.aquiferParent.headToPotential(pylayer,self.head)
        for i in range(self.Ncp):
            row = []
            for e in elementList:
                rowpot = e.getMatrixCoefficients(self.aquiferParent,pylayer,self.xcp[i],self.ycp[i],\
                                                lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                rowqxqy = e.getMatrixCoefficients(self.aquiferParent,pylayer,self.xcp[i],self.ycp[i],\
                                lambda el,aq,pylayer,x,y:el.dischargeInfluenceInLayer(aq,pylayer,x,y))
                if size(rowqxqy) > 0:
                    rowdis = rowqxqy[0] * cos( self.thetaNormIn ) + rowqxqy[1] * sin( self.thetaNormIn )
                rowpart = rowpot + rowdis * self.aquiferParent.Tcomp * self.res / self.width
                row = row + rowpart.tolist()
            row = row + [ pot - \
                          self.modelParent.potentialInLayer(self.aquiferParent,pylayer,self.xcp[i],self.ycp[i]) - \
                          self.modelParent.dischargeNormInLayer(self.aquiferParent,pylayer,self.xcp[i],self.ycp[i],self.thetaNormIn) *\
                          self.aquiferParent.Tcomp * self.res / self.width ]
            rows = rows + [ row ]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        self.parameters[:,0] = self.parameters[:,0] + xsol[icount:icount+self.Ncp]
        icount = icount+self.Ncp
        return icount
    def check(self):
        print 'CompResBoundaryElement from '+str(self.xy1)+' to '+str(self.xy2)
        for i in range(self.Ncp):
            dis = self.modelParent.dischargeNormVector(self.xcp[i],self.ycp[i],self.thetaNormIn)[self.pylayers[0]] *\
                self.aquiferParent.Tcomp / self.aquiferParent.T[self.pylayers[0]]
            print 'Parameter '+str(i)+' Strength from head difference: '+\
                  str( ( self.modelParent.head(self.layers[0],self.xcp[i],self.ycp[i]) - self.head) * self.width / self.res ) +\
                  ' Discharge to element: '+str(-dis)
        return None

class CompDischargeBoundaryElement(LineSinkHoLap):
    '''Class for comprehensive discharge specified (on left side) boundary element'''
    def __init__(self,modelParent,x1,y1,x2,y2,discomp,order,label=None):
        LineSinkHoLap.__init__(self,modelParent,x1,y1,x2,y2,zeros(order+1),order,layers=[1],addToModel=0)
        self.discomp = float(discomp)
        self.label = label
        self.modelParent.addElement(self)
    def __repr__(self):
	return 'CompDischargeBoundaryElement z1,z2,discomp,layer,order: ' + \
               str((self.x1,self.y1,self.x2,self.y2,self.discomp,self.order))
    def getMatrixRows(self,elementList):
        rows=[]
        pylayer = self.pylayers[0]
        discharge = self.discomp / self.aquiferParent.Tcomp * self.aquiferParent.T[pylayer]
        for i in range(self.Ncp):
            row = []
            for e in elementList:
                rowqxqy = e.getMatrixCoefficients(self.aquiferParent,pylayer,self.xcp[i],self.ycp[i],\
                                lambda el,aq,pylayer,x,y:el.dischargeInfluenceInLayer(aq,pylayer,x,y))
                if size(rowqxqy) > 0:
                    rowpart = rowqxqy[0] * cos( self.thetaNormIn ) + rowqxqy[1] * sin( self.thetaNormIn )
                    row = row + rowpart.tolist()
            row = row + [ discharge -
                self.modelParent.dischargeNormInLayer(self.aquiferParent,pylayer,self.xcp[i],self.ycp[i],self.thetaNormIn) ]
            rows = rows + [ row ]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        self.parameters[:,0] = self.parameters[:,0] + xsol[icount:icount+self.Ncp]
        icount = icount+self.Ncp
        return icount
    def check(self):
        print 'CompDischargeBoundaryElement from '+str(self.xy1)+' to '+str(self.xy2)
        for i in range(self.Ncp):
            dis = self.modelParent.dischargeNormVector(self.xcp[i],self.ycp[i],self.thetaNormIn)[self.pylayers[0]] *\
                self.aquiferParent.Tcomp / self.aquiferParent.T[self.pylayers[0]]
            print 'Control point '+str((self.xcp[i],self.ycp[i]))+' Specified discomp: '+str(self.discomp)+\
                  ' Computed discharge: '+ \
                  str(dis)+\
                  ' Strength: '+str(self.parameters[i,0])
        return None

class LineSinkHoSemi(Element):
    '''Prototype element for higher order line-sink in single semi-confined aquifer'''
    def __init__(self,modelParent,x1,y1,x2,y2,sigma=[1.0],order=0,layers=[1],label=None,aquiferParentFixed=None,addToModel=1):
	Element.__init__(self,modelParent)
        self.x1 = float(x1); self.y1 = float(y1); self.x2 = float(x2); self.y2 = float(y2);
        self.sigma = array( sigma, 'd' )
        self.order = order
        self.layers = array(layers)
        assert len(self.layers) == 1 and self.layers[0] == 1, 'TimML input error, LineSinkHoSemi must be in single semi-confined aquifer'
        assert len(self.sigma) == self.order + 1, 'TimML input error, sigma of LineSinkHoSemi must be length order+1'
        self.aquiferParentFixed = aquiferParentFixed
        self.label = label
        self.type = 'LineSinkHoBes'
        self.setCoefs()
        if addToModel:
            self.modelParent.addElement(self)
    def __repr__(self):
	return 'LineSink z1,z2,sigma,order,layers: ' + str((self.x1,self.y1,self.x2,self.y2,self.sigma,self.order,list(self.layers)))
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
        if self.aquiferParentFixed == None:
            self.aquiferParent = self.modelParent.aq.findAquiferData(self.xc,self.yc)  # Determined at xc,yc
        else:
            self.aquiferParent = self.aquiferParentFixed
        assert self.aquiferParent.type == self.aquiferParent.semi, 'TimML input error, LineSinkHoSemi must be screened in single semi-confined aqufier system'
        self.Naq = 2  # This is weird, but needed to call FORTRAN extension
        self.labpzero = self.aquiferParent.zeropluslab
      	self.pylayers = self.layers-1; self.NscreenedLayers = len(self.layers)
        self.parameters = self.sigma[:,newaxis]
        self.orderpone = self.order + 1
        self.rangeorderpone = range(self.orderpone)
        self.Ncp = self.order + 1
        thetacp = arange(pi,0,-pi/self.Ncp) - 0.5 * pi/self.Ncp
        Zcp = zeros( self.Ncp, 'D' )
        Zcp.real = cos(thetacp)
        Zcp.imag = 1e-6  # control point just on positive site (this is handy later on)
        zcp = Zcp * (self.z2 - self.z1) / 2.0 + 0.5 * (self.z1 + self.z2)
        self.xcp = zcp.real; self.ycp = zcp.imag
        self.potInf = zeros((self.orderpone,self.aquiferParent.Naquifers),'d')
        self.potBes = zeros(2,'d')
        self.disxInf = zeros((self.orderpone,self.aquiferParent.Naquifers),'d')
        self.disyInf = zeros((self.orderpone,self.aquiferParent.Naquifers),'d')
        self.disxBes = zeros(2,'d'); self.disyBes = zeros(2,'d')
    def potentialInfluence(self,aq,x,y):
        '''Returns array of (NscreenedLayers,aq.Naquifers)'''
        if self.aquiferParent == aq:
            for ior in range(self.orderpone):
                potbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,self.Naq,self.labpzero,ior,self.potBes)
                self.potInf[ior,0] = self.potBes[1]
        else:
            self.potInf[:] = 0.0
        return self.potInf
    def dischargeInfluence(self,aq,x,y):
        '''Returns two arrays of (order,aq.Naquifers)'''
        if self.aquiferParent == aq:
            for ior in range(self.orderpone):
                disbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,self.Naq,self.labpzero,ior,self.disxBes,self.disyBes)
                self.disxInf[ior,0] = self.disxBes[1]
                self.disyInf[ior,0] = self.disyBes[1]
        else:
            self.disxInf[:] = 0.0; self.disyInf[:] = 0.0
        return [self.disxInf, self.disyInf]
    def layout(self):
        rv = [ 2,[],[] ]
        rv[1] = [ self.x1,self.x2 ]
        rv[2] = [ self.y1,self.y2 ]
        return rv
    def check(self):
        print 'Linesink from '+str(self.xy1)+' to '+str(self.xy2)+' Layer: '+str(self.layers[0])+\
              ' Strength: '+str(self.parameters[0,0])+' has no unknown parameters'
        return None
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
            elif abs( sum(sum(self.parameters)) ) < 1e-10:  # Screened in multiple layers, but zero total discharge
                z1 = complex(x1,y1); z2 = complex(x2,y2)
                Z1 = (2.0*z1 - (self.z1+self.z2))/(self.z2-self.z1)
                Z2 = (2.0*z2 - (self.z1+self.z2))/(self.z2-self.z1)
                # If point 1 on one side and point 2 on other side, find out if segment intesects line-sink
                if (Z1.imag > 0 and Z2.imag < 0) or (Z1.imag < 0 and Z2.imag > 0):
                    Xintersec = Z1.real + (0.0-Z1.imag) / (Z2.imag-Z1.imag) * (Z2.real-Z1.real)
                    if abs(Xintersec) <= 1:  # There is an intersection point
                        changed = 1
                        theta = self.thetaNormOut
                        if Z1.imag < 0:  # old point on outside, Znew1 just out, Znew2 just in
                            Znew1 = complex(Xintersec,-1e-6); Znew2 = complex(Xintersec,1e-6)
                            if idir > 0: theta = self.thetaNormOut + pi  # Such that theta points in supposed direction of flow
                        else:
                            Znew1 = complex(Xintersec,1e-6); Znew2 = complex(Xintersec,-1e-6)
                            if idir < 0: theta = self.thetaNormOut + pi
                        znew1 = ( (self.z2-self.z1)*Znew1 + (self.z1+self.z2) ) / 2.0  # New complex coordinate
                        znew2 = ( (self.z2-self.z1)*Znew2 + (self.z1+self.z2) ) / 2.0
                        xnew1 = znew1.real; ynew1 = znew1.imag
                        xnew2 = znew2.real; ynew2 = znew2.imag
                        horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
                        horstepnew = sqrt( (xyz1[0]-xnew1)**2 + (xyz1[1]-ynew1)**2 )
                        zvertnew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )  # Just to be confusing, this is z (vertical coordinate) 
                        disnorm = self.modelParent.dischargeNormVector(xnew1,ynew1,theta)[pyLayer]  # Normal discharge at xnew1,ynew1
                        if disnorm < 0:
                            # Normal discharge at boundary points the other way; can only be if step size too large
                            # Half step to boundary
                            xyznew = 0.5 * ( xyz1 + array([xnew1,ynew1,zvertnew]) )
                        else:
                            dis = self.modelParent.dischargeNormBelowZ(xnew1,ynew1,zvertnew,theta)
                            [xnew,ynew,znew] = self.modelParent.zForGivenNormalDischarge\
                                            (zvertnew,xnew1,ynew1,xnew2,ynew2,theta,dis)
                            xyznew = array([xnew,ynew,znew])
            else:
                stop = 1  # Seems to be only logical step. Stop traceline, there is nothing else we know
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
    def distanceSquaredToElementOld(self,x,y):
        Z = ( 2.0 * complex(x,y) - (self.z1 + self.z2) ) / (self.z2 - self.z1)
        if abs(Z.real) <= 1.0:
            dissq = ( Z.imag * self.L / 2 )**2
        elif Z.real < -1.0:
            dissq = ( x - self.x1 )**2 + ( y - self.y1 ) **2
        else:
            dissq = ( x - self.x2 )**2 + ( y - self.y2 ) **2
        return dissq

class LeftDischargeSemiElement(LineSinkHoSemi):
    '''Class for comprehensive discharge specified (on left side) boundary element'''
    def __init__(self,modelParent,x1,y1,x2,y2,disleft,order,label=None):
        LineSinkHoSemi.__init__(self,modelParent,x1,y1,x2,y2,zeros(order+1),order,layers=[1],addToModel=0)
        self.disleft = float(disleft)
        self.label = label
        self.modelParent.addElement(self)
    def __repr__(self):
	return 'LeftDischargeSemiElement z1,z2,disleft,layer,order: ' + \
               str((self.x1,self.y1,self.x2,self.y2,self.disleft,self.order))
    def getMatrixRows(self,elementList):
        rows=[]
        pylayer = self.pylayers[0]
        for i in range(self.Ncp):
            row = []
            for e in elementList:
                rowqxqy = e.getMatrixCoefficients(self.aquiferParent,pylayer,self.xcp[i],self.ycp[i],\
                                lambda el,aq,pylayer,x,y:el.dischargeInfluenceInLayer(aq,pylayer,x,y))
                if size(rowqxqy) > 0:
                    rowpart = rowqxqy[0] * cos( self.thetaNormIn ) + rowqxqy[1] * sin( self.thetaNormIn )
                    row = row + rowpart.tolist()
            row = row + [ self.disleft -
                self.modelParent.dischargeNormInLayer(self.aquiferParent,pylayer,self.xcp[i],self.ycp[i],self.thetaNormIn) ]
            rows = rows + [ row ]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        self.parameters[:,0] = self.parameters[:,0] + xsol[icount:icount+self.Ncp]
        icount = icount+self.Ncp
        return icount
    def check(self):
        print 'LeftDischargeSemiElement from '+str(self.xy1)+' to '+str(self.xy2)
        for i in range(self.Ncp):
            dis = self.modelParent.dischargeNormVector(self.xcp[i],self.ycp[i],self.thetaNormIn)[self.pylayers[0]]
            print 'Control point '+str((self.xcp[i],self.ycp[i]))+' Specified disleft: '+str(self.disleft)+\
                  ' Computed disleft: '+ str(dis) + ' Strength: '+str(self.parameters[i,0])
        return None

from numpy import *
import numpy.linalg as linalg
## import popen2  # Needed if pipeline code is used
from mlelement import *
from besselaes import *

class LineSink(Element):
    def __init__(self,modelParent,x1,y1,x2,y2,sigma,layers=0,aquiferParentFixed=None,addToModel=1,label=None):
	Element.__init__(self,modelParent)
        self.x1 = float(x1); self.y1 = float(y1); self.x2 = float(x2); self.y2 = float(y2);
        self.sigma = float(sigma)
        self.layers = atleast_1d(layers)
        self.aquiferParentFixed = aquiferParentFixed
        self.label = label
        self.type = 'linesink'
        self.setCoefs()
	if addToModel:
            self.modelParent.addElement(self)
    def __repr__(self):
	return 'LineSink z1,z2,sigma,layers: ' + str((self.x1,self.y1,self.x2,self.y2,self.sigma,list(self.layers)))
    def setCoefs(self):
        self.xy1 = (self.x1,self.y1); self.xy2 = (self.x2,self.y2)
        self.xmin = min(self.x1,self.x2); self.xmax = max(self.x1,self.x2)
        self.ymin = min(self.y1,self.y2); self.ymax = max(self.y1,self.y2)
        self.thetaNormOut = arctan2(self.y2-self.y1,self.x2-self.x1) - pi/2.0
        self.alpha = arctan2(self.y2-self.y1,self.x2-self.x1)
        self.cosalpha = cos(self.alpha); self.sinalpha = sin(self.alpha)
        self.z1 = complex(self.x1,self.y1); self.z2 = complex(self.x2,self.y2);
        self.xc = 0.5*(self.x1+self.x2); self.yc = 0.5*(self.y1+self.y2)
        self.L = abs(self.z2-self.z1); self.Lover4pi = self.L/(4.0*pi); self.Lover2 = self.L / 2.0
        if self.aquiferParentFixed is None:
            self.aquiferParent = self.modelParent.aq.findAquiferData(self.xc,self.yc)  # Determined at xc,yc
        else:
            self.aquiferParent = self.aquiferParentFixed
      	self.pylayers = self.layers  # fixed to base zero
	self.NscreenedLayers = len(self.layers)
	if self.NscreenedLayers == 1:  # Screened in only one layer, no unknowns
            self.parameters = array([[self.sigma]])
        else:
            self.parameters = zeros((self.NscreenedLayers,1),'d')
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
                for i in range(self.NscreenedLayers):
                    pylayer = self.pylayers[i]
                    ramat = self.aquiferParent.eigvec
                    rb = zeros(self.aquiferParent.Naquifers,'d')
                    rb[pylayer] = 1.0
                    self.coef[i,:] = linalg.solve(ramat,rb)
        self.paramxcoef = sum( self.parameters * self.coef, 0 )  # Parameters times coefficients
        self.potInf = zeros(self.aquiferParent.Naquifers,'d')
        self.potLap = zeros(1,'d')
        self.disxInf = zeros(1,'d'); self.disyInf = zeros(1,'d')
        self.connected = True # Used to keep track of changes between connected and disconnected
    def potentialInfluence(self,aq,x,y):
        '''Returns array of (NscreenedLayers,aq.Naquifers)'''
        # Order of if-s is slightly different than for wells, because Laplace part is not computed separately
        if self.aquiferParent.type == self.aquiferParent.conf:
            if self.aquiferParent == aq:  # Same confined aquifer
                # lab = array( [0.]+list(self.aquiferParent.lab),'d' );
                lab = self.aquiferParent.zeropluslab
                potbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers,lab,0,self.potInf)  # Call FORTRAN extension
                rv = self.coef * self.potInf
            elif aq.type == aq.conf:  # Different confined aquifer, Laplace part only
                rv = zeros((self.NscreenedLayers,aq.Naquifers),'d')
                lab = array([0.0]); potInf = array([0.0])
                potbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,1,lab,0,potInf)  # Call FORTRAN extension
                rv[:,0] = potInf[0]
        elif self.aquiferParent.type == self.aquiferParent.semi:
            rv = zeros((self.NscreenedLayers,aq.Naquifers),'d')
            if self.aquiferParent == aq:  # Same semi-confined aquifer
                lab = self.aquiferParent.zeropluslab  # Still add zero, because that is what FORTRAN extension expects
                potInf = zeros(self.aquiferParent.Naquifers+1,'d')  # Should be changed if I change FORTRAN extension
                potbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers+1,lab,0,potInf)  # Call FORTRAN extension
                rv =  self.coef * potInf[1:]
                # Else we return all zeros
        return rv
    def dischargeInfluence(self,aq,x,y):
        '''Returns two arrays of (NscreenedLayers,aq.Naquifers)'''
        rvx = zeros((self.NscreenedLayers,aq.Naquifers),'d')
        rvy = zeros((self.NscreenedLayers,aq.Naquifers),'d')
        if self.aquiferParent.type == self.aquiferParent.conf:
            if self.aquiferParent == aq:  # Same confined aquifer
                disxInf = zeros(self.aquiferParent.Naquifers,'d')
                disyInf = zeros(self.aquiferParent.Naquifers,'d')
                lab = self.aquiferParent.zeropluslab
                disbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers,lab,0,disxInf,disyInf)  # Call FORTRAN extension
                rvx = self.coef * disxInf; rvy = self.coef * disyInf
            elif aq.type == aq.conf:  # Different confined aquifer, Laplace part only
                lab = zeros(1,'d')
                dislaplsho(x,y,self.x1,self.y1,self.x2,self.y2,0,self.disxInf,self.disyInf)  # Call FORTRAN extension
                rvx[:,0] = self.coef[:,0] * self.disxInf[0]; rvy[:,0] = self.coef[:,0] * self.disyInf[0]
        elif self.aquiferParent.type == self.aquiferParent.semi:
            if self.aquiferParent == aq:  # Same semi-confined aquifer
                lab = self.aquiferParent.zeropluslab  # Still add zero, because that is what FORTRAN extension expects
                disxInf = zeros(self.aquiferParent.Naquifers+1,'d')  # Should be changed if I change FORTRAN extension
                disyInf = zeros(self.aquiferParent.Naquifers+1,'d')
                disbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers+1,lab,0,disxInf,disyInf)  # Call FORTRAN extension
                rvx =  self.coef * disxInf[1:]; rvy =  self.coef * disyInf[1:] 
        return [rvx,rvy]

    def potentialCollection(self,potsum,potadd,elementList,aq,x,y):

        # Not tested for semi-confined aquifer (improvement can and should be made there; no zeros, for example)
        
        for el in elementList:
            if el.aquiferParent.type == el.aquiferParent.conf:
                if el.aquiferParent == aq:  # Same confined aquifer
                    lab = el.aquiferParent.zeropluslab
                    potbeslsho(x,y,el.x1,el.y1,el.x2,el.y2,el.aquiferParent.Naquifers,lab,0,potadd)  # Call FORTRAN extension
                    potadd = el.paramxcoef * potadd
                elif aq.type == aq.conf:  # Different confined aquifer, Laplace part only
                    potlaplsho(x,y,el.x1,el.y1,el.x2,el.y2,0,self.potLap)  # Call FORTRAN extension
                    potadd[0] = el.paramxcoef[0] * self.potLap[0]
                    potadd[1:] = 0.0
            elif el.aquiferParent.type == el.aquiferParent.semi:
                if el.aquiferParent == aq:  # Same semi-confined aquifer
                    lab = el.aquiferParent.zeropluslab  # Still add zero, because that is what FORTRAN extension expects
                    potInf = zeros(el.aquiferParent.Naquifers+1,'d')  # Should be changed if I change FORTRAN extension
                    potbeslsho(x,y,el.x1,el.y1,el.x2,el.y2,el.aquiferParent.Naquifers+1,lab,0,potInf)  # Call FORTRAN extension
                    potadd =  el.paramxcoef * potInf[1:]
                else:
                    potadd[:] = 0.0
            potsum = potsum + potadd  
        return potsum

    def dischargeCollection(self,dissum,disadd,elementList,aq,x,y):

        # Needs testing in other aquifers

        for el in elementList:
            if el.aquiferParent.type == el.aquiferParent.conf:
                if el.aquiferParent == aq:  # Same confined aquifer
                    lab = el.aquiferParent.zeropluslab
                    disbeslsho(x,y,el.x1,el.y1,el.x2,el.y2,el.aquiferParent.Naquifers,lab,0,disadd[0,:],disadd[1,:])  # Call FORTRAN extension
                    disadd = el.paramxcoef * disadd
                elif aq.type == aq.conf:  # Different confined aquifer, Laplace part only
                    dislaplsho(x,y,el.x1,el.y1,el.x2,el.y2,0,self.disxInf,self.disyInf)  # Call FORTRAN extension
                    disadd[0,0] = el.paramxcoef[0] * self.disxInf[0]
                    disadd[1,0] = el.paramxcoef[0] * self.disyInf[0]
                    disadd[:,1:] = 0.0
            elif el.aquiferParent.type == el.aquiferParent.semi:
                if el.aquiferParent == aq:  # Same semi-confined aquifer
                    lab = el.aquiferParent.zeropluslab  # Still add zero, because that is what FORTRAN extension expects
                    disInf = zeros((2,el.aquiferParent.Naquifers+1),'d')  # Should be changed if I change FORTRAN extension
                    disbeslsho(x,y,el.x1,el.y1,el.x2,el.y2,el.aquiferParent.Naquifers+1,lab,0,disInf[0,:],disInf[1,:])  # Call FORTRAN extension
                    disadd =  el.paramxcoef * disInf[:,1:]
            dissum = dissum + disadd
        return dissum
    def totalDischargeInfluence(self,aq,pylayer,x,y):
        rv = ones(self.NscreenedLayers,'d') * self.L
        return rv   
    def layout(self):
        rv = [ 2,[],[] ]
        rv[1] = [ self.x1,self.x2 ]
        rv[2] = [ self.y1,self.y2 ]
        return rv
    def getMatrixRows(self,elementList):
        rows=[]
        if self.NscreenedLayers > 1:  # Otherwise no unknowns !
            rowlast = zeros((0),'d')
            lastpylayer = self.pylayers[self.NscreenedLayers-1]; Tlast = self.aquiferParent.T[lastpylayer]
            for e in elementList:
                rowpart = e.getMatrixCoefficients(self.aquiferParent,lastpylayer,self.xc,self.yc,\
                                lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                rowlast = hstack(( rowlast, rowpart ))
            potlast = self.modelParent.potentialInLayer(self.aquiferParent,lastpylayer,self.xc,self.yc)
            for i in range(self.NscreenedLayers - 1):                
                row = zeros(0,'d'); T = self.aquiferParent.T[self.pylayers[i]]
                for e in elementList:
                    rowpart = e.getMatrixCoefficients(self.aquiferParent,self.pylayers[i],self.xc,self.yc,\
                                    lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                    row = hstack(( row, rowpart ))
                row = Tlast * row - T * rowlast    # Multiply with T's and subtract last row
                row = hstack(( row,  T * potlast - Tlast * \
                    self.modelParent.potentialInLayer(self.aquiferParent,self.pylayers[i],self.xc,self.yc) ))  # Add right hand side as last element
                rows = rows + [row.tolist()]
            row = zeros(0,'d')
            # Last equation is sum of strength parameters = sigma
            for e in elementList:
                if e == self:
                    row = hstack(( row, ones(self.NscreenedLayers,'d') ))
                else:
                    row = hstack(( row, e.getMatrixCoefficients(self.aqdum,self.ldum,self.xdum,self.ydum,\
                                             lambda el,aqdum,ldum,xdum,ydum:el.zeroFunction(aqdum,ldum,xdum,ydum)) ))
            #row = hstack(( row, self.sigma - sum(self.parameters,0)[0] ))
            row = hstack(( row, self.sigma - sum(self.parameters[:,0]) ))
            rows = rows + [row.tolist()]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,xc,yc,func):
        rv = zeros(0,'d')
        if self.NscreenedLayers > 1:
            rv = func(self,aq,pylayer,xc,yc)
        return rv
    def takeParameters(self,xsol,icount):
        if self.NscreenedLayers > 1:
            for i in range(self.NscreenedLayers):
                self.parameters[i,0] = self.parameters[i,0] + xsol[icount]
                icount = icount+1
            self.paramxcoef = sum( self.parameters * self.coef, 0 )
        return icount
    def check(self):
        if self.NscreenedLayers == 1:
            print 'Linesink from '+str(self.xy1)+' to '+str(self.xy2)+' Layer: '+str(self.layers[0])+\
                  ' Strength: '+str(self.parameters[0,0])+' has no unknown parameters'
        else:
            print 'Linesink from '+str(self.xy1)+' to '+str(self.xy2)
            print 'Specified strength: '+str(self.sigma)+' Computed strength: '+str(sum(self.parameters,0)[0])
            for i in range(self.NscreenedLayers):
                print 'Layer '+str(self.layers[i])+' Head: '+\
                      str(self.modelParent.head(self.layers[i],self.xc,self.yc,self.aquiferParent))+\
                      ' Strength: '+str(self.parameters[i,0])
        return None
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):
        changed = 0; stop = 0; xyznew = 0.0
        if ( pyLayer == self.pylayers ).any():  # In layer line-sink is screened. This is true if self.pylayers contains pyLayer
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
                        theta = self.thetaNormOut
                        if Z1.imag < 0:  # old point on outside, Znew1 just below, Znew2 just above
                            Znew1 = complex(Xintersec,-1e-6); Znew2 = complex(Xintersec,1e-6)
                            if idir > 0: theta = self.thetaNormOut + pi  # Such that theta points in supposed direction of flow
                        else:
                            Znew1 = complex(Xintersec,1e-6); Znew2 = complex(Xintersec,-1e-6)
                            if idir < 0: theta = self.thetaNormOut + pi
                        znew1 = ( (self.z2 - self.z1)*Znew1 + (self.z1 + self.z2) ) / 2.0  # New complex coordinate
                        znew2 = ( (self.z2 - self.z1)*Znew2 + (self.z1 + self.z2) ) / 2.0
                        xnew1 = znew1.real; ynew1 = znew1.imag
                        xnew2 = znew2.real; ynew2 = znew2.imag
                        horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
                        horstepnew = sqrt( (xyz1[0]-xnew1)**2 + (xyz1[1]-ynew1)**2 )
                        zvertnew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )  # Vertical coordinate at new1
                        disnorm1 = self.modelParent.dischargeNormVector(xnew1,ynew1,theta)[pyLayer]  # Normal discharge at new1
                        disnorm2 = self.modelParent.dischargeNormVector(xnew2,ynew2,theta)[pyLayer]
                        if ( self.parameters[0,0] > 0 and idir > 0 ) or ( self.parameters[0,0] < 0 and idir < 0 ):
                            # Takes water out ,check if we should stop trace.
                            # Allows for zero to allow for line-sinks used elsewhere
                            if sign(disnorm1) == sign(disnorm2):  # Some water leaves on other side
                                flowbelow = (zvertnew - self.aquiferParent.zb[pyLayer]) / self.aquiferParent.H[pyLayer] * abs(disnorm1)
                                if abs(disnorm2) > flowbelow:
                                    zvertnew = self.aquiferParent.zb[pyLayer] + flowbelow / abs(disnorm2) * self.aquiferParent.H[pyLayer] 
                                    print 'xnew2,ynew2,zvertnew ',xnew2,ynew2,zvertnew
                                    xyznew = array([xnew2,ynew2,zvertnew]) 
                                else:  # Taken out !
                                    xyznew = array([xnew1,ynew1,self.aquiferParent.zt[pyLayer]])
                                    stop = 1
                            else:  # All water taken out
                                Znew = complex(Xintersec,0.0)
                                znew = ( (self.z2-self.z1)*Znew + (self.z1+self.z2) ) / 2.0
                                xnew = znew.real; ynew = znew.imag
                                horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
                                horstepnew = sqrt( (xyz1[0]-xnew)**2 + (xyz1[1]-ynew)**2 )
                                #xyznew = array([xnew,ynew,self.aquiferParent.zt[pyLayer]]) # Changed so that elevation at exit is correct
                                znew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2])
                                xyznew = array([xnew,ynew,znew])
                                stop = 1
                        else: # Line-sink puts water in (should jump streamline down)
                            flowbelow = (zvertnew - self.aquiferParent.zb[pyLayer]) / self.aquiferParent.H[pyLayer] * abs(disnorm1)
                            zvertnew = self.aquiferParent.zb[pyLayer] + flowbelow / abs(disnorm2) * self.aquiferParent.H[pyLayer] 
                            xyznew = array([xnew2,ynew2,zvertnew]) 
                elif abs(Z2.imag) < 2e-6 * step / self.L and abs(Z2.real) <= 1:  # Point 2 is on linesink
                    stop = 1
#            elif abs( sum(sum(self.parameters)) ) < 1e-10:  # Screened in multiple layers, but zero total discharge
            elif abs( sum(self.parameters) ) < 1e-10:  # Screened in multiple layers, but zero total discharge
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
                
    def PYpotentialInfluence(self,x,y,zvert=0,t=0):
        ac = zeros((9),'d'); bc = zeros((9),'d')
        ac[0] = -0.500004211065677e0;   bc[0] = 0.115956920789028e0
        ac[1] = -0.124989431448700e0;   bc[1] = 0.278919134951974e0
        ac[2] = -0.781685695640975e-2;  bc[2] = 0.252752008621167e-1
        ac[3] = -0.216324415010288e-3;  bc[3] = 0.841879407543506e-3
        ac[4] = -0.344525393452639e-5;  bc[4] = 0.152425102734818e-4
        ac[5] = -0.315133836828774e-7;  bc[5] = 0.148292488399579e-6
        ac[6] = -0.296636186265427e-9;  bc[6] = 0.157622547107156e-8
        ac[7] = -0.313689942474032e-12; bc[7] = 0.117975437124933e-11
        ac[8] = -0.112031912249579e-13; bc[8] = 0.655345107753534e-13
        NBes = 8; Rconv = 7; Rconvsq = Rconv * Rconv; tiny = 1e-12;
        zin = complex(x,y)
	z = ( 2.0*zin - (self.z1+self.z2) ) / (self.z2-self.z1)
	zplus1 = z + 1.0; zmin1 = z - 1.0
	if abs(zplus1) < self.tiny: zplus1 = zplus1 + self.tiny
	if abs(zmin1) < self.tiny: zmin1 = zmin1 + self.tiny 
	omega = self.Lover4pi * ( zplus1*cmath.log(zplus1) - zmin1*cmath.log(zmin1) )
        potInf = [omega.real]
        for i in range(self.aquiferParent.Naquifers-1):
            pot = 0.0
            NLS = math.ceil(self.L/self.aquiferParent.lab[i])   # Compute whether linesink must be broken up

            for j in range(NLS):
                z1 = self.z1 + float(j)/NLS * (self.z2-self.z1)
                z2 = z1 + (self.z2-self.z1)/NLS
                L = abs(z2-z1)

                N=self.NBes
                # Compute local z coordinate
                biglab = 2.0 * self.aquiferParent.lab[i] / L;
                z = ( 2.0*zin - (z1+z2) ) / ( z2 - z1 ) / biglab; zbar = conjugate(z);
                if abs(z) > self.Rconv:
                    pot = pot + 0.0
                    continue
                        
                # Coefficients gamma(n,m) of (d-zbar)^n=g_n,1+g_n,2*(d-z)+g_n,3*(d-z)^2+...+g_n,n+1*(d-z)^n
                # Store coefficents in matrix. Note that arange had to go to n+1, which means it stops at n
                binom = zeros((N+1,N+1),'d')
                zzbarp = zeros((N+1),'D')
                for n in range(0,N+1):
                    zzbarp[n]=(z-zbar)**n
                    for m in range(0,n+1):
                        binom[n][m]=product(arange(m+1,n+1,1,'d'))/product(arange(1,n-m+1,1,'d'))
                   
                gamma=zeros((N+1,N+1),'D');
                for n in range(0,N+1):
                    for m in range(0,n+1):
                        gamma[n][m]=binom[n][m]*zzbarp[n-m]

                # Coefficients of a-tilde and b-tilde
                alpha = zeros((2*N+1),'D'); beta = zeros((2*N+1),'D');
                for n in range(0,2*N+1):
                    alpha[n] = 0.0; beta[n] = 0.0;
                    for m in range(max(0,n-N),n/2+1):
                        alpha[n] = alpha[n] + self.ac[n-m] * gamma[n-m][m];
                        beta[n] = beta[n] + self.bc[n-m] * gamma[n-m][m];
                 
                # Evaluation of integral
                B = 0;
                d1minz = -1.0/biglab - z; d2minz = 1.0/biglab - z;
                if abs(d1minz) < self.tiny: d1minz = d1minz + self.tiny  # Check if at singularity
                if abs(d2minz) < self.tiny: d2minz = d2minz + self.tiny  
                ln1 = cmath.log(d1minz); ln2 = cmath.log(d2minz);
                for n in range(0,2*N+1):
                    B = B + (2.0*alpha[n]*ln2 - 2.0*alpha[n]/(n+1) + beta[n]) * d2minz**(n+1)/(n+1);
                    B = B - (2.0*alpha[n]*ln1 - 2.0*alpha[n]/(n+1) + beta[n]) * d1minz**(n+1)/(n+1);

                pot = pot - B.real * self.aquiferParent.lab[i] / (2*pi);

            potInf = potInf + [pot]

        rv = []
        for coef in self.coef:
            rv = rv + [list(potInf * coef)]  # Returns list of lists
        return rv

# Added for impulse response testing
    def IRpotentialInfluence(self,x,y,zvert=0,t=0):
        ac = zeros((2),'d'); bc = zeros((2),'d')
        #ac[0] = 0.0;   bc[0] = 0.0
        ac[1] = 1.0;   bc[1] = -2.0
        NBes = 1; Rconv = 7; Rconvsq = Rconv * Rconv; tiny = 1e-12;
        zin = complex(x,y)
	z = ( 2.0*zin - (self.z1+self.z2) ) / (self.z2-self.z1)
	zplus1 = z + 1.0; zmin1 = z - 1.0
	if abs(zplus1) < tiny: zplus1 = zplus1 + tiny
	if abs(zmin1) < tiny: zmin1 = zmin1 + tiny 

        z1 = self.z1 
        z2 = self.z2
        L = abs(z2-z1)

        N = NBes
        # Compute local z coordinate
        biglab = 1.0;
        z = ( 2.0*zin - (self.z1+self.z2) ) / (self.z2-self.z1); zbar = conjugate(z);
                
        # Coefficients gamma(n,m) of (d-zbar)^n=g_n,1+g_n,2*(d-z)+g_n,3*(d-z)^2+...+g_n,n+1*(d-z)^n
        # Store coefficents in matrix. Note that arange had to go to n+1, which means it stops at n
        binom = zeros((N+1,N+1),'d')
        zzbarp = zeros((N+1),'D')
        for n in range(0,N+1):
            zzbarp[n]=(z-zbar)**n
            for m in range(0,n+1):
                binom[n][m]=product(arange(m+1,n+1,1,'d'))/product(arange(1,n-m+1,1,'d'))
           
        gamma=zeros((N+1,N+1),'D');
        for n in range(0,N+1):
            for m in range(0,n+1):
                gamma[n][m]=binom[n][m]*zzbarp[n-m]

        # Coefficients of a-tilde and b-tilde
        alpha = zeros((2*N+1),'D'); beta = zeros((2*N+1),'D');
        for n in range(0,2*N+1):
            alpha[n] = 0.0; beta[n] = 0.0;
            for m in range(max(0,n-N),n/2+1):
                alpha[n] = alpha[n] + ac[n-m] * gamma[n-m][m];
                beta[n] = beta[n] + bc[n-m] * gamma[n-m][m];
         
        # Evaluation of integral
        B = complex(0,0);
        d1minz = -1.0/biglab - z; d2minz = 1.0/biglab - z;
        if abs(d1minz) < tiny: d1minz = d1minz + tiny  # Check if at singularity
        if abs(d2minz) < tiny: d2minz = d2minz + tiny  
        ln1 = log(d1minz); ln2 = log(d2minz);
        #alpha[1] = alpha[1] / 2.0
        #beta[1] = beta[1] / 2.0
        print 'alpha ',alpha
        print 'beta ',beta
        for n in range(0,2*N+1):
            B = B + (2.0*alpha[n]*ln2 - 2.0*alpha[n]/(n+1) + beta[n]) * d2minz**(n+1)/(n+1);
            B = B - (2.0*alpha[n]*ln1 - 2.0*alpha[n]/(n+1) + beta[n]) * d1minz**(n+1)/(n+1);

        print 'B.real ',B.real

        pot = - B.real / (2*pi);

        return pot

class HeadLineSink(LineSink):
    def __init__(self,modelParent,x1,y1,x2,y2,head,layers=0,aquiferParentFixed=None,label=None):
        LineSink.__init__(self,modelParent,x1,y1,x2,y2,0.0,layers,aquiferParentFixed)
        self.head = float(head)
        self.label = label
    def __repr__(self):
	return 'HeadLineSink z1,z2,head,strength,layer: ' + str((self.x1,self.y1,self.x2,self.y2,self.head,self.parameters,list(self.layers)))
    def getMatrixRows(self,elementList):
        rows = []
        for i in range(self.NscreenedLayers):
            pot = self.aquiferParent.headToPotential(self.pylayers[i],self.head)     
            row = []
            for e in elementList:
                rowpart = e.getMatrixCoefficients(self.aquiferParent,self.pylayers[i],self.xc,self.yc,\
                                                  lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                row = row + rowpart.tolist()
            row = row + [pot - self.modelParent.potentialInLayer(self.aquiferParent,self.pylayers[i],self.xc,self.yc)]
            rows = rows + [ row ]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        for i in range(self.NscreenedLayers):
                self.parameters[i,0] = self.parameters[i,0] + xsol[icount]
                icount = icount+1
        self.paramxcoef = sum( self.parameters * self.coef, 0 )
        return icount
    def check(self):
        print 'HeadLineSink from '+str(self.xy1)+' to '+str(self.xy2)
        for i in range(self.NscreenedLayers):
            print 'Layer '+str(self.layers[i])+' Specified head: '+str(self.head)+\
                  ' Computed head: '+str(self.modelParent.headVector(self.xc,self.yc)[self.pylayers[i]])+\
                  ' Strength: '+str(self.parameters[i,0])
        return None

class ResLineSink(LineSink):
    def __init__(self,modelParent,x1,y1,x2,y2,head,res,width,layers=0,label=None,bottomelev=None):  
        LineSink.__init__(self,modelParent,x1,y1,x2,y2,0.0,layers)
        self.head = float(head); self.res = float(res); self.width = float(width)
        self.label = label
        if bottomelev != None:
            self.bottomelev = float(bottomelev) # Only used for determining percolating line-sinks
        else:
            self.bottomelev = None
	#zcp = complex(self.xc,self.yc) + complex(0.0, 0.5 * self.width) * exp( complex( 0.0, self.alpha ) )
        # Put the control point back at the center of the line-sink, else strange things may happen when width is large
        zcp = complex(self.xc,self.yc)
	self.xcp, self.ycp = zcp.real, zcp.imag
    def __repr__(self):
	return 'ResLineSink z1,z2,head,res,width,strength,layer: ' + \
               str((self.x1,self.y1,self.x2,self.y2,self.head,self.res,self.width,self.parameters,list(self.layers)))
    def getMatrixRows(self,elementList):  
        rows = []
        for i in range(self.NscreenedLayers):  
            pot = self.aquiferParent.headToPotential(self.pylayers[i],self.head)     
            row = []
            for e in elementList:
                rowpart = e.getMatrixCoefficients(self.aquiferParent,self.pylayers[i],self.xcp,self.ycp,\
                                lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                if e == self:
                    rowpart[i] = rowpart[i] - self.aquiferParent.T[self.pylayers[i]] * self.res / self.width
                row = row + rowpart.tolist()
            row = row + [ pot - \
                          self.modelParent.potentialInLayer(self.aquiferParent,self.pylayers[i],self.xcp,self.ycp) +\
                          self.parameters[i,0] * self.aquiferParent.T[self.pylayers[i]] * self.res / self.width ]
            rows = rows + [ row ]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        for i in range(self.NscreenedLayers):
                self.parameters[i,0] = self.parameters[i,0] + xsol[icount]
                icount = icount+1
        self.paramxcoef = sum( self.parameters * self.coef, 0 )
        return icount
    def check(self):
        print 'ResLineSink from '+str(self.xy1)+' to '+str(self.xy2)
        for i in range(self.NscreenedLayers):
            if self.connected:
                print 'Layer '+str(self.layers[i])+' Strength from head difference: '+\
                  str( ( self.modelParent.head(self.layers[i],self.xcp,self.ycp) - self.head) * self.width / self.res ) +\
                  ' Strength: '+str(self.parameters[i,0])
            else:
                print 'Layer '+str(self.layers[i])+' Percolating. '+\
                    'Bottom elev., head in aquifer, strength: '+\
                    str( ( self.bottomelev, self.modelParent.head(self.layers[i],self.xcp,self.ycp), self.parameters[i,0] ) )
        return None
    def getMatrixRows_nonlinear(self,elementList):
        '''Only checked if screened in only one layer'''
        if self.NscreenedLayers > 1: return None
        rows = []
        for i in range(self.NscreenedLayers):
            head_in_aquifer = self.modelParent.head( self.pylayers[i]+1, self.xcp, self.ycp )
            if head_in_aquifer >= self.bottomelev:  # There is a hydraulic connection
                if self.connected:
                    return None
                else:
                    self.connected = True
                    rows = self.getMatrixRows( elementList )
            else:  # There is no hydraulic connection, set strength accordingly
                if self.connected:
                    self.connected = False
                    strength = -( self.head - self.bottomelev ) * self.width / self.res
                    row = []
                    for e in elementList:
                        rowpart = e.getMatrixCoefficients(self.aqdum,self.ldum,self.xdum,self.ydum,\
                                        lambda el,aq,pylayer,x,y:el.zeroFunction(aq,pylayer,x,y))
                        if e == self:
                            rowpart[i] = 1.0
                        row.extend(rowpart)
                    row = row + [ strength - self.parameters[i,0] ]
                    rows = rows + [ row ]
                else:
                    return None
        return rows

### Alternative implementations of the potentialInfluence function
### Using a pipe that is opened to an executable file
### This is rather fast and easier to implement than a true Python extension,
### but is also, pick your favorite, a rather ugly/slick solution
### May still be useful for platforms where a Python extension is not available

##    lsin,lsout = popen2.popen2('/TimML/bessells')
##    def potentialInfluencePipe(self,aq,x,y,zvert=0,t=0):
##        output = str(x)+' '+str(y)+' '+str(self.xy1)+' '+str(self.xy2)
##        if aq == self.aquiferParent:
##            self.lsout.write((str(self.aquiferParent.Naquifers)+'\n'))      
##            for la in self.aquiferParent.lab:
##                output=output+' '+str(la)
##            self.lsout.write((output+'\n'))
##            potInf = list( eval(self.lsin.readline()) )
##            rv = []
##            for coef in self.coef:
##                rv = rv + [list(potInf * coef)]  # Returns list of lists
##        else:
##            self.lsout.write('1\n')
##            self.lsout.write((output+'\n'))
##            potInf = eval(self.lsin.readline()) + list( zeros(aq.Naquifers-1,'d') )
##            rv = []
##            for i in range(self.NscreenedLayers):
##                rv = rv + [potInf]
##        return rv

class LineSinkDitch(Element):
    def __init__(self,modelParent,xylist,Q,res,width,layers=0,addToModel=1,label=None,aquiferParentFixed=None):
	Element.__init__(self,modelParent)
	if addToModel:
            self.modelParent.addElement(self)
        self.NLS = len(xylist) - 1
        self.xylist = xylist
        self.lsList = []
        for i in range(self.NLS):
            ls = LineSink(modelParent,xylist[i][0],xylist[i][1],xylist[i+1][0],xylist[i+1][1],0,layers,aquiferParentFixed,0)
            self.lsList = self.lsList + [ls]
        self.Q = float(Q)
        self.res = float(res)
        self.width = float(width)
        self.layers = atleast_1d(layers)
        assert len(self.layers)==1, "TimML Input error: LineSinkDitch screened in more than 1 layer"
        if aquiferParentFixed == None:
            self.aquiferParent = self.modelParent.aq.findAquiferData(xylist[0][0],xylist[0][1])  # Determined at first node
        else:
            self.aquiferParent = aquiferParentFixed
        self.label = label
        self.type = 'linesink'
        self.setCoefs()
    def __repr__(self):
	return 'LineSinkDitch Q,res,width,layers ' + str((self.Q,self.res,self.width,list(self.layers))) +\
               ' consisting of the following ' + str(self.NLS) + ' linesinks ' + str(self.lsList)
    def setCoefs(self):
        self.x1 = self.xylist[0][0];  self.y1 = self.xylist[0][1]
        self.x2 = self.xylist[-1][0]; self.y2 = self.xylist[-1][1]
        self.z1 = complex(self.x1,self.y1); self.z2 = complex(self.x2,self.y2)
        self.Ltot = sqrt( (self.x2-self.x1)**2 + (self.y2-self.y1)**2 )
        self.parameters = zeros( (self.NLS,1), 'd' )
        self.pylayers = self.layers  # fixed to zero base
        self.xc = zeros(self.NLS,'d'); self.yc = zeros(self.NLS,'d')
        self.xcp = zeros(self.NLS,'d'); self.ycp = zeros(self.NLS,'d')
        self.L = zeros(self.NLS,'d')
        for i in range(self.NLS):
            self.xc[i] = self.lsList[i].xc
            self.yc[i] = self.lsList[i].yc
            self.L[i] = self.lsList[i].L
            zcp = complex(self.xc[i],self.yc[i]) + \
                  complex(0.0, 0.5 * self.width) * exp( complex( 0.0, self.lsList[i].alpha ) )
            self.xcp[i] = zcp.real
            self.ycp[i] = zcp.imag
        self.delhead = zeros(self.NLS-1,'d')  # Should be dependent on flow in future
        for ls in self.lsList:
            ls.setCoefs()
        return
    def addElementToCollection(self,ml):
        for ls in self.lsList:
            ml.addElementToCollection(ls)
    def potentialInfluence(self,aq,x,y,zvert=0,t=0):
        rv = zeros( (self.NLS,aq.Naquifers), 'd' )
        for i in range(self.NLS):
            rv[i,:] = self.lsList[i].potentialInfluence(aq,x,y)[0,:]  # This is an inconsistency in Numeric
            # The part [0,:] should not have been necessary, but alas
        return rv
    def dischargeInfluence(self,aq,x,y,zvert=0,t=0):
        rvx = zeros( (self.NLS,aq.Naquifers), 'd' )
        rvy = zeros( (self.NLS,aq.Naquifers), 'd' )
        for i in range(self.NLS):
            [disx,disy] = self.lsList[i].dischargeInfluence(aq,x,y)
            rvx[i,:] = disx[0,:]  # This is an inconsistency in Numeric
            rvy[i,:] = disy[0,:]
        return [rvx,rvy]
    def layout(self):
        rv = [ self.NLS+1,[],[] ]
        rv[1] = []; rv[2] = []
        for ls in self.lsList:
            rv[1] = rv[1] + [ls.x1]
            rv[2] = rv[2] + [ls.y1]
        rv[1] = rv[1] + [ls.x2]
        rv[2] = rv[2] + [ls.y2]
        return rv
    def getMatrixRowsOld(self,elementList):
        rows=[]
        pylayer = self.pylayers[0]  # For now screened in only one layer
        for i in range(self.NLS-1):
            row = zeros(0,'d')
            aq1 = self.lsList[i+1].aquiferParent; T1 = aq1.T[pylayer]
            aq0 = self.lsList[i].aquiferParent; T0 = aq0.T[pylayer]
            for e in elementList:
                rowpart = e.getMatrixCoefficients( aq1,pylayer,self.xc[i+1],self.yc[i+1],\
                    lambda el,aq,pylayer,x,y: T0*el.potentialInfluenceInLayer(aq,pylayer,x,y) - \
                    T1*el.potentialInfluenceInLayer(aq0,pylayer,self.xc[i],self.yc[i]) )
                if e == self:
                    rowpart[i+1] = rowpart[i+1] - T0 * T1 * self.res / self.width
                    rowpart[i] = rowpart[i] + T0 * T1 * self.res / self.width
                row = hstack(( row, rowpart ))
            delpot = T0 * T1 * ( self.delhead[i] - aq1.hstar + aq0.hstar )
            delpotold = T0 * self.modelParent.potentialInLayer(self.aquiferParent,pylayer,self.xc[i+1],self.yc[i+1]) -\
                        T1 * self.modelParent.potentialInLayer(self.aquiferParent,pylayer,self.xc[i],self.yc[i])
            delsigma = ( self.parameters[i+1,0] - self.parameters[i,0] ) * T0 * T1 * self.res / self.width
            row = hstack(( row, delpot - delpotold + delsigma ))
            rows = rows + [row.tolist()]
        # Last equation is total discharge = Q
        row = zeros(0,'d')
        for e in elementList:
            if e == self:
                row = hstack(( row, self.L ))
            else:
                row = hstack(( row, e.getMatrixCoefficients(self.aqdum,self.ldum,self.xdum,self.ydum,\
                                lambda el,aqdum,ldum,xdum,ydum:el.zeroFunction(aqdum,ldum,xdum,ydum)) ))
        row = hstack(( row, self.Q - sum(self.parameters[:,0] * self.L) ))
        rows = rows + [row.tolist()]
        return rows
    # Modified to evaluate at xcp,ycp
    def getMatrixRows(self,elementList):
        rows=[]
        pylayer = self.pylayers[0]  # For now screened in only one layer
        for i in range(self.NLS-1):
            row = zeros(0,'d')
            aq1 = self.lsList[i+1].aquiferParent; T1 = aq1.T[pylayer]
            aq0 = self.lsList[i].aquiferParent; T0 = aq0.T[pylayer]
            for e in elementList:
                rowpart = e.getMatrixCoefficients( aq1,pylayer,self.xcp[i+1],self.ycp[i+1],\
                    lambda el,aq,pylayer,x,y: T0*el.potentialInfluenceInLayer(aq,pylayer,x,y) - \
                    T1*el.potentialInfluenceInLayer(aq0,pylayer,self.xcp[i],self.ycp[i]) )
                if e == self:
                    rowpart[i+1] = rowpart[i+1] - T0 * T1 * self.res / self.width
                    rowpart[i] = rowpart[i] + T0 * T1 * self.res / self.width
                row = hstack(( row, rowpart ))
            delpot = T0 * T1 * ( self.delhead[i] - aq1.hstar + aq0.hstar )
            delpotold = T0 * self.modelParent.potentialInLayer(self.aquiferParent,pylayer,self.xcp[i+1],self.ycp[i+1]) -\
                        T1 * self.modelParent.potentialInLayer(self.aquiferParent,pylayer,self.xcp[i],self.ycp[i])
            delsigma = ( self.parameters[i+1,0] - self.parameters[i,0] ) * T0 * T1 * self.res / self.width
            row = hstack(( row, delpot - delpotold + delsigma ))
            rows = rows + [row.tolist()]
        # Last equation is total discharge = Q
        row = zeros(0,'d')
        for e in elementList:
            if e == self:
                row = hstack(( row, self.L ))
            else:
                row = hstack(( row, e.getMatrixCoefficients(self.aqdum,self.ldum,self.xdum,self.ydum,\
                                lambda el,aqdum,ldum,xdum,ydum:el.zeroFunction(aqdum,ldum,xdum,ydum)) ))
        row = hstack(( row, self.Q - sum(self.parameters[:,0] * self.L) ))
        rows = rows + [row.tolist()]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,xc,yc,func):
        return func(self,aq,pylayer,xc,yc)
##    def takeParameters(self,xsol,icount):
##        for i in range(self.NLS):
##            self.parameters[i,0] = self.parameters[i,0] + xsol[icount]
##            icount = icount+1
##        return icount
    def takeParameters(self,xsol,icount):
        for i in range(self.NLS):
            self.parameters[i,0] = self.parameters[i,0] + xsol[icount]
            self.lsList[i].parameters[0,0] = self.lsList[i].parameters[0,0] + xsol[icount]
            self.lsList[i].paramxcoef = sum( self.lsList[i].parameters * self.lsList[i].coef, 0 )
            icount = icount+1
        return icount
    def checkOld(self):
        print 'Linesink Ditch '+str(self.xylist[0])+' to '+str(self.xylist[-1])
        print 'Specified discharge: '+str(self.Q)+' Computed discharge: '+str( sum(self.parameters[:,0] * self.L) )
        h0 = self.modelParent.head(self.layers[0],self.xc[0],self.yc[0]) - self.parameters[0,0]*self.res/self.width
        print 'Head in first line-sink: ' + str(h0)
        for i in range(self.NLS-1):
            h1 = self.modelParent.head(self.layers[0],self.xc[i+1],self.yc[i+1]) - self.parameters[i+1,0]*self.res/self.width
            h0 = self.modelParent.head(self.layers[0],self.xc[i],self.yc[i]) - self.parameters[i,0]*self.res/self.width
            print 'Linesink '+ str(i+1) + ' head ' + str(h0) + \
                  ' Specified delh: ' + str(self.delhead[i]) + ' Computed delh ' + \
                  str(h1-h0) + ' Strength: '+str(self.parameters[i,0])
        return None
    # Modified for xcp,ypc
    def check(self):
        print 'Linesink Ditch '+str(self.xylist[0])+' to '+str(self.xylist[-1])
        print 'Specified discharge: '+str(self.Q)+' Computed discharge: '+str( sum(self.parameters[:,0] * self.L) )
        h0 = self.modelParent.head(self.layers[0],self.xcp[0],self.ycp[0]) - self.parameters[0,0]*self.res/self.width
        print 'Head in first line-sink: ' + str(h0)
        for i in range(self.NLS-1):
            h1 = self.modelParent.head(self.layers[0],self.xcp[i+1],self.ycp[i+1]) - self.parameters[i+1,0]*self.res/self.width
            h0 = self.modelParent.head(self.layers[0],self.xcp[i],self.ycp[i]) - self.parameters[i,0]*self.res/self.width
            print 'Linesink '+ str(i+1) + ' head ' + str(h0) + \
                  ' Specified delh: ' + str(self.delhead[i]) + ' Computed delh ' + \
                  str(h1-h0) + ' Strength: '+str(self.parameters[i,0])
        return None
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):
        changed = 0; stop = 0; xyznew = zeros(3,'d')
        for ls in self.lsList:
            [changed, stop, xyznew] = ls.nearElement(pyLayer,xyz1,xyz2,step,idir)
            if changed or stop: break
        return [changed, stop, xyznew]  # Didn't adjust time
    def distanceSquaredToElement(self,x,y):
        Z = ( 2.0 * complex(x,y) - (self.z1 + self.z2) ) / (self.z2 - self.z1)
        if abs(Z.real) <= 1.0:
            dissq = ( Z.imag * self.Ltot / 2.0 )**2
        elif Z.real < -1.0:
            dissq = ( x - self.x1 )**2 + ( y - self.y1 ) **2
        else:
            dissq = ( x - self.x2 )**2 + ( y - self.y2 ) **2
        return dissq

class LineSinkDitchNew(Element):
    '''LineSinkDitch with everything variable, each segment is screened in just one layer
    line-sinks may be in different layers, but there is still an aquiferParent.....'''
    def __init__(self,modelParent,xylist,Q,res,width,layers=0,addToModel=1,label=None,aquiferParentFixed=None,xoffset = 0.0):
	Element.__init__(self,modelParent)
	if addToModel:
            self.modelParent.addElement(self)
        self.NLS = len(xylist) - 1
        self.xylist = xylist
        self.lsList = []
        if iterable(layers):
            assert len(layers) == self.NLS, "TimML Input error: layers must be length of numer of line-sinks"
            self.layers = array(layers,'i')
        else:
            self.layers = layers * ones( self.NLS, 'i' )
        #print 'xoffset ',xoffset
        for i in range(self.NLS):
            ls = LineSink(modelParent,xylist[i][0]-xoffset,xylist[i][1],xylist[i+1][0]+xoffset,xylist[i+1][1],0,[self.layers[i]],aquiferParentFixed,0)
            self.lsList = self.lsList + [ls]
        self.Q = float(Q)
        self.res = float(res)
        self.width = float(width)
##        if iterable(res):
##            assert len(res) == self.NLS, "TimML Input error: res must be length of numer of line-sinks"
##            self.res = array(res,'d')
##        else:
##            self.res = res * ones( self.NLS )
##        if iterable(width):
##            assert len(width) == self.NLS, "TimML Input error: width must be length of numer of line-sinks"
##            self.width = array(width,'d')
##        else:
##            self.width = width * ones( self.NLS )        
        if aquiferParentFixed is None:
            self.aquiferParent = self.modelParent.aq.findAquiferData(xylist[0][0],xylist[0][1])  # Determined at first node
        else:
            self.aquiferParent = aquiferParentFixed
        self.label = label
        self.type = 'linesink'
        self.setCoefs()
    def __repr__(self):
	return 'LineSinkDitch Q,res,width,layers ' + str((self.Q,self.res,self.width,list(self.layers))) +\
               ' consisting of the following ' + str(self.NLS) + ' linesinks ' + str(self.lsList)
    def setCoefs(self):
        self.x1 = self.xylist[0][0];  self.y1 = self.xylist[0][1]
        self.x2 = self.xylist[-1][0]; self.y2 = self.xylist[-1][1]
        self.z1 = complex(self.x1,self.y1); self.z2 = complex(self.x2,self.y2)
        self.Ltot = sqrt( (self.x2-self.x1)**2 + (self.y2-self.y1)**2 )
        self.parameters = zeros( (self.NLS,1), 'd' )
        self.pylayers = self.layers  # fixed to base zero
        self.xc = zeros(self.NLS,'d'); self.yc = zeros(self.NLS,'d')
        self.xcp = zeros(self.NLS,'d'); self.ycp = zeros(self.NLS,'d')
        self.L = zeros(self.NLS,'d')
        for i in range(self.NLS):
            self.xc[i] = self.lsList[i].xc
            self.yc[i] = self.lsList[i].yc
            self.L[i] = self.lsList[i].L
            zcp = complex(self.xc[i],self.yc[i]) + \
                  complex(0.0, 0.5 * self.width) * exp( complex( 0.0, self.lsList[i].alpha ) )
            self.xcp[i] = zcp.real
            self.ycp[i] = zcp.imag
        self.delhead = zeros(self.NLS-1,'d')  # Should be dependent on flow in future
        for ls in self.lsList:
            ls.setCoefs()
        return
    def addElementToCollection(self,ml):
        for ls in self.lsList:
            ml.addElementToCollection(ls)
    def potentialInfluence(self,aq,x,y,zvert=0,t=0):
        rv = zeros( (self.NLS,aq.Naquifers), 'd' )
        for i in range(self.NLS):
            rv[i,:] = self.lsList[i].potentialInfluence(aq,x,y)[0,:]  # This is an inconsistency in Numeric
            # The part [0,:] should not have been necessary, but alas
        return rv
    def dischargeInfluence(self,aq,x,y,zvert=0,t=0):
        rvx = zeros( (self.NLS,aq.Naquifers), 'd' )
        rvy = zeros( (self.NLS,aq.Naquifers), 'd' )
        for i in range(self.NLS):
            [disx,disy] = self.lsList[i].dischargeInfluence(aq,x,y)
            rvx[i,:] = disx[0,:]  # This is an inconsistency in Numeric
            rvy[i,:] = disy[0,:]
        return [rvx,rvy]
    def layout(self):
        rv = [ self.NLS+1,[],[] ]
        rv[1] = []; rv[2] = []
        for ls in self.lsList:
            rv[1] = rv[1] + [ls.x1]
            rv[2] = rv[2] + [ls.y1]
        rv[1] = rv[1] + [ls.x2]
        rv[2] = rv[2] + [ls.y2]
        return rv
    def getMatrixRows(self,elementList):
        rows=[]
        for i in range(self.NLS-1):
            row = zeros(0,'d')
            aq1 = self.lsList[i+1].aquiferParent
	    pylayer1 = self.pylayers[i+1]
	    if aq1.fakesemi: pylayer1 += 1
	    T1 = aq1.T[pylayer1]; 
            aq0 = self.lsList[i].aquiferParent
	    pylayer0 = self.pylayers[i]
	    if aq0.fakesemi: pylayer0 += 1
	    T0 = aq0.T[pylayer0]
            for e in elementList:
                rowpart = e.getMatrixCoefficients( aq1,pylayer1,self.xcp[i+1],self.ycp[i+1],\
                    lambda el,aq,pylayer,x,y: T0*el.potentialInfluenceInLayer(aq,pylayer,x,y) - \
                    T1*el.potentialInfluenceInLayer(aq0,pylayer0,self.xcp[i],self.ycp[i]) )
                if e == self:
                    rowpart[i+1] = rowpart[i+1] - T0 * T1 * self.res / self.width
                    rowpart[i] = rowpart[i] + T0 * T1 * self.res / self.width
                row = hstack(( row, rowpart ))
            delpot = T0 * T1 * self.delhead[i]
            delpotold = T0 * self.modelParent.potentialInLayer(aq1,pylayer1,self.xcp[i+1],self.ycp[i+1]) -\
                        T1 * self.modelParent.potentialInLayer(aq0,pylayer0,self.xcp[i],self.ycp[i])
            delsigma = ( self.parameters[i+1,0] - self.parameters[i,0] ) * T0 * T1 * self.res / self.width
            row = hstack(( row, delpot - delpotold + delsigma ))
            rows = rows + [row.tolist()]
        # Last equation is total discharge = Q
        row = zeros(0,'d')
        for e in elementList:
            if e == self:
                row = hstack(( row, self.L ))
            else:
                row = hstack(( row, e.getMatrixCoefficients(self.aqdum,self.ldum,self.xdum,self.ydum,\
                                lambda el,aqdum,ldum,xdum,ydum:el.zeroFunction(aqdum,ldum,xdum,ydum)) ))
        row = hstack(( row, self.Q - sum(self.parameters[:,0] * self.L) ))
        rows = rows + [row.tolist()]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,xc,yc,func):
        return func(self,aq,pylayer,xc,yc)
    def takeParameters(self,xsol,icount):
        for i in range(self.NLS):
            self.parameters[i,0] = self.parameters[i,0] + xsol[icount]
            self.lsList[i].parameters[0,0] = self.lsList[i].parameters[0,0] + xsol[icount]
            self.lsList[i].paramxcoef = sum( self.lsList[i].parameters * self.lsList[i].coef, 0 )
            icount = icount+1
        return icount
    def check(self):
        print 'Linesink Ditch '+str(self.xylist[0])+' to '+str(self.xylist[-1])
        print 'Specified discharge: '+str(self.Q)+' Computed discharge: '+str( sum(self.parameters[:,0] * self.L) )
        h0 = self.modelParent.head(self.layers[0],self.xcp[0],self.ycp[0]) - self.parameters[0,0]*self.res/self.width
        print 'Head in first line-sink: ' + str(h0)
        for i in range(self.NLS-1):
            h1 = self.modelParent.head(self.layers[i+1],self.xcp[i+1],self.ycp[i+1]) - self.parameters[i+1,0]*self.res/self.width
            h0 = self.modelParent.head(self.layers[i],self.xcp[i],self.ycp[i]) - self.parameters[i,0]*self.res/self.width
            print 'Linesink '+ str(i+1) + ' head ' + str(h0) + \
                  ' Specified delh: ' + str(self.delhead[i]) + ' Computed delh ' + \
                  str(h1-h0) + ' Strength: '+str(self.parameters[i,0])
        return None
    def headAtDitch(self):
        return self.modelParent.head(self.layers[0],self.xcp[0],self.ycp[0])
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):
        changed = 0; stop = 0; xyznew = zeros(3,'d')
        for ls in self.lsList:
            [changed, stop, xyznew] = ls.nearElement(pyLayer,xyz1,xyz2,step,idir)
            if changed or stop: break
        return [changed, stop, xyznew]  # Didn't adjust time
    def distanceSquaredToElement(self,x,y):
        Z = ( 2.0 * complex(x,y) - (self.z1 + self.z2) ) / (self.z2 - self.z1)
        if abs(Z.real) <= 1.0:
            dissq = ( Z.imag * self.Ltot / 2.0 )**2
        elif Z.real < -1.0:
            dissq = ( x - self.x1 )**2 + ( y - self.y1 ) **2
        else:
            dissq = ( x - self.x2 )**2 + ( y - self.y2 ) **2
        return dissq

class DoubleLineSink(LineSink):
    '''Class for double linesinks one for outside, one for inside
    Prototype class; designed for one layer model for now'''
    def __init__(self,modelParent,x1,y1,x2,y2,layers,aqin,aqout,labelin=None):
        # DoubleLineSink is linesink, but we overload bunch of functions
        LineSink.__init__(self,modelParent,x1,y1,x2,y2,0,layers,aquiferParentFixed=None,addToModel=0,label=labelin)
        self.aqin = aqin; self.aqout = aqout
	self.lsin  = LineSink(modelParent,x1,y1,x2,y2,0,layers,aqin,0)
        self.lsout = LineSink(modelParent,x1,y1,x2,y2,0,layers,aqout,0)
        Zout = complex(0,-1e-6); zcout = Zout * (self.z2 - self.z1) / 2.0 + 0.5 * (self.z1 + self.z2)
        Zin = complex(0,1e-6); zcin = Zin * (self.z2 - self.z1) / 2.0 + 0.5 * (self.z1 + self.z2)
        self.xcout = zcout.real; self.ycout = zcout.imag
        self.xcin  = zcin.real;  self.ycin = zcin.imag
        self.Nparamin = self.aqin.Naquifers; self.Nparamout = self.aqout.Naquifers
        self.modelParent.addElement(self)
    def __repr__(self):
	return 'DoubleLineSink z1,z2,sigma,layers: ' + str((self.x1,self.y1,self.x2,self.y2,self.sigma,list(self.layers)))
    def potentialInfluence(self,aq,x,y):
        '''Returns array of (iorder,aq.Naquifers)'''
        if self.aqin == aq:
            rv = self.lsin.potentialInfluence(aq,x,y)
        elif self.aqout == aq:
            rv = self.lsout.potentialInfluence(aq,x,y)
        else:  # For now; this has to be changed later of course
            rv = self.lsout.potentialInfluence(aq,x,y)  # I think this does everything right
        return (rv[0,:])[newaxis,:]  # Just take first row; this is an ugly fix, but we are only interested in Fi_L and Fi_m's, not the coef part
    def potentialInfluenceInLayer(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in pylayer as array (1 value per parameter)
        Needs to be overloaded because there is inside and outside'''
        potInf = self.potentialInfluence(aq,x,y)
        rv = zeros(0,'d')
        for p in potInf:
            rv = hstack(( rv, p * aq.eigvec[pylayer,:] ))
        if aq == self.aqin:
            rv = hstack(( rv, zeros( self.Nparamout, 'd' ) ))
        else:
            rv = hstack(( zeros( self.Nparamin, 'd' ), rv ))
        return rv
    def potentialInfluenceAllLayers(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in all layers as an array'''
        potInf = self.potentialInfluence(aq,x,y)
        rv = zeros((aq.Naquifers,0),'d')
        for p in potInf:
            rv = hstack( ( rv, p*aq.eigvec ) )
        if aq == self.aqin:
            rv = hstack(( rv, zeros( ( aq.Naquifers, self.Nparamout ), 'd' ) ))
        else:
            rv = hstack(( zeros( (aq.Naquifers,self.Nparamin), 'd' ), rv ))
        return rv
    def potentialContribution(self,aq,x,y):
        '''Returns array of potentialContribution. Needs to be overloaded cause there is inside and outside'''
        if aq == self.aqin:
            potInf = self.lsin.parameters[:,0] * self.potentialInfluence(aq,x,y)
        elif aq == self.aqout:
            potInf = self.lsout.parameters[:,0] * self.potentialInfluence(aq,x,y)
        return potInf
    def dischargeInfluence(self,aq,x,y):
        '''Returns two arrays of (1,aq.Naquifers)'''
        if self.aqin == aq:
            [rvx,rvy] = self.lsin.dischargeInfluence(aq,x,y)
        elif self.aqout == aq:
            [rvx,rvy] = self.lsout.dischargeInfluence(aq,x,y)
        else:
            [rvx,rvy] = [zeros((1,aq.Naquifers),'d'),zeros((1,aq.Naquifers),'d')]
        return [(rvx[0,:])[newaxis,:],(rvy[0,:])[newaxis,:]]  # Same ugly fix as for potential
    def dischargeInfluenceAllLayers(self,aq,dumlayer,x,y):
        '''Returns dischargeInfluenceAllLayers function in aquifer aq as an array
        Needs to be overloaded because there is no parameter outside the functions
        Needs to be modified for inside and outside'''
        [disx,disy] = self.dischargeInfluence(aq,x,y)
        rvx = zeros((aq.Naquifers,0),'d')
        rvy = zeros((aq.Naquifers,0),'d')
        for d in disx:
            rvx = hstack( ( rvx, d*aq.eigvec ) )
        for d in disy:
            rvy = hstack( ( rvy, d*aq.eigvec ) )
        if aq == self.aqin:
            rvx = hstack(( rvx, zeros( ( aq.Naquifers, self.Nparamout ), 'd' ) ))
            rvy = hstack(( rvy, zeros( ( aq.Naquifers, self.Nparamout ), 'd' ) ))
        else:
            rvx = hstack(( zeros( (aq.Naquifers,self.Nparamin), 'd' ), rvx ))
            rvy = hstack(( zeros( (aq.Naquifers,self.Nparamin), 'd' ), rvy ))
        return [rvx,rvy]
    def dischargeContribution(self,aq,x,y):
        '''Returns matrix with two rowvectors of dischargeContributions Qx and Qy. Needs to be overloaded cause there is inside and outside'''
        disInf = self.dischargeInfluence(aq,x,y)
        if aq == self.aqin:
            param = self.lsin.parameters[:,0]  # Convert to row
        elif aq == self.aqout:
            param = self.lsout.parameters[:,0]  # Convert to row
        disxInf = sum( param * disInf[0], 0 )
        disyInf = sum( param * disInf[1], 0 )
        return array([disxInf,disyInf])
    def getMatrixRows(self,elementList):
        rows=[]
        pylayer = self.pylayers[0]
        rowin = zeros((self.aqin.Naquifers,0),'d'); rowout = zeros((self.aqout.Naquifers,0),'d')
        for e in elementList:
            rowinpart = e.getMatrixCoefficients(self.aqin,pylayer,self.xcin,self.ycin,\
                                            lambda el,aq,pylayer,x,y:el.potentialInfluenceAllLayers(aq,pylayer,x,y))
            if size(rowinpart) > 0:
                rowin = hstack(( rowin, rowinpart ))
            rowoutpart = e.getMatrixCoefficients(self.aqout,pylayer,self.xcout,self.ycout,\
                                            lambda el,aq,pylayer,x,y:el.potentialInfluenceAllLayers(aq,pylayer,x,y))
            if size(rowoutpart) > 0:
                rowout = hstack(( rowout, rowoutpart ))
        row = self.aqout.Tcol * rowin - self.aqin.Tcol * rowout
        row = hstack( ( row,
            self.aqin.Tcol * self.aqout.Tcol * ( self.aqout.hstar - self.aqin.hstar ) + \
            self.aqin.Tcol * self.modelParent.potentialVectorCol(self.xcout,self.ycout,self.aqout) - \
            self.aqout.Tcol * self.modelParent.potentialVectorCol(self.xcin,self.ycin,self.aqin) ) )
        rows = rows + row.tolist()
        rowin = zeros((self.aqin.Naquifers,0),'d'); rowout = zeros((self.aqout.Naquifers,0),'d')
        for e in elementList:
            rowqxqyin = e.getMatrixCoefficients(self.aqin,pylayer,self.xcin,self.ycin,\
                            lambda el,aq,pylayer,x,y:el.dischargeInfluenceAllLayers(aq,pylayer,x,y))
            if size(rowqxqyin) > 0:
                rowpart = rowqxqyin[0] * cos( self.thetaNormOut ) + rowqxqyin[1] * sin( self.thetaNormOut )
                rowin = hstack(( rowin, rowpart ))
            rowqxqyout = e.getMatrixCoefficients(self.aqout,pylayer,self.xcout,self.ycout,\
                            lambda el,aq,pylayer,x,y:el.dischargeInfluenceAllLayers(aq,pylayer,x,y))
            if size(rowqxqyout) > 0:
                rowpart = rowqxqyout[0] * cos( self.thetaNormOut ) + rowqxqyout[1] * sin( self.thetaNormOut )
                rowout = hstack(( rowout, rowpart ))
        row = rowin - rowout
        row = hstack(( row,
            self.modelParent.dischargeNormVectorCol(self.xcout,self.ycout,self.thetaNormOut,self.aqout) -\
            self.modelParent.dischargeNormVectorCol(self.xcin,self.ycin,self.thetaNormOut,self.aqin) ))
        rows = rows + row.tolist()
        return rows
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        self.lsin.parameters[:self.aqin.Naquifers,0] = self.lsin.parameters[:self.aqin.Naquifers,0] + xsol[icount:icount+self.aqin.Naquifers]
        icount = icount + self.aqin.Naquifers
        self.lsout.parameters[:self.aqout.Naquifers,0] = self.lsout.parameters[:self.aqout.Naquifers,0] + xsol[icount:icount+self.aqout.Naquifers]
        icount = icount + self.aqout.Naquifers
        return icount
    def check(self):
        print 'DoubleLineSink from '+str(self.z1)+' to '+str(self.z2)+' in layer '+str(self.pylayers[0]+1)
        print 'Control point: '+str((self.xc,self.yc))+\
              ' Dis inside : '+str(self.modelParent.dischargeNormInLayer(self.aqin,self.pylayers[0],self.xcin,self.ycin,self.thetaNormOut))+\
              ' Dis outside: '+str(self.modelParent.dischargeNormInLayer(self.aqout,self.pylayers[0],self.xcout,self.ycout,self.thetaNormOut))
        print 'Control point: '+str((self.xc,self.yc))+\
              ' head inside : '+str(self.modelParent.head(self.layers[0],self.xcin,self.ycin,self.aqin))+\
              ' Head outside: '+str(self.modelParent.head(self.layers[0],self.xcout,self.ycout,self.aqout))
        return None


class DoubleLineSinkConf(DoubleLineSink):
    '''Class for double linesinks one for outside, one for inside
    Prototype class; designed for one layer model for now'''
    def __init__(self,modelParent,x1,y1,x2,y2,layers,aqin,aqout,labelin=None):
        # DoubleLineSink is linesink, but we overload bunch of functions
        DoubleLineSink.__init__(self,modelParent,x1,y1,x2,y2,layers,aqin,aqout,labelin)
    def __repr__(self):
	return 'DoubleLineSinkConf z1,z2,sigma,layers: ' + str((self.x1,self.y1,self.x2,self.y2,self.sigma,list(self.layers)))
    def potentialInfluence(self,aq,x,y):
        '''Returns array of (iorder,aq.Naquifers)'''
        if self.aqin == aq:
            rv = self.lsin.potentialInfluence(aq,x,y)
        elif self.aqout == aq:
            rv = self.lsout.potentialInfluence(aq,x,y)
        else:  # For now; this has to be changed later of course
            rv = self.lsout.potentialInfluence(aq,x,y)  # I think this does everything right
        return rv  # Just take first row; this is an ugly fix, but we are only interested in Fi_L and Fi_m's, not the coef part
    def potentialInfluenceInLayer(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in pylayer as 1D array (1 value per parameter)'''
        potInf = self.potentialInfluence(aq,x,y)
        return sum( potInf * aq.eigvec[pylayer,:], 1 )
    def potentialInfluenceInLayer(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in pylayer as array (1 value per parameter)
        Needs to be overloaded because there is inside and outside'''
        potInf = self.potentialInfluence(aq,x,y)
        rv = zeros(0,'d')
        if aq == self.aqin:
            rv = hstack(( sum( potInf * aq.eigvec[pylayer,:], 1 ), zeros( self.Nparamout, 'd' ) ))
        else:
            rv = hstack(( zeros( self.Nparamin, 'd' ), sum( potInf * aq.eigvec[pylayer,:], 1 ) ))
        return rv
    def potentialInfluenceAllLayers(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in all layers as an array'''
        potInf = self.potentialInfluence(aq,x,y)
        rv = zeros((aq.Naquifers,0),'d')
        for p in potInf:  # Now different so that parameters are indeed strengths of elements
            pot = sum(p*aq.eigvec,1)  # Now a row vector; This assumes parameter is outside PotL and Potm-s
            rv = hstack( ( rv, pot[:,newaxis] ) )  
        if aq == self.aqin:
            rv = hstack(( rv, zeros( ( aq.Naquifers, self.Nparamout ), 'd' ) ))
        else:
            rv = hstack(( zeros( (aq.Naquifers,self.Nparamin), 'd' ), rv ))
        return rv
    def potentialContribution(self,aq,x,y):
        '''Returns array of potentialContribution. Needs to be overloaded cause there is inside and outside'''
        if aq == self.aqin:
            param = self.lsin.parameters
        elif aq == self.aqout:
            param = self.lsout.parameters
        return sum( param * self.potentialInfluence(aq,x,y), 0 )
    def dischargeInfluence(self,aq,x,y):
        '''Returns two arrays of (1,aq.Naquifers)'''
        if self.aqin == aq:
            [rvx,rvy] = self.lsin.dischargeInfluence(aq,x,y)
        elif self.aqout == aq:
            [rvx,rvy] = self.lsout.dischargeInfluence(aq,x,y)
        else:
            [rvx,rvy] = [zeros((1,aq.Naquifers),'d'),zeros((1,aq.Naquifers),'d')]
        return [rvx,rvy]  # Same ugly fix as for potential
    def dischargeInfluenceAllLayers(self,aq,dumlayer,x,y):
        '''Returns dischargeInfluenceAllLayers function in aquifer aq as an array
        Needs to be overloaded because there is no parameter outside the functions
        Needs to be modified for inside and outside'''
        [disx,disy] = self.dischargeInfluence(aq,x,y)
        rvx = zeros((aq.Naquifers,0),'d')
        rvy = zeros((aq.Naquifers,0),'d')
        for d in disx:
            dis = sum( d*aq.eigvec, 1 )
            rvx = hstack( ( rvx, dis[:,newaxis] ) )
        for d in disy:
            dis = sum( d*aq.eigvec, 1 )
            rvy = hstack( ( rvy, dis[:,newaxis] ) )
        if aq == self.aqin:
            rvx = hstack(( rvx, zeros( ( aq.Naquifers, self.Nparamout ), 'd' ) ))
            rvy = hstack(( rvy, zeros( ( aq.Naquifers, self.Nparamout ), 'd' ) ))
        else:
            rvx = hstack(( zeros( (aq.Naquifers,self.Nparamin), 'd' ), rvx ))
            rvy = hstack(( zeros( (aq.Naquifers,self.Nparamin), 'd' ), rvy ))
        return [rvx,rvy]
    def dischargeContribution(self,aq,x,y):
        '''Returns matrix with two rowvectors of dischargeContributions Qx and Qy. Needs to be overloaded cause there is inside and outside'''
        disInf = self.dischargeInfluence(aq,x,y)
        if aq == self.aqin:
            param = self.lsin.parameters
        elif aq == self.aqout:
            param = self.lsout.parameters
        return array([ sum( param * disInf[0], 0 ), sum( param * disInf[1], 0 ) ])
    def zeroFunction(self,aqdum,ldum,xdum,ydum):
        '''Returns list of zeros of length number of parameters'''
        return list( zeros( 2*len(self.parameters),'d' ) )
    def getMatrixRows(self,elementList):
        rows=[]
        pylayer = self.pylayers[0]
        rowin = zeros((self.aqin.Naquifers,0),'d'); rowout = zeros((self.aqout.Naquifers,0),'d')
        for e in elementList:
            rowinpart = e.getMatrixCoefficients(self.aqin,pylayer,self.xcin,self.ycin,\
                                            lambda el,aq,pylayer,x,y:el.potentialInfluenceAllLayers(aq,pylayer,x,y))
            if size(rowinpart) > 0:
                rowin = hstack(( rowin, rowinpart ))
            rowoutpart = e.getMatrixCoefficients(self.aqout,pylayer,self.xcout,self.ycout,\
                                            lambda el,aq,pylayer,x,y:el.potentialInfluenceAllLayers(aq,pylayer,x,y))
            if size(rowoutpart) > 0:
                rowout = hstack(( rowout, rowoutpart ))
        row = self.aqout.Tcol * rowin - self.aqin.Tcol * rowout
        row = hstack( ( row,
            self.aqin.Tcol * self.aqout.Tcol * ( self.aqout.hstar - self.aqin.hstar ) + \
            self.aqin.Tcol * self.modelParent.potentialVectorCol(self.xcout,self.ycout,self.aqout) - \
            self.aqout.Tcol * self.modelParent.potentialVectorCol(self.xcin,self.ycin,self.aqin) ) )
        rows = rows + row[1:,:].tolist()  # Took off first equation (that is where line-doublet is taking care)
        rowin = zeros((self.aqin.Naquifers,0),'d'); rowout = zeros((self.aqout.Naquifers,0),'d')
        for e in elementList:
            rowqxqyin = e.getMatrixCoefficients(self.aqin,pylayer,self.xcin,self.ycin,\
                            lambda el,aq,pylayer,x,y:el.dischargeInfluenceAllLayers(aq,pylayer,x,y))
            if size(rowqxqyin) > 0:
                rowpart = rowqxqyin[0] * cos( self.thetaNormOut ) + rowqxqyin[1] * sin( self.thetaNormOut )
                rowin = hstack(( rowin, rowpart ))
            rowqxqyout = e.getMatrixCoefficients(self.aqout,pylayer,self.xcout,self.ycout,\
                            lambda el,aq,pylayer,x,y:el.dischargeInfluenceAllLayers(aq,pylayer,x,y))
            if size(rowqxqyout) > 0:
                rowpart = rowqxqyout[0] * cos( self.thetaNormOut ) + rowqxqyout[1] * sin( self.thetaNormOut )
                rowout = hstack(( rowout, rowpart ))
        row = rowin - rowout
        row = hstack(( row,
            self.modelParent.dischargeNormVectorCol(self.xcout,self.ycout,self.thetaNormOut,self.aqout) -\
            self.modelParent.dischargeNormVectorCol(self.xcin,self.ycin,self.thetaNormOut,self.aqin) ))
        rows = rows + row.tolist()
        # Last equation is sum of strength parameters on inside = 0
        row = zeros(0,'d')
        for e in elementList:
            if e == self:
                row = hstack(( row, ones(self.NscreenedLayers,'d'), zeros(self.NscreenedLayers,'d') ))
            else:
                row = hstack(( row, e.getMatrixCoefficients(self.aqdum,self.ldum,self.xdum,self.ydum,\
                                         lambda el,aqdum,ldum,xdum,ydum:el.zeroFunction(aqdum,ldum,xdum,ydum)) ))
        #row = hstack(( row, - sum(self.lsin.parameters,0)[0] ))
        row = hstack(( row, - sum(self.lsin.parameters[:,0]) ))
        rows = rows + [row.tolist()]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        self.lsin.parameters[:self.aqin.Naquifers,0] = self.lsin.parameters[:self.aqin.Naquifers,0] + xsol[icount:icount+self.aqin.Naquifers]
        icount = icount + self.aqin.Naquifers
        self.lsout.parameters[:self.aqout.Naquifers,0] = self.lsout.parameters[:self.aqout.Naquifers,0] + xsol[icount:icount+self.aqout.Naquifers]
        icount = icount + self.aqout.Naquifers
        return icount
    def check(self):
        print 'DoubleLineSink from '+str(self.z1)+' to '+str(self.z2)+' in layer '+str(self.pylayers[0]+1)
        print 'Control point: '+str((self.xc,self.yc))+\
              '\n Dis inside : '+str(self.modelParent.dischargeNormVector(self.xcin,self.ycin,self.thetaNormOut,self.aqin))+\
              '\n Dis outside: '+str(self.modelParent.dischargeNormVector(self.xcout,self.ycout,self.thetaNormOut,self.aqout))
        print 'Control point: '+str((self.xc,self.yc))+\
              '\n Head inside : '+str(self.modelParent.headVector(self.xcin,self.ycin,self.aqin))+\
              '\n Head outside: '+str(self.modelParent.headVector(self.xcout,self.ycout,self.aqout))
        return None

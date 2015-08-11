'''
mllinedoubletimp.py contains the LineDoubletImp class.
This file is part of the TimML library and is distributed under
the GNU LPGL. See the TimML.py file for more details.
(c) Mark Bakker, 2002-2008
'''

from numpy import *
import numpy.linalg as linalg
from mlelement import *
from besselaes import *

class LineDoubletImp(Element):
    def __init__(self,modelParent,x1,y1,x2,y2,order=0,layers=[0],label=None,lab=[0.0],aquiferParentFixed=None,addToModel=1):
	Element.__init__(self,modelParent)
        self.x1 = float(x1); self.y1 = float(y1); self.x2 = float(x2); self.y2 = float(y2);
        self.order = order
        self.layers = array(layers)
        if len(lab) < self.order+1:
            self.lab = zeros(self.order+1)
        else:
            self.lab = array(lab) 
        self.aquiferParentFixed = aquiferParentFixed
        self.label = label
        self.type = 'linedoubletimp'
        self.setCoefs()
        if addToModel:
            self.modelParent.addElement(self)
    def __repr__(self):
	return 'LineDoubletImp z1,z2,lab,order,layers: ' + str((self.x1,self.y1,self.x2,self.y2,self.lab,self.order,list(self.layers)))
    def setCoefs(self):
        self.xy1 = (self.x1,self.y1); self.xy2 = (self.x2,self.y2)
        self.xmin = min(self.x1,self.x2); self.xmax = max(self.x1,self.x2)
        self.ymin = min(self.y1,self.y2); self.ymax = max(self.y1,self.y2)
        self.thetaNormOut = arctan2(self.y2-self.y1,self.x2-self.x1) - pi/2.0
        self.cosout = cos( self.thetaNormOut ); self.sinout = sin( self.thetaNormOut )
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
      	self.pylayers = self.layers  # fixed to base zero
	self.NscreenedLayers = len(self.layers)
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
                print 'semi-confined no longer supported'
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
        self.potLapLdInf = zeros(self.Ndegree,'d'); self.potInf = zeros(self.aquiferParent.Naquifers,'d')
        self.disxLapLdInf = zeros(self.Ndegree,'d'); self.disyLapLdInf = zeros(self.Ndegree,'d');
        self.disxInf = zeros(self.aquiferParent.Naquifers,'d'); self.disyInf = zeros(self.aquiferParent.Naquifers,'d')
    def potentialInfluence(self,aq,x,y):
        '''Returns array of (NscreenedLayers*Order,aq.Naquifers)
        Sequence: first all order 0 linedoublets, than all order 1, etc.'''
        rv = zeros((self.NscreenedLayers * self.Ndegree,aq.Naquifers),'d')
        # Order of if-s is slightly different than for wells, because Laplace part is not computed separately
        if self.aquiferParent.type == self.aquiferParent.conf:
            if self.aquiferParent == aq:  # Same confined aquifer
                lab = self.aquiferParent.zeropluslab
                for iorder in range(self.Ndegree):
                    potbesldho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers,lab,iorder,self.potInf)
                    # Call FORTRAN extension
                    rv[iorder*self.NscreenedLayers:iorder*self.NscreenedLayers+self.NscreenedLayers,:] = \
                        self.coef * self.potInf
            else:  # Different confined aquifer, Laplace part only
                potlapldho(x,y,self.x1,self.y1,self.x2,self.y2,self.order,self.potLapLdInf)
                for iorder in range(self.Ndegree):
                    rv[iorder*self.NscreenedLayers:iorder*self.NscreenedLayers+self.NscreenedLayers,0] = \
                        self.potLapLdInf[iorder]
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
                    disbesldho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers,lab,ior,self.disxInf,self.disyInf)
                    # Call FORTRAN extension
                    rvx[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,:] = \
                        self.coef * self.disxInf
                    rvy[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,:] = \
                        self.coef * self.disyInf                
            else:  # Different confined aquifer, Laplace part only
                dislapldho(x,y,self.x1,self.y1,self.x2,self.y2,self.order,self.disxLapLdInf,self.disyLapLdInf)
                for ior in range(self.Ndegree):
                    rvx[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,0] = \
                        self.disxLapLdInf[ior]
                    rvy[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,0] = \
                        self.disyLapLdInf[ior]
        elif self.aquiferParent.type == self.aquiferParent.semi:
            print 'not yet implemented'
        return [rvx,rvy]
    def potentialCollection(self,potsum,potadd,elementList,aq,x,y):
        for el in elementList:
            if el.aquiferParent.type == el.aquiferParent.conf:
                if el.aquiferParent == aq:  # Same confined aquifer
                    lab = el.aquiferParent.zeropluslab
                    for ior in range(el.Ndegree):
                        potbesldho(x,y,el.x1,el.y1,el.x2,el.y2,el.aquiferParent.Naquifers,lab,ior,el.potInf)
                        potsum = potsum + el.paramxcoef[ior,:] * el.potInf
                elif aq.type == aq.conf:  # Different confined aquifer, Laplace part only
                    potlapldho(x,y,el.x1,el.y1,el.x2,el.y2,el.order,el.potLapLdInf)  # Call FORTRAN extension
                    for ior in range(el.Ndegree):
                        potsum[0] = potsum[0] + el.paramxcoef[ior,0] * el.potLapLdInf[ior]
            elif el.aquiferParent.type == el.aquiferParent.semi:
                print 'not yet implemented'
        return potsum
    def dischargeCollection(self,dissum,disadd,elementList,aq,x,y):
        for el in elementList:
            if el.aquiferParent.type == el.aquiferParent.conf:
                if el.aquiferParent == aq:  # Same confined aquifer
                    lab = el.aquiferParent.zeropluslab
                    for ior in range(el.Ndegree):
                        disbesldho(x,y,el.x1,el.y1,el.x2,el.y2,el.aquiferParent.Naquifers,lab,ior,disadd[0,:],disadd[1,:])
                        dissum = dissum + el.paramxcoef[ior,:] * disadd
                elif aq.type == aq.conf:  # Different confined aquifer, Laplace part only
                    dislapldho(x,y,el.x1,el.y1,el.x2,el.y2,el.order,el.disxLapLdInf,el.disyLapLdInf)
                    for ior in range(el.Ndegree):
                        dissum[0,0] = dissum[0,0] + el.paramxcoef[ior,0] * el.disxLapLdInf[ior]
                        dissum[1,0] = dissum[1,0] + el.paramxcoef[ior,0] * el.disyLapLdInf[ior]
            elif el.aquiferParent.type == el.aquiferParent.semi:
                print 'not yet implemented'
        return dissum
    def layout(self):
        rv = [ 2,[],[] ]
        rv[1] = [ self.x1,self.x2 ]
        rv[2] = [ self.y1,self.y2 ]
        return rv
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):
        # I don't think there is anything to do
        changed = 0; stop = 0; xyznew = 0.0
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
        for pylayer in self.pylayers:
            for icp in range(self.Ncp):
                xc = self.xcp[icp]; yc = self.ycp[icp]
                row = []
                for e in elementList:
                    rowqxqy = e.getMatrixCoefficients(self.aquiferParent,pylayer,xc,yc,\
                                lambda el,aq,pylayer,x,y:el.dischargeInfluenceInLayer(aq,pylayer,x,y))                    
                    if size(rowqxqy) > 0:
                        rowpart = rowqxqy[0] * self.cosout + rowqxqy[1] * self.sinout
                        row = row + rowpart.tolist()
                rhs = -self.modelParent.dischargeNormInLayer(self.aquiferParent,pylayer,xc,yc,self.thetaNormOut)
                row = row + [rhs] # Add rhs as last term
                rows = rows + [row]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,xc,yc,func):
        rv = func(self,aq,pylayer,xc,yc)
        return rv
    def takeParameters(self,xsol,icount):
        for i in range(self.Ncp*self.NscreenedLayers):
            self.parameters[i,0] = self.parameters[i,0] + xsol[icount]
            icount = icount+1
        for ior in range(self.Ndegree):
            self.paramxcoef[ior,:] = \
            sum( self.parameters[ior*self.NscreenedLayers:ior*self.NscreenedLayers+self.NscreenedLayers,:] * self.coef, 0 )  # Parameters times coefficients
        return icount
    def check(self):
        print 'Linedoublet from '+str(self.xy1)+' to '+str(self.xy2)
        for icp in range(self.Ncp):
            print 'Control point '+str(icp)
            for i in range(self.NscreenedLayers):
                print 'Layer '+str(self.layers[i])+' Dis: '+\
                  str( self.modelParent.dischargeNormInLayer\
                      (self.aquiferParent,self.pylayers[i],self.xcp[icp],self.ycp[icp],self.thetaNormOut) )


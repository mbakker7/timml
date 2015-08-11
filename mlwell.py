from numpy import *
import numpy.linalg as linalg
import scipy.special
from mlelement import *

class Well(Element):
    '''Well class; layers is list of aquifers in which well is screened
    Well class. Derived from Element
    Attributes provided on input:
    - modelParent: model it is added to
    - xw: x location of well
    - yw: y location of well
    - discharge: discharge of well
    - rw: radius of well
    - layers: row vector of aquifers in which well is screened (entered as list); doesn't have to be consecutive
    Attributes computed:
    - aquiferParent: AquiferData parent
    - rwsq: square of well radius
    - pylayers: list of Python layer numbers
    - NscreenedLayers: number of aquifers well is screend in
    - coef: coefficients
    All attributes from element.
    '''
    Rconvsq = 7 * 7
    def __init__(self,modelParent,xw,yw,Qw,rw,layers=0,label=None):
	Element.__init__(self,modelParent)
	self.xw = float(xw); self.yw = float(yw)
	self.discharge = float(Qw); self.rw = float(rw)
	self.layers = atleast_1d(layers)
	self.label = label
	self.type = 'well'
	self.setCoefs()
        self.modelParent.addElement(self)
    def __repr__(self):
	return 'Well xw,yw,Qw,rw,layers: ' + str((self.xw,self.yw,self.discharge,self.rw,list(self.layers)))
    def setCoefs(self):
	self.aquiferParent = self.modelParent.aq.findAquiferData(self.xw,self.yw)
        self.rwsq = self.rw**2;
        self.pylayers = self.layers  # fixed to base zero
	self.NscreenedLayers = len(self.layers)
	if self.NscreenedLayers == 1:  # Screened in only one layer, no unknown parameters
            self.parameters = array([[self.discharge]])
        else:
            self.parameters = zeros((self.NscreenedLayers,1),'d')
        self.coef = ones((self.NscreenedLayers,self.aquiferParent.Naquifers),'d')
	if self.aquiferParent.Naquifers > 1:   # Multiple aquifers, must compute coefficients
            if self.aquiferParent.type == self.aquiferParent.conf:
                for i in range(self.NscreenedLayers):
                    pylayer = self.pylayers[i]
                    ramat = self.aquiferParent.eigvec[:,1:]  # All eigenvectors, but not first normalized transmissivity
                    ramat = vstack(( ramat[0:pylayer,:], ramat[pylayer+1:,:] ))  # Remove row pylayer
                    rb = self.aquiferParent.eigvec[:,0] / (2.0*pi)   # Store Tn vector in rb
                    rb = hstack(( rb[0:pylayer], rb[pylayer+1:] )) # solve takes rowvector
                    self.coef[i,1:self.aquiferParent.Naquifers] = linalg.solve(ramat,rb)
            elif self.aquiferParent.type == self.aquiferParent.semi:
                for i in range(self.NscreenedLayers):
                    pylayer = self.pylayers[i]
                    ramat = self.aquiferParent.eigvec
                    rb = zeros(self.aquiferParent.Naquifers,'d')
                    rb[pylayer] = - 1.0 / (2.0 * pi)
                    self.coef[i,:] = linalg.solve(ramat,rb)
        if self.aquiferParent.Naquifers == 1 and self.aquiferParent.type == self.aquiferParent.semi:
            self.coef[0,0] = -1.0 / (2.0 * pi)
        self.paramxcoef = sum( self.parameters * self.coef, 0 )  # Parameters times coefficients
#        if self.NscreenedLayers == 1:
#            self.paramxcoef = self.parameters[0,0] * self.coef[0,:]
    def potentialInfluence(self,aq,x,y,z=0,t=0):
        # Makes use of scipy.special.K0 to compute Bessel function
        rv = zeros( (self.NscreenedLayers,aq.Naquifers), 'd' )
        if self.aquiferParent.type == self.aquiferParent.conf:
            if aq.type == aq.conf:
                rsq = (x-self.xw)**2 + (y-self.yw)**2
                rsq = max( rsq, self.rwsq )
                rv[:,0] = log(rsq/self.rwsq)/(4*pi)
                if self.aquiferParent == aq:            # Compute Bessel functions
                    for i in range(aq.Naquifers - 1):
                        rsqolam = rsq / aq.lab[i]**2
                        if rsqolam <= self.Rconvsq:
                            rolam = sqrt(rsqolam)
                            rv[:,i+1] = scipy.special.k0(rolam)
                    rv = self.coef * rv     # Only multiply with coefficients if Bessel part is computed, cause in fakesemi there are more aquifers, and first coef is zero
        elif self.aquiferParent.type == self.aquiferParent.semi:
            if self.aquiferParent == aq:                # In same semi-confined aquifer
                rsq = (x-self.xw)**2 + (y-self.yw)**2
                rsq = max( rsq, self.rwsq )
                for i in range(aq.Naquifers):
                    rsqolam = rsq / aq.lab[i]**2
                    if rsqolam <= self.Rconvsq:
                        rolam = sqrt(rsqolam)
                        rv[:,i] = scipy.special.k0(rolam)
                rv = self.coef * rv
                # Else we return zeros
        return rv
    def dischargeInfluence(self,aq,x,y):
        rvx = zeros((self.NscreenedLayers,aq.Naquifers),'d'); rvy = zeros((self.NscreenedLayers,aq.Naquifers),'d')
        if self.aquiferParent.type == self.aquiferParent.conf:
            if aq.type == aq.conf:
                rsq = (x-self.xw)**2 + (y-self.yw)**2
                xminxw = x-self.xw ; yminyw = y-self.yw
                if rsq < self.rwsq:
                    rsq = self.rwsq
                    xminxw = self.rw; yminyw = 0.0
                rvx[:,0] = -1.0/(2.0*pi) * xminxw / rsq
                rvy[:,0] = -1.0/(2.0*pi) * yminyw / rsq 
                if aq == self.aquiferParent:            # Compute Bessel functions
                    r = sqrt(rsq)
                    for i in range(aq.Naquifers - 1):
                        rsqolam = rsq / aq.lab[i]**2
                        if rsqolam <= self.Rconvsq:
                            rolam = sqrt(rsqolam)
                            kone = scipy.special.k1(rolam)
                            rvx[:,i+1] = kone * xminxw / (r * aq.lab[i])
                            rvy[:,i+1] = kone * yminyw / (r * aq.lab[i])
                    rvx = self.coef * rvx; rvy = self.coef * rvy  # Only multiply with coefficients if Bessel part is computed, cause in fakesemi there are more aquifers
##        elif self.aquiferParent.type == self.aquiferParent.semi:
##            if self.aquiferParent == aq:                # In same semi-confined aquifer
##                rsq = (x-self.xw)**2 + (y-self.yw)**2
##                xminxw = x-self.xw ; yminyw = y-self.yw
##                if rsq < self.rwsq:
##                    rsq = self.rwsq
##                    xminxw = self.rw; yminyw = 0.0
##        	r = sqrt(rsq)
##                for i in range(aq.Naquifers):
##                    rsqolam = rsq / aq.lab[i]**2
##                    if rsqolam <= self.Rconvsq:
##                        rolam = sqrt(rsqolam)
##                        kone = scipy.special.k1(rolam)
##                        rvx[i] = kone * xminxw / (r * aq.lab[i])
##                        rvy[i] = kone * yminyw / (r * aq.lab[i])
##                    rvx = self.coef * rvx; rvy = self.coef * rvy
        return [rvx,rvy]

    def potentialCollection(self,potsum,potadd,elementList,aq,x,y):

        # Not tested for semi-confined aquifer (improvement can and should be made there; no zeros, for example)
        # Could be nicer in evaluation of bessel functions (all at the same time), but may be troublesome to figure3
        # out which ones are far away
        
        for el in elementList:

            potadd[:] = 0.0
            
            if el.aquiferParent.type == el.aquiferParent.conf:
                if aq.type == aq.conf:
                    rsq = (x-el.xw)*(x-el.xw) + (y-el.yw)*(y-el.yw)
                    if rsq < el.rwsq: rsq = el.rwsq
                    potadd[0] = log(rsq/el.rwsq)/(4*pi)
                    if el.aquiferParent == aq:            # Compute Bessel functions
                        for i in range(aq.Naquifers - 1):
                            rsqolam = rsq / aq.lab[i]**2
                            if rsqolam <= el.Rconvsq:
                                rolam = sqrt(rsqolam)
                                potadd[i+1] = scipy.special.k0(rolam)
                        potadd = el.paramxcoef * potadd
                    else:
                        potadd[0] = el.paramxcoef[0] * potadd[0]
            elif el.aquiferParent.type == el.aquiferParent.semi:
                if el.aquiferParent == aq:                # In same semi-confined aquifer
                    rsq = (x-el.xw)*(x-el.xw) + (y-el.yw)*(y-el.yw)
                    if rsq < el.rwsq: rsq = el.rwsq
                    for i in range(aq.Naquifers):
                        rsqolam = rsq / aq.lab[i]**2
                        if rsqolam <= el.Rconvsq:
                            rolam = sqrt(rsqolam)
                            potadd[i] = scipy.special.k0(rolam)
                    potadd = el.paramxcoef * potadd
                    # Else we return zeros

            potsum = potsum + potadd
        return potsum

    def dischargeCollection(self,dissum,disadd,elementList,aq,x,y):
        # Not tested in other aquifers than the one it is screened in
        for el in elementList:
            disadd[:] = 0.0

            if el.aquiferParent.type == el.aquiferParent.conf:
                if aq.type == aq.conf:
                    rsq = (x-el.xw)**2 + (y-el.yw)**2
                    xminxw = x-el.xw ; yminyw = y-el.yw
                    if rsq < el.rwsq:
                        rsq = el.rwsq
                        xminxw = el.rw; yminyw = 0.0
                    xyminxwyw = array([xminxw,yminyw])
                    disadd[:,0] = el.paramxcoef[0] * (-1.0)/(2.0*pi) * xyminxwyw / rsq
                    if aq == el.aquiferParent:            # Compute Bessel functions
                        r = sqrt(rsq)
                        for i in range(aq.Naquifers - 1):
                            rsqolam = rsq / aq.lab[i]**2
                            if rsqolam <= el.Rconvsq:
                                rolam = sqrt(rsqolam)
                                kone = scipy.special.k1(rolam)
                                disadd[:,i+1] = el.paramxcoef[i+1] * kone * xyminxwyw / (r * aq.lab[i])
            elif el.aquiferParent.type == el.aquiferParent.semi:
                if el.aquiferParent == aq:                # In same semi-confined aquifer
                    rsq = (x-el.xw)**2 + (y-el.yw)**2
                    xminxw = x-el.xw ; yminyw = y-el.yw
                    if rsq < el.rwsq:
                        rsq = el.rwsq
                        xminxw = el.rw; yminyw = 0.0
                    r = sqrt(rsq)
                    for i in range(aq.Naquifers):
                        rsqolam = rsq / aq.lab[i]**2
                        if rsqolam <= el.Rconvsq:
                            rolam = sqrt(rsqolam)
                            kone = scipy.special.k1(rolam)
                            disadd[0,i] = kone * xminxw / (r * aq.lab[i])
                            disadd[1,i] = kone * yminyw / (r * aq.lab[i])
                    disadd = el.paramxcoef * disadd
            dissum = dissum + disadd
        return dissum
    def totalDischargeInfluence(self,aq,pylayer,x,y):
        rv = ones(self.NscreenedLayers,'d')
        return rv    
    def layout(self):
        return [ 1,[self.xw],[self.yw] ]
    def getMatrixRows(self,elementList):
        rows = []
        if self.NscreenedLayers > 1:
            rowlast = zeros(0,'d')
            lastpylayer = self.pylayers[self.NscreenedLayers-1]; Tlast = self.aquiferParent.T[lastpylayer]
            for e in elementList:
                rowpart = e.getMatrixCoefficients(self.aquiferParent,lastpylayer,self.xw,self.yw,\
                                                  lambda el,aq,layer,x,y:el.potentialInfluenceInLayer(aq,layer,x,y))
                rowlast = hstack(( rowlast, rowpart ))
            potlast = self.modelParent.potentialInLayer(self.aquiferParent,lastpylayer,self.xw,self.yw)
            for i in range(self.NscreenedLayers - 1):                
                row = zeros(0,'d'); T = self.aquiferParent.T[self.pylayers[i]]
                for e in elementList:
                    rowpart = e.getMatrixCoefficients(self.aquiferParent,self.pylayers[i],self.xw,self.yw,\
                                                      lambda el,aq,layer,x,y:el.potentialInfluenceInLayer(aq,layer,x,y))
                    row = hstack(( row, rowpart ))
                row = Tlast * row - T * rowlast
                row = hstack(( row, T * potlast - Tlast *\
                    self.modelParent.potentialInLayer(self.aquiferParent,self.pylayers[i],self.xw,self.yw) ))
                rows = rows + [row.tolist()]
            row = zeros(0,'d')
            for e in elementList:
                if e == self:
                    row = hstack(( row, ones(self.NscreenedLayers,'d') ))
                else:
                    row = hstack(( row, e.getMatrixCoefficients(self.aqdum,self.ldum,self.xdum,self.ydum,\
                                            lambda el,aqdum,ldum,xdum,ydum:el.zeroFunction(aqdum,ldum,xdum,ydum)) ))
            row = hstack(( row, self.discharge - sum(self.parameters) ))
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
            print 'Well at '+str((self.xw,self.yw))+' Layer: '+str(self.layers[0])+\
                  ' Discharge '+str(self.parameters[0,0])+' has no unknown parameters'
        else:
            print 'Well at '+str((self.xw,self.yw))
            print 'Specified discharge: '+str(self.discharge)+' Computed discharge: '+str(sum(self.parameters))
            for i in range(self.NscreenedLayers):
                print 'Layer '+str(self.layers[i])+' Head: '+str(self.modelParent.head(self.layers[i],self.xw,self.yw))+\
                      ' Discharge: '+str(self.parameters[i,0])
        return None
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):
        changed = 0; stop = 0; xyznew = 0.0
        if ( pyLayer == self.pylayers ).any():  # In layer that well is screened. This syntax is allowed! I checked
            # If point 1 within step of well and discharge of well positive, then terminate trace at well
            if ( ( self.parameters[0,0] > 0 and idir > 0 ) or ( self.parameters[0,0] < 0 and idir < 0 ) ):
                dissq = ( xyz1[0] - self.xw )*( xyz1[0] - self.xw ) + ( xyz1[1] - self.yw )*( xyz1[1] - self.yw )
                if dissq < (step+self.rw)*(step+self.rw):
                    horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
                    horstepnew = sqrt( (xyz1[0]-self.xw)**2 + (xyz1[1]-self.yw)**2 )
                    znew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )
                    if self.aquiferParent.inWhichPyLayer(znew) != pyLayer: znew = xyz1[2]  # If previous calculation gives point in other aquifer, take old point
                    xyznew = array([ self.xw, self.yw, znew ])  
                    changed = 1; stop = 1
        return [changed, stop, xyznew]
    def distanceSquaredToElement(self,x,y):
        '''Returns distance squared to element. Used for deciding tracing step'''
        dissq = ( x - self.xw )*( x - self.xw ) + ( y - self.yw )*( y - self.yw )
        dissq = dissq - sqrt( dissq ) * 2.0 * self.rw + self.rwsq # Believe me, this matters
        return dissq
    def PYpotentialInfluence(self,aq,x,y,z=0,t=0):
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
        NBes = 8; 
	rsq = (x-self.xw)**2 + (y-self.yw)**2
	rsq = max( rsq, self.rwsq )
	potInf = [ log(rsq/self.rwsq)/(4.0*pi) ]
	if aq == self.aquiferParent:
            for i in range(self.aquiferParent.Naquifers - 1):
                rsqolam = rsq / self.aquiferParent.lab[i]**2
                if rsqolam > self.Rconvsq:
                    potInf = potInf + [0.0]
                    continue
                logrsqolam = log(rsqolam)
                pot = ac[0] * logrsqolam + bc[0]
                rsqdum = 1.0
                for j in range(1,NBes+1):
                    rsqdum = rsqdum * rsqolam
                    pot = pot + ( ac[j] * logrsqolam + bc[j] ) * rsqdum
                potInf = potInf + [pot]
            rv = []
            for coef in self.coef:
                rv = rv + [list(potInf * coef)]  # Returns list of lists
        else:
            pot = potInf + list( zeros(aq.Naquifers-1,'d') )
            rv = []
            for i in range(self.NscreenedLayers):
                rv = rv + [pot]
        return rv


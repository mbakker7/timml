from numpy import *
import numpy.linalg as linalg
import scipy.special
from mlelement import *

# Not fully adjusted for semi-confined conditions
class CircAreaSink(Element):
    def __init__(self,modelParent,xp,yp,Rp,infil,layer=0,label=None):
	Element.__init__(self,modelParent)
        self.xp = float(xp); self.yp = float(yp); self.Rp = float(Rp); self.infil = float(infil)
        self.layer = layer
        if iterable(layer):
	    assert len(self.layer)==1, "TimML Input error: CircAreasink can only be defined in one layer"
	    self.layer = self.layer[0]        
	self.pylayer = self.layer  # fixed base to zero
	self.label = label
        self.type = 'circareasink'
        self.setCoefs()
        self.modelParent.addElement(self)
    def __repr__(self):
	return 'CircAreaSink xp,yp,Rp,infil,layer: ' + str((self.xp,self.yp,self.Rp,self.infil,self.pylayer))
    def setCoefs(self):
        self.Rpsq = self.Rp**2
        self.parameters = array([[self.infil]])
	self.aquiferParent = self.modelParent.aq.findAquiferData(self.xp,self.yp)
        assert (self.aquiferParent.type == self.aquiferParent.conf or self.aquiferParent.Naquifers ==1), \
               "Input error: CircAreaSink centered in semi-confined system with more than one layer"
	self.coef = ones((1,self.aquiferParent.Naquifers),'d')
	self.Rlarge = 500.0  # If R/lab > Rlarge, then we use asymptotic approximation to compute potential
        if self.aquiferParent.Naquifers > 1:
            self.isLargeCircle = zeros(self.aquiferParent.Naquifers-1)
            self.I1Rolab = zeros(self.aquiferParent.Naquifers-1,'d')
            self.K1Rolab = zeros(self.aquiferParent.Naquifers-1,'d')
            for i in range(self.aquiferParent.Naquifers-1):
                if self.Rp/self.aquiferParent.lab[i] > self.Rlarge:
                    self.isLargeCircle[i] = 1
                else:                    
                    self.I1Rolab[i] = scipy.special.i1( self.Rp / self.aquiferParent.lab[i] )
                    self.K1Rolab[i] = scipy.special.k1( self.Rp / self.aquiferParent.lab[i] )
            if self.aquiferParent.type == self.aquiferParent.conf:
                ramat = list( transpose(self.aquiferParent.eigvec) )  # Take transpose and make into a list
                rb = -ramat.pop(0)   # Remove Tn vector from ramat and store in rb with minus sign
                ramat = transpose( array(ramat) ) # Make into matrix and transpose back
                ramat = ramat / self.aquiferParent.lab / self.Rp # Divide vectors by lambda_m and Rp
                #ramat = take( ramat,range(0,self.pylayer)+range(self.pylayer+1,self.aquiferParent.Naquifers) ) # Remove screened layer
                ramat = vstack(( ramat[0:self.pylayer,:], ramat[self.pylayer+1:,:] ))  # Remove row pylayer
                rb = take( rb,range(0,self.pylayer)+range(self.pylayer+1,self.aquiferParent.Naquifers) )  # solve takes rowvector, so no need to convert
                self.coef[0,1:] = linalg.solve(ramat,rb)
            elif self.aquiferParent.type == self.aquiferParent.semi:
                # do something different
                print 'center of CircAreaSink in semi-confined aquifer'
        if self.aquiferParent.type == self.aquiferParent.semi:  # For now only with one layer
            self.I1Rolab = array(1,'d'); self.K1Rolab = array(1,'d')
            self.I1Rolab[0] = scipy.special.i1( self.Rp / self.aquiferParent.lab[0] )
            self.K1Rolab[0] = scipy.special.k1( self.Rp / self.aquiferParent.lab[0] )
    def potentialInfluence(self,aq,x,y,z=0,t=0):
        # In and Kn approximations from (9.7.1) and (9.7.2), Abramowitz and Stegun
        # K1(R)*I0(r)
        # fk1i0 = sqrt( 1./(4.0*r*R) ) * exp( r-R ) * (1. + 3./(8*R) - 15./(128*R**2) + 315./(3072*R**3) ) \
        #                                           * (1. + 1./(8*r) +  9./(128*r**2) + 225./(3072*r**3) )
        # I1(R)*K0(r)
        # fi1k0 = sqrt( 1./(4.0*r*R) ) * exp( R-r ) * (1. - 3./(8*R) - 15./(128*R**2) - 315./(3072*R**3) ) \
        #                                           * (1. - 1./(8*r) +  9./(128*r**2) - 225./(3072*r**3) )
        rv = zeros( (1,aq.Naquifers), 'd' )
        if self.aquiferParent.type == self.aquiferParent.conf:
            if aq.type == aq.conf:
                rsq = (x-self.xp)**2 + (y-self.yp)**2
                r = sqrt(rsq)
                if r <= self.Rp:
                    rv[0,0] = (self.Rp**2 - rsq)/4.0
                else:
                    rv[0,0] = -self.Rp**2 / 2.0 * log( r/self.Rp )
                if self.aquiferParent == aq:
                    for i in range(self.aquiferParent.Naquifers - 1):
                        rolab = r / self.aquiferParent.lab[i]; Rolab = self.Rp / self.aquiferParent.lab[i]
                        if r <= self.Rp:
                            if self.isLargeCircle[i] == 0:
                                pot = self.aquiferParent.lab[i] /self.Rp - self.K1Rolab[i] * scipy.special.i0(rolab)
                            else:
                                pot = self.aquiferParent.lab[i] /self.Rp
                                if (Rolab - rolab) < 10.0: # zero after 10 lambda
                                    pot = pot - \
                                        sqrt( 1./(4.0 * rolab * Rolab) ) * exp( rolab - Rolab ) * \
                                        (1. + 3./(8*Rolab) - 15./(128*Rolab**2) + 315./(3072*Rolab**3) ) * \
                                        (1. + 1./(8*rolab) +  9./(128*rolab**2) + 225./(3072*rolab**3) )
                        else:
                            if self.isLargeCircle[i] == 0:                    
                                pot = self.I1Rolab[i] * scipy.special.k0(rolab)
                            else:
                                pot = 0.0
                                if (rolab - Rolab) < 10.0: # zero after 10 lambda
                                    pot = sqrt( 1./(4.0 * rolab * Rolab) ) * exp( Rolab - rolab ) * \
                                        (1. - 3./(8*Rolab) - 15./(128*Rolab**2) - 315./(3072*Rolab**3) ) * \
                                        (1. - 1./(8*rolab) +  9./(128*rolab**2) - 225./(3072*rolab**3) )
                        rv[0,i+1] = pot
                    rv = self.coef * rv
                else:  # Need to add difference in heads between layers
                    if r <= self.Rp:
                        rv[0,:] = rv[0,:] + aq.constVecPot
        elif self.aquiferParent.type == self.aquiferParent.semi:  # For now only works in 1 semi-confined aquifer
            if self.aquiferParent == aq:                # In same semi-confined aquifer
                rsq = (x-self.xp)**2 + (y-self.yp)**2
                r = sqrt(rsq)
                lab = self.aquiferParent.lab[0]
                if r <= self.Rp:  # Inside pond
                    A = -lab * self.Rp * self.K1Rolab[0]
                    rv[0,0] = lab**2 + A * scipy.special.i0( r/lab )
                else:
                    B = lab * self.Rp * self.I1Rolab[0]
                    rv[0,0] = B * scipy.special.k0( r/lab )
                rv = self.coef * rv
        return rv
    def potentialContribution(self,aq,x,y):
        return self.parameters[0,0] * self.potentialInfluence(aq,x,y)[0,:]
    def potentialCollection(self,potsum,potadd,elementList,aq,x,y):
        for e in elementList:
            potsum = potsum + e.parameters[0,0] * e.potentialInfluence(aq,x,y)[0,:]
        return potsum
    def dischargeInfluence(self,aq,x,y,z=0,t=0):
        rvx = zeros((1,aq.Naquifers),'d'); rvy = zeros((1,aq.Naquifers),'d')
	rsq = (x-self.xp)**2 + (y-self.yp)**2
	r = sqrt(rsq)
        if r < 1e-6:  # Else evaluation at center blows up
            r = 1e-6
            rsq = r**2
	if r <= self.Rp:
            rvx[0,0] = (x - self.xp) / 2.0
            rvy[0,0] = (y - self.yp) / 2.0
        else:
            rvx[0,0] = self.Rp**2 * (x-self.xp) / (2.0 * rsq)
            rvy[0,0] = self.Rp**2 * (y-self.yp) / (2.0 * rsq)
        if self.aquiferParent == aq:
            for i in range(self.aquiferParent.Naquifers - 1):
                rolab = r / self.aquiferParent.lab[i]; Rolab = self.Rp / self.aquiferParent.lab[i]
                if r <= self.Rp:
                    if self.isLargeCircle[i] == 0:
                        dis = self.K1Rolab[i] * scipy.special.i1(rolab) / self.aquiferParent.lab[i]
                        rvx[0,i+1] = dis * (x-self.xp) / r
                        rvy[0,i+1] = dis * (y-self.yp) / r
                    else:
                        if (Rolab - rolab) < 10.0: # zero after 10 lambda
                            K1RI1r = sqrt( 1./(4.0 * rolab * Rolab) ) * exp( rolab - Rolab ) * \
                                (1. + 3./(8*Rolab) - 15./(128*Rolab**2) + 315./(3072*Rolab**3) ) * \
                                (1. - 3./(8*rolab) - 15./(128*rolab**2) - 315/(3072*rolab**3) )
                            rvx[0,i+1] = K1RI1r * (x-self.xp) / ( r * self.aquiferParent.lab[i] )
                            rvy[0,i+1] = K1RI1r * (y-self.yp) / ( r * self.aquiferParent.lab[i] )
                else:
                    if self.isLargeCircle[i] == 0:                    
                        dis = self.I1Rolab[i] * scipy.special.k1(rolab) / self.aquiferParent.lab[i]
                        rvx[0,i+1] = dis * (x-self.xp) / r
                        rvy[0,i+1] = dis * (y-self.yp) / r
                    else:
                        if (rolab - Rolab) < 10.0: # zero after 10 lambda
                            I1RK1r = sqrt( 1./(4.0 * rolab * Rolab) ) * exp( Rolab - rolab ) * \
                                (1. - 3./(8*Rolab) - 15./(128*Rolab**2) - 315./(3072*Rolab**3) ) * \
                                (1. + 3./(8*rolab) - 15./(128*rolab**2) + 315./(3072*rolab**3) )
                            rvx[0,i+1] = I1RK1r * (x-self.xp) / ( r * self.aquiferParent.lab[i] )
                            rvy[0,i+1] = I1RK1r * (y-self.yp) / ( r * self.aquiferParent.lab[i] )
            # Only if aq == aquiferParent, otherwise not
            rvx = self.coef * rvx
            rvy = self.coef * rvy
        return [rvx,rvy]
    def dischargeCollection(self,dissum,disadd,elementList,aq,x,y):
        for e in elementList:
            disInf = e.dischargeInfluence(aq,x,y)
            dissum[0,:] = dissum[0,:] + e.parameters[0,0] * disInf[0][0,:]
            dissum[1,:] = dissum[1,:] + e.parameters[0,0] * disInf[1][0,:]
        return dissum
    def dischargeContribution(self,aq,x,y):
        '''Returns matrix with two rowvectors of dischargeContributions Qx and Qy'''
        #  Should be improved for speed. Take eigvec out here!
        disInf = self.dischargeInfluence(aq,x,y)
        disxInf = self.parameters[0,0] * disInf[0][0,:]
        disyInf = self.parameters[0,0] * disInf[1][0,:]
        return array([disxInf,disyInf])
    def totalDischargeInfluence(self,aq,pylayer,x,y):
        rv = - ones(1,'d') * self.Rpsq * pi
        return rv  
    def layout(self):
        rv = [ 100,[],[] ]
        theta = arange(0,2*pi+pi/50,pi/50)
        rv[1] = list( self.xp + self.Rp * cos(theta) )
        rv[2] = list( self.yp + self.Rp * sin(theta) )
        return rv
    def check(self):
        print 'CircAreaSink at '+str((self.xp,self.yp))+' Layer: '+str(self.layer)+\
                  ' Infiltration '+str(self.parameters[0,0])+' has no unknown parameters'
        return None
    def qzTop(self,x,y):
        '''Returns qz at top of aquifer at (x,y) due to element
        Only non-zero for areal recharge'''
        rsq = (x-self.xp)**2 + (y-self.yp)**2
        qz = 0.0
        if rsq <= self.Rpsq:
            qz = - self.parameters[0,0]  # Positive parameter is infiltration (qz<0)
        return qz
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):
        changed = 0; stop = 0; xyznew = zeros(3,'d')
        r1sq = ( xyz1[0] - self.xp )**2 + ( xyz1[1] - self.yp ) **2
        r2sq = ( xyz2[0] - self.xp )**2 + ( xyz2[1] - self.yp ) **2
        if ( r1sq < self.Rpsq and r2sq > self.Rpsq  ) or ( r1sq > self.Rpsq and r2sq < self.Rpsq  ):
            (x1,y1) = xyz1[0:2]; (x2,y2) = xyz2[0:2]
            a = (x2-x1)**2 + (y2-y1)**2
            b = 2.0 * ( (x2-x1) * (x1-self.xp) + (y2-y1) * (y1-self.yp) )
            c = self.xp**2 + self.yp**2 + x1**2 + y1**2 - 2.0 * (self.xp*x1 +self.yp*y1) - self.Rpsq
            u1 = ( -b - sqrt(b**2 - 4.0*a*c) ) / (2.0*a)
            u2 = ( -b + sqrt(b**2 - 4.0*a*c) ) / (2.0*a)
            if u1 > 0:
                u = u1 * 1.000001 # Go just beyond circle
            else:
                u = u2 * 1.000001 # Go just beyond circle
            xn = x1 + u * (x2-x1); yn = y1 + u * (y2-y1)
            zn = xyz1[2] + u * ( xyz2[2] - xyz1[2] )
            xyznew = array([xn,yn,zn],'d')
            changed = 1
        return [changed, stop, xyznew]  # Didn't adjust time
    def distanceSquaredToElement(self,x,y):
        '''Returns distance squared to element. Used for deciding tracing step'''
        dis = sqrt( ( x - self.xp )**2 + ( y - self.yp ) **2 )
        dissq = ( dis - self.Rp ) **2
        return dissq

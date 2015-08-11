from numpy import *
from mlelement import *
from mllinesink import *  # Can go?l
from besselaes import *

class LineDoublet(Element):
    '''Linedoublets for inhomogeneity boundary; single aquifer only, can
    only connect to background aquifer'''
    # Prototype version for variable order
    oneOverTwoPii = 1.0 / ( 2.0 * pi * complex(0.0,1.0) )
    tiny = 1e-8; tinyimag = tiny * complex(0.0,1.0)  # 1e-6 is needed, as I got problems with the gnu compiler with 1e-8
    def __init__(self,modelParent,x1,y1,x2,y2,aqin,aqout,Ndegree=2,overspec=1,addToModel=1):
        Element.__init__(self,modelParent)
        self.x1 = float(x1); self.y1 = float(y1); self.x2 = float(x2); self.y2 = float(y2);
        self.aqin = aqin; self.aqout = aqout
        assert Ndegree < 9, "TimML Input error: Maximum order is 8"
        self.Ndegree = Ndegree
        self.overspec = overspec # Overspec doesn't work yet.
        # self.ldLeft = None; self.ldRight = None  # Note: these have to be set appropriately. Use Makeinhomogeneity utility.
        # Note that ldRight is never used at this time. ldLeft not anymore either !
        self.type = 'linedoublet'
        self.setCoefs()
        if addToModel: self.modelParent.addElement(self)
    def __repr__(self):
	return 'LineDoublet z1,z2,strengths: ' + str((self.x1,self.y1,self.x2,self.y2,self.parameters))
    def setCoefs(self):
        self.aquiferParent = self.aqout
        self.z1 = complex(self.x1,self.y1); self.z2 = complex(self.x2,self.y2)
        self.xc = 0.5*(self.x1+self.x2); self.yc = 0.5*(self.y1+self.y2)
        self.L = abs(self.z2-self.z1); self.Lover2 = self.L / 2.0
        self.xmin = min(self.x1,self.x2); self.xmax = max(self.x1,self.x2)
        self.ymin = min(self.y1,self.y2); self.ymax = max(self.y1,self.y2)
        self.thetaNormOut = arctan2(self.y2-self.y1,self.x2-self.x1) - pi/2.0
        self.alpha = arctan2(self.y2-self.y1,self.x2-self.x1)
        self.cosalpha = cos(self.alpha); self.sinalpha = sin(self.alpha)
        self.parameters = zeros((self.Ndegree+1,1),'d')
#        zcp1 = self.z1 + .25 * (self.z2 - self.z1)
#        zcp2 = self.z1 + .75 * (self.z2 - self.z1)
        self.Ncp = int( ( self.Ndegree+1 ) * self.overspec )
        an = arange( pi / (2.0*self.Ncp), pi, pi/self.Ncp )
        Zcp = zeros( self.Ncp, 'D' )
        Zcp.real = cos(an)
        Zcp.imag = 1e-6  # control point just on inside
        zcp = Zcp * (self.z2 - self.z1) / 2.0 + 0.5 * (self.z1 + self.z2)
        self.xcpin = zcp.real; self.ycpin = zcp.imag
        Zcp.imag = -1e-6  # control point just on inside
        zcp = Zcp * (self.z2 - self.z1) / 2.0 + 0.5 * (self.z1 + self.z2)
        self.xcpout = zcp.real; self.ycpout = zcp.imag
        #
        self.potInf = zeros(self.Ndegree+1,'d')
        self.disxInf = zeros(self.Ndegree+1,'d'); self.disyInf = zeros(self.Ndegree+1,'d')
##        assert self.aqout.Naquifers == 1 or self.aqin.Naquifers == 1, \
##               "TimML Input error:  either inside or outside aquifer must have one layer only"        
    def potentialInfluence(self,aq,x,y):
        zin = complex(x,y)
        bigZ = ( 2.0*zin - (self.z1 + self.z2) )/ (self.z2 - self.z1)
        if abs(bigZ.imag) < self.tiny and abs(bigZ.real) < 1.0-self.tiny:  # point is on boundary; this could be nicer by setting the log
                if self.aqin == aq:
                    bigZ = complex(bigZ.real,self.tiny)
                elif self.aqout == aq:
                    bigZ = complex(bigZ.real,-self.tiny)
                zin = ( (self.z2-self.z1)*bigZ + (self.z1+self.z2) ) / 2.0
        bigZmin1 = bigZ - 1.0; bigZplus1 = bigZ + 1.0
##        if abs(bigZmin1) < self.tiny:
##            # point is at right corner; move back a little along element; may fail if corner is very pointy
##            # This needs to be moved by the PolygonInhom class; after all, that class decided it was inside.
##            if self.aqin == aq:
##                zin = self.z2 + complex(4.0,-1.0) * self.tiny * (self.z1 - self.z2)
##            else:
##                zin = self.z2 + complex(4.0,1.0) * self.tiny * (self.z1 - self.z2)
##        elif abs(bigZplus1) < self.tiny:  # point is at left corner, move up a little along ldLeft
##            if self.aqin == aq:
##                zin = self.ldLeft.z2 + complex(4.0,-1.0) * self.tiny * (self.ldLeft.z1 - self.ldLeft.z2)
##            else:
##                zin = self.ldLeft.z2 + complex(4.0,1.0) * self.tiny * (self.ldLeft.z1 - self.ldLeft.z2)
        if abs(bigZmin1) < self.tiny or abs(bigZplus1) < self.tiny: # gotta be inside; move point
            # zin = self.aqin.movePoint(zin)
            # doesn't have to be inside when inhoms butt-up
            zin = aq.movePoint(zin)
        #  Only works for line-doublet in one layer, so that is hardcoded (the zero lists)
        rv = zeros((self.Ndegree+1,aq.Naquifers),'d')
        x = zin.real; y = zin.imag
        potlapldho(x,y,self.x1,self.y1,self.x2,self.y2,self.Ndegree,self.potInf)
        rv[:,0] = self.potInf
        return rv
    def potentialCollection(self,potsum,potadd,elementList,aq,xin,yin):
        for el in elementList:
            # Check first if distance to element is small; this saves a little time
            x = xin; y = yin
            X = (x-el.xc) * el.cosalpha + (y-el.yc) * el.sinalpha
            #print 'X ',X
            if X < -el.Lover2:
                #print 'x,y ',x,y
                #print 'x1,y1 ',el.x1,el.y1
                dissq = ( x - el.x1 )*( x - el.x1 ) + ( y - el.y1 )*( y - el.y1 )
            elif X > el.Lover2:
                dissq = ( x - el.x2 )*( x - el.x2 ) + ( y - el.y2 )*( y - el.y2 )
                #print 'x,y ',x,y
                #print 'x2,y2 ',el.x2,el.y2
            else:
                Y = -(x-el.xc) * el.sinalpha + (y-el.yc) * el.cosalpha
                #print 'Y ',Y
                dissq = Y * Y
            #print 'X,dissq ',X,dissq
            if dissq < 0.0001:
                zin = complex(x,y)
                bigZ = ( 2.0*zin - (el.z1 + el.z2) )/ (el.z2 - el.z1)
                if abs(bigZ.imag) < el.tiny and abs(bigZ.real) < 1.0-el.tiny:  # point is on boundary; this could be nicer by setting the log
                        if el.aqin == aq:
                            bigZ = complex(bigZ.real,el.tiny)
                        elif el.aqout == aq:
                            bigZ = complex(bigZ.real,-el.tiny)
                        zin = ( (el.z2-el.z1)*bigZ + (el.z1+el.z2) ) / 2.0
                bigZmin1 = bigZ - 1.0; bigZplus1 = bigZ + 1.0
##                if abs(bigZmin1) < el.tiny:  # point is at right corner; move back a little along element; may fail if corner is very pointy
##                    #print 'hello 1'
##                    if el.aqin == aq:
##                        zin = el.z2 + complex(4.0,-1.0) * el.tiny * (el.z1 - el.z2)
##                    else:
##                        zin = el.z2 + complex(4.0,1.0) * el.tiny * (el.z1 - el.z2)
##                elif abs(bigZplus1) < el.tiny:  # point is at left corner, move up a little along ldLeft
##                    #print 'hello 2'
##                    if el.aqin == aq:
##                        zin = el.ldLeft.z2 + complex(4.0,-1.0) * el.tiny * (el.ldLeft.z1 - el.ldLeft.z2)
##                    else:
##                        zin = el.ldLeft.z2 + complex(4.0,1.0) * el.tiny * (el.ldLeft.z1 - el.ldLeft.z2)
                #  Only works for line-doublet in one layer, so that is hardcoded (the zero lists)
                if abs(bigZmin1) < self.tiny or abs(bigZplus1) < self.tiny: # gotta be inside; move point
                    # zin = el.aqin.movePoint(zin)
                    # doesn't have to be inside when inhoms butt-up
                    zin = aq.movePoint(zin)
                x = zin.real; y = zin.imag
            potlapldho(x,y,el.x1,el.y1,el.x2,el.y2,el.Ndegree,el.potInf)
            for i in range(el.Ndegree+1):
                potsum[0] = potsum[0] + el.parameters[i,0] * el.potInf[i]
        return potsum
    def dischargeInfluence(self,aq,x,y):
        zin = complex(x,y)
        rv = []
        bigZ = ( 2.0*zin - (self.z1 + self.z2) )/ (self.z2 - self.z1)
        if abs(bigZ.imag) < self.tiny and abs(bigZ.real) < 1.0-self.tiny:  # point is on boundary; this could be nicer by setting the log
                if self.aqin == aq:
                    bigZ = complex(bigZ.real,self.tiny)
                elif self.aqout == aq:
                    bigZ = complex(bigZ.real,-self.tiny)
                zin = ( (self.z2-self.z1)*bigZ + (self.z1+self.z2) ) / 2.0
        bigZmin1 = bigZ - 1.0; bigZplus1 = bigZ + 1.0
##        if abs(bigZmin1) < self.tiny:  # point is at right corner; move back a little along element; may fail if corner is very pointy
##            if self.aqin == aq:
##                zin = self.z2 + complex(4.0,-1.0) * self.tiny * (self.z1 - self.z2)
##            else:
##                zin = self.z2 + complex(4.0,1.0) * self.tiny * (self.z1 - self.z2)
##        elif abs(bigZplus1) < self.tiny:  # point is at left corner, move up a little along ldLeft
##            if self.aqin == aq:
##                zin = self.ldLeft.z2 + complex(4.0,-1.0) * self.tiny * (self.ldLeft.z1 - self.ldLeft.z2)
##            else:
##                zin = self.ldLeft.z2 + complex(4.0,1.0) * self.tiny * (self.ldLeft.z1 - self.ldLeft.z2)
        if abs(bigZmin1) < self.tiny or abs(bigZplus1) < self.tiny: # gotta be inside; move point
            zin = self.aqin.movePoint(zin)
        #  Only works for line-doublet in one layer, so that is hardcoded (the zero lists)
        rvx = zeros((self.Ndegree+1,aq.Naquifers),'d')
        rvy = zeros((self.Ndegree+1,aq.Naquifers),'d')
        x = zin.real; y = zin.imag
##        lab = zeros(1,'d')
##        disx = zeros(1,'d'); disy = zeros(1,'d')
##        disbesldho(x,y,self.x1,self.y1,self.x2,self.y2,1,lab,0,disx,disy)
##        rvx[0,0] = disx[0]; rvy[0,0] = disy[0]
##        disbesldho(x,y,self.x1,self.y1,self.x2,self.y2,1,lab,1,disx,disy)
##        rvx[1,0] = disx[0]; rvy[1,0] = disy[0]
##        disbesldho(x,y,self.x1,self.y1,self.x2,self.y2,1,lab,2,disx,disy)
##        rvx[2,0] = disx[0]; rvy[2,0] = disy[0]
        dislapldho(x,y,self.x1,self.y1,self.x2,self.y2,self.Ndegree,self.disxInf,self.disyInf)
        rvx[:,0] = self.disxInf[:]; rvy[:,0] = self.disyInf[:]
        return [rvx,rvy]
    def dischargeCollection(self,dissum,disadd,elementList,aq,x,y):
        for el in elementList:
            # Check first if distance to element is small; this saves a little time
            X = (x-el.xc) * el.cosalpha + (y-el.yc) * el.sinalpha
            if X < -el.Lover2:
                dissq = ( x - el.x1 )*( x - el.x1 ) + ( y - el.y1 )*( y - el.y1 )
            elif X > el.Lover2:
                dissq = ( x - el.x2 )*( x - el.x2 ) + ( y - el.y2 )*( y - el.y2 )
            else:
                Y = -(x-el.xc) * el.sinalpha + (y-el.yc) * el.cosalpha
                dissq = Y * Y
            if dissq < 0.0001:
                zin = complex(x,y)
                bigZ = ( 2.0*zin - (el.z1 + el.z2) )/ (el.z2 - el.z1)
                if abs(bigZ.imag) < el.tiny and abs(bigZ.real) < 1.0-el.tiny:  # point is on boundary; this could be nicer by setting the log
                        if el.aqin == aq:
                            bigZ = complex(bigZ.real,el.tiny)
                        elif el.aqout == aq:
                            bigZ = complex(bigZ.real,-el.tiny)
                        zin = ( (el.z2-el.z1)*bigZ + (el.z1+el.z2) ) / 2.0
                bigZmin1 = bigZ - 1.0; bigZplus1 = bigZ + 1.0
##                if abs(bigZmin1) < el.tiny:  # point is at right corner; move back a little along element; may fail if corner is very pointy
##                    if el.aqin == aq:
##                        zin = el.z2 + complex(4.0,-1.0) * el.tiny * (el.z1 - el.z2)
##                    else:
##                        zin = el.z2 + complex(4.0,1.0) * el.tiny * (el.z1 - el.z2)
##                elif abs(bigZplus1) < el.tiny:  # point is at left corner, move up a little along ldLeft
##                    if el.aqin == aq:
##                        zin = el.ldLeft.z2 + complex(4.0,-1.0) * el.tiny * (el.ldLeft.z1 - el.ldLeft.z2)
##                    else:
##                        zin = el.ldLeft.z2 + complex(4.0,1.0) * el.tiny * (el.ldLeft.z1 - el.ldLeft.z2)
##                #  Only works for line-doublet in one layer, so that is hardcoded (the zero lists)
                if abs(bigZmin1) < self.tiny or abs(bigZplus1) < self.tiny: # gotta be inside; move point
                    zin = el.aqin.movePoint(zin)
                x = zin.real; y = zin.imag
            dislapldho(x,y,el.x1,el.y1,el.x2,el.y2,el.Ndegree,el.disxInf,el.disyInf)
            for i in range(el.Ndegree+1):
                dissum[0,0] = dissum[0,0] + el.parameters[i,0] * el.disxInf[i]
                dissum[1,0] = dissum[1,0] + el.parameters[i,0] * el.disyInf[i]     
        return dissum
    def getMatrixRows(self,elementList):
        rows=[]
        for i in range(self.Ncp):
            # Hardcoded for pylayer 0 
            pylayer = 0
            row = zeros(0,'d'); rowin = zeros(0,'d'); rowout = zeros(0,'d')
            for e in elementList:
                # Modification for fakesemi
                if not self.aqin.fakesemi:
                    rowinpart = e.getMatrixCoefficients(self.aqin,pylayer,self.xcpin[i],self.ycpin[i],\
                                                lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                else:
                    rowinpart = e.getMatrixCoefficients(self.aqin,pylayer+1,self.xcpin[i],self.ycpin[i],\
                                                lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                rowin = hstack(( rowin, rowinpart ))
                rowoutpart = e.getMatrixCoefficients(self.aqout,pylayer,self.xcpout[i],self.ycpout[i],\
                                                lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                rowout = hstack(( rowout, rowoutpart ))
                if not self.aqin.fakesemi:
                    row = self.aqout.T[pylayer] * array(rowin) - self.aqin.T[pylayer] * array(rowout)
                else:
                    row = self.aqout.T[pylayer] * array(rowin) - self.aqin.T[pylayer+1] * array(rowout)
            if not self.aqin.fakesemi:
                row = hstack(( row,
                        self.aqin.T[pylayer] * self.modelParent.potentialInLayer(self.aqout,pylayer,self.xcpout[i],self.ycpout[i]) -\
                        self.aqout.T[pylayer] * self.modelParent.potentialInLayer(self.aqin,pylayer,self.xcpin[i],self.ycpin[i]) ))
            else:
                row = hstack(( row,
                        self.aqin.T[pylayer+1] * self.modelParent.potentialInLayer(self.aqout,pylayer,self.xcpout[i],self.ycpout[i]) -\
                        self.aqout.T[pylayer] * self.modelParent.potentialInLayer(self.aqin,pylayer+1,self.xcpin[i],self.ycpin[i]) ))
            rows = rows + [row.tolist()]
        return rows
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        for i in range(self.Ndegree+1):
            self.parameters[i,0] = self.parameters[i,0] + xsol[icount]
            icount = icount + 1
        return icount
    def layout(self):
        return [0]
    def check(self):
        print 'LineDoublet from '+str(self.z1)+' to '+str(self.z2)+' Top Layer '
        for i in range(self.Ndegree+1):
            print 'Control point: '+str(i)+\
                  ' Head inside : '+str(self.modelParent.head(1,self.xcpin[i],self.ycpin[i],self.aqin))+\
                  ' Head outside: '+str(self.modelParent.head(1,self.xcpout[i],self.ycpout[i],self.aqout))
        return None
    def check_normalflow(self):
        print 'LineDoublet from '+str(self.z1)+' to '+str(self.z2)+' Top Layer '
        for i in range(self.Ndegree+1):
            print 'Control point: '+str(i)+\
                  ' Qnorm inside : '+str(self.modelParent.dischargeNormInLayer(self.aqin,0,self.xcpin[i],self.ycpin[i],self.thetaNormOut))+\
                  ' Qnorm outside: '+str(self.modelParent.dischargeNormInLayer(self.aqout,0,self.xcpout[i],self.ycpout[i],self.thetaNormOut))
        return None
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):
        # Not doing anything right now. Should really only do something if 1 aquifer both in and out. Else line-sink is doing the work
        changed = 0; stop = 0; xyznew = 0.0
##        [x1,y1,z1] = xyz1; [x2,y2,z2] = xyz2
##        if x1 < self.xmin - step or x1 > self.xmax + step or y1 < self.ymin - step or y1 > self.ymax + step:
##           return [changed, stop, xyznew, extrat]  # Away from line-doublet
##        z1 = complex(x1,y1); z2 = complex(x2,y2)
##        Z1 = (2.0*z1 - (self.z1+self.z2))/(self.z2-self.z1)
##        Z2 = (2.0*z2 - (self.z1+self.z2))/(self.z2-self.z1)
##        # If point 1 on one side and point 2 on other side, find out if segment intesects line-doublet
##        if (Z1.imag > 0 and Z2.imag < 0) or (Z1.imag < 0 and Z2.imag > 0):
##            Xintersec = Z1.real + (0.0-Z1.imag) / (Z2.imag-Z1.imag) * (Z2.real-Z1.real)
##            if abs(Xintersec) <= 1:  # There is an intersection point
##                changed = 1
##                Znew = complex(Xintersec,0.0)
##                znew = ( (self.z2-self.z1)*Znew + (self.z1+self.z2) ) / 2.0  # New complex coordinate
##                xnew = znew.real; ynew = znew.imag
##                horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
##                horstepnew = sqrt( (xyz1[0]-xnew)**2 + (xyz1[1]-ynew)**2 )
##                zvertnew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )  # Just to be confusing, this is z (vertical coordinate)                 
##                if Z1.imag < 0:  # old point on outside, so new point on inside
##                    dis = self.modelParent.dischargeNormBelowZ(self.aqout,xnew,ynew,zvertnew,self.thetaNormOut)
##                    zvertnew2 = self.modelParent.zForGivenNormalDischarge\
##                                (self.aqin,zvertnew,xnew,ynew,self.thetaNormOut,dis)
##                else:
##                    dis = self.modelParent.dischargeNormBelowZ(self.aqin,xnew,ynew,zvertnew,self.thetaNormOut)
##                    zvertnew2 = self.modelParent.zForGivenNormalDischarge\
##                                (self.aqout,zvertnew,xnew,ynew,self.thetaNormOut,dis)
##                xyznew[:] = [xnew,ynew,zvertnew2]  # Should modify z
        return [changed, stop, xyznew]  # Didn't adjust time
    def distanceSquaredToElementOld(self,x,y):
        Z = ( 2.0 * complex(x,y) - (self.z1 + self.z2) ) / (self.z2 - self.z1)
        if abs(Z.real) <= 1.0:
            dissq = ( Z.imag * self.L / 2 )**2
        elif Z.real < -1.0:
            dissq = ( x - self.x1 )**2 + ( y - self.y1 ) **2
        else:
            dissq = ( x - self.x2 )**2 + ( y - self.y2 ) **2
        return dissq
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

class LineDoubletString(Element):
    '''Laplace Linedoublet string for inhomogeneity boundary; single aquifer only, hardcoded to degree=2
    Coded such that there are only 2 unknowns and the potential is continuous across nodes
    Each line-doublet has unknowns of strength at right corner and quadratic strength'''
    # Prototype version with variable order or line-doublet
    oneOverTwoPii = 1.0 / ( 2.0 * pi * complex(0.0,1.0) )
    tiny = 1e-8; tinyimag = tiny * complex(0.0,1.0)
    def __init__(self,modelParent,xylist,aqin,aqout,Ndegree=2,overspec=1,closed=1):
        Element.__init__(self,modelParent)
        self.modelParent.addElement(self)
        self.xylist = xylist
        self.aqin = aqin; self.aqout = aqout
        self.type = 'linedoubletstring'
        self.Ndegree = Ndegree
        self.overspec = overspec
        if closed != 1: print 'LineDoubletString must be closed (for now)'
        self.closed = closed
        self.setCoefs()
    def __repr__(self):
	return 'LineDoubletString z1,z2,strengths: ' + str((self.xylist))
    def setCoefs(self):
        # If not closed, then this should change (less sides) !
        self.Nsides = len(self.xylist)
        self.ldList = []
        xy = self.xylist
        self.xlist = [xy[0][0]]; self.ylist = [xy[0][1]]
        self.z1 = []; self.z2 = []
        for i in range(self.Nsides-1):
            ld = LineDoublet(self.modelParent,xy[i][0],xy[i][1],xy[i+1][0],xy[i+1][1],\
                             self.aqin,self.aqout,self.Ndegree,self.overspec,addToModel=0)
            self.ldList = self.ldList + [ld]
            self.xlist = self.xlist + [xy[i+1][0]]; self.ylist = self.ylist + [xy[i+1][1]]
            self.z1 = self.z1 + [ complex( xy[i][0], xy[i][1] ) ]
            self.z2 = self.z2 + [ complex( xy[i+1][0], xy[i+1][1] ) ]
        if self.closed:
            ld = LineDoublet(self.modelParent,xy[-1][0],xy[-1][1],xy[0][0],xy[0][1],\
                             self.aqin,self.aqout,self.Ndegree,self.overspec,addToModel=0)
            self.ldList = self.ldList + [ld]
            self.xlist = self.xlist + [xy[0][0]]; self.ylist = self.ylist + [xy[0][1]]
            self.z1 = self.z1 + [ complex( xy[-1][0], xy[-1][1] ) ]
            self.z2 = self.z2 + [ complex( xy[0][0], xy[0][1] ) ]
        self.z1 = array(self.z1); self.z2 = array(self.z2)
        self.Nparam = 0
        for ld in self.ldList:
            self.Nparam = self.Nparam + ld.Ndegree  # Cause the constant and linear portions are joined
        self.parameters = zeros( (self.Nparam,1), 'd' )
    def potentialInfluence(self,aq,x,y):
        # Coded such that strength is continuous across nodes; assumes closed string
        # Gotta make sure that if (x,y) is on element it gets moved but is still in correct aquifer.
        # Gotta make sure that if (x,y) is at corner point it gets moved but is still in correct aquifer.
        zin = complex(x,y)
        bigZ = ( 2.0*zin - (self.z1 + self.z2) )/ (self.z2 - self.z1)  # This is an array
        bigZplus1 = bigZ + 1.0
        # print 'zin incoming ',zin
        for i in range(self.Nsides):
            if abs(bigZ.imag[i]) < self.tiny and abs(bigZ[i].real) < 1.0-self.tiny:  # point is on boundary; this could be nicer by setting the log
                if self.aqin == aq:
                    bigZnew = complex(bigZ[i].real,self.tiny)
                elif self.aqout == aq:
                    bigZnew = complex(bigZ[i].real,-self.tiny)
                zin = ( (self.z2[i]-self.z1[i])*bigZnew + (self.z1[i]+self.z2[i]) ) / 2.0
                # print 'on line, new zin ',zin
                break
            elif abs(bigZplus1[i]) < self.tiny:  # point is at left corner; move back a little along element; may fail if corner is very pointy
                if self.aqin == aq:
                    bigZnew = complex( -1.0+10*self.tiny, self.tiny )  # Assumes tan between angles is larger than 5.71 degrees
                else:
                    bigZnew = complex( -1.0+10*self.tiny, -self.tiny )
                zin = ( (self.z2[i]-self.z1[i])*bigZnew + (self.z1[i]+self.z2[i]) ) / 2.0
                # print 'at corner, new zin ',zin
                break
        xin = zin.real; yin = zin.imag
        pot = []
        for ld in self.ldList:
            potlapldho(xin,yin,ld.x1,ld.y1,ld.x2,ld.y2,ld.Ndegree,ld.potInf)
            pot = pot + [ ld.potInf ]
        pot = pot + [pot[0]]  #  Add last term at end
        for i in range(self.Nsides):
            for j in range( 2, self.ldList[i].Ndegree+1, 2 ):
                pot[i][j] = pot[i][j] - pot[i][0]
            for j in range( 3, self.ldList[i].Ndegree+1, 2 ):
                pot[i][j] = pot[i][j] - pot[i][2]
        rv = zeros( (self.Nparam,aq.Naquifers), 'd' )
        icount = 0
        for i in range(self.Nsides):
            rv[icount,0] = ( pot[i][0] + pot[i][1] ) + ( pot[i+1][0] - pot[i+1][1] )
            rv[icount+1:icount+self.ldList[i].Ndegree,0] = pot[i][2:]
            icount = icount + self.ldList[i].Ndegree
        return rv
    def dischargeInfluence(self,aq,x,y):
        zin = complex(x,y)
        bigZ = ( 2.0*zin - (self.z1 + self.z2) )/ (self.z2 - self.z1)  # This is an array
        bigZplus1 = bigZ + 1.0
        # print 'zin incoming ',zin
        for i in range(self.Nsides):
            if abs(bigZ.imag[i]) < self.tiny and abs(bigZ[i].real) < 1.0-self.tiny:  # point is on boundary; this could be nicer by setting the log
                if self.aqin == aq:
                    bigZnew = complex(bigZ[i].real,self.tiny)
                elif self.aqout == aq:
                    bigZnew = complex(bigZ[i].real,-self.tiny)
                zin = ( (self.z2[i]-self.z1[i])*bigZnew + (self.z1[i]+self.z2[i]) ) / 2.0
                # print 'on line, new zin ',zin
                break
            elif abs(bigZplus1[i]) < self.tiny:  # point is at left corner; move back a little along element; may fail if corner is very pointy
                if self.aqin == aq:
                    bigZnew = complex( -1.0+10*self.tiny, self.tiny )  # Assumes tan between angles is larger than 5.71 degrees
                else:
                    bigZnew = complex( -1.0+10*self.tiny, -self.tiny )
                zin = ( (self.z2[i]-self.z1[i])*bigZnew + (self.z1[i]+self.z2[i]) ) / 2.0
                # print 'at corner, new zin ',zin
                break
        xin = zin.real; yin = zin.imag
        rvx = zeros( (2*self.Nsides,aq.Naquifers), 'd' )
        rvy = zeros( (2*self.Nsides,aq.Naquifers), 'd' )
        disx = []; disy = []
        for ld in self.ldList:
            dislapldho(xin,yin,ld.x1,ld.y1,ld.x2,ld.y2,ld.Ndegree,ld.disxInf,ld.disyInf)
            disx = disx + [ ld.disxInf ]
            disy = disy + [ ld.disyInf ]
        for i in range(self.Nsides-1):        
            rvx[2*i,0] = ( disx[i][0] + disx[i][1] ) + ( disx[i+1][0] - disx[i+1][1] )
            rvx[2*i+1,0] = disx[i][2] - disx[i][0]
            rvy[2*i,0] = ( disy[i][0] + disy[i][1] ) + ( disy[i+1][0] - disy[i+1][1] )
            rvy[2*i+1,0] = disy[i][2] - disy[i][0]
        rvx[2*self.Nsides-2,0] = ( disx[-1][0] + disx[-1][1] ) + ( disx[0][0] - disx[0][1] )
        rvx[2*self.Nsides-1,0] = disx[-1][2] - disx[-1][0]
        rvy[2*self.Nsides-2,0] = ( disy[-1][0] + disy[-1][1] ) + ( disy[0][0] - disy[0][1] )
        rvy[2*self.Nsides-1,0] = disy[-1][2] - disy[-1][0]
        return [rvx,rvy]
    def getMatrixRows(self,elementList):
        rows = []; rhs = []; matleast = []
        # Hardcoded for pylayer 0 
        pylayer = 0
        for i in range(self.Nsides):
            ld = self.ldList[i]
            for icp in range(ld.Ncp):  # Modified for arbitrary number of control points
                row = []; rowin = []; rowout = []
                #xcp = ld.xcp[icp]; ycp = ld.ycp[icp]
                for e in elementList:
                    rowinpart = e.getMatrixCoefficients(self.aqin,pylayer,ld.xcpin[icp],ld.ycpin[icp],\
                                                    lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                    rowin = rowin + rowinpart.tolist()
                    rowoutpart = e.getMatrixCoefficients(self.aqout,pylayer,ld.xcpout[icp],ld.ycpout[icp],\
                                                    lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                    rowout = rowout + rowoutpart.tolist()
                    row = self.aqout.T[pylayer] * array(rowin) - self.aqin.T[pylayer] * array(rowout)
                    # Build least-squares matrix to multiply
                    if e == self:
                        matrow = self.aqout.T[pylayer] * rowinpart - self.aqin.T[pylayer] * rowoutpart
                        matleast = matleast + [matrow.tolist()]
                rhs = rhs + [\
                    self.aqin.T[pylayer] * self.modelParent.potentialInLayer(self.aqout,pylayer,ld.xcpout[icp],ld.ycpout[icp]) -\
                    self.aqout.T[pylayer] * self.modelParent.potentialInLayer(self.aqin,pylayer,ld.xcpin[icp],ld.ycpin[icp]) ]
                rows = rows + [row]
        # Now convert with least squares procedure
        matleastT = transpose(matleast)
##        print 'matleastT ',shape(matleastT)
##        print 'rows ',shape(array(rows))
##        rows = matrixmultiply( matleastT, array(rows) )
##        rhs = matrixmultiply( matleastT, array(rhs) )
        rows = dot( matleastT, array(rows) )
        rhs = dot( matleastT, array(rhs) )
        rows = hstack(( rows, rhs[:,newaxis] ))
        return rows.tolist()
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        self.parameters[:,0] = self.parameters[:,0] + xsol[icount:icount+self.Nparam]
        icount = icount + self.Nparam
        return icount
    def layout(self):
        return [0]
    def check(self):
        for ld in self.ldList:
            print 'LineDoublet from '+str(ld.z1)+' to '+str(ld.z2)+' Top Layer '
            for i in range(ld.Ndegree+1):
                if not ld.aqin.fakesemi:
                    print 'Control point: '+str(i)+\
                          ' Head inside : '+str(ld.modelParent.head(1,ld.xcpin[i],ld.ycpin[i],ld.aqin))+\
                          ' Head outside: '+str(ld.modelParent.head(1,ld.xcpout[i],ld.ycpout[i],ld.aqout))
                else:
                    print 'Control point: '+str(i)+\
                          ' Head inside : '+str(ld.modelParent.head(2,ld.xcpin[i],ld.ycpin[i],ld.aqin))+\
                          ' Head outside: '+str(ld.modelParent.head(1,ld.xcpout[i],ld.ycpout[i],ld.aqout))
        return None
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):
        # Not doing anything right now. Should really only do something if 1 aquifer both in and out. Else line-sink is doing the work
        changed = 0; stop = 0; xyznew = 0.0
##        [x1,y1,z1] = xyz1; [x2,y2,z2] = xyz2
##        if x1 < self.xmin - step or x1 > self.xmax + step or y1 < self.ymin - step or y1 > self.ymax + step:
##           return [changed, stop, xyznew, extrat]  # Away from line-doublet
##        z1 = complex(x1,y1); z2 = complex(x2,y2)
##        Z1 = (2.0*z1 - (self.z1+self.z2))/(self.z2-self.z1)
##        Z2 = (2.0*z2 - (self.z1+self.z2))/(self.z2-self.z1)
##        # If point 1 on one side and point 2 on other side, find out if segment intesects line-doublet
##        if (Z1.imag > 0 and Z2.imag < 0) or (Z1.imag < 0 and Z2.imag > 0):
##            Xintersec = Z1.real + (0.0-Z1.imag) / (Z2.imag-Z1.imag) * (Z2.real-Z1.real)
##            if abs(Xintersec) <= 1:  # There is an intersection point
##                changed = 1
##                Znew = complex(Xintersec,0.0)
##                znew = ( (self.z2-self.z1)*Znew + (self.z1+self.z2) ) / 2.0  # New complex coordinate
##                xnew = znew.real; ynew = znew.imag
##                horstepold = sqrt( (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 )
##                horstepnew = sqrt( (xyz1[0]-xnew)**2 + (xyz1[1]-ynew)**2 )
##                zvertnew = xyz1[2] + horstepnew / horstepold * ( xyz2[2] - xyz1[2] )  # Just to be confusing, this is z (vertical coordinate)                 
##                if Z1.imag < 0:  # old point on outside, so new point on inside
##                    dis = self.modelParent.dischargeNormBelowZ(self.aqout,xnew,ynew,zvertnew,self.thetaNormOut)
##                    zvertnew2 = self.modelParent.zForGivenNormalDischarge\
##                                (self.aqin,zvertnew,xnew,ynew,self.thetaNormOut,dis)
##                else:
##                    dis = self.modelParent.dischargeNormBelowZ(self.aqin,xnew,ynew,zvertnew,self.thetaNormOut)
##                    zvertnew2 = self.modelParent.zForGivenNormalDischarge\
##                                (self.aqout,zvertnew,xnew,ynew,self.thetaNormOut,dis)
##                xyznew[:] = [xnew,ynew,zvertnew2]  # Should modify z
        return [changed, stop, xyznew]  # Didn't adjust time
    def distanceSquaredToElement(self,x,y):
        dissq = 1e20
        for ld in self.ldList:
            dissqld = ld.distanceSquaredToElement(x,y)
            if dissqld < dissq:
                dissq = dissqld
        return dissq

            

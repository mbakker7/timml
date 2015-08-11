'''
mluflow.py contains the Uniform flow class
This file is part of the TimML library and is distributed under
the GNU LPGL. See the TimML.py file for more details.
(c) Mark Bakker, 2002-2007
'''

from numpy import *
from mlelement import *

class Uflow(Element):
    '''Uflow class, angle in degrees
    Uflow class. Derived from Element
    Attributes provided on input:
    - modelParent: model it is added to
    - grad: gradient
    - angle: angle in degrees with positive x
    Attributes computed:
    - aquiferParent: AquiferData parent
    - Qx: comprehensive discharge in x direction
    - Qy: comprehensive discharge in y direction
    All attributes from element.
    '''
    def __init__(self,modelParent,grad,angle,label=None):
	Element.__init__(self,modelParent)
	self.grad = float(grad); self.angle = float(angle)
	self.label = label
	self.type = 'uflow'
	self.setCoefs()
	self.modelParent.addElement(self)
    def __repr__(self):
	return 'Uflow grad, angle: ' + str((self.grad,self.angle))
    def setCoefs(self):
	self.aquiferParent = self.modelParent.aq
	assert self.aquiferParent.type == self.aquiferParent.conf, 'TimML input error: Model aquifer most be confined'
        self.Qx = self.grad * cos(self.angle*pi/180) * self.aquiferParent.Tcomp
        self.Qy = self.grad * sin(self.angle*pi/180) * self.aquiferParent.Tcomp
        self.parameters = array([[self.Qx],[self.Qy]])
    def potentialInfluence(self,aq,x,y,z=0,t=0):
        rv = zeros( (2,aq.Naquifers), 'd' )
        if aq.type == aq.conf:
            rv[0,0] = -x
            rv[1,0] = -y
        return rv
    def dischargeInfluence(self,aq,x,y):
        rvx = zeros( (2,aq.Naquifers), 'd' ); rvy = zeros( (2,aq.Naquifers), 'd' )
        if aq.type == aq.conf:
            rvx[0,0] = 1.0
            rvy[1,0] = 1.0
        return [rvx,rvy]
    def check(self):
        print 'Uflow with gradient '+str(self.grad)+' and angle '+str(self.angle)+\
              ' has no unknown parameters'
        return None

class Triple(Element):
    '''Triple class for three reference points giving a constant and uniform flow
    '''
    def __init__(self,modelParent,x1,y1,h1,layer1,x2,y2,h2,layer2,x3,y3,h3,layer3,label=None):
	Element.__init__(self,modelParent)
	self.x = array([x1,x2,x3],'d'); self.y = array([y1,y2,y3],'d')
	self.h = array([h1,h2,h3],'d'); self.layers = array([layer1,layer2,layer3])
	self.label = label
	self.type = 'triple'
	self.setCoefs()
	self.modelParent.addElement(self)
    def __repr__(self):
	return 'Triple points x,y,head,layer ' + str((self.x,self.y,self.h,self.layers)) 
    def setCoefs(self):
	self.aquiferParent = self.modelParent.aq
	assert self.aquiferParent.type == self.aquiferParent.conf, 'TimML input error: Model aquifer most be confined'
        self.pylayers = self.layers  # fixed to zero base
        self.parameters = zeros((3,1),'d')
    def potentialInfluence(self,aq,x,y,z=0,t=0):
        rv = zeros( (3,aq.Naquifers), 'd' )
        if aq.type == aq.conf:
            rv[0,0] = 1.0
            rv[1,0] = -x
            rv[2,0] = -y
        return rv
    def dischargeInfluence(self,aq,x,y):
        rvx = zeros( (3,aq.Naquifers), 'd' ); rvy = zeros( (3,aq.Naquifers), 'd' )
        if aq.type == aq.conf:
            rvx[1,0] = 1.0
            rvy[2,0] = 1.0
        return [rvx,rvy]
    def getMatrixRows(self,elementList):
        rows = []
        for i in range(3):
            row = zeros(0,'d')
            for e in elementList:
                rowpart = e.getMatrixCoefficients(self.aquiferParent,self.pylayers[i],self.x[i],self.y[i],\
                                                  lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
                row = hstack(( row, rowpart ))
            row = hstack(( row, self.aquiferParent.headToPotential(self.pylayers[i],self.h[i]) - \
                           self.modelParent.potentialInLayer(self.aquiferParent,self.pylayers[i],self.x[i],self.y[i]) ))
            rows.append(row.tolist())
        return rows
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        self.parameters[:,0] = self.parameters[:,0] + xsol[icount:icount+3]
        return icount+3
    def check(self):
        print 'Triple point, Constant, Qx, Qy '+str(self.parameters[:,0])
        for i in range(3):
            print 'i,x,y,h,layer '+str((self.layers[i],self.x[i],self.y[i],self.h[i],self.layers[i]))+\
                  ', Computed head '+str(self.modelParent.head(self.layers[i],self.x[i],self.y[i]))
        return None

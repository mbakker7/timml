'''
mlconstant.py contains the Constant class
This file is part of the TimML library and is distributed under
the GNU LPGL. See the TimML.py file for more details.
(c) Mark Bakker, 2002-2007
'''

from mlelement import *

class Constant(Element):
    '''Reference point; layer is list of aquifer in which ref pt is placed
    Base class for constant function; derived from Element
    Attributes provided on input:
    - modelParent: model it belongs to
    - xr: x location of reference point
    - yr: y location of reference point
    - head: head at reference point
    - layer: layer number in which reference point is placed
    Attributes computed:
    - pylayer: Python layer number in which reference point is placed
    - pot: potential equivalent of specified head
    All attributes from Element
    '''
    def __init__(self,modelParent,xr,yr,head,layer=0,label=None):
	Element.__init__(self,modelParent)
	self.aquiferParent = self.modelParent.aq.findAquiferData(xr,yr)
        self.xr = float(xr); self.yr = float(yr); self.head = float(head);
        self.layer = layer; 
        if iterable(layer):
	    assert len(self.layer)==1, "TimML Input error: Constant can only be defined in one layer"
	    self.layer = self.layer[0]
        assert self.aquiferParent.type == self.aquiferParent.conf, "TimML Input error: Constant can only be defined in confined aquifer"
        self.label = label
        self.type = 'constant'
        self.setCoefs()
        self.modelParent.addElement(self)
    def __repr__(self):
	return 'Constant xr,yr,head,C,layer: ' + str((self.xr,self.yr,self.head,self.parameters[0,0],self.pylayer))
    def setCoefs(self):
        self.pylayer = self.layer;  # Fixed to zero base
        self.pot = self.aquiferParent.headToPotential(self.pylayer,self.head)
        self.parameters = array([[0.0]])        
    def potentialInfluence(self,aq,x,y,z=0,t=0):
        pot = zeros((1,aq.Naquifers),'d')
        if aq.type == aq.conf: pot[0,0] = 1.0
        return pot
    def potentialContribution(self,aq,x,y):
        '''Returns VECTOR of potentialContribution; doesn't multiply with eigenvectors'''
        return self.parameters[0,0] * self.potentialInfluence(aq,x,y)[0,:]
    def potentialCollection(self,potsum,potadd,elementList,aq,x,y):
        if aq.type == aq.conf:
            potsum[0] = potsum[0] + self.parameters[0,0]
        return potsum
    def dischargeInfluence(self,aq,x,y):
        dis = zeros((1,aq.Naquifers),'d')
        return [ dis,dis ]
    def dischargeCollection(self,dissum,disadd,elementList,aq,x,y):
        return dissum
    def totalDischargeInfluence(self,aq,pylayer,x,y):
        rv = zeros(1,'d')
        return rv    
    def getMatrixRows(self,elementList):
        row = zeros(0,'d')
        for e in elementList:
            rowpart = e.getMatrixCoefficients(self.aquiferParent,self.pylayer,self.xr,self.yr,\
                                              lambda el,aq,pylayer,x,y:el.potentialInfluenceInLayer(aq,pylayer,x,y))
            row = hstack(( row, rowpart ))
        row = hstack(( row, self.pot - self.modelParent.potentialInLayer(self.aquiferParent,self.pylayer,self.xr,self.yr) ))
        return [row.tolist()]
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        return func(self,aq,pylayer,x,y)
    def takeParameters(self,xsol,icount):
        self.parameters[0,0] = self.parameters[0,0] + xsol[icount]
        return icount+1
    def check(self):
        print 'Constant at '+str((self.xr,self.yr))+' in layer '+str(self.layer)
        print 'Specified head: '+str(self.head)+' Computed head: '+str(self.modelParent.head(self.layer,self.xr,self.yr))+\
              ' Constant: '+str(self.parameters[0,0])
        return None
        
class NoFlowFromInfinity(Constant):
    def __init__(self,modelParent,label=None):
	Constant.__init__(self,modelParent,0,0,0,[1],label)
    def __repr__(self):
	return 'NoFlowFromInfinity '
    def getMatrixRows(self,elementList):
        row = zeros(0,'d')
        for e in elementList:
            rowpart = e.getMatrixCoefficients(self.aquiferParent,self.pylayer,self.xr,self.yr,\
                                              lambda el,aq,pylayer,x,y:el.totalDischargeInfluence(aq,pylayer,x,y))
            row = hstack(( row, rowpart ))
        row = hstack(( row, -self.modelParent.totalDischargeFromInf() ))
        return [row.tolist()]
    def check(self):
        print 'NoFlorFromInfinity. Constant: '+str(self.parameters[0,0])
        print 'Computed flow from infinity: '+str(self.modelParent.totalDischargeFromInf())
        return None

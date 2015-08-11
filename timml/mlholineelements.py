'''
mlholineelements.py contains the higherorder LineSinkHO and LineDoubletHO classes.
Elements are screened in all layers and strength must be specified for all layers.
This file is part of the TimML library and is distributed under
the GNU LPGL. See the TimML.py file for more details.
(c) Mark Bakker, 2003-2007
'''

from numpy import *
import numpy.linalg as linalg
from mlelement import *
from besselaes import *

class LineSinkHo(Element):
    '''sigma is array of (order,Naquifers)
    This is a special line-sink, for which a strength is specified in each layer, for each order'''
    def __init__(self,modelParent,x1,y1,x2,y2,order,sigma,aquiferParentFixed=None,addToModel=1,Bessel=1):
	Element.__init__(self,modelParent)
        self.x1 = float(x1); self.y1 = float(y1); self.x2 = float(x2); self.y2 = float(y2);
        self.sigma = array(sigma,'d')
        self.order = order
        self.aquiferParentFixed = aquiferParentFixed
        self.Bessel = Bessel
        self.setCoefs()
	if addToModel:
            self.modelParent.addElement(self)
    def __repr__(self):
	return 'LineSinkHo z1,z2,order,sigma: ' + str((self.x1,self.y1,self.x2,self.y2,self.order,list(self.sigma) ))
    def setCoefs(self):
        self.xy1 = (self.x1,self.y1); self.xy2 = (self.x2,self.y2)
        self.z1 = complex(self.x1,self.y1); self.z2 = complex(self.x2,self.y2);
        self.xc = 0.5*(self.x1+self.x2); self.yc = 0.5*(self.y1+self.y2)
        self.L = abs(self.z2-self.z1)
        if self.aquiferParentFixed == None:
            self.aquiferParent = self.modelParent.aq.findAquiferData(self.xc,self.yc)  # Determined at xc,yc
        else:
            self.aquiferParent = self.aquiferParentFixed
        self.parameters = ones((self.order+1,1),'d')  # Somewhat unsatisfactory as this is not the strength! Implemented this way because it is hard to get a paramter to be zero
        if self.Bessel:
            assert shape(self.sigma) == (self.order+1,self.aquiferParent.Naquifers) , \
                'TimML Input error: sigma must be of size '+str((self.order+1,self.aquiferParent.Naquifers))
            self.coef = zeros(shape(self.sigma),'d')
            for i in range(self.order+1):
                ramat = self.aquiferParent.eigvec
                rb = self.sigma[i,:]
                self.coef[i,:] = linalg.solve(ramat,rb)
        else:
            assert shape(self.sigma) == (self.order+1,) , \
               'TimML Input error: sigma must be of size '+str(self.order+1)
            self.coef = self.sigma[:,newaxis]
    def potentialInfluence(self,aq,x,y):
        pot = zeros((self.order+1,aq.Naquifers),'d')
        if aq == self.aquiferParent and self.Bessel:
            lab = self.aquiferParent.zeropluslab
            potInf = zeros((self.aquiferParent.Naquifers),'d')
            for i in range(self.order+1):
                potbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers,lab,i,potInf)  # Call FORTRAN extension
                pot[i,:] = potInf
            rv = self.coef * pot
        else:
            lab = zeros(1,'d')
            potInf = zeros(1,'d')
            for i in range(self.order+1):
                potbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,1,lab,i,potInf)  # Call FORTRAN extension
                pot[i,0] = potInf[0]
            pot[:,0] = pot[:,0] * self.coef[:,0]
            rv = pot
        return rv
    def dischargeInfluence(self,aq,x,y):
        disx = zeros((self.order+1,self.aquiferParent.Naquifers),'d')
        disy = zeros((self.order+1,self.aquiferParent.Naquifers),'d')
        if aq == self.aquiferParent and self.Bessel:
            lab = array( [0.]+list(self.aquiferParent.lab),'d' )
            disxInf = zeros((self.aquiferParent.Naquifers),'d')
            disyInf = zeros((self.aquiferParent.Naquifers),'d')
            for i in range(self.order+1):
                disbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers,lab,i,disxInf,disyInf)  # Call FORTRAN extension
                disx[i,:] = disxInf; disy[i,:] = disyInf
            rvx = self.coef * disx; rvy = self.coef * disy
        else:
            lab = zeros((1),'d')
            disxInf = zeros((1),'d'); disyInf = zeros((1),'d')
            for i in range(self.order+1):
                disbeslsho(x,y,self.x1,self.y1,self.x2,self.y2,1,lab,i,disxInf,disyInf)  # Call FORTRAN extension
                disx[i,0] = disxInf[0]; disy[i,0] = disyInf[0]
            disx[:,0] = disx[:,0] * self.coef[:,0]
            disy[:,0] = disy[:,0] * self.coef[:,0]
            rvx = disx; rvy = disy
        return [rvx,rvy]
    def layout(self):
        rv = [ 2,[],[] ]
        rv[1] = [ self.x1,self.x2 ]
        rv[2] = [ self.y1,self.y2 ]
        return rv
    def check(self):
        print 'LinesinkHo from '+str(self.xy1)+' to '+str(self.xy2)+' has no unknown parameters'
        return None

class LineDoubletHo(Element):
    '''potjump is array of (order,Naquifers)'''
    def __init__(self,modelParent,x1,y1,x2,y2,order,potjump,aquiferParentFixed=None,addToModel=1,Bessel=1):
	Element.__init__(self,modelParent)
        self.x1 = float(x1); self.y1 = float(y1); self.x2 = float(x2); self.y2 = float(y2);
        self.potjump = array(potjump,'d')
        self.order = order
        self.aquiferParentFixed = aquiferParentFixed
        self.Bessel = Bessel
        self.setCoefs()
	if addToModel:
            self.modelParent.addElement(self)
    def __repr__(self):
	return 'LineDoubletHo z1,z2,order,potjump: ' + str((self.x1,self.y1,self.x2,self.y2,self.order,list(self.potjump) ))
    def setCoefs(self):
        self.xy1 = (self.x1,self.y1); self.xy2 = (self.x2,self.y2)
        self.z1 = complex(self.x1,self.y1); self.z2 = complex(self.x2,self.y2);
        self.xc = 0.5*(self.x1+self.x2); self.yc = 0.5*(self.y1+self.y2)
        self.L = abs(self.z2-self.z1)
        if self.aquiferParentFixed == None:
            self.aquiferParent = self.modelParent.aq.findAquiferData(self.xc,self.yc)  # Determined at xc,yc
        else:
            self.aquiferParent = self.aquiferParentFixed
        self.parameters = ones((self.order+1,1),'d')  # Somewhat unsatisfactory as this is not the strength! Implemented this way because it is hard to get a paramter to be zero
        if self.Bessel:
            assert shape(self.potjump) == (self.order+1,self.aquiferParent.Naquifers) , \
                'TimML Input error: potjump must be of size '+str((self.order+1,self.aquiferParent.Naquifers))
            self.coef = zeros(shape(self.potjump),'d')
            for i in range(self.order+1):
                ramat = self.aquiferParent.eigvec
                rb = self.potjump[i,:]
                self.coef[i,:] = linalg.solve(ramat,rb)
        else:
            assert shape(self.potjump) == (self.order+1,) , \
               'TimML Input error: potjump must be of size '+str(self.order+1)
            self.coef = self.potjump[:,newaxis]
    def potentialInfluence(self,aq,x,y):
        pot = zeros((self.order+1,aq.Naquifers),'d')
        if aq == self.aquiferParent and self.Bessel:
            lab = array( [0.]+list(self.aquiferParent.lab),'d' )
            potInf = zeros((self.aquiferParent.Naquifers),'d')
            for i in range(self.order+1):
                potbesldho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers,lab,i,potInf)  # Call FORTRAN extension
                pot[i,:] = potInf
            rv = self.coef * pot
        else:  # Only Laplace part
            potInf = zeros(self.order+1,'d')
            potlapldho(x,y,self.x1,self.y1,self.x2,self.y2,self.order,potInf)  # Call FORTRAN extension
            pot[:,0] = potInf * self.coef[:,0]
            rv = pot
        return rv
    def dischargeInfluence(self,aq,x,y):
        disx = zeros((self.order+1,self.aquiferParent.Naquifers),'d')
        disy = zeros((self.order+1,self.aquiferParent.Naquifers),'d')
        if aq == self.aquiferParent and self.Bessel:
            lab = array( [0.]+list(self.aquiferParent.lab),'d' )
            disxInf = zeros((self.aquiferParent.Naquifers),'d')
            disyInf = zeros((self.aquiferParent.Naquifers),'d')
            for i in range(self.order+1):
                disbesldho(x,y,self.x1,self.y1,self.x2,self.y2,self.aquiferParent.Naquifers,lab,i,disxInf,disyInf)  # Call FORTRAN extension
                disx[i,:] = disxInf; disy[i,:] = disyInf
            rvx = self.coef * disx; rvy = self.coef * disy
        else:  # Only Laplace part
##            lab = zeros((1),'d')
##            disxInf = zeros((1),'d'); disyInf = zeros((1),'d')
##            for i in range(self.order+1):
##                disbesldho(x,y,self.x1,self.y1,self.x2,self.y2,1,lab,i,disxInf,disyInf)  # Call FORTRAN extension
##                disx[i,0] = disxInf[0]; disy[i,0] = disyInf[0]
##            disx[:,0] = disx[:,0] * self.coef[:,0]
##            disy[:,0] = disy[:,0] * self.coef[:,0]
##            rvx = disx; rvy = disy

            disxInf = zeros((self.order+1),'d'); disyInf = zeros((self.order+1),'d')
            dislapldho(x,y,self.x1,self.y1,self.x2,self.y2,self.order,disxInf,disyInf)  # Call FORTRAN extension
            disx[:,0] = disxInf * self.coef[:,0]
            disy[:,0] = disyInf * self.coef[:,0]
            rvx = disx; rvy = disy
        return [rvx,rvy]
    def layout(self):
        rv = [ 2,[],[] ]
        rv[1] = [ self.x1,self.x2 ]
        rv[2] = [ self.y1,self.y2 ]
        return rv
    def check(self):
        print 'LineDoubletHo from '+str(self.xy1)+' to '+str(self.xy2)+' has no unknown parameters'
        return None

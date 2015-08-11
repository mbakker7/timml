'''
mlelement.py contains the Element base class
This file is part of the TimML library and is distributed under
the GNU LPGL. See the TimML.py file for more details.
(c) Mark Bakker, 2002-2007
'''

from numpy import *

class Element:
    def __init__(self,modelParent):
        self.modelParent = modelParent
        self.aquiferParent = None
        self.parameters = array([[]],'d')
	self.elementList = [ self ]
	self.label = None
	self.type = None  # Type of element is needed for fast evaluation functions
	self.aqdum=None; self.ldum=0; self.xdum=0.0; self.ydum=0.0
	self.pylayers = [0] # In case it is not specified
    def setCoefs(self):
        '''Sets coefficients for the element. Must be overloaded
        if it is supposed to do something. This method is called upon
        construction, and by solve if all elements are re-initialized'''
        return
    def potentialInfluence(self,aq,x,y):
        '''Returns array with Phi_L in first spot of each row and
        Phi_m in remaining spots; one row for each parameter.
        When aq is not aquifer in which element is situated, only returns Phi_L
        in first spot, with zeros in remaining spots, according to number of aquifers'''
	raise 'Must overload Element.potentialInfluence()'
    def potentialInfluenceInLayer(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in pylayer as 1D array (1 value per parameter)'''
        potInf = self.potentialInfluence(aq,x,y)
        return sum( potInf * aq.eigvec[pylayer,:], 1 )
    def potentialInfluenceSpecLayers(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in specified layers'''
        # Faster version, the old one was DOG slow
        potInf = self.potentialInfluence(aq,x,y)
        rv = []
        for i in pylayer:
            rv.append( sum( potInf * aq.eigvec[i,:], 1 ) )
        return array(rv)
    def potentialContribution(self,aq,x,y):
        '''Returns VECTOR of potentialContribution; doesn't multiply with eigenvectors'''
        return sum( self.parameters * self.potentialInfluence(aq,x,y), 0 )
    def potentialCollection(self,potsum,potadd,elementList,aq,x,y):
        for e in elementList:
            potsum = potsum + e.potentialContribution(aq,x,y)
        return potsum
    def potentialInfluenceAllLayers(self,aq,pylayer,x,y):
        '''Returns PotentialInfluence function in aquifer aq in all layers'''
        # Faster version, the old one was DOG slow
        potInf = self.potentialInfluence(aq,x,y)
        rv = []
        for i in range(aq.Naquifers):
            rv.append( sum( potInf * aq.eigvec[i,:], 1 ) )
        return array(rv)
    def dischargeCollection(self,dissum,disadd,elementList,aq,x,y):
        for e in elementList:
            dissum = dissum + e.dischargeContribution(aq,x,y)
        return dissum
    def dischargeInfluence(self,aq,x,y):
        '''Returns an array containing two lists, one for Qx, one for Qy.
        Each list is a list of lists with Laplacian solution in first spot and
        Bessel solutions in remaining spots; one list for each parameter.
        When aq is not aquifer in which element is situated, only returns Laplace portion
        with zeros in remaining spots'''
	raise 'Must overload Element.dischargeInfluence()'
    def dischargeContribution(self,aq,x,y):
        '''Returns matrix with two rowvectors of dischargeContributions Qx and Qy'''
        #  Should be improved for speed. Take eigvec out here!
        disInf = self.dischargeInfluence(aq,x,y)
        disxInf = sum( self.parameters * disInf[0], 0 )
        disyInf = sum( self.parameters * disInf[1], 0 )
        return array([disxInf,disyInf])
    def dischargeInfluenceAllLayers(self,aq,dumlayer,x,y):
        '''Returns dischargeInfluenceAllLayers function in aquifer aq as an array
        Needs to be overloaded because there is no parameter outside the functions
        Needs to be modified for inside and outside'''
        [disx,disy] = self.dischargeInfluence(aq,x,y)
        rvx = []; rvy = []
        for i in range(aq.Naquifers):
            rvx.append( sum( disx * aq.eigvec[i,:], 1 ) )
            rvy.append( sum( disy * aq.eigvec[i,:], 1 ) )
        return [ array(rvx), array(rvy) ]
    def dischargeInfluenceSpecLayers(self,aq,pylayer,x,y):
        '''Returns dischargeInfluenceAllLayers function in aquifer aq as an array
        Needs to be overloaded because there is no parameter outside the functions
        Needs to be modified for inside and outside'''
        [disx,disy] = self.dischargeInfluence(aq,x,y)
        rvx = []; rvy = []
        for i in pylayer:
            rvx.append( sum( disx * aq.eigvec[i,:], 1 ) )
            rvy.append( sum( disy * aq.eigvec[i,:], 1 ) )
        return [ array(rvx), array(rvy) ]
    def dischargeInfluenceInLayer(self,aq,pylayer,x,y):
        '''Returns dischargeInfluenceInLayer function in aquifer aq as an array'''
        [disx,disy] = self.dischargeInfluence(aq,x,y)
        rvx = sum( disx * aq.eigvec[pylayer,:], 1 )
        rvy = sum( disy * aq.eigvec[pylayer,:], 1 )
        return [rvx,rvy]
    def totalDischargeInfluence(self,aq,pylayer,x,y):
        '''Returns array with Qtot_L in first spot of each row and
        Qtot_m in remaining spots; one row for each parameter.'''
	raise 'Must overload totalDischargeInfluence()'
    def totalDischarge(self):
        return sum( self.parameters[:,0] * self.totalDischargeInfluence(0,0,0,0) )
    def getMatrixRows(self,elementList):
        '''This method is passed the elementList and returns a list of M lists
        (one list for each unknown parameter). Each list represents a full row
        of the matrix with  the right-hand side of the equation in the last spot.
        Returns an empty list where the element does not have any unknown
        parameters. The method loops through all the elements in elementList
        and calls its getMatrixCoefficients method with the appropriate (combination of)
        influence function as func to build its rows of the matrix.'''
        return []    
    def getMatrixCoefficients(self,aq,pylayer,x,y,func):
        '''This method is passed the AquiferData instance aq, pylayer, an (x,y)
        pair, and a temporary function func. It will evaluate the function func
        at (x,y) for aq in layer pylayer if the element has unknown parameters;
        if the element has no unknown parameters, it returns an empty list. func
        can be any of the influence functions that have the form
        func(self,aq,pylayer,x,y), or a combination there-off.'''
        return zeros(0,'d')
    def takeParameters(self,xsol,icount):
        ''' This method is passed the solution vector and a counter. It sets its
        unknown parameters equal to the values of the solution vector, starting at
        icount, and increases icount by its number of unknown parameters. icount
        is then passed back. If the element does not have any unknown parameters
        it returns icount without changing it.'''
        return icount
    def layout(self):
        '''Returns list consisting of 3 items.
        The first item is an integer indicating the number of (x,y) pairs.
        The second item is a list of x coordinates
        The third item is a list of y coordinates'''
        return [0]
    def zeroFunction(self,aqdum,ldum,xdum,ydum):
        '''Returns list of zeros of length number of parameters'''
        return list( zeros( len(self.parameters),'d' ) )
    def check(self):
        '''Presents data from which may be verified whether correct solution is obtained'''
        print 'check not implemented for this element'
        return None
    def parSum(self):
        return sum( sum(abs(self.parameters)) )
    def qzTop(self,x,y):
        '''Returns qz at top of aquifer at (x,y) due to element
        Only non-zero for areal recharge'''
        return 0.0
    def nearElement(self,pyLayer,xyz1,xyz2,step,idir):
        '''Peforms something special when tracing near/over/into element
        Input: xyz1: last point of traceline, xyz2: candidate next point of trace line
        Returns: list with
        - change: logical set to 1 if new xyz has been modified
	- stop: logical set to 1 if trace must be terminated at new xyz
	- xyznew: XYZ object of new point of trace line
	- extrat: extra time it takes to reach point that cannot be
	computed from standard stepsize (for example leakage through resistance layer)'''
        changed = 0; stop = 0; xyznew = zeros(3,'d')
        return [changed, stop, xyznew]
    def distanceSquaredToElement(self,x,y):
        '''Returns distance squared to element. Used for deciding tracing step'''
        dissq = 1e30
        return dissq
    def iterationStep(self,alpha):
        '''Does whatever adjustment the element should do during an iterative solve'''
        return
    def addElementToCollection(self,ml):
        ml.addElementToCollection(self)
    def getMatrixRows_nonlinear(self,elementList):
        '''Function that returns new equations when a non-linear relationshipe results
        in a different condition. Returns None if no such condition exists'''
        return None

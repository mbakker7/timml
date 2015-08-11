'''
ml.py contains the Model class
This file is part of the TimML library and is distributed under
the GNU LPGL. See the TimML.py file for more details.
(c) Mark Bakker, 2002-2007
'''

from sys import stdout
from numpy import *
import numpy.linalg as linalg
from mlaquifer import *
from mlelement import *
#from TimMLgui import setActiveModel
from mltrace import *
from scipy.integrate import romberg

class Model:
    '''Model class; all parameters from top down as lists
    Attributes:
    - elementList: List of elements in model
    - aq: instance of (background) AquiferData class
    For parameters provided on input see AquiferData class
    '''
    def __init__(self,k=[1],zb=[0],zt=[1],c=[],n=[],nll=[],semi=False,type='conf',hstar=0.0):
        self.elementList = []
        self.elementDict = {}
        self.collectionDict = {}
        self.lakeList = []
        self.aq = Aquifer(k,zb,zt,c,n,nll,semi=False,type=type,hstar=hstar)
#        setActiveModel(self)
    def __repr__(self):
	return 'elementList ' + str(self.elementList)
    def addElement(self,el):
        '''Adds Element instance el to elementList'''
        self.elementList.append(el)
        if el.label != None:
            self.elementDict[el.label] = el
    def addElementToCollection(self,el):
        if el.type == None:
            print 'element with no type given'
        elif self.collectionDict.has_key(el.type):
            self.collectionDict[el.type].append(el)
        else:
            self.collectionDict.__setitem__(el.type,[el])
    def removeElement(self,elementLabel):  ### Doesn't remove element from dictionary !!!!
        el = self.elementDict[elementLabel]
        for i in range(len(self.elementList)):
            if el == self.elementList[i]:
                del self.elementList[i]
                print 'The following element was removed: ',el
                return
    def potentialVector(self,x,y,aq=None):
        '''Returns row vector of potentials at (x,y) for AquiferData aq; if aq=None, program will find it.
        Performs multiplication with eigenvectors here'''
        if aq == None: aq = self.aq.findAquiferData(x,y)
        pot = 0.0
        for el in self.elementList:
            pot = pot + el.potentialContribution(aq,x,y)
        return sum( pot * aq.eigvec , 1 )
    def potentialCollection(self,x,y,aq=None):
        if aq == None: aq = self.aq.findAquiferData(x,y)
        potsum = zeros(aq.Naquifers,'d')
        potadd = zeros(aq.Naquifers,'d')
        for col in self.collectionDict:
            potsum = self.collectionDict[col][0].potentialCollection(potsum,potadd,self.collectionDict[col],aq,x,y)
        return sum( potsum * aq.eigvec , 1 )        
    def potentialVectorCol(self,x,y,aq=None):
        '''Returns column vector of potentials at (x,y) for AquiferData aq; if aq=None, program will find it'''
        pot = self.potentialVector(x,y,aq)
        return pot[:,newaxis]
    def potential(self,layer,x,y,aq=None):
        '''Returns potential in layer at (x,y) for AquiferData aq; if aq=None, program will find it''' 
        pot = self.potentialVector(x,y,aq)
	pylayer = layer # fixed to base zero
        return pot[pylayer]
    def potentialInLayer(self,aq,pylayer,x,y):
        '''Used for internal purposes only. Returns potential in aq and pylayer at (x,y)'''
        pot = self.potentialCollection(x,y,aq)
        return pot[pylayer]
    def headVectorOld(self,x,y,aq=None):
        '''Returns row vector of heads at (x,y) for AquiferData aq; if aq=None, program will find it'''
        if aq==None: aq = self.aq.findAquiferData(x,y)
        potVec = self.potentialVector(x,y,aq)
        return aq.potentialVectorToHeadVector(potVec)
    def headVector(self,x,y,aq=None):  ### Only works after a solve (I don't think it is needed otherwise)
        '''Returns row vector of heads at (x,y) for AquiferData aq; if aq=None, program will find it'''
        if aq==None: aq = self.aq.findAquiferData(x,y)
        potVec = self.potentialCollection(x,y,aq)
        return aq.potentialVectorToHeadVector(potVec)
    def dischargeVector(self,x,y,aq=None):
        '''Returns a list of Qx array and Qy array at (x,y) for AquiferData aq; if aq=None, program will find it'''
        if aq == None: aq = self.aq.findAquiferData(x,y)
        Qx = 0.0; Qy = 0.0
        for el in self.elementList:
            dis = el.dischargeContribution(aq,x,y)
            Qx = Qx + dis[0]; Qy = Qy + dis[1]
        Qx = sum( Qx * aq.eigvec , 1 )
        Qy = sum( Qy * aq.eigvec , 1 )
        return [Qx,Qy]
    def dischargeCollection(self,x,y,aq=None):
        if aq == None: aq = self.aq.findAquiferData(x,y)
        dissum = zeros((2,aq.Naquifers),'d')
        disadd = zeros((2,aq.Naquifers),'d')
        for col in self.collectionDict:
            dissum = self.collectionDict[col][0].dischargeCollection(dissum,disadd,self.collectionDict[col],aq,x,y)
        disx = sum( dissum[0,:] * aq.eigvec , 1 )
        disy = sum( dissum[1,:] * aq.eigvec , 1 )
        return [disx, disy]
    def head(self,layer,x,y,aq=None):
        '''Returns head in layer at (x,y) for AquiferData aq; if aq=None, program will find it''' 
        if aq == None: aq = self.aq.findAquiferData(x,y)
        if aq.fakesemi: layer = layer+1
	pylayer = layer # fixed to base zero
        pot = self.potential(pylayer,x,y,aq)
        return aq.potentialToHead(pylayer,pot)
    def head3D(self,x,y,z):
        '''Returns head at (x,y,z), if z in aquifer; otherwise returns 0 ''' 
        aq = self.aq.findAquiferData(x,y)
        pylayer = aq.inWhichPyLayer(z)
        if pylayer == 9999: return 0.0
        if pylayer >= 0:
            pot = self.potential(pylayer,x,y,aq)  # fixed to base zero
            return aq.potentialToHead(pylayer,pot)
        return 0.0
    def head3Dinterp(self,x,y,z):
        '''Returns head at (x,y,z), if z in aquifer; otherwise returns 0, linearly interpolated in the vertical ''' 
        aq = self.aq.findAquiferData(x,y)
        pylayer = aq.inWhichPyLayer(z)
        h = 0.0
        if pylayer == 9999: return 0.0
        if pylayer >= 0:
            headvec = self.headVector(x,y,aq)
            if z >= aq.zcenter[0]:
                h = headvec[0]
            elif z <= aq.zcenter[-1]:
                h = headvec[-1]
            else:
                if z > aq.zcenter[pylayer]:
                    h = headvec[pylayer] + \
                        (z - aq.zcenter[pylayer]) / (aq.zcenter[pylayer-1] - aq.zcenter[pylayer]) *\
                        (headvec[pylayer-1]-headvec[pylayer])
                else:
                    h = headvec[pylayer+1] + \
                        (z - aq.zcenter[pylayer+1]) / (aq.zcenter[pylayer] - aq.zcenter[pylayer+1]) *\
                        (headvec[pylayer]-headvec[pylayer+1])
        return h
    def dischargeNormVector(self,x,y,theta,aq=None):
        '''Returns normal dischargevector in each layer in direction theta '''
        if aq == None: aq = self.aq.findAquiferData(x,y)
        [Qx,Qy] = self.dischargeCollection(x,y,aq)
        Qr = Qx * cos(theta) + Qy * sin(theta);
        return Qr
    def integratedLeftDischarge(self,xylist,offset=1e-6,tol=1e-3,divmax=10):
        rv = []
        N = len(xylist) - 1
        X1 = -1.0 + offset
        X2 = 1.0 - offset
        for i in range(N):
            Qn = []
            z1 = xylist[i][0] + xylist[i][1]*1j
            z2 = xylist[i+1][0] + xylist[i+1][1]*1j
            L = abs(z2-z1)
            theta = arctan2( z2.imag-z1.imag, z2.real-z1.real ) + pi/2
            aq = self.aq.findAquiferData(z1.real,z1.imag)
            for n in range(aq.Naquifers):
                def func(X):
                    Z = X + offset*1j
                    z = 0.5 * ( Z * (z2-z1) + (z1+z2) )
                    x,y = z.real,z.imag
                    return self.dischargeNormVector(x,y,theta)[n]
                Qn.append(romberg(func,X1,X2,tol=tol,divmax=divmax)*L/2)
            rv.append(array(Qn))
        return array(rv)
    def dischargeNormVectorOld(self,x,y,theta,aq=None):
        '''Returns normal dischargevector in each layer in direction theta '''
        if aq == None: aq = self.aq.findAquiferData(x,y)
        [Qx,Qy] = self.dischargeVector(x,y,aq)
        Qr = Qx * cos(theta) + Qy * sin(theta);
        return Qr
    def dischargeNormVectorCol(self,x,y,theta,aq=None):
        '''Returns normal dischargevector in each layer in direction theta as column vector'''
        if aq == None: aq = self.aq.findAquiferData(x,y)
        [Qx,Qy] = self.dischargeVector(x,y,aq)
        Qr = Qx * cos(theta) + Qy * sin(theta);
        # return transpose(Qr[newaxis,:])
        return Qr[:,newaxis]
    def dischargeNormInLayer(self,aq,pylayer,x,y,theta):
        '''Returns discharge in direction theta in aquifer aq and layer pylayer'''
        [Qx,Qy] = self.dischargeVector(x,y,aq)
        Qr = Qx[pylayer] * cos(theta) + Qy[pylayer] * sin(theta);
        return Qr
    def totalDischargeFromInf(self):
        rv = 0.0
        for e in self.elementList:
            rv = rv + e.totalDischarge()
        return rv
    def dischargeNormBelowZ(self,x,y,z,theta,aq=None):
        if aq == None: aq = self.aq.findAquiferData(x,y)        
        Qn = self.dischargeNormVector(x,y,theta,aq)
        pyLayer = aq.inWhichPyLayer(z)
        if pyLayer < 0: print 'z in leaky layer in dischargeNormBelowZ'
        Qntot = 0.0
        #  if pyLayer < aq.Naquifers - 1: Qntot = cumsum(Qn[pyLayer+1:])[0]  # could have done differently, but want float not array
        # I think this is bug and should be
        if pyLayer < aq.Naquifers - 1: Qntot = sum(Qn[pyLayer+1:])
        Qntot = Qntot + (z-aq.zb[pyLayer]) * Qn[pyLayer] / aq.H[pyLayer]
        return Qntot
    def zForGivenNormalDischarge(self,zbegin,xmin,ymin,xplus,yplus,theta,Qngiven,aq=None):
        '''Returns z below which normal discharge is Qngiven in direction theta; search is started at zbegin '''
        x = xplus; y = yplus
        if aq == None: aq = self.aq.findAquiferData(x,y)
        Qn = self.dischargeNormVector(x,y,theta,aq)
        pyLayer = aq.inWhichPyLayer(zbegin)
        if pyLayer < 0:  # in resistance layer
            pyLayer = -pyLayer
            zbegin = aq.zt[pyLayer-1]
        if pyLayer == -9999:  # above top
            pyLayer = 0
            zbegin = aq.zt[0]
        if pyLayer == 9999:  # below bottom
            pyLayer = aq.Naquifers - 1
            zbegin = aq.zb[0]
        # Compute discharge below zbegin
        Qnbegin = self.dischargeNormBelowZ(x,y,zbegin,theta)
        # Compute vertical flow according to Strack, 1995, WRR 31(12), p.3016
        Qvert = Qngiven - Qnbegin
        # Cumulative discharge at x,y
        Qncum = zeros(aq.Naquifers+1,'d')
        Qncum[1:] = cumsum(Qn)  # Summing from the bottom up
        Qncum = sum(Qn) - Qncum  # 0 is begin top of aquifer 0, 1 is bot of aquifer 0 top of aquifer 1, etc.
        # Set ifound and isign
        ifound = 0; isign = 1.0
        # Find new elevation
        if Qvert > 0 :  # Search up
            print 'search up'
            if Qngiven < Qnbegin: isign = -1.0
            for i in range(pyLayer,-1,-1):
                if isign * Qncum[i] > isign * Qngiven:  # Found layer
                    znew = aq.zb[i] + abs( ( Qngiven - Qncum[i+1] ) / Qn[i] ) * aq.H[i]
                    ifound = 1
                    break
        else:  # Search down
            print 'search down'
            if Qngiven < Qnbegin: isign = -1.0
            for i in range(pyLayer,aq.Naquifers):
                if isign * Qncum[i+1] > isign * Qngiven:  # Found layer
                    znew = aq.zb[i] + abs( ( Qngiven - Qncum[i+1] ) / Qn[i] ) * aq.H[i]
                    ifound = 1
                    break
        if ifound == 0:  # Does not flow to xplus,yplus, but stays in xmin,ymin
            isign = 1.0
            x = xmin; y = ymin
            aq = self.aq.findAquiferData(x,y)
            Qn = self.dischargeNormVector(x,y,theta,aq)
            pyLayer = aq.inWhichPyLayer(zbegin)
            # Cumulative discharge at x,y
            Qncum = zeros(aq.Naquifers+1,'d')
            Qncum[1:] = cumsum(Qn)  # Summing from the bottom up
            Qncum = sum(Qn) - Qncum  # 0 is begin top of aquifer 0, 1 is bot of aquifer 0 top of aquifer 1, etc.
            if Qvert > 0:  # Search up
                print 'search up on minus side'
                if Qngiven < Qncum[pyLayer]: isign = -1.0
                if pyLayer == 0: print 'Error: cannot search up'
                for i in range(pyLayer-1,-1,-1):
                    if isign * Qncum[i] > isign * Qngiven:  # Found layer
                        znew = aq.zb[i] + abs( ( Qngiven - Qncum[i+1] ) / Qn[i] ) * aq.H[i]
                        ifound = 1
                        break
            else:  # search down
                print 'search down on minus side'
                if Qngiven < Qncum[pyLayer+1]: isign = -1.0
                if pyLayer == aq.Naquifers-1: print 'Error: cannot search down'
                for i in range(pyLayer,aq.Naquifers): # Start layer below
                    if isign * Qncum[i+1] > isign * Qngiven:  # Found layer
                        znew = aq.zb[i] + abs( ( Qngiven - Qncum[i+1] ) / Qn[i] ) * aq.H[i]
                        ifound = 1
                        break
        if ifound == 0: print 'Could not find new z in zForGivenNormalDischarge'
        return [x,y,znew]
##    def solve(self,reInitializeAllElements=1):
##        '''Compute solution. If elements don't have to be initialized, set parameter to 0'''
##        print 'starting solve'
##        if reInitializeAllElements:
##            self.aq.setCoefs()
##            for p in self.aq.inhomList:
##                p.setCoefs()
##            for e in self.elementList:
##                e.setCoefs()
##        # Rebuild element collection
##        self.collectionDict = {}
##        for e in self.elementList:
##            e.addElementToCollection(self)
##        matrix = []
##        Nel = len(self.elementList)
##        print 'Number of elements: ',Nel
##        print 'Percent progress: ',
##        imilestone = (Nel-1)*arange(0,11,1,'i')/10
##        icount = 0; iprog = 0
##        eq_list = []
##        for e in self.elementList:
##            rows = e.getMatrixRows(self.elementList)
##            eq_list.append( len(rows) )
##            matrix.extend( rows )
##            if icount == imilestone[iprog]:
##                print int(10.0*iprog),
##                iprog = iprog + 1
##            icount = icount + 1
##        print eq_list
##        matrix = array(matrix)
##        size = shape(matrix)
##        if size[0] == 0:
##            print 'No unknown parameters'
##            return
##        rhs = take( matrix, range(size[1]-1,size[1]),1)
##        rhs = transpose(rhs)[0]   # Transfer to row vector
##        matrix = take( matrix, range(0,size[1]-1), 1)
##        size = shape(matrix)
##        print 'size of matrix '+str(size)
##        # print matrix
##        if size[0] > size[1]:
##            print 'more rows than columns: attempting least squares solution'
##            xsol = linalg.linear_least_squares(matrix,rhs)
##            print 'rank ',xsol[2]
##            xsol = xsol[0]
##        elif size[0] == size[1]:
##            xsol = transpose( linalg.solve(matrix,rhs) )
##        else:
##            print 'Error: more unknowns than equations'
##            return
##        icount = 0
##        for e in self.elementList:
##            icount = e.takeParameters(xsol,icount)
##        print 'solution complete'
    def solveIter(self,Niter,alpha=1.0,reInitializeAllElements=1):
        if reInitializeAllElements:
            self.solve()
        for i in range(Niter):
            print 'iteration ',i
            for e in self.elementList:
                e.iterationStep(alpha)
            self.solve(0)
        for e in self.elementList:
            e.iterationStep(alpha)
        print 'Iterations complete'
##    def solveIter(self,tol=1e-8):
##        '''Compute solution iteratively. STILL EXPERIMENTAL !!!
##        Designed to solve every element seperately'''
####        if reInitializeAllElements:
####            self.aq.setCoefs()
####            for p in self.aq.inhomList:
####                p.setCoefs()
####            for e in self.elementList:
####                e.setCoefs()
##        sumOld = 0.0; error = 1.0
##        while error > tol:
##            sumNew = 0.0
##            for e in self.elementList:
##                matrix = e.getMatrixRows(e.elementList)
##                if shape(matrix)[0] == 0: continue  # No unknowns; also not taken into account for convergence (nice!)
##                matrix = array(matrix); size = shape(matrix)
##                rhs = take( matrix, range(size[1]-1,size[1]),1)
##                rhs = transpose(rhs)[0]   # Transfer to row vector
##                matrix = take( matrix, range(0,size[1]-1), 1)
##                xsol = transpose( linalg.solve_linear_equations(matrix,rhs) )
##                icount = 0
##                icount = e.takeParameters(xsol,icount)
##                sumNew = sumNew + e.parSum()
##            error = abs( (sumNew-sumOld)/sumNew ); sumOld = sumNew
##            print 'error: '+str(error)
##        print 'iteration complete'
    def check(self):
        '''Prints data to screen to check succes of solution'''
        for e in self.elementList:
            e.check()
    def inWhichPyLayer(self,x,y,z,aq=None):
        if aq == None: aq = self.aq.findAquiferData(x,y)
        return aq.inWhichPyLayer(z)
    def inWhichPyLayers(self,x,y,z,aq=None):
        pylayers = []
        for zz in z:
            if aq == None: aq = self.aq.findAquiferData(x,y)
            pylayers.append( aq.inWhichPyLayer(zz) )
        return pylayers
    def qzTop(self,x,y):
        qz = 0.0
        for e in self.elementList:
            qz = qz + e.qzTop(x,y)
        return qz
    def velocity(self,x,y,z):
        head = self.headVector(x,y)
        [disx, disy] = self.dischargeCollection(x,y)
        aqdata = self.aq.findAquiferData(x,y)
        pyLayer = self.inWhichPyLayer(x,y,z,aqdata)
        assert pyLayer != -9999 and pyLayer != 9999, 'TimML error: (x,y,z) outside aquifer '+str((x,y,z))
        if pyLayer >= 0:  # In aquifer
            vx = disx[pyLayer] / ( aqdata.H[pyLayer] * aqdata.n[pyLayer] )
            vy = disy[pyLayer] / ( aqdata.H[pyLayer] * aqdata.n[pyLayer] )
            if pyLayer > 0:
                vztop = ( head[pyLayer] - head[pyLayer-1] ) / ( aqdata.c[pyLayer] * aqdata.n[pyLayer] )
            else:
                if aqdata.type == aqdata.conf:
                    vztop = self.qzTop(x,y) / aqdata.n[pyLayer]
                elif aqdata.type == aqdata.semi:
                    vztop = ( head[0] - aqdata.hstar ) / ( aqdata.c[0] * aqdata.n[0] )
            if pyLayer < aqdata.Naquifers-1:
                vzbot = ( head[pyLayer+1] - head[pyLayer] ) / ( aqdata.c[pyLayer+1] * aqdata.n[pyLayer] )
            else:
                vzbot = 0.0
            vz = (z - aqdata.zb[pyLayer]) * (vztop - vzbot) / aqdata.H[pyLayer] + vzbot
        else:  # In leaky layer
            vx = 0.0
            vy = 0.0
            vz = ( head[-pyLayer] - head[-pyLayer-1] ) / ( aqdata.c[-pyLayer] * aqdata.nll[-pyLayer] ) 
        return array([vx,vy,vz])
    def velocity2(self,x,y,z,aqdata,pyLayer):
        head = self.headVector(x,y,aqdata)
        if pyLayer >= 0:  # In aquifer
            [disx, disy] = self.dischargeCollection(x,y,aqdata)
            vx = disx[pyLayer] / ( aqdata.H[pyLayer] * aqdata.n[pyLayer] )
            vy = disy[pyLayer] / ( aqdata.H[pyLayer] * aqdata.n[pyLayer] )
            if pyLayer > 0:
                vztop = ( head[pyLayer] - head[pyLayer-1] ) / ( aqdata.c[pyLayer] * aqdata.n[pyLayer] )
            else:
                if aqdata.type == aqdata.conf:
                    vztop = self.qzTop(x,y) / aqdata.n[pyLayer]
                elif aqdata.type == aqdata.semi:
                    vztop = ( head[0] - aqdata.hstar ) / ( aqdata.c[0] * aqdata.n[0] )
            if pyLayer < aqdata.Naquifers-1:
                vzbot = ( head[pyLayer+1] - head[pyLayer] ) / ( aqdata.c[pyLayer+1] * aqdata.n[pyLayer] )
            else:
                vzbot = 0.0
            vz = (z - aqdata.zb[pyLayer]) * (vztop - vzbot) / aqdata.H[pyLayer] + vzbot
        else:  # In leaky layer
            vx = 0.0
            vy = 0.0
            vz = ( head[-pyLayer] - head[-pyLayer-1] ) / ( aqdata.c[-pyLayer] * aqdata.nll[-pyLayer] ) 
        return array([vx,vy,vz])
    # Function needed to call traceline from object
    def traceline(self,xstart,ystart,zstart,stepin,tmax=1e30,maxsteps=10,tstart=0.0,window=[-1e30,-1e30,1e30,1e30],\
              labfrac = 2.0, Hfrac = 2.0):
        a,b,c,d = traceline(self,xstart,ystart,zstart,stepin,tmax,maxsteps,tstart,window,labfrac,Hfrac)
        return a,b,c,d

    def solve(self,reInitializeAllElements=1,doIterations=False,maxIter=5,storematrix=False,conditionnumber=False):
        '''Compute solution
        reInitializeAllElements: Initializes all elements and aquifers (default is 1)
        doIterations: Iterates for non-linear conditions if present
        maxIter: Maximum number of iterations. If reached before convergence a message is written to the screen
        storematrix: Logical to to indicate whether the matrix should remain stored. Only useful during development
        '''
        self.matrix = 0; self.rhs = 0; self.xsol = 0; self.eqlist = 0; self.eqcumlist = 0  # If not created (no unknowns) they need to exist to be deleted
        print 'Starting solve'
        newsolution = self.solveNonLinear( 1, reInitializeAllElements )
        lakechange = False  # Need to start with false in case there are no lakes
        if doIterations:
            i = 0
            while (newsolution or lakechange) and i < maxIter:
                print 'Iteration ',i+1
                if len(self.lakeList) == 0:
                    newsolution = self.solveNonLinear( start=0 ) # Only reinitializes if start = 1
                else:
                    # Make changes for lakes that percolate
                    lakechange = False
                    for lake in self.lakeList:
                        change = lake.change_lake()
                        if change: lakechange = True
                    if lakechange: solution = self.solveNonLinear( 1, reInitializeAllElements ) #Cannot modify existing matrix; need to recompile
                    # Check for percolating line-sinks
                    newsolution = self.solveNonLinear( start=0 )
                i+=1
            if (newsolution or lakechange) and i == maxIter:
                print 'Warning: Maximum number of iterations reached before convergence'
            else:
                print 'Iterations complete'
        if conditionnumber:
            u,s,vh = linalg.svd(self.matrix)
            print 'Condition number: ',s[0]/s[-1]
        if not storematrix:
            del self.matrix; del self.rhs; del self.xsol; del self.eqlist; del self.eqcumlist
        print 'Solution complete'
        return
        
    def solveNonLinear(self,start=1,reInitializeAllElements=1):
        '''Do one iteration step
        When start=1 it solves in the old fashioned way
        When start=0 is will use the stored matrix and only adjust the rows that have changed
        '''
        newsolution_computed = False
        if start == 1:  # Generate matrix
            self.collectionDict = {}
            for e in self.elementList: e.addElementToCollection(self)
            if reInitializeAllElements == 1:
                # Initialize elements and rebuild collection
                self.aq.setCoefs()
                for p in self.aq.inhomList:
                    p.setCoefs()
                for e in self.elementList:
                    e.setCoefs()
            matrix = []
            Nel = len(self.elementList)
            print 'Number of elements: ',Nel
            print 'Percent progress: ',
	    stdout.flush()
            imilestone = (Nel-1)*arange(0,11,1,'i')/10
            icount = 0; iprog = 0
            eq_list = []
            # Compile matrix
            for e in self.elementList:
                rows = e.getMatrixRows(self.elementList)
                eq_list.append( len(rows) )
                matrix.extend( rows )
                if icount == imilestone[iprog]:
                    print int(10.0*iprog),
		    stdout.flush()
                    iprog = iprog + 1
                icount = icount + 1
            print ' ' # Just to end the printing in a row
            eqcumlist = cumsum(eq_list)
            self.eqlist = eq_list
            self.eqcumlist = eqcumlist - 1
            #print 'matrix ',matrix
            matrix = array(matrix)
            size = shape(matrix)
            if size[0] == 0:
                print 'No unknown parameters'
                return newsolution_computed
            rhs = take( matrix, range(size[1]-1,size[1]),1)
            rhs = transpose(rhs)[0]   # Transfer to row vector
            matrix = take( matrix, range(0,size[1]-1), 1)
        else:  # Matrix is already compiled, just need correction
            changed = False; ichange = 0
            matrix = self.matrix; rhs = self.rhs; xsol = self.xsol; eqcumlist = self.eqcumlist
            rhs = rhs - dot( matrix, xsol )
            for e,iel in zip( self.elementList, range(len(self.elementList)) ):
                newrows = e.getMatrixRows_nonlinear(self.elementList)
                if newrows != None:
                    ichange += 1
                    if not changed: changed = True
                    if len( newrows ) > 1:
                        print 'Non-linear iteration for element with multiple unknown strengths not implemented'
                        return
                    else:
                        matrix[ eqcumlist[iel],: ] = array( newrows[0][:-1] )
                        rhs[ eqcumlist[iel] ] = newrows[0][-1]
            print 'Number of elements changed (such as percolating line-sinks): ',ichange
            if changed == False: return newsolution_computed
        size = shape(matrix)
        print 'size of matrix '+str(size)
        if size[0] == size[1]:
            xsol = transpose( linalg.solve(matrix,rhs) )
        else:
            print 'Error: different number of unknowns than equations'
            return
        icount = 0
        for e in self.elementList:
            icount = e.takeParameters(xsol,icount)
        self.matrix = matrix
        self.rhs = rhs
        self.xsol = xsol       
        newsolution_computed = True
        return newsolution_computed
    
    def press(self,event):
        print 'Hello ml'
        #if event.inaxes is None: return
        #print 'event.button ',event.button
        #if event.button != 3: return
        #ax = gcf().axes[0]
        #x1,x2,y1,y2 = ax.axis()
        #step = (x2 - x1) / 100.0
        #if not self.ml.trace.forward: step = -step
        #tmax = self.ml.trace.tmax
        #zbegin = self.ml.trace.zbegin
        #xsec = self.ml.trace.xsec
        #Npoints = len(zbegin)
        #timtracelines( self.ml, Npoints*[event.xdata], Npoints*[event.ydata], zbegin, step, tmax=tmax, Nmax=200, window=(x1,y1,x2,y2), xsec=xsec )

def Model3D(z=[1,0.5,0],kh=1.0,kzoverkh=1.0):
    z=array(z,'d')
    Naq = len(z) - 1
    zt = z[0:-1]
    zb = z[1:]
    H = zt - zb
    if iterable( kh ):
        assert len(kh) == Naq, "TimML Input error: length of kh must equal number of layers" 
        k = array(kh,'d')
    else:
        k = kh * ones(Naq)
    if iterable( kzoverkh ):
        assert len(kzoverkh) == Naq, "TimML Input error: length of kh must equal number of layers" 
        kzoverkh = array(kzoverkh,'d')
    else:
        kzoverkh = kzoverkh * ones(Naq)
    c = 0.5 * H[:-1] / ( kzoverkh[:-1] * k[:-1] ) + 0.5 * H[1:] / ( kzoverkh[1:] * k[1:] )     
    ml = Model( k=list(k), zb=list(zb), zt=list(zt), c=list(c) )
    return ml

    

 

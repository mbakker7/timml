'''
mlcircinhom.py contains the CircleInhom class
This file is part of the TimML library and is distributed under
the GNU LPGL. See the TimML.py file for more details.
(c) Mark Bakker, 2002-2007
'''

import scipy.special
from element import *
from numpy import *


class CircleInhom(Element):
    '''CircleInhom class
    Note that aquiferparent doesn't have a meaning for this element
    All attributes from element.
    '''
    Rconvsq = 7 * 7

    def __init__(self, modelParent, order, aqin, aqout, label=None):
        Element.__init__(self, modelParent)
        self.order = order
        self.aqin = aqin
        self.aqout = aqout
        assert aqin.Naquifers == aqout.Naquifers, "TimML Input error: Number of aquifers inside and outside must be equal"
        self.label = label
        self.type = 'circleinhom'
        self.setCoefs()
        self.modelParent.addElement(self)

    def __repr__(self):
        return 'CircInhom xc,yc,order ' + str((self.xc, self.yc, self.order)) + ' with inside aquifer data ' + str(self.aqin)

    def setCoefs(self):
        self.xc = self.aqin.xc
        self.yc = self.aqin.yc
        self.R = self.aqin.R
        self.Rsq = self.aqin.Rsq
        self.paramin = zeros((2 * self.order + 1, self.aqin.Naquifers), 'd')
        self.paramout = zeros((2 * self.order + 1, self.aqin.Naquifers), 'd')
        self.Nparamin = (2 * self.order + 1) * self.aqin.Naquifers
        self.Nparamout = self.Nparamin
        # Compute control points
        self.Ncp = 2 * self.order + 1
        self.thetacp = arange(0, 2 * pi - 0.01, 2 * pi / self.Ncp)
        self.xcp = self.xc + self.R * cos(self.thetacp)
        self.ycp = self.yc + self.R * sin(self.thetacp)
        # Compute values on the edge (to scale functions)
        # values on inside edge
        self.matRin = zeros((self.order + 1, self.aqin.Naquifers), 'd')
        rolab = self.R / self.aqin.lab  # vector with all R/labin values
        if self.aqin.type == self.aqin.conf:
            for p in range(self.order + 1):
                self.matRin[p, 0] = self.R**p  # first value in row is R**p
                self.matRin[p, 1:] = scipy.special.iv(
                    p, rolab)  # other values are Ip(R/lab)
        elif self.aqin.type == self.aqin.semi:
            for p in range(self.order + 1):
                self.matRin[p, :] = scipy.special.iv(
                    p, rolab)  # all values are Ip(R/lab)
        self.matRin = 1.0 / self.matRin
        # values on outside edge
        self.matRout = zeros((self.order + 1, self.aqout.Naquifers), 'd')
        rolab = self.R / self.aqout.lab  # vector with all R/labout values
        if self.aqout.type == self.aqin.conf:
            for p in range(self.order + 1):
                # first value in row is 1/r**p
                self.matRout[p, 0] = 1. / self.R**(p)
                self.matRout[p, 1:] = scipy.special.kn(
                    p, rolab)  # other values are Kp
            # first term is logarithm, so scale by log(self.R+1)to make sure that self.R can be 1
            self.matRout[0, 0] = log(self.R + 1.0) / (2.0 * pi)
        elif self.aqout.type == self.aqout.semi:
            for p in range(self.order + 1):
                self.matRout[p, :] = scipy.special.kn(
                    p, rolab)  # other values are Kp
        self.matRout = 1.0 / self.matRout

    def potentialInfluence(self, aq, x, y):
        r = sqrt((x - self.xc)**2 + (y - self.yc)**2)
        theta = arctan2((y - self.yc), (x - self.xc))
        mat = zeros((self.order + 1, aq.Naquifers),
                    'd')  # matrix to store values
     # Inside Cylinder
        if aq == self.aqin:
            rolab = r / aq.lab  # vector with all r/lab values
            if aq.type == aq.conf:  # Confined aquifer
                for p in range(self.order + 1):
                    mat[p, 0] = r**p  # first value in row is r**p
                    mat[p, 1:] = scipy.special.iv(
                        p, rolab)  # other values are Ip
            elif aq.type == aq.semi:  # Semi-confined aquifer
                for p in range(self.order + 1):
                    mat[p, :] = scipy.special.iv(p, rolab)  # all values are Ip
            mat = mat * self.matRin
     # In aquifer outside cylinder
        elif aq == self.aqout:
            rolab = r / aq.lab  # vector with all r/lab values
            if aq.type == aq.conf:  # Confined aquifer
                for p in range(self.order + 1):
                    mat[p, 0] = 1. / r**(p)  # first value in row is 1/r**p
                    mat[p, 1:] = scipy.special.kn(
                        p, rolab)  # other values are Kp
                mat[0, 0] = log(r) / (2.0 * pi)  # logarithm on outside
            if aq.type == aq.semi:  # Semi-confined aquifer
                for p in range(self.order + 1):
                    mat[p, :] = scipy.special.kn(p, rolab)  # all values are Kp
            mat = mat * self.matRout
     # Somewhere else, assume outside
        else:
            if self.aqout.type == self.aqout.conf and aq.type == aq.conf:
                for p in range(self.order + 1):
                    mat[p, 0] = 1. / r**(p)  # first value in row is 1/r**p
                mat[0, 0] = log(r) / (2.0 * pi)
                mat = mat * self.matRout
     # Store values in return value
        rv = zeros((2 * self.order + 1, aq.Naquifers), 'd')
        rv[0, :] = mat[0, :]  # First row only has cosine part; sin(0)=0
        for p in range(1, self.order + 1):
            rv[2 * p - 1, :] = mat[p, :] * cos(p * theta)
            rv[2 * p, :] = mat[p, :] * sin(p * theta)
        return rv

    def potentialInfluenceInLayer(self, aq, pylayer, x, y):
        '''Returns PotentialInfluence function in aquifer aq in pylayer as array (1 value per parameter)
        Needs to be overloaded because there is no parameter outside the functions
        Needs to be modified for inside and outside'''
        potInf = self.potentialInfluence(aq, x, y)
        rv = zeros(0, 'd')
        for p in potInf:
            rv = hstack((rv, p * aq.eigvec[pylayer, :]))
        if aq == self.aqin:
            rv = hstack((rv, zeros(self.Nparamout, 'd')))
        else:
            rv = hstack((zeros(self.Nparamin, 'd'), rv))
        return rv

    def potentialInfluenceAllLayers(self, aq, pylayer, x, y):
        '''Returns PotentialInfluence function in aquifer aq in all layers as an array'''
        potInf = self.potentialInfluence(aq, x, y)
        rv = zeros((aq.Naquifers, 0), 'd')
        for p in potInf:
            rv = hstack((rv, p * aq.eigvec))
        if aq == self.aqin:
            rv = hstack((rv, zeros((aq.Naquifers, self.Nparamout), 'd')))
        else:
            rv = hstack((zeros((aq.Naquifers, self.Nparamin), 'd'), rv))
        return rv

    def potentialInfluenceSpecLayers(self, aq, pylayer, x, y):
        '''Returns PotentialInfluence function in aquifer aq in all layers as an array'''
        potInf = self.potentialInfluence(aq, x, y)
        pylen = len(pylayer)
        rv = zeros((pylen, 0), 'd')
        eigvec = take(aq.eigvec, pylayer, 0)
        for p in potInf:
            rv = hstack((rv, p * eigvec))
        if aq == self.aqin:
            rv = hstack((rv, zeros((pylen, self.Nparamout), 'd')))
        else:
            rv = hstack((zeros((pylen, self.Nparamin), 'd'), rv))
        return rv

    def potentialContribution(self, aq, x, y):
        '''Returns array of potentialContribution. Needs to be overloaded cause there is inside and outside'''
        if aq == self.aqin:
            potInf = sum(self.paramin * self.potentialInfluence(aq, x, y), 0)
        else:
            potInf = sum(self.paramout * self.potentialInfluence(aq, x, y), 0)
        return potInf

    def dischargeInfluence(self, aq, x, y):
        r = sqrt((x - self.xc)**2 + (y - self.yc)**2)
        theta = arctan2((y - self.yc), (x - self.xc))
        mat = zeros((self.order + 1, aq.Naquifers),
                    'd')  # matrix to store values
        if aq == self.aqin:  # Inside cylinder
            rolab = r / aq.lab  # vector with all r/lab values
            if aq.type == aq.conf:  # Confined aquifer
                for p in range(self.order + 1):
                    if p == 0:
                        mat[0, 0] = 0.0
                    else:
                        # first value in row is pr**(p-1)
                        mat[p, 0] = - float(p) * r**(p - 1)
                    mat[p, 1:] = - (scipy.special.iv(p - 1, rolab) +
                                    scipy.special.iv(p + 1, rolab)) / (2 * aq.lab)
            elif aq.type == aq.semi:
                for p in range(self.order + 1):
                    mat[p, :] = - (scipy.special.iv(p - 1, rolab) +
                                   scipy.special.iv(p + 1, rolab)) / (2 * aq.lab)
            mat = mat * self.matRin
        elif aq == self.aqout:  # In aquifer outside cylinder
            rolab = r / aq.lab  # vector with all r/lab values
            if aq.type == aq.conf:
                for p in range(self.order + 1):
                    mat[p, 0] = p / r**(p + 1)  # first value in row is 1/r**p
                    mat[p, 1:] = (scipy.special.kn(p - 1, rolab) +
                                  scipy.special.kn(p + 1, rolab)) / (2 * aq.lab)
                mat[0, 0] = - 1.0 / (2.0 * pi * r)  # log
            elif aq.type == aq.semi:
                for p in range(self.order + 1):
                    mat[p, :] = (scipy.special.kn(p - 1, rolab) +
                                 scipy.special.kn(p + 1, rolab)) / (2 * aq.lab)
            mat = mat * self.matRout
        else:  # Somewhere else
            if self.aqout.type == self.aqout.conf and aq.type == aq.conf:
                mat[0, 0] = - 1.0 / (2.0 * pi * r)
                for p in range(1, self.order + 1):
                    mat[p, 0] = p / r**(p + 1)  # first value in row is 1/r**p
                mat = mat * self.matRout
        # Store values in return value
        rvrad = zeros((2 * self.order + 1, aq.Naquifers), 'd')
        rvrad[0, :] = mat[0, :]  # First row only has cosine part; sin(0)=0
        for p in range(1, self.order + 1):
            rvrad[2 * p - 1, :] = mat[p, :] * cos(p * theta)
            rvrad[2 * p, :] = mat[p, :] * sin(p * theta)
    # Tangential part
        mat = zeros((self.order + 1, aq.Naquifers),
                    'd')  # matrix to store values
        if aq == self.aqin:  # Inside cylinder
            rolab = r / aq.lab  # vector with all r/lab values
            if aq.type == aq.conf:
                # all values for p=0 are zero (no dependence on theta)
                for p in range(1, self.order + 1):
                    # first value in row is pr**(p-1)
                    mat[p, 0] = float(p) * r**(p - 1)
                    mat[p, 1:] = float(p) * scipy.special.iv(p, rolab) / r
            elif aq.type == aq.semi:
                # all values for p=0 are zero (no dependence on theta)
                for p in range(1, self.order + 1):
                    mat[p, :] = float(p) * scipy.special.iv(p, rolab) / r
            mat = mat * self.matRin
        elif aq == self.aqout:  # In aquifer outside cylinder
            rolab = r / aq.lab  # vector with all r/lab values
            if aq.type == aq.conf:
                # all values for p=0 are zero (no dependence on theta)
                for p in range(1, self.order + 1):
                    # first value in row is 1/r**p
                    mat[p, 0] = float(p) / r**(p + 1)
                    mat[p, 1:] = float(p) * scipy.special.kn(p, rolab) / r
            elif aq.type == aq.semi:
                # all values for p=0 are zero (no dependence on theta)
                for p in range(1, self.order + 1):
                    mat[p, :] = float(p) * scipy.special.kn(p, rolab) / r
            mat = mat * self.matRout
        else:  # Somewhere else
            if self.aqout.type == self.aqout.conf and aq.type == aq.conf:
                for p in range(1, self.order + 1):
                    mat[p, 0] = float(p) / r**(p + 1)
            mat = mat * self.matRout
        # Store values in return value
        rvtan = zeros((2 * self.order + 1, aq.Naquifers), 'd')
        for p in range(1, self.order + 1):
            rvtan[2 * p - 1, :] = mat[p, :] * sin(p * theta)
            rvtan[2 * p, :] = - mat[p, :] * cos(p * theta)
        qx = rvrad * cos(theta) - rvtan * sin(theta)
        qy = rvrad * sin(theta) + rvtan * cos(theta)
        return [qx, qy]

    def dischargeInfluenceRad(self, aq, x, y):
        r = sqrt((x - self.xc)**2 + (y - self.yc)**2)
        theta = arctan2((y - self.yc), (x - self.xc))
        mat = zeros((self.order + 1, aq.Naquifers),
                    'd')  # matrix to store values
        if aq == self.aqin:  # Inside cylinder
            rolab = r / aq.lab  # vector with all r/lab values
            for p in range(self.order + 1):
                if p == 0:
                    mat[0, 0] = 0.0
                else:
                    # first value in row is pr**(p-1)
                    mat[p, 0] = - float(p) * r**(p - 1)
                mat[p, 1:] = - (scipy.special.iv(p - 1, rolab) +
                                scipy.special.iv(p + 1, rolab)) / (2 * aq.lab)
            mat = mat * self.matRin
        elif aq == self.aqout:  # In aquifer outside cylinder
            rolab = r / aq.lab  # vector with all r/lab values
            for p in range(self.order + 1):
                mat[p, 0] = p / r**(p + 1)  # first value in row is 1/r**p
                mat[p, 1:] = (scipy.special.kn(p - 1, rolab) +
                              scipy.special.kn(p + 1, rolab)) / (2 * aq.lab)
            mat[0, 0] = 0.0  # no constant on outside
            mat = mat * self.matRout
        else:  # Somewhere else
            for p in range(self.order + 1):
                mat[p, 0] = p / r**(p + 1)  # first value in row is 1/r**p
            mat[0, 0] = 0.0
            mat = mat * self.matRout
        # Store values in return value
        rvrad = zeros((2 * self.order + 1, aq.Naquifers), 'd')
        rvrad[0, :] = mat[0, :]  # First row only has cosine part; sin(0)=0
        for p in range(1, self.order + 1):
            rvrad[2 * p - 1, :] = mat[p, :] * cos(p * theta)
            rvrad[2 * p, :] = mat[p, :] * sin(p * theta)
        return rvrad

    def dischargeInfluenceRadInLayer(self, aq, pylayer, x, y):
        '''Returns dischargeInfluenceRadInLayer function in aquifer aq in pylayer as list (1 value per parameter)
        Needs to be overloaded because there is no parameter outside the functions
        Needs to be modified for inside and outside'''
        disInf = self.dischargeInfluenceRad(aq, x, y)
        rv = []
        for d in disInf:
            rv = rv + list(d * aq.eigvec[pylayer, :])
        if aq != self.aqin:
            dump = rv.pop(0)  # Get rid of very first term (always zero)
        return rv

    def dischargeInfluenceRadAllLayers(self, aq, dumlayer, x, y):
        '''Returns dischargeInfluenceRadAllLayers function in aquifer aq as an array
        Needs to be overloaded because there is no parameter outside the functions
        Needs to be modified for inside and outside'''
        disInf = self.dischargeInfluenceRad(aq, x, y)
        rv = zeros((aq.Naquifers, 0), 'd')
        for d in disInf:
            rv = hstack((rv, d * aq.eigvec))
        if aq != self.aqin:
            rv = rv[:, 1:]  # Take off first column (parameter always zero)
        return rv

    def dischargeInfluenceAllLayers(self, aq, dumlayer, x, y):
        '''Returns dischargeInfluenceAllLayers function in aquifer aq as an array
        Needs to be overloaded because there is no parameter outside the functions
        Needs to be modified for inside and outside'''
        [disx, disy] = self.dischargeInfluence(aq, x, y)
        rvx = zeros((aq.Naquifers, 0), 'd')
        rvy = zeros((aq.Naquifers, 0), 'd')
        for d in disx:
            rvx = hstack((rvx, d * aq.eigvec))
        for d in disy:
            rvy = hstack((rvy, d * aq.eigvec))
        if aq == self.aqin:
            rvx = hstack((rvx, zeros((aq.Naquifers, self.Nparamout), 'd')))
            rvy = hstack((rvy, zeros((aq.Naquifers, self.Nparamout), 'd')))
        else:
            rvx = hstack((zeros((aq.Naquifers, self.Nparamin), 'd'), rvx))
            rvy = hstack((zeros((aq.Naquifers, self.Nparamin), 'd'), rvy))
        return [rvx, rvy]

    def dischargeInfluenceInLayer(self, aq, pylayer, x, y):
        '''Returns dischargeInfluence in pylayer, modified for paramin and paramout'''
        [disx, disy] = self.dischargeInfluence(aq, x, y)
        rvx = disx * aq.eigvec[pylayer, :]
        rvy = disy * aq.eigvec[pylayer, :]
        if aq == self.aqin:
            rvx = hstack((rvx.ravel(), zeros(self.Nparamout, 'd')))
            rvy = hstack((rvy.ravel(), zeros(self.Nparamout, 'd')))
        else:
            rvx = hstack((zeros(self.Nparamin, 'd'), rvx.ravel()))
            rvy = hstack((zeros(self.Nparamin, 'd'), rvy.ravel()))
        return [rvx, rvy]

    def dischargeInfluenceSpecLayers(self, aq, pylayer, x, y):
        '''Returns dischargeInfluenceAllLayers function in aquifer aq as an array
        Needs to be overloaded because there is no parameter outside the functions
        Needs to be modified for inside and outside'''
        [disx, disy] = self.dischargeInfluence(aq, x, y)
        pylen = len(pylayer)
        rvx = zeros((pylen, 0), 'd')
        rvy = zeros((pylen, 0), 'd')
        eigvec = take(aq.eigvec, pylayer, 0)
        for d in disx:
            rvx = hstack((rvx, d * eigvec))
        for d in disy:
            rvy = hstack((rvy, d * eigvec))
        if aq == self.aqin:
            rvx = hstack((rvx, zeros((pylen, self.Nparamout), 'd')))
            rvy = hstack((rvy, zeros((pylen, self.Nparamout), 'd')))
        else:
            rvx = hstack((zeros((pylen, self.Nparamin), 'd'), rvx))
            rvy = hstack((zeros((pylen, self.Nparamin), 'd'), rvy))
        return [rvx, rvy]

    def dischargeContribution(self, aq, x, y):
        '''Returns matrix with two rowvectors of dischargeContributions Qx and Qy. Needs to be overloaded cause there is inside and outside'''
        disInf = self.dischargeInfluence(aq, x, y)
        if aq == self.aqin:
            disxInf = sum(self.paramin * disInf[0], 0)
            disyInf = sum(self.paramin * disInf[1], 0)
        else:
            disxInf = sum(self.paramout * disInf[0], 0)
            disyInf = sum(self.paramout * disInf[1], 0)
        return array([disxInf, disyInf])

    def zeroFunction(self, aqdum, ldum, xdum, ydum):
        '''Returns list of zeros of length number of parameters'''
        return list(zeros(self.Nparamin + self.Nparamout, 'd'))

    def getMatrixRows(self, elementList):
        rows = []
        # Jump in potential
        for i in range(self.Ncp):  # Hardcoded for same number of aquifers in and out
            rowin = zeros((self.aqin.Naquifers, 0), 'd')
            rowout = zeros((self.aqin.Naquifers, 0), 'd')  # Zero columns!
            for e in elementList:
                rowpart = e.getMatrixCoefficients(self.aqin, self.ldum, self.xcp[i], self.ycp[i],
                                                  lambda el, aq, ldum, x, y: el.potentialInfluenceAllLayers(aq, ldum, x, y))
                if size(rowpart) > 0:
                    rowin = hstack((rowin, rowpart))
            for e in elementList:
                rowpart = e.getMatrixCoefficients(self.aqout, self.ldum, self.xcp[i], self.ycp[i],
                                                  lambda el, aq, ldum, x, y: el.potentialInfluenceAllLayers(aq, ldum, x, y))
                if size(rowpart) > 0:
                    rowout = hstack((rowout, rowpart))
            row = self.aqout.Tcol * rowin - self.aqin.Tcol * rowout
            row = hstack((row,
                          self.aqin.Tcol * self.aqout.Tcol * (self.aqout.hstar - self.aqin.hstar) +
                          self.aqin.Tcol * self.modelParent.potentialVectorCol(self.xcp[i], self.ycp[i], self.aqout) -
                          self.aqout.Tcol * self.modelParent.potentialVectorCol(self.xcp[i], self.ycp[i], self.aqin)))
            rows = rows + row.tolist()
        # Jump in radial discharge
        for i in range(self.Ncp):  # Hardcoded for same number of aquifers in and out
            rowin = zeros((self.aqin.Naquifers, 0), 'd')
            rowout = zeros((self.aqin.Naquifers, 0), 'd')  # Zero columns!
            for e in elementList:
                rowqxqy = e.getMatrixCoefficients(self.aqin, self.ldum, self.xcp[i], self.ycp[i],
                                                  lambda el, aq, ldum, x, y: el.dischargeInfluenceAllLayers(aq, ldum, x, y))
                if size(rowqxqy) > 0:
                    rowpart = rowqxqy[0] * cos(self.thetacp[i]) + \
                        rowqxqy[1] * sin(self.thetacp[i])
                    rowin = hstack((rowin, rowpart))
            for e in elementList:
                rowqxqy = e.getMatrixCoefficients(self.aqout, self.ldum, self.xcp[i], self.ycp[i],
                                                  lambda el, aq, ldum, x, y: el.dischargeInfluenceAllLayers(aq, ldum, x, y))
                if size(rowqxqy) > 0:
                    rowpart = rowqxqy[0] * cos(self.thetacp[i]) + \
                        rowqxqy[1] * sin(self.thetacp[i])
                    rowout = hstack((rowout, rowpart))
            row = rowin - rowout
            row = hstack((row,
                          self.modelParent.dischargeNormVectorCol(self.xcp[i], self.ycp[i], self.thetacp[i], self.aqout) -
                          self.modelParent.dischargeNormVectorCol(self.xcp[i], self.ycp[i], self.thetacp[i], self.aqin)))
            rows = rows + row.tolist()
        return rows

    def getMatrixCoefficients(self, aq, pylayer, x, y, func):
        return func(self, aq, pylayer, x, y)

    def takeParameters(self, xsol, icount):
        par = xsol[icount: icount + self.Nparamin]
        self.paramin = self.paramin + \
            reshape(par, (2 * self.order + 1, self.aqin.Naquifers))
        icount = icount + self.Nparamin
        par = xsol[icount: icount + self.Nparamout]
        self.paramout = self.paramout + \
            reshape(par, (2 * self.order + 1, self.aqin.Naquifers))
        icount = icount + self.Nparamout
        return icount

    def check(self):
        for i in range(self.Ncp):
            print('Control point ' + str(i))
            print('head inside:  ' + str(self.modelParent.headVector(self.xcp[i], self.ycp[i], self.aqin)))
            print('head outside: ' + str(self.modelParent.headVector(self.xcp[i], self.ycp[i], self.aqout)))
            print('Qrad inside:  ' + str(self.modelParent.dischargeNormVector(self.xcp[i], self.ycp[i], self.thetacp[i], self.aqin)))
            print('Qrad outside: ' + str(self.modelParent.dischargeNormVector(self.xcp[i], self.ycp[i], self.thetacp[i], self.aqout)))

    def layout(self):
        return [0]

    def parSum(self):
        return sum(sum(abs(self.paramin))) + sum(sum(abs(self.paramout)))

    def nearElement(self, pyLayer, xyz1, xyz2, step, idir):
        changed = 0
        stop = 0
        xyznew = zeros(3, 'd')
        r1sq = (xyz1[0] - self.xc)**2 + (xyz1[1] - self.yc) ** 2
        r2sq = (xyz2[0] - self.xc)**2 + (xyz2[1] - self.yc) ** 2
        if (r1sq < self.Rsq and r2sq > self.Rsq) or (r1sq > self.Rsq and r2sq < self.Rsq):
            (x1, y1) = xyz1[0:2]
            (x2, y2) = xyz2[0:2]
            a = (x2 - x1)**2 + (y2 - y1)**2
            b = 2.0 * ((x2 - x1) * (x1 - self.xc) + (y2 - y1) * (y1 - self.yc))
            c = self.xc**2 + self.yc**2 + x1**2 + y1**2 - \
                2.0 * (self.xc * x1 + self.yc * y1) - self.Rsq
            u1 = (-b - sqrt(b**2 - 4.0 * a * c)) / (2.0 * a)
            u2 = (-b + sqrt(b**2 - 4.0 * a * c)) / (2.0 * a)
            if u1 > 0:
                u = u1 * 1.000001  # Go just beyond circle
            else:
                u = u2 * 1.000001  # Go just beyond circle
            xn = x1 + u * (x2 - x1)
            yn = y1 + u * (y2 - y1)
            zn = xyz1[2] + u * (xyz2[2] - xyz1[2])
            xyznew = array([xn, yn, zn], 'd')
            changed = 1
        return [changed, stop, xyznew]

    def distanceSquaredToElement(self, x, y):
        '''Returns distance squared to element. Used for deciding tracing step'''
        dis = sqrt((x - self.xc)**2 + (y - self.yc) ** 2)
        dissq = (dis - self.R) ** 2
        return dissq

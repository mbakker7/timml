import numpy as np

class HeadEquation:
    def equation(self):
        '''Mix-in class that returns matrix rows for head-specified conditions. (really written as constant potential element)
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            rhs[istart:istart+self.Nlayers] = self.pc
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    e.potinflayers(self.xc[icp], self.yc[icp], self.pylayers)
                    ieq += e.Nunknowns
                else:
                    rhs[istart:istart+self.Nlayers] -= \
                    e.potentiallayers(self.xc[icp], self.yc[icp], self.pylayers)  # Pretty cool that this works, really
        return mat, rhs
    
class DisvecEquation:
    def equation(self):
        '''Mix-in class that returns matrix rows for head-specified conditions. (really written as constant potential element)
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    qx, qy = e.disinflayers(self.xc[icp], self.yc[icp], self.pylayers)
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
                    ieq += e.Nunknowns
                else:
                    qx, qy = e.disveclayers(self.xc[icp], self.yc[icp], self.pylayers)
                    rhs[istart:istart+self.Nlayers] -= qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
        return mat, rhs
    
class DisvecEquationOut:
    def equation(self):
        '''Mix-in class that returns matrix rows for head-specified conditions. (really written as constant potential element)
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    qx, qy = e.disinflayers(self.xcout[icp], self.ycout[icp], self.pylayers)
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
                    ieq += e.Nunknowns
                else:
                    qx, qy = e.disveclayers(self.xcout[icp], self.ycout[icp], self.pylayers)
                    rhs[istart:istart+self.Nlayers] -= qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
        return mat, rhs
    
class HeadDiffEquation:
    def equation(self):
        '''Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    e.potinflayers(self.xcin[icp], self.ycin[icp], self.pylayers, aq=self.aqin) / self.aqin.Tcol - \
                    e.potinflayers(self.xcout[icp], self.ycout[icp], self.pylayers, aq=self.aqout) / self.aqout.Tcol
                    ieq += e.Nunknowns
                else:
                    rhs[istart:istart+self.Nlayers] -= \
                    e.potentiallayers(self.xcin[icp], self.ycin[icp], self.pylayers, aq=self.aqin) / self.aqin.T - \
                    e.potentiallayers(self.xcout[icp], self.ycout[icp], self.pylayers, aq=self.aqout) / self.aqout.T
        return mat, rhs
    
class HeadDiffEquation2:
    # Integrated head inside and outside are equal
    def equation(self):
        '''Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    headin = self.intpot(e.potinflayers, self.xcin[icp], self.ycin[icp], \
                                        self.xcin[icp+1], self.ycin[icp+1], self.pylayers, \
                                        ndeg=self.ndeg, aq=self.aqin) / self.aqin.Tcol
                    headout = self.intpot(e.potinflayers, self.xcout[icp], self.ycout[icp], \
                                         self.xcout[icp+1], self.ycout[icp+1], self.pylayers, \
                                         ndeg=self.ndeg, aq=self.aqout) / self.aqout.Tcol
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = headin - headout
                    ieq += e.Nunknowns
                else:
                    headin = self.intpot(e.potentiallayers, self.xcin[icp], self.ycin[icp], \
                                         self.xcin[icp+1], self.ycin[icp+1], self.pylayers, \
                                         ndeg=self.ndeg, aq=self.aqin) / self.aqin.T
                    headout = self.intpot(e.potentiallayers, self.xcout[icp], self.ycout[icp], \
                                          self.xcout[icp+1], self.ycout[icp+1], self.pylayers, \
                                          ndeg=self.ndeg, aq=self.aqout) / self.aqout.T
                    rhs[istart:istart+self.Nlayers] -= headin - headout
        return mat, rhs
    
class DisvecDiffEquation:
    def equation(self):
        '''Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    qxin, qyin = e.disinflayers(self.xcin[icp], self.ycin[icp], self.pylayers, aq=self.aqin)
                    qxout, qyout = e.disinflayers(self.xcout[icp], self.ycout[icp], self.pylayers, aq=self.aqout)
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = \
                    (qxin - qxout) * self.cosnorm[icp] + (qyin - qyout) * self.sinnorm[icp]
                    ieq += e.Nunknowns
                else:
                    qxin, qyin = e.disveclayers(self.xcin[icp], self.ycin[icp], self.pylayers, aq=self.aqin)
                    qxout, qyout = e.disveclayers(self.xcout[icp], self.ycout[icp], self.pylayers, aq=self.aqout)
                    rhs[istart:istart+self.Nlayers] -= (qxin - qxout) * self.cosnorm[icp] + (qyin - qyout) * self.sinnorm[icp]
        return mat, rhs
    
class DisvecDiffEquation2:
    def equation(self):
        '''Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (Nunknowns,Neq)
        Returns rhs part Nunknowns
        '''
        mat = np.empty((self.Nunknowns, self.model.Neq))
        rhs = np.zeros(self.Nunknowns)  # Needs to be initialized to zero
        for icp in range(self.Ncp):
            istart = icp * self.Nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.Nunknowns > 0:
                    fluxin = self.intflux(e.disinflayers, self.xcin[icp], self.ycin[icp], \
                                          self.xcin[icp+1], self.ycin[icp+1], self.pylayers, ndeg=self.ndeg, aq=self.aqin)
                    fluxout = self.intflux(e.disinflayers, self.xcout[icp], self.ycout[icp], \
                                          self.xcout[icp+1], self.ycout[icp+1], self.pylayers, ndeg=self.ndeg, aq=self.aqout)
                    mat[istart:istart+self.Nlayers, ieq:ieq+e.Nunknowns] = fluxin - fluxout
                    ieq += e.Nunknowns
                else:
                    fluxin = self.intflux(e.disveclayers, self.xcin[icp], self.ycin[icp], \
                                          self.xcin[icp+1], self.ycin[icp+1], self.pylayers, ndeg=self.ndeg, aq=self.aqin)
                    fluxout = self.intflux(e.disveclayers, self.xcout[icp], self.ycout[icp], \
                                          self.xcout[icp+1], self.ycout[icp+1], self.pylayers, ndeg=self.ndeg, aq=self.aqout)                    
                    rhs[istart:istart+self.Nlayers] -= fluxin - fluxout
        return mat, rhs
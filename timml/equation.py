import numpy as np


class PotentialEquation:
    def equation(self):
        """Mix-in class that returns matrix rows for potential-specified conditions.
        Returns matrix part (nunknowns,neq)
        Returns rhs part nunknowns
        """
        mat = np.empty((self.nunknowns, self.model.neq))
        # rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        rhs = self.pc.copy()
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            # rhs[istart:istart+self.nlayers] = self.pc[]
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    mat[
                        istart : istart + self.nlayers, ieq : ieq + e.nunknowns
                    ] = e.potinflayers(self.xc[icp], self.yc[icp], self.layers)
                    if e == self:
                        mat[
                            istart : istart + self.nlayers, ieq : ieq + e.nunknowns
                        ] -= self.resfac[icp]
                    ieq += e.nunknowns
                else:
                    rhs[istart : istart + self.nlayers] -= e.potentiallayers(
                        self.xc[icp], self.yc[icp], self.layers
                    )  # Pretty cool that this works, really
        return mat, rhs


class HeadEquation:
    def equation(self):
        """Mix-in class that returns matrix rows for head-specified conditions. (now written as heads)
        Returns matrix part (nunknowns,neq)
        Returns rhs part nunknowns
        """
        mat = np.empty((self.nunknowns, self.model.neq))
        # rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        rhs = self.hc.copy()
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            # rhs[istart:istart+self.nlayers] = self.pc[]
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    mat[istart : istart + self.nlayers, ieq : ieq + e.nunknowns] = (
                        e.potinflayers(self.xc[icp], self.yc[icp], self.layers)
                        / self.aq.Tcol[self.layers]
                    )
                    if e == self:
                        mat[
                            istart : istart + self.nlayers, ieq : ieq + e.nunknowns
                        ] -= self.resfac[icp]
                    ieq += e.nunknowns
                else:
                    rhs[istart : istart + self.nlayers] -= (
                        e.potentiallayers(self.xc[icp], self.yc[icp], self.layers)
                        / self.aq.T[self.layers]
                    )  # Pretty cool that this works, really
        return mat, rhs


class HeadEquationNoRes:  # This class can be deleted when HeadEquation works with zero resistance
    def equation(self):
        """Mix-in class that returns matrix rows for head-specified conditions. (really written as constant potential element)
        Returns matrix part (nunknowns,neq)
        Returns rhs part nunknowns
        """
        mat = np.empty((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            rhs[istart : istart + self.nlayers] = self.pc
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    mat[
                        istart : istart + self.nlayers, ieq : ieq + e.nunknowns
                    ] = e.potinflayers(self.xc[icp], self.yc[icp], self.layers)
                    ieq += e.nunknowns
                else:
                    rhs[istart : istart + self.nlayers] -= e.potentiallayers(
                        self.xc[icp], self.yc[icp], self.layers
                    )  # Pretty cool that this works, really
        return mat, rhs


class MscreenWellEquation:
    def equation(self):
        """Mix-in class that returns matrix rows for mscreen condition.
        Mscreen condition applied at each control point separately (so not like in ditch).
        Returns matrix part (nunknowns,neq)
        Returns rhs part (nunknowns)
        """
        mat = np.zeros((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        rhs[0 : self.nlayers - 1] = 0.0
        rhs[self.nlayers - 1] = self.Qc
        ieq = 0
        for e in self.model.elementlist:
            if e.nunknowns > 0:
                head = (
                    e.potinflayers(self.xc, self.yc, self.layers)
                    / self.aq.Tcol[self.layers, :]
                )
                mat[0 : self.nlayers - 1, ieq : ieq + e.nunknowns] = (
                    head[:-1] - head[1:]
                )
                if e == self:
                    for i in range(self.nlayers - 1):
                        mat[i, ieq + i] -= self.resfac[i]
                        mat[i, ieq + i + 1] += self.resfac[i + 1]
                    mat[self.nlayers - 1, ieq : ieq + self.nlayers] = 1.0
                ieq += e.nunknowns
            else:
                head = (
                    e.potentiallayers(self.xc, self.yc, self.layers)
                    / self.aq.T[self.layers]
                )
                rhs[0 : self.nlayers - 1] -= head[:-1] - head[1:]
        return mat, rhs


class DisvecEquation:
    def equation(self):
        """Mix-in class that returns matrix rows for zero normal flux conditions.
        Returns matrix part (nunknowns,neq)
        Returns rhs part nunknowns
        """
        mat = np.empty((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    qx, qy = e.disvecinflayers(self.xc[icp], self.yc[icp], self.layers)
                    mat[istart : istart + self.nlayers, ieq : ieq + e.nunknowns] = (
                        qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
                    )
                    ieq += e.nunknowns
                else:
                    qx, qy = e.disveclayers(self.xc[icp], self.yc[icp], self.layers)
                    rhs[istart : istart + self.nlayers] -= (
                        qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
                    )
        return mat, rhs


class DisvecEquationOut:
    def equation(self):
        """Mix-in class that returns matrix rows for zero normal flux condition.
        Using the control point on the outside
        Returns matrix part (nunknowns,neq)
        Returns rhs part nunknowns
        """
        mat = np.empty((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    qx, qy = e.disvecinflayers(
                        self.xcout[icp], self.ycout[icp], self.layers
                    )
                    mat[istart : istart + self.nlayers, ieq : ieq + e.nunknowns] = (
                        qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
                    )
                    ieq += e.nunknowns
                else:
                    qx, qy = e.disveclayers(
                        self.xcout[icp], self.ycout[icp], self.layers
                    )
                    rhs[istart : istart + self.nlayers] -= (
                        qx * self.cosnorm[icp] + qy * self.sinnorm[icp]
                    )
        return mat, rhs


class LeakyWallEquation:
    def equation(self):
        """Mix-in class that returns matrix rows for leaky wall condition.
        Qnormal = resfac * (headin - headout)
        Returns matrix part (nunknowns,neq)
        Returns rhs part nunknowns
        """
        mat = np.empty((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    qx, qy = e.disvecinflayers(self.xc[icp], self.yc[icp], self.layers)
                    mat[istart : istart + self.nlayers, ieq : ieq + e.nunknowns] = (
                        qx * self.cosnorm[icp]
                        + qy * self.sinnorm[icp]
                        - self.resfac[:, np.newaxis]
                        * (
                            e.potinflayers(
                                self.xcin[icp], self.ycin[icp], self.layers, aq=self.aq
                            )
                            / self.aq.Tcol[self.layers]
                            - e.potinflayers(
                                self.xcout[icp],
                                self.ycout[icp],
                                self.layers,
                                aq=self.aq,
                            )
                            / self.aq.Tcol[self.layers]
                        )
                    )
                    ieq += e.nunknowns
                else:
                    qx, qy = e.disveclayers(self.xc[icp], self.yc[icp], self.layers)
                    rhs[istart : istart + self.nlayers] -= (
                        qx * self.cosnorm[icp]
                        + qy * self.sinnorm[icp]
                        + self.resfac
                        * (
                            e.potentiallayers(
                                self.xcin[icp], self.ycin[icp], self.layers, aq=self.aq
                            )
                            / self.aq.T[self.layers]
                            - e.potentiallayers(
                                self.xcout[icp],
                                self.ycout[icp],
                                self.layers,
                                aq=self.aq,
                            )
                            / self.aq.T[self.layers]
                        )
                    )
        return mat, rhs


class HeadDiffEquation:
    def equation(self):
        """Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (nunknowns,neq)
        Returns rhs part nunknowns
        """
        mat = np.empty((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    mat[istart : istart + self.nlayers, ieq : ieq + e.nunknowns] = (
                        e.potinflayers(
                            self.xcin[icp], self.ycin[icp], self.layers, aq=self.aqin
                        )
                        / self.aqin.Tcol
                        - e.potinflayers(
                            self.xcout[icp], self.ycout[icp], self.layers, aq=self.aqout
                        )
                        / self.aqout.Tcol
                    )
                    ieq += e.nunknowns
                else:
                    rhs[istart : istart + self.nlayers] -= (
                        e.potentiallayers(
                            self.xcin[icp], self.ycin[icp], self.layers, aq=self.aqin
                        )
                        / self.aqin.T
                        - e.potentiallayers(
                            self.xcout[icp], self.ycout[icp], self.layers, aq=self.aqout
                        )
                        / self.aqout.T
                    )
        return mat, rhs


class HeadDiffEquation2:
    # Integrated head inside and outside are equal
    def equation(self):
        """Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (nunknowns,neq)
        Returns rhs part nunknowns
        """
        mat = np.empty((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    headin = (
                        self.intpot(
                            e.potinflayers,
                            self.xcin[icp],
                            self.ycin[icp],
                            self.xcin[icp + 1],
                            self.ycin[icp + 1],
                            self.layers,
                            aq=self.aqin,
                        )
                        / self.aqin.Tcol[self.layers]
                    )
                    headout = (
                        self.intpot(
                            e.potinflayers,
                            self.xcout[icp],
                            self.ycout[icp],
                            self.xcout[icp + 1],
                            self.ycout[icp + 1],
                            self.layers,
                            aq=self.aqout,
                        )
                        / self.aqout.Tcol[self.layers]
                    )
                    mat[istart : istart + self.nlayers, ieq : ieq + e.nunknowns] = (
                        headin - headout
                    )
                    ieq += e.nunknowns
                else:
                    headin = (
                        self.intpot(
                            e.potentiallayers,
                            self.xcin[icp],
                            self.ycin[icp],
                            self.xcin[icp + 1],
                            self.ycin[icp + 1],
                            self.layers,
                            aq=self.aqin,
                        )
                        / self.aqin.T[self.layers]
                    )
                    headout = (
                        self.intpot(
                            e.potentiallayers,
                            self.xcout[icp],
                            self.ycout[icp],
                            self.xcout[icp + 1],
                            self.ycout[icp + 1],
                            self.layers,
                            aq=self.aqout,
                        )
                        / self.aqout.T[self.layers]
                    )
                    rhs[istart : istart + self.nlayers] -= headin - headout
        return mat, rhs


class DisvecDiffEquation:
    def equation(self):
        """Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (nunknowns,neq)
        Returns rhs part nunknowns
        """
        mat = np.empty((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    qxin, qyin = e.disvecinflayers(
                        self.xcin[icp], self.ycin[icp], self.layers, aq=self.aqin
                    )
                    qxout, qyout = e.disvecinflayers(
                        self.xcout[icp], self.ycout[icp], self.layers, aq=self.aqout
                    )
                    mat[istart : istart + self.nlayers, ieq : ieq + e.nunknowns] = (
                        qxin - qxout
                    ) * self.cosnorm[icp] + (qyin - qyout) * self.sinnorm[icp]
                    ieq += e.nunknowns
                else:
                    qxin, qyin = e.disveclayers(
                        self.xcin[icp], self.ycin[icp], self.layers, aq=self.aqin
                    )
                    qxout, qyout = e.disveclayers(
                        self.xcout[icp], self.ycout[icp], self.layers, aq=self.aqout
                    )
                    rhs[istart : istart + self.nlayers] -= (
                        qxin - qxout
                    ) * self.cosnorm[icp] + (qyin - qyout) * self.sinnorm[icp]
        return mat, rhs


class DisvecDiffEquation2:
    def equation(self):
        """Mix-in class that returns matrix rows for difference in head between inside and
        outside equals zeros
        Returns matrix part (nunknowns,neq)
        Returns rhs part nunknowns
        """
        mat = np.empty((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    fluxin = self.intflux(
                        e.disvecinflayers,
                        self.xcin[icp],
                        self.ycin[icp],
                        self.xcin[icp + 1],
                        self.ycin[icp + 1],
                        self.layers,
                        aq=self.aqin,
                    )
                    fluxout = self.intflux(
                        e.disvecinflayers,
                        self.xcout[icp],
                        self.ycout[icp],
                        self.xcout[icp + 1],
                        self.ycout[icp + 1],
                        self.layers,
                        aq=self.aqout,
                    )
                    mat[istart : istart + self.nlayers, ieq : ieq + e.nunknowns] = (
                        fluxin - fluxout
                    )
                    ieq += e.nunknowns
                else:
                    fluxin = self.intflux(
                        e.disveclayers,
                        self.xcin[icp],
                        self.ycin[icp],
                        self.xcin[icp + 1],
                        self.ycin[icp + 1],
                        self.layers,
                        aq=self.aqin,
                    )
                    fluxout = self.intflux(
                        e.disveclayers,
                        self.xcout[icp],
                        self.ycout[icp],
                        self.xcout[icp + 1],
                        self.ycout[icp + 1],
                        self.layers,
                        aq=self.aqout,
                    )
                    rhs[istart : istart + self.nlayers] -= fluxin - fluxout
        return mat, rhs


class IntDisVecEquation:
    def equation(self):
        mat = np.empty((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    flux = self.intflux(
                        e.disvecinflayers,
                        self.xc[icp],
                        self.yc[icp],
                        self.xc[icp + 1],
                        self.yc[icp + 1],
                        self.layers,
                        aq=self.aq,
                    )
                    mat[istart : istart + self.nlayers, ieq : ieq + e.nunknowns] = flux
                    ieq += e.nunknowns
                else:
                    flux = self.intflux(
                        e.disveclayers,
                        self.xc[icp],
                        self.yc[icp],
                        self.xc[icp + 1],
                        self.yc[icp + 1],
                        self.layers,
                        aq=self.aq,
                    )
                    rhs[istart : istart + self.nlayers] -= flux

        return mat, rhs


class IntLeakyWallEquation:
    def equation(self):
        mat = np.empty((self.nunknowns, self.model.neq))
        rhs = np.zeros(self.nunknowns)  # Needs to be initialized to zero
        for icp in range(self.ncp):
            istart = icp * self.nlayers
            ieq = 0
            for e in self.model.elementlist:
                if e.nunknowns > 0:
                    flux = self.intflux(
                        e.disvecinflayers,
                        self.xc[icp],
                        self.yc[icp],
                        self.xc[icp + 1],
                        self.yc[icp + 1],
                        self.layers,
                        aq=self.aq,
                    )
                    headin = (
                        self.intpot(
                            e.potinflayers,
                            self.xcin[icp],
                            self.ycin[icp],
                            self.xcin[icp + 1],
                            self.ycin[icp + 1],
                            self.layers,
                            aq=self.aqin,
                        )
                        / self.aqin.Tcol[self.layers]
                    )
                    headout = (
                        self.intpot(
                            e.potinflayers,
                            self.xcout[icp],
                            self.ycout[icp],
                            self.xcout[icp + 1],
                            self.ycout[icp + 1],
                            self.layers,
                            aq=self.aqout,
                        )
                        / self.aqout.Tcol[self.layers]
                    )

                    mat[istart : istart + self.nlayers, ieq : ieq + e.nunknowns] = (
                        flux - self.resfac * (headin - headout)
                    )
                    ieq += e.nunknowns
                else:
                    flux = self.intflux(
                        e.disveclayers,
                        self.xc[icp],
                        self.yc[icp],
                        self.xc[icp + 1],
                        self.yc[icp + 1],
                        self.layers,
                        aq=self.aq,
                    )
                    headin = (
                        self.intpot(
                            e.potentiallayers,
                            self.xcin[icp],
                            self.ycin[icp],
                            self.xcin[icp + 1],
                            self.ycin[icp + 1],
                            self.layers,
                            aq=self.aqin,
                        )
                        / self.aqin.T[self.layers]
                    )
                    headout = (
                        self.intpot(
                            e.potentiallayers,
                            self.xcout[icp],
                            self.ycout[icp],
                            self.xcout[icp + 1],
                            self.ycout[icp + 1],
                            self.layers,
                            aq=self.aqout,
                        )
                        / self.aqout.T[self.layers]
                    )

                    rhs[istart : istart + self.nlayers] += (
                        -flux + self.resfac.squeeze() * (headin - headout)
                    )
        return mat, rhs

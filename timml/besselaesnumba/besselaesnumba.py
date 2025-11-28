import numba
import numpy as np


def initialize():
    pass


@numba.njit(nogil=True, cache=True)
def make_rbinom():
    rrange = np.arange(9.0)
    rbinom = np.zeros((9, 9))
    for n in range(9):
        for m in range(9):
            if m > n:
                rbinom[n, m] = 1
            else:
                rbinom[n, m] = np.prod(rrange[m + 1 : n + 1]) / np.prod(
                    rrange[1 : n - m + 1]
                )
    return rbinom


RBINOM = make_rbinom()
AC = np.zeros(9)
BC = np.zeros(9)
AC1 = np.zeros(9)
BC1 = np.zeros(9)
NTERMS = 8

# Coefficients of Table 1
AC[0] = -0.500004211065677e0
BC[0] = 0.115956920789028e0
AC[1] = -0.124989431448700e0
BC[1] = 0.278919134951974e0
AC[2] = -0.781685695640975e-2
BC[2] = 0.252752008621167e-1
AC[3] = -0.216324415010288e-3
BC[3] = 0.841879407543506e-3
AC[4] = -0.344525393452639e-5
BC[4] = 0.152425102734818e-4
AC[5] = -0.315133836828774e-7
BC[5] = 0.148292488399579e-6
AC[6] = -0.296636186265427e-9
BC[6] = 0.157622547107156e-8
AC[7] = -0.313689942474032e-12
BC[7] = 0.117975437124933e-11
AC[8] = -0.112031912249579e-13
BC[8] = 0.655345107753534e-13

# Coefficients of K1
AC1[0] = 0.250000197208863e0
BC1[0] = -0.307966963840932e0
AC1[1] = 0.312495047266388e-1
BC1[1] = -0.853676915840295e-1
AC1[2] = 0.130228768540005e-2
BC1[2] = -0.464343185899275e-2
AC1[3] = 0.270943632576982e-4
BC1[3] = -0.112338260301993e-3
AC1[4] = 0.341642640876988e-6
BC1[4] = -0.157491472277549e-5
AC1[5] = 0.271286480571077e-8
BC1[5] = -0.133414321160211e-7
AC1[6] = 0.197096143802491e-10
BC1[6] = -0.106342159633141e-9
AC1[7] = 0.329351103353054e-13
BC1[7] = -0.159531523882074e-12
AC1[8] = 0.573031034976631e-15
BC1[8] = -0.340195779923156e-14


@numba.njit(nogil=True, cache=True)
def prepare_z(x, y, z1, z2):
    zin = complex(x, y)
    z1in = z1
    z2in = z2
    Lin = abs(z2in - z1in)
    z = (2.0 * zin - (z1in + z2in)) / (z2in - z1in)
    zplus1 = z + 1.0
    zmin1 = z - 1.0
    # If at cornerpoint, move slightly
    if abs(zplus1) < 1.0e-8 * 2.0 / Lin:
        zplus1 = zplus1 + 1.0e-8
    if abs(zmin1) < 1.0e-8 * 2.0 / Lin:
        zmin1 = zmin1 + 1.0e-8
    return zin, z1in, z2in, Lin, z, zplus1, zmin1


@numba.njit(nogil=True, cache=True)
def potbeslsho(x, y, z1, z2, labda, order, ilap, naq):
    """Potbeslsho.

    Parameters
    ----------
       x,y: Point where potential is computed
       z1: Complex begin point of line-sink
       z2: Complex end point of line-sink
       labda(naq): labda's (zero for first labda if Laplace)
       order: Order of the line-sink
       ilap: equals 1 when first value is Laplace line-sink and first labda equals zero
       naq: Number of aquifers
       rv(naq): Array to store return value (must be pre-allocated)

    Returns
    -------
       rv(naq): Potentials. First spot is Laplace value if ilap=1
    """
    rv = np.zeros(naq)

    # lstype = 1 means line-sink
    lstype = 1

    # Radius of convergence
    if order > 5:
        Rconv = 5.0
    else:
        Rconv = 7.0

    # if (ilap==1) :
    #    istart = 1
    # else:
    #    istart = 0
    zin, z1in, z2in, Lin, z, zplus1, zmin1 = prepare_z(x, y, z1, z2)
    # Laplace linesink
    if ilap == 1:
        power = order + 1
        pcor = complex(0.0, 0.0)
        for n in range(1, int((power + 1) / 2) + 1):
            pcor = pcor + z ** (power - 2 * n + 1) / (2 * n - 1)
        pcor = 2.0 * pcor
        comega = (
            z**power * np.log((zmin1) / (zplus1))
            + pcor
            - np.log(zmin1)
            + (-1.0) ** power * np.log(zplus1)
        )
        comega = -comega * Lin / (4.0 * np.pi * power)
        rv[0] = np.real(comega)
    # N-1 leakage factors
    for i in range(ilap, naq):
        # Check whether entire linesink is outside radius of convergence
        # Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L,
        # or |zeta|>1/biglab+7 (zeta is called z here)
        biglab = 2.0 * labda[i] / Lin
        z = (2.0 * zin - (z1in + z2in)) / (z2in - z1in) / biglab

        if abs(z) < (Rconv + 1.0 / biglab):
            pot = IntegralF(zin, z1in, z2in, Lin, labda[i], order, Rconv, lstype)
            rv[i] = -Lin / 2.0 * pot
        else:
            rv[i] = 0.0
    return rv


@numba.njit(nogil=True, cache=True)
def potbeslsv(x, y, z1, z2, lab, order, ilap, naq):
    # Check if endpoints need to be adjusted using the largest labda (the first one)
    pot = np.zeros((order + 1, naq))
    for n in range(0, order + 1):
        pot[n, 0 : naq + 1] = potbeslsho(x, y, z1, z2, lab, n, ilap, naq)
    return pot


@numba.njit(nogil=True, cache=True)
def disbeslsho(x, y, z1, z2, labda, order, ilap, naq):
    # Input:
    #   x,y: Point where discharge is computed
    #   z1: Complex begin point of line-sink
    #   z2: Complex end point of line-sink
    #   labdain(Naquifers): Array with zero in first spot and labda's in remaining spots
    #   order: Order of the line-sink
    #   Naquifers: Number of aquifers
    #   rvx(Naquifers),rvy(Naquifers): Arrays to store return values
    # Output:
    #   rvx(Naquifers),rvy(Naquifers): Values of Qx and Qy with Laplace value in first
    # spot and mod.Helmholtz potentials in remaining spots

    rv = np.zeros((2, naq))
    # Radius of convergence
    if order > 5:
        Rconv = 5.0
    else:
        Rconv = 7.0

    # lstype = 1 means line-sink
    lstype = 1

    zin, z1in, z2in, Lin, z, zplus1, zmin1 = prepare_z(x, y, z1, z2)

    pcor = 0.0
    # Laplace linesink
    if ilap == 1:
        pcor = complex(0.0, 0.0)
        for n in range(1, int((order + 1) / 2) + 1):
            pcor = pcor + float(order - 2 * n + 2) * z ** (order + 1 - 2 * n) / float(
                2 * n - 1
            )

        pcor = 2.0 * pcor

        cdum = 1.0 / (
            order + 1
        )  # Without this intermediate statement it didn't seem to work
        wdis = float(order + 1) * z**order * np.log((zmin1) / (zplus1)) + pcor

        wdis = (
            wdis
            + (z ** (order + 1) - 1.0) / zmin1
            - (z ** (order + 1) - (-1.0) ** (order + 1)) / zplus1
        )

        wdis = wdis * Lin / 2.0 / (z2in - z1in) / np.pi * cdum

        rv[0, 0] = np.real(wdis)
        rv[1, 0] = -np.imag(wdis)

    for i in range(ilap, naq):
        wdis = complex(0.0, 0.0)

        # Check whether entire linesink is outside radius of convergence
        # Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L,
        # or |zeta|>1/biglab+7 (zeta is called z here)
        biglab = 2.0 * labda[i] / Lin
        z = (2.0 * zin - (z1in + z2in)) / (z2in - z1in) / biglab

        if abs(z) < (Rconv + 1.0 / biglab):
            wdis = IntegralG(zin, z1in, z2in, Lin, labda[i], order, Rconv, lstype)
            wdis = 2.0 * Lin / (z2in - z1in) / biglab * wdis

            rv[0, i] = np.real(wdis)
            rv[1, i] = -np.imag(wdis)
        else:
            rv[0, i] = 0.0
            rv[1, i] = 0.0

    return rv


@numba.njit(nogil=True, cache=True)
def disbeslsv(x, y, z1, z2, lab, order, ilap, naq):
    # locals
    qxqy = np.zeros((2 * (order + 1), naq))
    # Check if endpoints need to be adjusted using the largest labda (the first one)
    for n in range(0, order + 1):
        rv = disbeslsho(x, y, z1, z2, lab, n, ilap, naq)
        qxqy[n, 0 : naq + 1] = rv[0, 0 : naq + 1]
        qxqy[n + order + 1, 0 : naq + 1] = rv[1, 0 : naq + 1]
    return qxqy


@numba.njit(nogil=True, cache=True)
def potbesldho(x, y, z1, z2, labda, order, ilap, naq):
    # Input:
    #   x,y: Point where potential is computed
    #   z1: Complex begin point of line-doublet
    #   z2: Complex end point of line-doublet
    #   labda(naq): labda's (zero for first labda if Laplace)
    #   order: Order of the line-doublet
    #   ilap: equals 1 when first value is Laplace line-doublet and first labda is zero
    #   naq: Number of aquifers
    #   rv(naq): Array to store return value (must be pre-allocated)
    # Output:
    #   rv(naq): Potentials. First spot is Laplace value if ilap=1

    rv = np.zeros(naq)

    # Radius of convergence
    if order > 5:
        Rconv = 5.0
    else:
        Rconv = 7.0

    # lstype=2 means line-doublet
    lstype = 2
    zin, z1in, z2in, Lin, z, zplus1, zmin1 = prepare_z(x, y, z1, z2)

    # Laplace line-doublet
    if ilap == 1:
        comega = z**order * np.log(zmin1 / zplus1)
        qm = complex(0.0, 0.0)
        for n in range(1, int((order + 1) / 2) + 1):
            qm = qm + z ** (order - 2.0 * float(n) + 1.0) / (2.0 * float(n) - 1.0)

        comega = 1.0 / (2.0 * np.pi * complex(0.0, 1.0)) * (comega + 2.0 * qm)
        rv[0] = np.real(comega)

    # N-1 leakage factors
    for i in range(ilap, naq):
        pot = 0.0
        # Check whether entire linedoublet is outside radius of convergence
        # Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L,
        # or |zeta|>1/biglab+7 (zeta is called z here)
        biglab = 2.0 * labda[i] / Lin
        z = (2.0 * zin - (z1in + z2in)) / (z2in - z1in) / biglab

        if abs(z) < (Rconv + 1.0 / biglab):
            m1, m2, NLS = findm1m2(zin, z1in, z2in, Lin, labda[i], Rconv)
            comega = complex(0.0, 0.0)
            if m1 > 0:  # Otherwise outside radius of convergence
                z1new = z1in + float(m1 - 1) / float(NLS) * (z2in - z1in)
                z2new = z1in + float(m2) / float(NLS) * (z2in - z1in)
                del0 = float(1 - m1 - m2 + NLS) / float(1 - m1 + m2)
                ra = float(NLS) / float(1 + m2 - m1)
                comega = IntegralLapLineDipole(zin, z1new, z2new, del0, ra, order)

            pot = IntegralF(zin, z1in, z2in, Lin, labda[i], order, Rconv, lstype)
            rv[i] = (
                np.real(comega / complex(0.0, 1.0)) + np.imag(z) / biglab * pot
            )  # Note that z is really zeta in analysis
        else:
            rv[i] = 0.0

    return rv


@numba.njit(nogil=True, cache=True)
def potbesldv(x, y, z1, z2, lab, order, ilap, naq):
    pot = np.zeros((order + 1, naq))
    # Check if endpoints need to be adjusted using the largest labda (the first one)
    for n in range(0, order + 1):
        pot[n, 0 : naq + 1] = potbesldho(x, y, z1, z2, lab, n, ilap, naq)
    return pot


@numba.njit(nogil=True, cache=True)
def disbesldho(x, y, z1, z2, labda, order, ilap, naq):
    # Input:
    #   x,y: Point where discharge is computed
    #   z1: Complex begin point of line-sink
    #   z2: Complex end point of line-sink
    #   labdain(Naquifers): Array with zero in first spot and labda's in remaining spots
    #   order: Order of the line-sink
    #   naq: Number of aquifers
    # Output:
    #   rv(2, Naquifers),rvy(Naquifers): Values of Qx and Qy with Laplace value in
    # first spot and mod.Helmholtz potentials in remaining spots
    rv = np.zeros((2, naq))
    # Radius of convergence
    if order > 5:
        Rconv = 5.0
    else:
        Rconv = 7.0

    # lstype=2 means line-doublet
    lstype = 2

    zin, z1in, z2in, Lin, z, zplus1, zmin1 = prepare_z(x, y, z1, z2)

    # Laplace line-doublet
    qm = complex(0.0, 0.0)
    if ilap == 1:
        if order == 0:
            wdis = -(1.0 / zmin1 - 1.0 / zplus1) / (
                np.pi * complex(0.0, 1.0) * (z2in - z1in)
            )
        else:
            wdis = float(order) * z ** (order - 1) * np.log(zmin1 / zplus1)
            wdis = wdis + z**order * (1.0 / zmin1 - 1.0 / zplus1)
            qm = complex(0.0, 0.0)
            if order > 1:  # To avoid a possible problem of 0 * 0^(-1)
                for n in range(1, int(order / 2) + 1):
                    qm = qm + float(order - 2 * n + 1) * z ** (order - 2 * n) / float(
                        2 * n - 1
                    )

            wdis = -(wdis + 2.0 * qm) / (np.pi * complex(0.0, 1.0) * (z2in - z1in))

        rv[0, 0] = np.real(wdis)
        rv[1, 0] = -np.imag(wdis)

    # N-1 or N leakage factors
    for i in range(ilap, naq):
        wdis = complex(0.0, 0.0)

        # Check whether entire line-doublet is outside radius of convergence
        # Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L,
        # or |zeta|>1/biglab+7 (zeta is called z here)
        biglab = 2.0 * labda[i] / Lin
        z = (2.0 * zin - (z1in + z2in)) / (z2in - z1in) / biglab

        if abs(z) < (Rconv + 1.0 / biglab):
            m1, m2, NLS = findm1m2(zin, z1in, z2in, Lin, labda[i], Rconv)
            wdis1 = complex(0.0, 0.0)
            if m1 > 0:
                z1new = z1in + float(m1 - 1) / float(NLS) * (z2in - z1in)
                z2new = z1in + float(m2) / float(NLS) * (z2in - z1in)
                del0 = float(1 - m1 - m2 + NLS) / float(1 - m1 + m2)
                ra = float(NLS) / float(1 + m2 - m1)
                wdis1 = IntegralLapLineDipoleDis(zin, z1new, z2new, del0, ra, order)
                wdis1 = -2.0 * wdis1 / (complex(0.0, 1.0) * (z2new - z1new))

            pot = IntegralF(zin, z1in, z2in, Lin, labda[i], order, Rconv, lstype)
            wdis2 = IntegralG(zin, z1in, z2in, Lin, labda[i], order, Rconv, lstype)
            wdis3 = pot / (2.0 * complex(0.0, 1.0)) + wdis2 * np.imag(z)

            wdis = wdis1 - 4.0 * wdis3 / (biglab**2 * (z2in - z1in))
        rv[0, i] = np.real(wdis)
        rv[1, i] = -1.0 * np.imag(wdis)
    return rv


@numba.njit(nogil=True, cache=True)
def disbesldv(x, y, z1, z2, lab, order, ilap, naq):
    qxqy = np.zeros((2 * (order + 1), naq))
    # Check if endpoints need to be adjusted using the largest labda (the first one)
    for n in range(0, order + 1):
        rv = disbesldho(x, y, z1, z2, lab, n, ilap, naq)
        qxqy[n, 0 : naq + 1] = rv[0, 0:naq]
        qxqy[n + order + 1, 0 : naq + 1] = rv[1, 0 : naq + 1]
    return qxqy


@numba.njit(nogil=True, cache=True)
def IntegralF(zin, z1in, z2in, Lin, labda, order, Rconv, lstype):
    czmzbarp = np.full(NTERMS + 1, complex(0.0, 0.0))
    cgamma = np.full((NTERMS + 1, NTERMS + 1), complex(0.0, 0.0))
    calphat = np.full(2 * NTERMS + 1, complex(0.0, 0.0))
    cbetat = np.full(2 * NTERMS + 1, complex(0.0, 0.0))
    cc = np.full(order + 2, complex(0.0, 0.0))
    calpha = np.full(2 * NTERMS + order + 1, complex(0.0, 0.0))
    cbeta = np.full(2 * NTERMS + order + 1, complex(0.0, 0.0))

    m1, m2, NLS = findm1m2(zin, z1in, z2in, Lin, labda, Rconv)
    if m1 == 0:
        # pot = 0.0
        return 0.0

    # Compute zeta (called z here). This is the regular value of the entire element
    L = abs(z2in - z1in)
    biglab = 2.0 * labda / L
    z = (2.0 * zin - (z1in + z2in)) / (z2in - z1in) / biglab
    zbar = np.conj(z)

    # Coefficients gamma(n,m), Eq. 21
    # Store coefficents in matrix.
    for n in range(0, NTERMS + 1):
        czmzbarp[n] = (z - zbar) ** n

    for n in range(0, NTERMS + 1):
        for m in range(0, n + 1):
            cgamma[n, m] = RBINOM[n, m] * czmzbarp[n - m]

    # Eq. 23 These coefficients should be modified for a higher order linesink
    for n in range(0, 2 * NTERMS + 1):
        calphat[n] = complex(0.0, 0.0)
        cbetat[n] = complex(0.0, 0.0)
        for m in range(max(0, n - NTERMS), int(n / 2) + 1):
            if lstype == 1:
                calphat[n] = calphat[n] + AC[n - m] * cgamma[n - m, m]
                cbetat[n] = cbetat[n] + BC[n - m] * cgamma[n - m, m]
            else:
                calphat[n] = calphat[n] + AC1[n - m] * cgamma[n - m, m]
                cbetat[n] = cbetat[n] + BC1[n - m] * cgamma[n - m, m]

    # Compute coefficients of delta^p
    for m in range(0, order + 1):
        cc[m] = RBINOM[order, m] * z ** (order - m) * biglab**order
    if order > 0:
        for n in range(0, 2 * NTERMS + order + 1):
            calpha[n] = complex(0.0, 0.0)
            cbeta[n] = complex(0.0, 0.0)
            for m in range(max(0, n - 2 * NTERMS), min(n, order) + 1):
                calpha[n] = calpha[n] + cc[m] * calphat[n - m]
                cbeta[n] = cbeta[n] + cc[m] * cbetat[n - m]
    else:
        calpha = calphat
        cbeta = cbetat

    # Evaluation of integral, Eq. 25
    cInt = complex(0.0, 0.0)
    del1 = -1.0 + 2.0 * (float(m1) - 1.0) / float(NLS)
    del2 = -1.0 + 2.0 * float(m2) / float(NLS)
    cd1minz = del1 / biglab - z
    cd2minz = del2 / biglab - z
    if abs(cd1minz) < 1.0e-8 / labda:
        cd1minz = cd1minz + 1.0e-8
    if abs(cd2minz) < 1.0e-8 / labda:
        cd2minz = cd2minz + 1.0e-8
    cln1 = np.log(cd1minz)
    cln2 = np.log(cd2minz)
    for n in range(0, 2 * NTERMS + order + 1):
        cInt = cInt + (
            2.0 * calpha[n] * cln2 - 2.0 * calpha[n] / (n + 1) + cbeta[n]
        ) * (cd2minz) ** (n + 1) / float(n + 1)
        cInt = cInt - (
            2.0 * calpha[n] * cln1 - 2.0 * calpha[n] / (n + 1) + cbeta[n]
        ) * (cd1minz) ** (n + 1) / float(n + 1)
    pot = np.real(cInt) * biglab / (2.0 * np.pi)
    return pot


@numba.njit(nogil=True, cache=True)
def IntegralG(zin, z1in, z2in, Lin, labda, order, Rconv, lstype):
    czmzbarp = np.full(NTERMS + 1, complex(0.0, 0.0))
    cgamma = np.full((NTERMS + 1, NTERMS + 1), complex(0.0, 0.0))
    cahat = np.full(2 * (NTERMS - 1) + 1, complex(0.0, 0.0))
    cbhat = np.full(2 * (NTERMS - 1) + 1, complex(0.0, 0.0))
    calphat = np.full(2 * NTERMS + 1, complex(0.0, 0.0))
    cbetat = np.full(2 * NTERMS + 1, complex(0.0, 0.0))
    cc = np.full(order + 2, complex(0.0, 0.0))
    calpha = np.full(2 * NTERMS + order + 1, complex(0.0, 0.0))
    cbeta = np.full(2 * NTERMS + order + 1, complex(0.0, 0.0))

    biglabin = 2.0 * labda / Lin

    m1, m2, NLS = findm1m2(zin, z1in, z2in, Lin, labda, Rconv)
    if m1 == 0:
        # wdis = complex(0.0, 0.0)
        return complex(0.0, 0.0)

    # Compute zeta (called z here). This is the regular value of the entire element
    L = abs(z2in - z1in)
    biglab = 2.0 * labda / L
    z = (2.0 * zin - (z1in + z2in)) / (z2in - z1in) / biglab
    zbar = np.conj(z)

    # Coefficients gamma(n,m), Eq. 21
    # Store coefficents in matrix.
    for n in range(0, NTERMS + 1):
        czmzbarp[n] = (z - zbar) ** n

    for n in range(0, NTERMS + 1):
        for m in range(0, n + 1):
            cgamma[n, m] = RBINOM[n, m] * czmzbarp[n - m]

    # Integral g1
    # Implemented with different z, rather than Delta1 and Delta2
    z1 = z1in + float(m1 - 1) / float(NLS) * (z2in - z1in)
    z2 = z1in + float(m2) / float(NLS) * (z2in - z1in)
    del0 = float(1 - m1 - m2 + NLS) / float(1 - m1 + m2)
    ra = float(NLS) / float(1 + m2 - m1)
    # comega = complex(0.0, 0.0)
    comega = IntegralLapLineDipole(zin, z1, z2, del0, ra, order)

    if lstype == 1:
        g1 = -AC[0] * biglabin * comega
        # Integral g2
        # Compute hat coefficients
        for n in range(0, NTERMS):
            cahat[n] = float(n + 1) * AC[n + 1]
            cbhat[n] = AC[n + 1] + float(n + 1) * BC[n + 1]
    else:
        g1 = -AC1[0] * biglabin * comega
        # Integral g2
        # Compute hat coefficients
        for n in range(0, NTERMS):
            cahat[n] = float(n + 1) * AC1[n + 1]
            cbhat[n] = AC1[n + 1] + float(n + 1) * BC1[n + 1]

    # Eq. 23
    for n in range(0, 2 * NTERMS):
        calphat[n] = complex(0.0, 0.0)
        cbetat[n] = complex(0.0, 0.0)
        for m in range(max(0, n - NTERMS + 1), int((n + 1) / 2) + 1):
            calphat[n] = calphat[n] + cahat[n - m] * cgamma[n - m + 1, m]
            cbetat[n] = cbetat[n] + cbhat[n - m] * cgamma[n - m + 1, m]

    # Compute coefficients of delta^p
    for m in range(0, order + 1):
        cc[m] = RBINOM[order, m] * z ** (order - m) * biglab**order
    if order > 0:
        for n in range(0, 2 * NTERMS + order):
            calpha[n] = complex(0.0, 0.0)
            cbeta[n] = complex(0.0, 0.0)
            for m in range(max(0, n - 2 * NTERMS + 1), min(n, order) + 1):
                calpha[n] = calpha[n] + cc[m] * calphat[n - m]
                cbeta[n] = cbeta[n] + cc[m] * cbetat[n - m]
    else:
        calpha = calphat
        cbeta = cbetat

    # Computation of integral
    g2 = complex(0.0, 0.0)
    del1 = -1.0 + 2.0 * (float(m1) - 1.0) / float(NLS)
    del2 = -1.0 + 2.0 * float(m2) / float(NLS)
    cd1minz = del1 / biglab - z
    cd2minz = del2 / biglab - z
    if abs(cd1minz) < 1.0e-8 / labda:
        cd1minz = cd1minz + 1.0e-8
    if abs(cd2minz) < 1.0e-8 / labda:
        cd2minz = cd2minz + 1.0e-8
    cln1 = np.log(cd1minz)
    cln2 = np.log(cd2minz)
    for n in range(0, 2 * NTERMS - 1 + order + 1):
        g2 = g2 - (calpha[n] * cln2 - calpha[n] / (n + 1) + cbeta[n]) * (cd2minz) ** (
            n + 1
        ) / float(n + 1)
        g2 = g2 + (calpha[n] * cln1 - calpha[n] / (n + 1) + cbeta[n]) * (cd1minz) ** (
            n + 1
        ) / float(n + 1)
    g2 = biglabin * g2 / (2 * np.pi)

    # Integral g3
    # Eq. 23
    calphat[0] = complex(0.0, 0.0)
    for n in range(1, 2 * NTERMS):  # Loop start at 1, because of bug in Digital Fortran
        calphat[n] = complex(0.0, 0.0)
        for m in range(max(0, n - NTERMS), int((n - 1) / 2) + 1):
            calphat[n] = calphat[n] + cahat[n - m - 1] * cgamma[n - m - 1, m] * (
                -1.0
            ) ** (n - 1 - 2 * m)

    # Compute coefficients of delta^p
    for m in range(0, order + 1):
        cc[m] = RBINOM[order, m] * zbar ** (order - m) * biglab**order

    if order > 0:
        for n in range(0, 2 * NTERMS + order):
            calpha[n] = complex(0.0, 0.0)
            for m in range(max(0, n - 2 * NTERMS + 1), min(n, order) + 1):
                calpha[n] = calpha[n] + cc[m] * calphat[n - m]
    else:
        calpha = calphat

    # Computation of integral
    g3 = complex(0.0, 0.0)
    # cd1minz = del1 / biglab - zbar  cd2minz = del2 / biglab - zbar
    # if ( abs(cd1minz) < 1.0e-8) cd1minz = cd1minz + 1.0e-8
    # if ( abs(cd2minz) < 1.0e-8) cd2minz = cd2minz + 1.0e-8
    # By definition log is conjugate of previous log this avoids problems with signs
    # along the line (and saves logs).
    cd1minz = np.conj(cd1minz)
    cd2minz = np.conj(cd2minz)
    cln1 = np.conj(cln1)
    cln2 = np.conj(cln2)
    for n in range(0, 2 * NTERMS + order):
        g3 = g3 - (calpha[n] * cln2 - calpha[n] / (n + 1)) * (cd2minz) ** (
            n + 1
        ) / float(n + 1)
        g3 = g3 + (calpha[n] * cln1 - calpha[n] / (n + 1)) * (cd1minz) ** (
            n + 1
        ) / float(n + 1)
    g3 = biglabin * g3 / (2.0 * np.pi)
    wdis = g1 + g2 + g3
    return wdis


@numba.njit(nogil=True, cache=True)
def IntegralLapLineDipole(zin, z1, z2, del0, ra, order):
    cg = np.full(order + 2, complex(0.0, 0.0))
    z = (2.0 * zin - (z1 + z2)) / (z2 - z1)
    zplus1 = z + 1.0
    zmin1 = z - 1.0

    # We don't always have to do this, so maybe put in condition?
    # Determine coefficients of powers of Delta
    for m in range(0, order + 1):
        cg[m] = RBINOM[order, m] * (-del0) ** (order - m) / ra**order

    zterm = complex(0.0, 0.0)
    for n in range(0, order + 1):
        zterm = zterm + cg[n] * z**n

    qmtot = complex(0.0, 0.0)
    for m in range(1, order + 1):
        qm = complex(0.0, 0.0)
        for n in range(1, int((m + 1) / 2) + 1):
            qm = qm + z ** (m - 2 * n + 1) / float(2 * n - 1)
        qmtot = qmtot + 2.0 * cg[m] * qm

    comega = (zterm * np.log(zmin1 / zplus1) + qmtot) / (2.0 * np.pi)
    return comega


@numba.njit(nogil=True, cache=True)
def IntegralLapLineDipoleDis(zin, z1, z2, del0, ra, order):
    cg = np.full(order + 2, complex(0.0, 0.0))

    z = (2.0 * zin - (z1 + z2)) / (z2 - z1)
    zplus1 = z + 1.0
    zmin1 = z - 1.0

    # Determine coefficients of powers of Delta for [ (Delta-Delta_0)/a ] ^ p
    for m in range(0, order + 1):
        cg[m] = RBINOM[order, m] * (-del0) ** (order - m) / ra**order

    zterm1 = complex(0.0, 0.0)
    zterm2 = complex(0.0, 0.0)
    for n in range(1, order + 1):
        zterm1 = zterm1 + cg[n] * float(n) * z ** (n - 1)
    for n in range(0, order + 1):
        zterm2 = zterm2 + cg[n] * z**n

    qmtot = complex(0.0, 0.0)
    for m in range(2, order + 1):
        qm = complex(0.0, 0.0)
        for n in range(1, int(m / 2) + 1):
            qm = qm + float(m - 2 * n + 1) * z ** (m - 2 * n) / float(2 * n - 1)
        qmtot = qmtot + 2.0 * cg[m] * qm

    wdis = (
        zterm1 * np.log(zmin1 / zplus1) + zterm2 * (1.0 / zmin1 - 1.0 / zplus1) + qmtot
    ) / (2.0 * np.pi)
    return wdis


@numba.njit(nogil=True, cache=True)
def findm1m2(zin, z1in, z2in, Lin, labda, Rconv):
    # Break integral up in sections of max one labda
    # and find first (m1) and last (m2) section within radius of convergence
    if labda == 0.0:
        NLS = 0
    else:
        NLS = int(np.ceil(Lin / labda))

    m1 = 0
    m2 = 0
    for j in range(1, NLS + 1):
        z1 = z1in + float(j - 1) / NLS * (z2in - z1in)
        z2 = z1 + (z2in - z1in) / NLS
        L = abs(z2 - z1)
        biglab = 2.0 * labda / L
        z = (2.0 * zin - (z1 + z2)) / (z2 - z1) / biglab
        if m1 == 0:
            if abs(z) < Rconv:
                m1 = j
        else:
            if abs(z) > Rconv:
                m2 = j - 1
                break
    if m2 == 0:
        m2 = NLS

    return m1, m2, NLS

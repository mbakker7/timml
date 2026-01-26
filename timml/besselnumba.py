import numba
import numpy as np

"""
real(kind=8) :: pi, tiny
real(kind=8), dimension(0:20) :: a, b, afar, a1, b1
real(kind=8), dimension(0:20) :: nrange
real(kind=8), dimension(0:20,0:20) :: gam
real(kind=8), dimension(8) :: xg, wg

initialize
----------
implicit none
real(kind=8) :: c, fac, twologhalf
real(kind=8), dimension(0:20) :: bot
real(kind=8), dimension(1:21) :: psi
integer :: n,m
"""

tiny = 1e-10
c = np.log(0.5) + 0.577215664901532860

fac = 1.0

nrange = np.arange(21, dtype=np.float64)

a = np.zeros(21, dtype=np.float64)
a[0] = 1.0
b = np.zeros(21, dtype=np.float64)

for n in range(1, 21):
    fac = n * fac
    a[n] = 1.0 / (4.0 ** nrange[n] * fac**2)
    b[n] = b[n - 1] + 1 / nrange[n]

b = (b - c) * a
a = -a / 2.0

gam = np.zeros((21, 21), dtype=np.float64)
for n in range(21):
    for m in range(n + 1):
        gam[n, m] = np.prod(nrange[m + 1 : n + 1]) / np.prod(nrange[1 : n - m + 1])

# gotta predefine these i.o. gam which is used in the old code
binom = np.zeros((21, 21), dtype=np.float64)
for n in range(21):
    for m in range(n + 1):
        binom[n, m] = np.prod(nrange[m + 1 : n + 1]) / np.prod(nrange[1 : n - m + 1])

afar = np.zeros(21, dtype=np.float64)
afar[0] = np.sqrt(np.pi / 2.0)

for n in range(1, 21):
    afar[n] = -((2.0 * n - 1.0) ** 2) / (n * 8) * afar[n - 1]

fac = 1.0
bot = np.zeros(21, dtype=np.float64)
bot[0] = 4.0
for n in range(1, 21):
    fac = n * fac
    bot[n] = fac * (n + 1) * fac * 4.0 ** (n + 1)

psi = np.zeros(21, dtype=np.float64)
for n in range(2, 22):
    psi[n - 1] = psi[n - 2] + 1 / (n - 1)
psi = psi - 0.577215664901532860

a1 = np.empty(21, dtype=np.float64)
b1 = np.empty(21, dtype=np.float64)
twologhalf = 2 * np.log(0.5)
for n in range(21):
    a1[n] = 1 / bot[n]
    b1[n] = (twologhalf - (2.0 * psi[n] + 1 / (n + 1))) / bot[n]


wg = np.zeros(8, dtype=np.float64)
xg = np.zeros(8, dtype=np.float64)

wg[0] = 0.101228536290378
wg[1] = 0.22238103445338
wg[2] = 0.31370664587789
wg[3] = 0.36268378337836
wg[4] = 0.36268378337836
wg[5] = 0.313706645877890
wg[6] = 0.22238103445338
wg[7] = 0.10122853629038

xg[0] = -0.960289856497536
xg[1] = -0.796666477413626
xg[2] = -0.525532409916329
xg[3] = -0.183434642495650
xg[4] = 0.183434642495650
xg[5] = 0.525532409916329
xg[6] = 0.796666477413626
xg[7] = 0.960289856497536


@numba.njit(nogil=True, cache=True)
def besselk0near(z, Nt):
    """besselk0near.

    implicit none
    complex(kind=8), intent(in) :: z
    integer, intent(in) :: Nt
    complex(kind=8) :: omega
    complex(kind=8) :: rsq, log1, term
    integer :: n
    """
    rsq = z**2
    term = 1.0 + 0.0j
    log1 = np.log(rsq)
    omega = a[0] * log1 + b[0]

    for n in range(1, Nt + 1):
        term = term * rsq
        omega = omega + (a[n] * log1 + b[n]) * term

    return omega


@numba.njit(nogil=True, cache=True)
def besselk0cheb(z, Nt):
    """besselk0cheb.

    implicit none
    complex(kind=8), intent(in) :: z
    integer, intent(in) :: Nt
    complex(kind=8) :: omega
    integer :: n, n2, ts
    real(kind=8) :: a, b, c, A3, u
    complex(kind=8) :: A1, A2, cn, cnp1, cnp2, cnp3
    complex(kind=8) :: z1, z2, S, T

    """
    cnp1 = complex(1.0, 0.0)
    cnp2 = complex(0.0, 0.0)
    cnp3 = complex(0.0, 0.0)
    a = 0.5
    c = 1.0
    b = 1.0 + a - c

    z1 = 2.0 * z
    z2 = 2.0 * z1
    ts = (-1) ** (Nt + 1)
    S = ts
    T = 1.0

    for n in range(Nt, -1, -1):
        u = (n + a) * (n + b)
        n2 = 2 * n
        A1 = 1.0 - (z2 + (n2 + 3.0) * (n + a + 1.0) * (n + b + 1.0) / (n2 + 4.0)) / u
        A2 = 1.0 - (n2 + 2.0) * (n2 + 3.0 - z2) / u
        A3 = -(n + 1.0) * (n + 3.0 - a) * (n + 3.0 - b) / (u * (n + 2.0))
        cn = (2.0 * n + 2.0) * A1 * cnp1 + A2 * cnp2 + A3 * cnp3
        ts = -ts
        S = S + ts * cn
        T = T + cn
        cnp3 = cnp2
        cnp2 = cnp1
        cnp1 = cn
    cn = cn / 2.0
    S = S - cn
    T = T - cn
    omega = 1.0 / np.sqrt(z1) * T / S
    omega = np.sqrt(np.pi) * np.exp(-z) * omega

    return omega


@numba.njit(nogil=True, cache=True)
def besselk0(x, y, lab):
    """besselk0.

    implicit none
    real(kind=8), intent(in) :: x,y
    complex(kind=8), intent(in) :: lab
    complex(kind=8) :: z, omega
    real(kind=8) :: cond
    """
    z = np.sqrt(x**2 + y**2) / lab
    cond = np.abs(z)

    if cond < 6:
        omega = besselk0near(z, 17)
    else:
        omega = besselk0cheb(z, 6)

    return omega


# zminzbar = np.zeros(21, dtype=np.complex_)
exprange = np.zeros(21, dtype=np.complex128)
anew = np.zeros(21, dtype=np.complex128)
bnew = np.zeros(21, dtype=np.complex128)


@numba.njit(nogil=True, cache=True)
def bessells_int(x, y, z1, z2, lab):
    """bessells_int.

    implicit none
    real(kind=8), intent(in) :: x,y
    complex(kind=8), intent(in) :: z1,z2,lab
    real(kind=8) :: biglab, biga, L, ang, tol
    complex(kind=8) :: zeta, zetabar, omega, log1, log2, term1, term2,
        d1minzeta, d2minzeta
    complex(kind=8), dimension(0:20) :: zminzbar, anew, bnew, exprange
    complex(kind=8), dimension(0:20,0:20) :: gamnew, gam2
    complex(kind=8), dimension(0:40) :: alpha, beta, alpha2
    integer :: n
    """
    zminzbar = np.zeros(21, dtype=np.complex128)

    L = np.abs(z2 - z1)
    biga = np.abs(lab)
    ang = np.arctan2(lab.imag, lab.real)
    biglab = 2 * biga / L

    tol = 1e-12

    exprange = np.exp(-complex(0, 2) * ang * nrange)
    anew = a * exprange
    bnew = (b - a * complex(0, 2) * ang) * exprange

    zeta = (2 * complex(x, y) - (z1 + z2)) / (z2 - z1) / biglab
    zetabar = np.conj(zeta)
    # #for n in range(21):
    # #    zminzbar[n] = (zeta-zetabar)**(20-n)  # Ordered from high power to low power
    zminzbar[20] = 1

    for n in range(1, 21):
        # Ordered from high power to low power
        zminzbar[20 - n] = zminzbar[21 - n] * (zeta - zetabar)

    gamnew = np.zeros((21, 21), dtype=np.complex128)
    gam2 = np.zeros((21, 21), dtype=np.complex128)
    for n in range(21):
        gamnew[n, 0 : n + 1] = gam[n, 0 : n + 1] * zminzbar[20 - n : 20 + 1]
        gam2[n, 0 : n + 1] = np.conj(gamnew[n, 0 : n + 1])

    alpha = np.zeros(41, dtype=np.complex128)
    beta = np.zeros(41, dtype=np.complex128)
    alpha2 = np.zeros(41, dtype=np.complex128)

    alpha[0] = anew[0]
    beta[0] = bnew[0]
    alpha2[0] = anew[0]

    for n in range(1, 21):
        alpha[n : 2 * n + 1] = alpha[n : 2 * n + 1] + anew[n] * gamnew[n, 0 : n + 1]
        beta[n : 2 * n + 1] = beta[n : 2 * n + 1] + bnew[n] * gamnew[n, 0 : n + 1]
        alpha2[n : 2 * n + 1] = alpha2[n : 2 * n + 1] + anew[n] * gam2[n, 0 : n + 1]

    omega = 0
    d1minzeta = -1 / biglab - zeta
    d2minzeta = 1 / biglab - zeta

    if np.abs(d1minzeta) < tol:
        d1minzeta = d1minzeta + complex(tol, 0)
    if np.abs(d2minzeta) < tol:
        d2minzeta = d2minzeta + complex(tol, 0)

    log1 = np.log(d1minzeta)
    log2 = np.log(d2minzeta)
    term1 = 1
    term2 = 1

    # I tried to serialize this, but it didn't speed things up
    for n in range(41):
        term1 = term1 * d1minzeta
        term2 = term2 * d2minzeta
        omega = omega + (alpha[n] * log2 - alpha[n] / (n + 1) + beta[n]) * term2 / (
            n + 1
        )
        omega = omega - (alpha[n] * log1 - alpha[n] / (n + 1) + beta[n]) * term1 / (
            n + 1
        )
        omega = omega + (alpha2[n] * np.conj(log2) - alpha2[n] / (n + 1)) * np.conj(
            term2
        ) / (n + 1)
        omega = omega - (alpha2[n] * np.conj(log1) - alpha2[n] / (n + 1)) * np.conj(
            term1
        ) / (n + 1)

    omega = -biga / (2 * np.pi) * omega

    return omega


@numba.njit(nogil=True, cache=True)
def bessells_gauss(x, y, z1, z2, lab):
    """bessells_gauss.

    implicit none
    real(kind=8), intent(in) :: x,y
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), intent(in) :: lab
    complex(kind=8) :: omega
    integer :: n
    real(kind=8) :: L, x0
    complex(kind=8) :: bigz, biglab
    """
    L = np.abs(z2 - z1)
    biglab = 2 * lab / L
    bigz = (2 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    omega = complex(0, 0)
    for n in range(1, 9):
        x0 = bigz.real - xg[n - 1]
        omega = omega + wg[n - 1] * besselk0(x0, bigz.imag, biglab)

    omega = -L / (4 * np.pi) * omega
    return omega


# @numba.njit(nogil=True, cache=True)
# def bessellsuni(x, y, z1, z2, lab):
#     """Bessellsuni.

#     # Uniform strength
#     implicit none
#     real(kind=8), intent(in) :: x,y
#     complex(kind=8), intent(in) :: z1,z2
#     complex(kind=8), intent(in) :: lab
#     complex(kind=8) :: omega

#     integer :: Nls, n
#     real(kind=8) :: Lnear, L
#     complex(kind=8) :: z, delz, za, zb
#     """
#     Lnear = 3.0
#     z = complex(x, y)
#     omega = complex(0.0, 0.0)
#     L = np.abs(z2 - z1)
#     if L < Lnear * np.abs(lab):  # No need to break integral up
#         if np.abs(z - 0.5 * (z1 + z2)) < 0.5 * Lnear * L:  # Do integration
#             omega = bessells_int(x, y, z1, z2, lab)
#         else:
#             omega = bessells_gauss(x, y, z1, z2, lab)
#     else:  # Break integral up in parts
#         Nls = int(np.ceil(L / (Lnear * np.abs(lab))))
#         delz = (z2 - z1) / Nls
#         L = np.abs(delz)
#         for n in range(1, Nls + 1):
#             za = z1 + (n - 1) * delz
#             zb = z1 + n * delz
#             if np.abs(z - 0.5 * (za + zb)) < 0.5 * Lnear * L:  # integration
#                 omega = omega + bessells_int(x, y, za, zb, lab)
#             else:
#                 omega = omega + bessells_gauss(x, y, za, zb, lab)
#     return omega


# @numba.njit(nogil=True, cache=True)
# def bessellsuniv(x, y, z1, z2, lab, rzero):
#     """Bessellsuniv.

#     # Uniform strength
#     implicit none
#     real(kind=8), intent(in) :: x,y
#     complex(kind=8), intent(in) :: z1,z2
#     integer, intent(in) :: nlab
#     complex(kind=8), dimension(nlab), intent(in) :: lab
#     complex(kind=8), dimension(nlab), intent(inout) :: omega
#     integer :: n
#     """
#     nlab = len(lab)
#     omega = np.zeros(nlab, dtype=np.complex128)
#     za, zb, N = circle_line_intersection(z1, z2, x + y * 1j, rzero * abs(lab[0]))
#     if N > 0:
#         for n in range(nlab):
#             omega[n] = bessellsuni(x, y, za, zb, lab[n])
#     return omega


@numba.njit(nogil=True, cache=True)
def circle_line_intersection(z1, z2, zc, R):
    """circle_line_intersection.

    implicit none
    complex(kind=8), intent(in) :: z1, z2, zc
    real(kind=8), intent(in) :: R
    real(kind=8), intent(inout) :: xouta, youta, xoutb, youtb
    integer, intent(inout) :: N
    real(kind=8) :: Lover2, d, xa, xb
    complex(kind=8) :: bigz, za, zb
    """
    N = 0
    za = complex(0, 0)
    zb = complex(0, 0)
    Lover2 = np.abs(z2 - z1) / 2
    bigz = (2 * zc - (z1 + z2)) * Lover2 / (z2 - z1)
    if abs(bigz.imag) < R:
        d = np.sqrt(R**2 - bigz.imag**2)
        xa = bigz.real - d
        xb = bigz.real + d
        if (xa < Lover2) and (xb > -Lover2):
            N = 2
            if xa < -Lover2:
                za = z1
            else:
                za = (xa * (z2 - z1) / Lover2 + (z1 + z2)) / 2.0
            if xb > Lover2:
                zb = z2
            else:
                zb = (xb * (z2 - z1) / Lover2 + (z1 + z2)) / 2.0
    return za, zb, N


@numba.njit(nogil=True, cache=True)
def bessellsv2(x, y, z1, z2, lab, order, R):
    """bessellsv2.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,R
    complex(kind=8), intent(in) :: z1,z2
    integer, intent(in) :: nlab
    real(kind=8) :: d1, d2
    complex(kind=8), dimension(nlab), intent(in) :: lab
    complex(kind=8), dimension(order+1,nlab) :: omega
    integer :: n, nterms
    """
    nlab = len(lab)
    nterms = order + 1
    omega = np.zeros((order + 1, nlab), dtype=np.complex128)
    # Check if endpoints need to be adjusted using the largest lambda (the first one)
    d1, d2 = find_d1d2(z1, z2, complex(x, y), R * np.abs(lab[0]))
    for n in range(nlab):
        omega[: nterms + 1, n] = bessells(x, y, z1, z2, lab[n], order, d1, d2)
    return omega


@numba.njit(nogil=True, cache=True)
def find_d1d2(z1, z2, zc, R):
    """find_d1d2.

    implicit none
    complex(kind=8), intent(in) :: z1, z2, zc
    real(kind=8), intent(in) :: R
    real(kind=8), intent(inout) :: d1, d2
    real(kind=8) :: Lover2, d, xa, xb
    complex(kind=8) :: bigz
    """
    d1 = -1.0
    d2 = 1.0
    Lover2 = np.abs(z2 - z1) / 2
    bigz = (2 * zc - (z1 + z2)) * Lover2 / (z2 - z1)
    if np.abs((bigz.imag)) < R:
        d = np.sqrt(R**2 - bigz.imag**2)
        xa = bigz.real - d
        xb = bigz.real + d
        if (xa < Lover2) and (xb > -Lover2):
            if xa < -Lover2:
                d1 = -1.0
            else:
                d1 = xa / Lover2
            if xb > Lover2:
                d2 = 1.0
            else:
                d2 = xb / Lover2
    return d1, d2


@numba.njit(nogil=True, cache=True)
def bessells(x, y, z1, z2, lab, order, d1in, d2in):
    """Bessells.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,d1in,d2in
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), intent(in) :: lab
    complex(kind=8), dimension(0:order) :: omega

    integer :: Nls, n
    real(kind=8) :: Lnear, L, d1, d2, delta
    complex(kind=8) :: z, delz, za, zb
    """
    omega = np.zeros(order + 1, dtype=np.complex128)
    Lnear = 3
    z = complex(x, y)
    L = np.abs(z2 - z1)
    if L < Lnear * np.abs(lab):  # No need to break integral up
        if np.abs(z - 0.5 * (z1 + z2)) < 0.5 * Lnear * L:  # Do integration
            omega = bessells_int_ho(x, y, z1, z2, lab, order, d1in, d2in)
        else:
            omega = bessells_gauss_ho_d1d2(x, y, z1, z2, lab, order, d1in, d2in)
    else:  # Break integral up in parts
        Nls = int(np.ceil(L / (Lnear * np.abs(lab))))
        delta = 2 / Nls
        delz = (z2 - z1) / Nls
        L = np.abs(delz)
        for n in range(1, Nls + 1):
            d1 = -1 + (n - 1) * delta
            d2 = -1 + n * delta
            if (d2 < d1in) or (d1 > d2in):
                continue
            d1 = np.max(np.array([d1, d1in]))
            d2 = np.min(np.array([d2, d2in]))
            za = z1 + (n - 1) * delz
            zb = z1 + n * delz
            if np.abs(z - 0.5 * (za + zb)) < 0.5 * Lnear * L:  # Do integration
                omega = omega + bessells_int_ho(x, y, z1, z2, lab, order, d1, d2)
            else:
                omega = omega + bessells_gauss_ho_d1d2(x, y, z1, z2, lab, order, d1, d2)
    return omega


@numba.njit(nogil=True, cache=True)
def bessells_gauss_ho(x, y, z1, z2, lab, order):
    """bessells_gauss_ho.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), intent(in) :: lab
    complex(kind=8), dimension(0:order) :: omega
    integer :: n, p
    real(kind=8) :: L, x0
    complex(kind=8) :: bigz, biglab
    complex(kind=8), dimension(8) :: k0
    """
    L = np.abs(z2 - z1)
    biglab = 2 * lab / L
    bigz = (2 * complex(x, y) - (z1 + z2)) / (z2 - z1)

    k0 = np.zeros(8, dtype=np.complex128)
    for n in range(8):
        x0 = bigz.real - xg[n]
        k0[n] = besselk0(x0, bigz.imag, biglab)

    omega = np.zeros(order + 1, dtype=np.complex128)
    for p in range(order + 1):
        omega[p] = complex(0, 0)
        for n in range(8):
            omega[p] = omega[p] + wg[n] * xg[n] ** p * k0[n]
        omega[p] = -L / (4 * np.pi) * omega[p]

    return omega


@numba.njit(nogil=True, cache=True)
def bessells_gauss_ho_d1d2(x, y, z1, z2, lab, order, d1, d2):
    """Returns integral from d1 to d2 along real axis.

    While strength is still Delta^order from -1 to +1.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,d1,d2
    complex(kind=8), intent(in) :: z1,z2,lab
    complex(kind=8), dimension(0:order) :: omega, omegac
    integer :: n, m
    real(kind=8) :: xp, yp, dc, fac
    complex(kind=8) :: z1p,z2p,bigz1,bigz2
    """
    omega = np.zeros(order + 1, dtype=np.complex128)
    bigz1 = complex(d1, 0)
    bigz2 = complex(d2, 0)
    z1p = 0.5 * (z2 - z1) * bigz1 + 0.5 * (z1 + z2)
    z2p = 0.5 * (z2 - z1) * bigz2 + 0.5 * (z1 + z2)
    omegac = bessells_gauss_ho(x, y, z1p, z2p, lab, order)
    dc = (d1 + d2) / (d2 - d1)
    for n in range(order + 1):
        for m in range(n + 1):
            omega[n] = omega[n] + gam[n, m] * dc ** (n - m) * omegac[m]
        omega[n] = (0.5 * (d2 - d1)) ** n * omega[n]
    return omega


@numba.njit(nogil=True, cache=True)
def isinside(z1, z2, zc, R):
    """Checks whether point zc is within oval with 'radius' R from line element.

    implicit none
    complex(kind=8), intent(in) :: z1, z2, zc
    real(kind=8), intent(in) :: R
    integer :: irv
    real(kind=8) :: Lover2, d, xa, xb
    complex(kind=8) :: bigz
    """
    irv = 0
    Lover2 = np.abs(z2 - z1) / 2
    bigz = (2 * zc - (z1 + z2)) * np.abs(z2 - z1) / (2 * (z2 - z1))
    if np.abs(bigz.imag) < R:
        d = np.sqrt(R**2 - bigz.imag**2)
        xa = bigz.real - d
        xb = bigz.real + d
        if (xa < Lover2) and (xb > -Lover2):
            irv = 1
    return irv


@numba.njit(nogil=True, cache=True)
def bessellsqxqyv2(x, y, z1, z2, lab, order, R):
    """bessellsqxqyv2.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,R
    complex(kind=8), intent(in) :: z1,z2
    integer, intent(in) :: nlab
    real(kind=8) :: d1, d2
    complex(kind=8), dimension(nlab), intent(in) :: lab
    complex(kind=8), dimension(2*(order+1),nlab) :: qxqy
    complex(kind=8), dimension(0:2*order+1) :: qxqylab
    integer :: n, nterms, nhalf
    """
    nlab = len(lab)
    qxqy = np.zeros((2 * (order + 1), nlab), dtype=np.complex128)
    nterms = order + 1
    # nhalf = nlab * (order + 1)
    d1, d2 = find_d1d2(z1, z2, complex(x, y), R * np.abs(lab[0]))
    for n in range(nlab):
        qxqylab = bessellsqxqy(x, y, z1, z2, lab[n], order, d1, d2)
        qxqy[:nterms, n] = qxqylab[0 : order + 1]
        qxqy[nterms : 2 * nterms, n] = qxqylab[order + 1 : 2 * (order + 1)]
    return qxqy


@numba.njit(nogil=True, cache=True)
def bessellsqxqy(x, y, z1, z2, lab, order, d1in, d2in):
    """Bessellsqxqy.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,d1in,d2in
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), intent(in) :: lab
    complex(kind=8), dimension(0:2*order+1) :: qxqy

    integer :: Nls, n
    real(kind=8) :: Lnear, L, d1, d2, delta
    complex(kind=8) :: z, delz, za, zb
    """
    Lnear = 3.0
    z = complex(x, y)
    qxqy = np.zeros(2 * order + 2, dtype=np.complex128)
    L = np.abs(z2 - z1)
    # print *,'Lnear*np.abs(lab) ',Lnear*np.abs(lab)
    if L < Lnear * np.abs(lab):  # No need to break integral up
        if np.abs(z - 0.5 * (z1 + z2)) < 0.5 * Lnear * L:  # Do integration
            qxqy = bessells_int_ho_qxqy(x, y, z1, z2, lab, order, d1in, d2in)
        else:
            qxqy = bessells_gauss_ho_qxqy_d1d2(x, y, z1, z2, lab, order, d1in, d2in)

    else:  # Break integral up in parts
        Nls = int(np.ceil(L / (Lnear * np.abs(lab))))
        # print *,'NLS ',Nls
        delta = 2.0 / Nls
        delz = (z2 - z1) / Nls
        L = np.abs(delz)
        for n in range(1, Nls + 1):
            d1 = -1.0 + (n - 1) * delta
            d2 = -1.0 + n * delta
            if (d2 < d1in) or (d1 > d2in):
                continue
            d1 = np.max(np.array([d1, d1in]))
            d2 = np.min(np.array([d2, d2in]))
            za = z1 + (n - 1) * delz
            zb = z1 + n * delz
            if np.abs(z - 0.5 * (za + zb)) < 0.5 * Lnear * L:  # Do integration
                qxqy = qxqy + bessells_int_ho_qxqy(x, y, z1, z2, lab, order, d1, d2)
            else:
                qxqy = qxqy + bessells_gauss_ho_qxqy_d1d2(
                    x, y, z1, z2, lab, order, d1, d2
                )
    return qxqy


@numba.njit(nogil=True, cache=True)
def bessells_gauss_ho_qxqy_d1d2(x, y, z1, z2, lab, order, d1, d2):
    """Returns integral from d1 to d2 along real axis.

    While strength is still Delta^order from -1 to +1.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,d1,d2
    complex(kind=8), intent(in) :: z1,z2,lab
    complex(kind=8), dimension(0:2*order+1) :: qxqy, qxqyc
    integer :: n, m
    real(kind=8) :: xp, yp, dc, fac
    complex(kind=8) :: z1p,z2p,bigz1,bigz2
    """
    qxqy = np.zeros(2 * order + 2, dtype=np.complex128)

    bigz1 = complex(d1, 0.0)
    bigz2 = complex(d2, 0.0)
    z1p = 0.5 * (z2 - z1) * bigz1 + 0.5 * (z1 + z2)
    z2p = 0.5 * (z2 - z1) * bigz2 + 0.5 * (z1 + z2)
    qxqyc = bessells_gauss_ho_qxqy(x, y, z1p, z2p, lab, order)
    dc = (d1 + d2) / (d2 - d1)
    for n in range(order + 1):
        for m in range(n + 1):
            qxqy[n] = qxqy[n] + gam[n, m] * dc ** (n - m) * qxqyc[m]
            qxqy[n + order + 1] = (
                qxqy[n + order + 1] + gam[n, m] * dc ** (n - m) * qxqyc[m + order + 1]
            )
        qxqy[n] = (0.5 * (d2 - d1)) ** n * qxqy[n]
        qxqy[n + order + 1] = (0.5 * (d2 - d1)) ** n * qxqy[n + order + 1]

    return qxqy


@numba.njit(nogil=True, cache=True)
def lapls_int_ho(x, y, z1, z2, order):
    """lapls_int_ho.

    ! Near field only
    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), dimension(0:order) :: omega, qm
    integer :: m, n
    real(kind=8) :: L
    complex(kind=8) :: z, zplus1, zmin1
    """
    omega = np.zeros(order + 1, dtype=np.complex128)
    L = np.abs(z2 - z1)
    z = (2.0 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    zplus1 = z + 1.0
    zmin1 = z - 1.0
    if np.abs(zplus1) < tiny:
        zplus1 = tiny
    if np.abs(zmin1) < tiny:
        zmin1 = tiny

    qm = np.zeros(order + 2, dtype=np.complex128)
    qm[1] = 2.0
    for m in range(3, order + 2, 2):
        qm[m] = qm[m - 2] * z * z + 2.0 / m
    for m in range(2, order + 2, 2):
        qm[m] = qm[m - 1] * z

    logterm = np.log(zmin1 / zplus1)
    logzmin1 = np.log(zmin1)
    logzplus1 = np.log(zplus1)
    for p in range(order + 1):
        omega[p] = (
            z ** (p + 1) * logterm + qm[p + 1] - logzmin1 + (-1) ** (p + 1) * logzplus1
        )
        omega[p] = -L / (4 * np.pi * (p + 1)) * omega[p]
    return omega.real


@numba.njit(nogil=True, cache=True)
def lapls_int_ho_wdis(x, y, z1, z2, order):
    """Note this is W andReturns Qx - iQy."""
    wdis = np.zeros(order + 1, dtype=np.complex128)
    L = np.abs(z2 - z1)
    z = (2.0 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    zplus1 = z + 1.0
    zmin1 = z - 1.0
    if np.abs(zplus1) < tiny:
        zplus1 = tiny
    if np.abs(zmin1) < tiny:
        zmin1 = tiny

    qm = np.zeros(order + 2, dtype=np.complex128)
    qm[0:1] = 0.0
    for m in range(2, order + 2):
        for n in range(1, m // 2 + 1):
            qm[m] = qm[m] + (m - 2 * n + 1) * z ** (m - 2 * n) / (2 * n - 1)
        qm[m] = 2 * qm[m]

    termzmin = 1.0 / zmin1
    termzplus = 1.0 / zplus1
    termlog = np.log(zmin1 / zplus1)
    for p in range(0, order + 1):
        wdis[p] = (p + 1) * z**p * termlog + z ** (p + 1) * (termzmin - termzplus)
        wdis[p] = wdis[p] + qm[p + 1] - termzmin + (-1) ** (p + 1) * termzplus
        wdis[p] = L / (2 * np.pi * (z2 - z1) * (p + 1)) * wdis[p]
    return wdis


@numba.njit(nogil=True, cache=True)
def lapld_int_ho_d1d2(x, y, z1, z2, order, d1, d2):
    """lapld_int_ho_d1d2.

    Near field only
    Returns integral from d1 to d2 along real axis while strength is still
    Delta^order from -1 to +1
    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,d1,d2
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), dimension(0:order) :: omega, omegac
    integer :: n, m
    real(kind=8) :: xp, yp, dc, fac
    complex(kind=8) :: z1p,z2p,bigz1,bigz2
    """
    omega = np.zeros(order + 1, dtype=np.complex128)

    bigz1 = complex(d1, 0.0)
    bigz2 = complex(d2, 0.0)
    z1p = 0.5 * (z2 - z1) * bigz1 + 0.5 * (z1 + z2)
    z2p = 0.5 * (z2 - z1) * bigz2 + 0.5 * (z1 + z2)
    omegac = lapld_int_ho(x, y, z1p, z2p, order)
    dc = (d1 + d2) / (d2 - d1)
    for n in range(order + 1):
        for m in range(n + 1):
            omega[n] = omega[n] + gam[n, m] * dc ** (n - m) * omegac[m]
        omega[n] = (0.5 * (d2 - d1)) ** n * omega[n]

    return omega


@numba.njit(nogil=True, cache=True)
def lapld_int_ho(x, y, z1, z2, order):
    """lapld_int_ho.

    ! Near field only
    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), dimension(0:order) :: omega, qm
    integer :: m, n
    real(kind=8) :: L
    complex(kind=8) :: z, zplus1, zmin1
    """
    omega = np.zeros(order + 1, dtype=np.complex128)
    qm = np.zeros(order + 1, dtype=np.complex128)

    # L = np.abs(z2 - z1)
    z = (2.0 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    zplus1 = z + 1.0
    zmin1 = z - 1.0
    # Not sure if this gives correct answer at corner point (z also appears in qm);
    # should really be caught in code that calls this function
    if np.abs(zplus1) < tiny:
        zplus1 = tiny
    if np.abs(zmin1) < tiny:
        zmin1 = tiny

    omega[0] = np.log(zmin1 / zplus1)
    for n in range(1, order + 1):
        omega[n] = z * omega[n - 1]

    if order > 0:
        qm[1] = 2.0
    for m in range(3, order + 1, 2):
        qm[m] = qm[m - 2] * z * z + 2.0 / m

    for m in range(2, order + 1, 2):
        qm[m] = qm[m - 1] * z

    omega = 1.0 / (complex(0.0, 2.0) * np.pi) * (omega + qm)
    return omega


@numba.njit(nogil=True, cache=True)
def bessells_gauss_ho_qxqy(x, y, z1, z2, lab, order):
    """bessells_gauss_ho_qxqy.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), intent(in) :: lab
    complex(kind=8), dimension(0:2*order+1) :: qxqy
    integer :: n, p
    real(kind=8) :: L, bigy, angz
    complex(kind=8) :: bigz, biglab
    real(kind=8), dimension(8) :: r, xmind
    complex(kind=8), dimension(8) :: k1
    complex(kind=8), dimension(0:order) :: qx,qy
    """
    qxqy = np.zeros(2 * order + 2, dtype=np.complex128)
    xmind = np.zeros(8, dtype=np.complex128)
    k1 = np.zeros(8, dtype=np.complex128)
    r = np.zeros(8, dtype=np.complex128)

    L = np.abs(z2 - z1)
    biglab = 2 * lab / L
    bigz = (2 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    bigy = bigz.imag
    for n in range(8):
        xmind[n] = bigz.real - xg[n]
        r[n] = np.sqrt(xmind[n] ** 2 + bigz.imag**2)
        k1[n] = besselk1(xmind[n], bigz.imag, biglab)

    qx = np.zeros(order + 1, dtype=np.complex128)
    qy = np.zeros(order + 1, dtype=np.complex128)
    for p in range(order + 1):
        for n in range(8):
            qx[p] = qx[p] + wg[n] * xg[n] ** p * xmind[n] * k1[n] / r[n]
            qy[p] = qy[p] + wg[n] * xg[n] ** p * bigy * k1[n] / r[n]

    qx = -qx * L / (4 * np.pi * biglab) * 2 / L
    qy = -qy * L / (4 * np.pi * biglab) * 2 / L

    angz = np.arctan2((z2 - z1).imag, (z2 - z1).real)
    qxqy[0 : order + 1] = qx * np.cos(angz) - qy * np.sin(angz)
    qxqy[order + 1 : 2 * order + 2] = qx * np.sin(angz) + qy * np.cos(angz)

    return qxqy


@numba.njit(nogil=True, cache=True)
def besselk1cheb(z, Nt):
    """besselk1cheb.

    implicit none
    complex(kind=8), intent(in) :: z
    integer, intent(in) :: Nt
    complex(kind=8) :: omega
    integer :: n, n2, ts
    real(kind=8) :: a, b, c, A3, u
    complex(kind=8) :: A1, A2, cn, cnp1, cnp2, cnp3
    complex(kind=8) :: z1, z2, S, T
    """
    cnp1 = 1.0 + 0.0j
    cnp2 = 0.0 + 0.0j
    cnp3 = 0.0 + 0.0j

    a = 1.5
    c = 3.0
    b = 1.0 + a - c

    z1 = 2 * z
    z2 = 2 * z1
    ts = (-1) ** (Nt + 1)
    S = ts
    T = 1.0

    for n in range(Nt, -1, -1):
        u = (n + a) * (n + b)
        n2 = 2 * n
        A1 = 1 - (z2 + (n2 + 3) * (n + a + 1) * (n + b + 1) / (n2 + 4)) / u
        A2 = 1 - (n2 + 2) * (n2 + 3 - z2) / u
        A3 = -(n + 1) * (n + 3 - a) * (n + 3 - b) / (u * (n + 2))
        cn = (2 * n + 2) * A1 * cnp1 + A2 * cnp2 + A3 * cnp3
        ts = -ts
        S = S + ts * cn
        T = T + cn
        cnp3 = cnp2
        cnp2 = cnp1
        cnp1 = cn

    cn = cn / 2
    S = S - cn
    T = T - cn
    omega = 1 / (np.sqrt(z1) * z1) * T / S
    omega = 2 * z * np.sqrt(np.pi) * np.exp(-z) * omega

    return omega


@numba.njit(nogil=True, cache=True)
def besselk1(x, y, lab):
    """besselk1.

    implicit none
    real(kind=8), intent(in) :: x,y
    complex(kind=8), intent(in) :: lab
    complex(kind=8) :: z, omega
    real(kind=8) :: cond
    """
    z = np.sqrt(x**2 + y**2) / lab
    cond = np.abs(z)

    if cond < 6:
        omega = besselk1near(z, 20)
    else:
        omega = besselk1cheb(z, 6)

    return omega


@numba.njit(nogil=True, cache=True)
def besselk1near(z, Nt):
    """besselk1near.

    implicit none
    complex(kind=8), intent(in) :: z
    integer, intent(in) :: Nt
    complex(kind=8) :: omega
    complex(kind=8) :: zsq, log1, term
    integer :: n
    """
    zsq = z**2
    term = z
    log1 = np.log(zsq)
    omega = 1.0 / z + (a1[0] * log1 + b1[0]) * z

    for n in range(1, Nt + 1):
        term = term * zsq
        omega = omega + (a1[n] * log1 + b1[n]) * term

    return omega


@numba.njit(nogil=True, cache=True)
def besselldv2(x, y, z1, z2, lab, order, R):
    """besselldv2.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,R
    complex(kind=8), intent(in) :: z1,z2
    integer, intent(in) :: nlab
    real(kind=8) :: d1, d2
    complex(kind=8), dimension(nlab), intent(in) :: lab
    complex(kind=8), dimension(order+1,nlab) :: omega
    integer :: n, nterms
    """
    nlab = len(lab)
    omega = np.zeros((order + 1, nlab), dtype=np.complex128)

    nterms = order + 1
    # Check if endpoints need to be adjusted using the largest lambda (the first one)
    d1, d2 = find_d1d2(z1, z2, complex(x, y), R * np.abs(lab[0]))
    for n in range(nlab):
        omega[: nterms + 1, n] = besselld(x, y, z1, z2, lab[n], order, d1, d2)

    return omega


@numba.njit(nogil=True, cache=True)
def besselld(x, y, z1, z2, lab, order, d1in, d2in):
    """Besselld.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,d1in,d2in
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), intent(in) :: lab
    complex(kind=8), dimension(0:order) :: omega

    integer :: Nls, n
    real(kind=8) :: Lnear, L, d1, d2, delta
    complex(kind=8) :: z, delz, za, zb
    """
    omega = np.zeros(order + 1, dtype=np.complex128)

    Lnear = 3.0
    z = complex(x, y)
    L = np.abs(z2 - z1)
    if L < Lnear * np.abs(lab):  # No need to break integral up
        if np.abs(z - 0.5 * (z1 + z2)) < 0.5 * Lnear * L:  # Do integration
            omega = besselld_int_ho(x, y, z1, z2, lab, order, d1in, d2in)
        else:
            omega = besselld_gauss_ho_d1d2(x, y, z1, z2, lab, order, d1in, d2in)
    else:  # Break integral up in parts
        Nls = int(np.ceil(L / (Lnear * np.abs(lab))))
        delta = 2.0 / Nls
        delz = (z2 - z1) / Nls
        L = np.abs(delz)
        for n in range(1, Nls + 1):
            d1 = -1.0 + (n - 1) * delta
            d2 = -1.0 + n * delta
            if (d2 < d1in) or (d1 > d2in):
                continue
            d1 = max(d1, d1in)
            d2 = min(d2, d2in)
            za = z1 + (n - 1) * delz
            zb = z1 + n * delz
            if np.abs(z - 0.5 * (za + zb)) < 0.5 * Lnear * L:  # Do integration
                omega = omega + besselld_int_ho(x, y, z1, z2, lab, order, d1, d2)
            else:
                omega = omega + besselld_gauss_ho_d1d2(x, y, z1, z2, lab, order, d1, d2)
    return omega


@numba.njit(nogil=True, cache=True)
def besselld_gauss_ho_d1d2(x, y, z1, z2, lab, order, d1, d2):
    """besselld_gauss_ho_d1d2.

    # Returns integral from d1 to d2 along real axis while strength is still
    # Delta^order from -1 to +1
    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,d1,d2
    complex(kind=8), intent(in) :: z1,z2,lab
    complex(kind=8), dimension(0:order) :: omega, omegac
    integer :: n, m
    real(kind=8) :: xp, yp, dc, fac
    complex(kind=8) :: z1p,z2p,bigz1,bigz2
    """
    omega = np.zeros(order + 1, dtype=np.complex128)

    bigz1 = complex(d1, 0.0)
    bigz2 = complex(d2, 0.0)
    z1p = 0.5 * (z2 - z1) * bigz1 + 0.5 * (z1 + z2)
    z2p = 0.5 * (z2 - z1) * bigz2 + 0.5 * (z1 + z2)
    omegac = besselld_gauss_ho(x, y, z1p, z2p, lab, order)
    dc = (d1 + d2) / (d2 - d1)
    omega[0 : order + 1] = 0.0
    for n in range(order + 1):
        for m in range(n + 1):
            omega[n] = omega[n] + gam[n, m] * dc ** (n - m) * omegac[m]
        omega[n] = (0.5 * (d2 - d1)) ** n * omega[n]
    return omega


@numba.njit(nogil=True, cache=True)
def besselld_gauss_ho(x, y, z1, z2, lab, order):
    """besselld_gauss_ho.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), intent(in) :: lab
    complex(kind=8), dimension(0:order) :: omega
    integer :: n, p
    real(kind=8) :: L, x0, r
    complex(kind=8) :: bigz, biglab
    complex(kind=8), dimension(8) :: k1overr
    """
    k1overr = np.zeros(8, dtype=np.complex128)
    omega = np.zeros(order + 1, dtype=np.complex128)

    L = np.abs(z2 - z1)
    biglab = 2.0 * lab / L
    bigz = (2.0 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    for n in range(8):
        x0 = bigz.real - xg[n]
        r = np.sqrt(x0**2 + bigz.imag**2)
        k1overr[n] = besselk1(x0, bigz.imag, biglab) / r
    for p in range(order + 1):
        omega[p] = complex(0.0, 0.0)
        for n in range(8):
            omega[p] = omega[p] + wg[n] * xg[n] ** p * k1overr[n]
        omega[p] = bigz.imag / (2.0 * np.pi * biglab) * omega[p]
    return omega


@numba.njit(nogil=True, cache=True)
def besselldqxqyv2(x, y, z1, z2, lab, order, R):
    """besselldqxqyv2.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,R
    complex(kind=8), intent(in) :: z1,z2
    integer, intent(in) :: nlab
    real(kind=8) :: d1, d2
    complex(kind=8), dimension(nlab), intent(in) :: lab
    complex(kind=8), dimension(2*(order+1),nlab) :: qxqy
    complex(kind=8), dimension(0:2*order+1) :: qxqylab
    integer :: n, nterms, nhalf
    """
    nlab = len(lab)
    qxqy = np.zeros((2 * (order + 1), nlab), dtype=np.complex128)
    nterms = order + 1
    # nhalf = nlab * (order + 1)
    d1, d2 = find_d1d2(z1, z2, complex(x, y), R * np.abs(lab[0]))
    for n in range(nlab):
        qxqylab = besselldqxqy(x, y, z1, z2, lab[n], order, d1, d2)
        qxqy[:nterms, n] = qxqylab[0 : order + 1]
        qxqy[nterms : 2 * nterms, n] = qxqylab[order + 1 : 2 * order + 1 + 1]
    return qxqy


@numba.njit(nogil=True, cache=True)
def besselldqxqy(x, y, z1, z2, lab, order, d1in, d2in):
    """Besselldqxqy.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,d1in,d2in
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), intent(in) :: lab
    complex(kind=8), dimension(0:2*order+1) :: qxqy

    integer :: Nls, n
    real(kind=8) :: Lnear, L, d1, d2, delta
    complex(kind=8) :: z, delz, za, zb

    """
    Lnear = 3
    z = complex(x, y)
    qxqy = np.zeros(2 * order + 2, dtype=np.complex128)

    L = np.abs(z2 - z1)

    # print *,'Lnear*np.abs(lab) ',Lnear*np.abs(lab)
    if L < Lnear * np.abs(lab):  # No need to break integral up
        if np.abs(z - 0.5 * (z1 + z2)) < 0.5 * Lnear * L:  # Do integration
            qxqy = besselld_int_ho_qxqy(x, y, z1, z2, lab, order, d1in, d2in)
        else:
            qxqy = besselld_gauss_ho_qxqy_d1d2(x, y, z1, z2, lab, order, d1in, d2in)

    else:  # Break integral up in parts
        Nls = int(np.ceil(L / (Lnear * np.abs(lab))))
        # print *,'NLS ',Nls
        delta = 2 / Nls
        delz = (z2 - z1) / Nls
        L = np.abs(delz)
        for n in range(1, Nls + 1):
            d1 = -1 + (n - 1) * delta
            d2 = -1 + n * delta
            if (d2 < d1in) or (d1 > d2in):
                continue
            d1 = np.max(np.array([d1, d1in]))
            d2 = np.min(np.array([d2, d2in]))
            za = z1 + (n - 1) * delz
            zb = z1 + n * delz
            if np.abs(z - 0.5 * (za + zb)) < 0.5 * Lnear * L:  # Do integration
                qxqy = qxqy + besselld_int_ho_qxqy(x, y, z1, z2, lab, order, d1, d2)
            else:
                qxqy = qxqy + besselld_gauss_ho_qxqy_d1d2(
                    x, y, z1, z2, lab, order, d1, d2
                )
    return qxqy


@numba.njit(nogil=True, cache=True)
def besselld_gauss_ho_qxqy_d1d2(x, y, z1, z2, lab, order, d1, d2):
    """Returns integral from d1 to d2 along real axis.

    While strength is still Delta^order from -1 to +1.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,d1,d2
    complex(kind=8), intent(in) :: z1,z2,lab
    complex(kind=8), dimension(0:2*order+1) :: qxqy, qxqyc
    integer :: n, m
    real(kind=8) :: xp, yp, dc, fac
    complex(kind=8) :: z1p,z2p,bigz1,bigz2
    """
    qxqy = np.zeros(2 * order + 2, dtype=np.complex128)

    bigz1 = complex(d1, 0)
    bigz2 = complex(d2, 0)
    z1p = 0.5 * (z2 - z1) * bigz1 + 0.5 * (z1 + z2)
    z2p = 0.5 * (z2 - z1) * bigz2 + 0.5 * (z1 + z2)
    qxqyc = besselld_gauss_ho_qxqy(x, y, z1p, z2p, lab, order)
    dc = (d1 + d2) / (d2 - d1)
    for n in range(order + 1):
        for m in range(n + 1):
            qxqy[n] = qxqy[n] + gam[n, m] * dc ** (n - m) * qxqyc[m]
            qxqy[n + order + 1] = (
                qxqy[n + order + 1] + gam[n, m] * dc ** (n - m) * qxqyc[m + order + 1]
            )

        qxqy[n] = (0.5 * (d2 - d1)) ** n * qxqy[n]
        qxqy[n + order + 1] = (0.5 * (d2 - d1)) ** n * qxqy[n + order + 1]

    return qxqy


@numba.njit(nogil=True, cache=True)
def besselld_gauss_ho_qxqy(x, y, z1, z2, lab, order):
    """besselld_gauss_ho_qxqy.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), intent(in) :: lab
    complex(kind=8), dimension(0:2*order+1) :: qxqy
    integer :: n, p
    real(kind=8) :: L, bigy, angz
    complex(kind=8) :: bigz, biglab
    real(kind=8), dimension(8) :: r, xmind
    complex(kind=8), dimension(8) :: k0,k1
    complex(kind=8), dimension(0:order) :: qx,qy
    """
    xmind = np.zeros(8, dtype=np.float64)
    r = np.zeros(8, dtype=np.float64)
    k0 = np.zeros(8, dtype=np.complex128)
    k1 = np.zeros(8, dtype=np.complex128)
    qxqy = np.zeros(2 * order + 2, dtype=np.complex128)

    L = np.abs(z2 - z1)
    biglab = 2.0 * lab / L
    bigz = (2.0 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    bigy = bigz.imag
    for n in range(8):
        xmind[n] = bigz.real - xg[n]
        r[n] = np.sqrt(xmind[n] ** 2 + bigz.imag**2)
        k0[n] = besselk0(xmind[n], bigz.imag, biglab)
        k1[n] = besselk1(xmind[n], bigz.imag, biglab)

    qx = np.zeros(order + 1, dtype=np.complex128)
    qy = np.zeros(order + 1, dtype=np.complex128)
    for p in range(order + 1):
        for n in range(8):
            qx[p] = qx[p] + wg[n] * xg[n] ** p * (-bigy) * xmind[n] / r[n] ** 3 * (
                r[n] * k0[n] / biglab + 2.0 * k1[n]
            )
            qy[p] = qy[p] + wg[n] * xg[n] ** p * (
                k1[n] / r[n]
                - bigy**2 / r[n] ** 3 * (r[n] * k0[n] / biglab + 2.0 * k1[n])
            )

    qx = -qx / (2 * np.pi * biglab) * 2 / L
    qy = -qy / (2 * np.pi * biglab) * 2 / L

    angz = np.arctan2((z2 - z1).imag, (z2 - z1).real)
    qxqy[0 : order + 1] = qx * np.cos(angz) - qy * np.sin(angz)
    qxqy[order + 1 : 2 * order + 1 + 1] = qx * np.sin(angz) + qy * np.cos(angz)

    return qxqy


@numba.njit(nogil=True, cache=True)
def besselldpart(x, y, z1, z2, lab, order, d1, d2):
    """Besselldpart.

    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,d1,d2
    complex(kind=8), intent(in) :: z1,z2,lab
    complex(kind=8), dimension(0:order) :: omega
    real(kind=8) :: biglab, biga, L, ang, tol, bigy
    complex(kind=8) :: zeta, zetabar, log1, log2, term1, term2, d1minzeta,
        d2minzeta, bigz
    complex(kind=8) :: cm, biglabcomplex
    complex(kind=8), dimension(0:20) :: zminzbar, anew, bnew, exprange
    complex(kind=8), dimension(0:20,0:20) :: gamnew, gam2
    complex(kind=8), dimension(0:40) :: alpha, beta, alpha2
    complex(kind=8), dimension(0:50) :: alphanew, betanew, alphanew2 ! Order fixed to 10
    integer :: m, n, p
    """
    zminzbar = np.zeros(21, dtype=np.complex128)

    L = np.abs(z2 - z1)
    # bigz = (2.0 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    # bigy = bigz.imag
    biga = np.abs(lab)
    ang = np.arctan2(lab.imag, lab.real)
    biglab = 2.0 * biga / L
    biglabcomplex = 2.0 * lab / L

    tol = 1e-12

    exprange = np.exp(-complex(0, 2) * ang * nrange)
    anew = a1 * exprange
    bnew = (b1 - a1 * complex(0, 2) * ang) * exprange

    zeta = (2.0 * complex(x, y) - (z1 + z2)) / (z2 - z1) / biglab
    zetabar = np.conj(zeta)
    zminzbar[-1] = 1.0
    for n in range(1, 21):
        # Ordered from high power to low power
        zminzbar[20 - n] = zminzbar[21 - n] * (zeta - zetabar)

    gamnew = np.zeros((21, 21), dtype=np.complex128)
    gam2 = np.zeros((21, 21), dtype=np.complex128)
    for n in range(21):
        gamnew[n, 0 : n + 1] = gam[n, 0 : n + 1] * zminzbar[20 - n : 20 + 1]
        gam2[n, 0 : n + 1] = np.conj(gamnew[n, 0 : n + 1])

    alpha = np.zeros(41, dtype=np.complex128)
    beta = np.zeros(41, dtype=np.complex128)
    alpha2 = np.zeros(41, dtype=np.complex128)
    alpha[0] = anew[0]
    beta[0] = bnew[0]
    alpha2[0] = anew[0]
    for n in range(1, 21):
        alpha[n : 2 * n + 1] = alpha[n : 2 * n + 1] + anew[n] * gamnew[n, 0 : n + 1]
        beta[n : 2 * n + 1] = beta[n : 2 * n + 1] + bnew[n] * gamnew[n, 0 : n + 1]
        alpha2[n : 2 * n + 1] = alpha2[n : 2 * n + 1] + anew[n] * gam2[n, 0 : n + 1]

    d1minzeta = d1 / biglab - zeta
    d2minzeta = d2 / biglab - zeta

    if np.abs(d1minzeta) < tol:
        d1minzeta = d1minzeta + complex(tol, 0.0)
    if np.abs(d2minzeta) < tol:
        d2minzeta = d2minzeta + complex(tol, 0.0)
    log1 = np.log(d1minzeta)
    log2 = np.log(d2minzeta)

    alphanew = np.zeros(51, dtype=np.complex128)
    alphanew2 = np.zeros(51, dtype=np.complex128)
    betanew = np.zeros(51, dtype=np.complex128)
    omega = np.zeros(order + 1, dtype=np.complex128)

    for p in range(order + 1):
        alphanew[0 : 40 + p + 1] = 0.0
        betanew[0 : 40 + p + 1] = 0.0
        alphanew2[0 : 40 + p + 1] = 0.0
        for m in range(p + 1):
            cm = biglab**p * gam[p, m] * zeta ** (p - m)
            alphanew[m : 40 + m + 1] = alphanew[m : 40 + m + 1] + cm * alpha[0 : 40 + 1]
            betanew[m : 40 + m + 1] = betanew[m : 40 + m + 1] + cm * beta[0 : 40 + 1]
            cm = biglab**p * gam[p, m] * zetabar ** (p - m)
            alphanew2[m : 40 + m + 1] = (
                alphanew2[m : 40 + m + 1] + cm * alpha2[0 : 40 + 1]
            )

        omega[p] = 0.0
        term1 = 1.0 + 0j
        term2 = 1.0 + 0j
        for n in range(40 + p + 1):
            term1 = term1 * d1minzeta
            term2 = term2 * d2minzeta
            omega[p] = omega[p] + (
                alphanew[n] * log2 - alphanew[n] / (n + 1) + betanew[n]
            ) * term2 / (n + 1)
            omega[p] = omega[p] - (
                alphanew[n] * log1 - alphanew[n] / (n + 1) + betanew[n]
            ) * term1 / (n + 1)
            omega[p] = omega[p] + (
                alphanew2[n] * np.conj(log2) - alphanew2[n] / (n + 1)
            ) * np.conj(term2) / (n + 1)
            omega[p] = omega[p] - (
                alphanew2[n] * np.conj(log1) - alphanew2[n] / (n + 1)
            ) * np.conj(term1) / (n + 1)

    # + real( lapld_int_ho(x,y,z1,z2,order) )
    omega = biglab / (2.0 * np.pi * biglabcomplex**2) * omega
    # omega = real( lapld_int_ho(x,y,z1,z2,order) )

    return omega


@numba.njit(nogil=True, cache=True)
def lapld_int_ho_wdis_d1d2(x, y, z1, z2, order, d1, d2):
    """lapld_int_ho_wdis_d1d2.

    # Near field only
    # Returns integral from d1 to d2 along real axis while strength is still
    # Delta^order from -1 to +1
    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y,d1,d2
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), dimension(0:order) :: wdis, wdisc
    integer :: n, m
    real(kind=8) :: xp, yp, dc, fac
    complex(kind=8) :: z1p,z2p,bigz1,bigz2
    """
    wdis = np.zeros(order + 1, dtype=np.complex128)

    bigz1 = complex(d1, 0.0)
    bigz2 = complex(d2, 0.0)
    z1p = 0.5 * (z2 - z1) * bigz1 + 0.5 * (z1 + z2)
    z2p = 0.5 * (z2 - z1) * bigz2 + 0.5 * (z1 + z2)
    wdisc = lapld_int_ho_wdis(x, y, z1p, z2p, order)
    dc = (d1 + d2) / (d2 - d1)
    wdis[0 : order + 1] = 0.0
    for n in range(order + 1):
        for m in range(n + 1):
            wdis[n] = wdis[n] + gam[n, m] * dc ** (n - m) * wdisc[m]
        wdis[n] = (0.5 * (d2 - d1)) ** n * wdis[n]
    return wdis


@numba.njit(nogil=True, cache=True)
def lapld_int_ho_wdis(x, y, z1, z2, order):
    """lapld_int_ho_wdis.

    # Near field only
    implicit none
    integer, intent(in) :: order
    real(kind=8), intent(in) :: x,y
    complex(kind=8), intent(in) :: z1,z2
    complex(kind=8), dimension(0:order) :: wdis
    complex(kind=8), dimension(0:10) :: qm  # Max order is 10
    integer :: m, n
    complex(kind=8) :: z, zplus1, zmin1, term1, term2, zterm
    """
    qm = np.zeros(11, dtype=np.complex128)
    wdis = np.zeros(order + 1, dtype=np.complex128)

    z = (2.0 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    zplus1 = z + 1.0
    zmin1 = z - 1.0
    # Not sure if this gives correct answer at corner point (z also appears in qm);
    # should really be caught in code that calls this function
    if np.abs(zplus1) < tiny:
        zplus1 = tiny
    if np.abs(zmin1) < tiny:
        zmin1 = tiny

    qm[0:1] = 0.0
    for m in range(2, order + 1):
        qm[m] = 0.0
        for n in range(1, m // 2 + 1):
            qm[m] = qm[m] + (m - 2 * n + 1) * z ** (m - 2 * n) / (2 * n - 1)

    term1 = 1.0 / zmin1 - 1.0 / zplus1
    term2 = np.log(zmin1 / zplus1)
    wdis[0] = term1
    zterm = complex(1.0, 0.0)
    for m in range(1, order + 1):
        wdis[m] = m * zterm * term2 + z * zterm * term1 + 2.0 * qm[m]
        zterm = zterm * z

    wdis = -wdis / (np.pi * complex(0.0, 1.0) * (z2 - z1))
    return wdis


# Fp function


@numba.njit(nogil=True, cache=True)
def Fp(x, y, z1, z2, biga, order, d1, d2, a, b, nt):
    tol = 1e-12
    zeta = (2 * complex(x, y) - (z1 + z2)) / (z2 - z1) / biga
    zetabar = np.conj(zeta)
    zminzbar = np.zeros(nt + 1, dtype=np.complex128)
    zminzbar[0] = 1
    for n in range(1, nt + 1):
        zminzbar[n] = zminzbar[n - 1] * (zeta - zetabar)

    eta = np.zeros((nt + 1, nt + 1), dtype=np.complex128)  # lower triangular
    etabar = np.zeros((nt + 1, nt + 1), dtype=np.complex128)
    for n in range(nt + 1):
        for m in range(0, n + 1):
            eta[n, m] = binom[n, m] * zminzbar[n - m]
            etabar[n, m] = np.conj(eta[n, m])

    atil = np.zeros(2 * nt + 1, dtype=np.complex128)
    btil = np.zeros(2 * nt + 1, dtype=np.complex128)
    ctil = np.zeros(2 * nt + 1, dtype=np.complex128)
    for n in range(2 * nt + 1):
        for m in range(max(0, n - nt), int(n / 2) + 1):
            atil[n] = atil[n] + a[n - m] * eta[n - m, m]
            btil[n] = btil[n] + b[n - m] * eta[n - m, m]
            ctil[n] = ctil[n] + a[n - m] * etabar[n - m, m]

    d1minzeta = d1 / biga - zeta
    d2minzeta = d2 / biga - zeta
    if np.abs(d1minzeta) < tol:
        d1minzeta = d1minzeta + complex(tol, 0)
    if np.abs(d2minzeta) < tol:
        d2minzeta = d2minzeta + complex(tol, 0)
    log1 = np.log(d1minzeta)
    log2 = np.log(d2minzeta)

    alpha = np.zeros(2 * nt + order + 1, dtype=np.complex128)
    beta = np.zeros(2 * nt + order + 1, dtype=np.complex128)
    gamma = np.zeros(2 * nt + order + 1, dtype=np.complex128)

    omega = np.zeros(order + 1, dtype=np.complex128)

    for p in range(order + 1):
        alpha[0 : 2 * nt + p + 1] = 0
        beta[0 : 2 * nt + p + 1] = 0
        gamma[0 : 2 * nt + p + 1] = 0

        d = np.zeros(p + 1, dtype=np.complex128)
        dbar = np.zeros(p + 1, dtype=np.complex128)
        for m in range(p + 1):
            d[m] = biga**p * binom[p, m] * zeta ** (p - m)
            dbar[m] = np.conj(d[m])
        for n in range(2 * nt + p + 1):
            for m in range(max(0, n - 2 * nt), min(p, n) + 1):
                alpha[n] = alpha[n] + d[m] * atil[n - m]
                beta[n] = beta[n] + d[m] * btil[n - m]
                gamma[n] = gamma[n] + dbar[m] * ctil[n - m]

        term1 = 1
        term2 = 1
        for n in range(2 * nt + p + 1):
            term1 = term1 * d1minzeta
            term2 = term2 * d2minzeta
            omega[p] = omega[p] + (
                alpha[n] * log2 - alpha[n] / (n + 1) + beta[n]
            ) * term2 / (n + 1)
            omega[p] = omega[p] - (
                alpha[n] * log1 - alpha[n] / (n + 1) + beta[n]
            ) * term1 / (n + 1)
            omega[p] = omega[p] + (
                gamma[n] * np.conj(log2) - gamma[n] / (n + 1)
            ) * np.conj(term2) / (n + 1)
            omega[p] = omega[p] - (
                gamma[n] * np.conj(log1) - gamma[n] / (n + 1)
            ) * np.conj(term1) / (n + 1)

    return biga * omega


@numba.njit(nogil=True, cache=True)
def bessells_int_ho(x, y, z1, z2, lab, order, d1, d2, nt=20):
    """
    Docs.

    To come here
    """
    L = np.abs(z2 - z1)
    ang = np.arctan2(lab.imag, lab.real)
    biga = 2 * np.abs(lab) / L

    exprange = np.exp(-complex(0, 2) * ang * nrange)
    ahat = a * exprange
    bhat = (b - a * complex(0, 2) * ang) * exprange

    omega = Fp(x, y, z1, z2, biga, order, d1, d2, ahat, bhat, nt)
    return -L / (4 * np.pi) * omega


@numba.njit(nogil=True, cache=True)
def bessells_int_ho_qxqy(x, y, z1, z2, lab, order, d1, d2):
    """
    Docs.

    To come here
    """
    nt = 20  # number of terms in series is nt + 1
    bigz = (2 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    bigx = bigz.real
    bigy = bigz.imag
    L = np.abs(z2 - z1)
    ang = np.arctan2(lab.imag, lab.real)
    angz = np.arctan2((z2 - z1).imag, (z2 - z1).real)
    biglab = 2 * lab / L
    biga = np.abs(biglab)

    exprange = np.exp(-complex(0, 2) * ang * nrange)
    ahat = a * exprange
    bhat = (b - a * complex(0, 2) * ang) * exprange

    atil = 2 * nrange[1:] * ahat[1:]
    btil = 2 * nrange[1:] * bhat[1:] + 2 * ahat[1:]

    omega = Fp(x, y, z1, z2, biga, order + 1, d1, d2, atil, btil, nt - 1)
    omegalap = lapld_int_ho_d1d2(x, y, z1, z2, order, d1, d2)
    term1 = 1 / (2 * np.pi * biga**2) * bigx * omega[:-1]
    term2 = -1 / (2 * np.pi * biga**2) * omega[1:]
    term3 = 2 * ahat[0] * omegalap.imag
    qx = term1 + term2 + term3
    term1 = 1 / (2 * np.pi * biga**2) * bigy * omega[:-1]
    term3 = 2 * ahat[0] * omegalap.real
    qy = term1 + term3

    qxqy = np.zeros(2 * order + 2, dtype=np.complex128)
    qxqy[: order + 1] = qx * np.cos(angz) - qy * np.sin(angz)
    qxqy[order + 1 :] = qx * np.sin(angz) + qy * np.cos(angz)
    return qxqy


@numba.njit(nogil=True, cache=True)
def besselld_int_ho(x, y, z1, z2, lab, order, d1, d2):
    """
    Docs.

    To come here
    """
    nt = 20  # number of terms in series is nt + 1
    bigz = (2 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    bigy = bigz.imag
    L = np.abs(z2 - z1)
    ang = np.arctan2(lab.imag, lab.real)
    biglab = 2 * lab / L
    biga = np.abs(biglab)

    exprange = np.exp(-complex(0, 2) * ang * nrange)
    ahat = a1 * exprange
    bhat = (b1 - a1 * complex(0, 2) * ang) * exprange

    omega = Fp(x, y, z1, z2, biga, order, d1, d2, ahat, bhat, nt)

    rv = (
        bigy / (2.0 * np.pi * biglab**2) * omega
        + lapld_int_ho_d1d2(x, y, z1, z2, order, d1, d2).real
    )

    return rv


@numba.njit(nogil=True, cache=True)
def besselld_int_ho_qxqy(x, y, z1, z2, lab, order, d1, d2):
    """
    Docs.

    To come here
    """
    nt = 20  # number of terms in series is nt + 1
    bigz = (2 * complex(x, y) - (z1 + z2)) / (z2 - z1)
    bigx = bigz.real
    bigy = bigz.imag
    L = np.abs(z2 - z1)
    ang = np.arctan2(lab.imag, lab.real)
    angz = np.arctan2((z2 - z1).imag, (z2 - z1).real)
    biglab = 2 * lab / L
    biga = np.abs(biglab)

    exprange = np.exp(-complex(0, 2) * ang * nrange)
    ahat = a1 * exprange
    bhat = (b1 - a1 * complex(0, 2) * ang) * exprange

    atil = 2 * nrange[1:] * ahat[1:]
    btil = 2 * nrange[1:] * bhat[1:] + 2 * ahat[1:]

    omega_pot = Fp(x, y, z1, z2, biga, order, d1, d2, ahat, bhat, nt)
    omega = Fp(x, y, z1, z2, biga, order + 1, d1, d2, atil, btil, nt - 1)
    omegalap = lapld_int_ho_d1d2(x, y, z1, z2, order, d1, d2)
    wlap = lapld_int_ho_wdis_d1d2(x, y, z1, z2, order, d1, d2)

    term1 = bigx / (2 * np.pi * biga**2) * omega[:-1]
    term2 = -1 / (2 * np.pi * biga**2) * omega[1:]
    term3 = 2 * ahat[0] * omegalap.imag
    qx = -2 * bigy / (L * biglab**2) * (term1 + term2 + term3)  # + wlap.real

    term1 = 1 / (2.0 * np.pi * biglab**2) * 2 / L * omega_pot
    term2 = bigy / (2 * np.pi * biga**2) * omega[:-1]
    term3 = 2 * ahat[0] * omegalap.real
    qy = -term1 - 2 * bigy / (L * biglab**2) * (term2 + term3)  # - wlap.imag

    qxqy = np.zeros(2 * order + 2, dtype=np.complex128)
    qxqy[: order + 1] = qx * np.cos(angz) - qy * np.sin(angz) + wlap.real
    qxqy[order + 1 :] = qx * np.sin(angz) + qy * np.cos(angz) - wlap.imag
    return qxqy


@numba.njit(nogil=True, cache=True)
def potbeslsv(x, y, z1, z2, lab, order, ilap, naq, R=8):
    """Potential of line-sink for use in timml."""
    z = x + y * 1j
    pot = np.zeros((order + 1, naq))
    if ilap:
        pot[:, 0] = lapls_int_ho(x, y, z1, z2, order).real
    for n in range(ilap, len(lab)):
        if isinside(z1, z2, z, R * lab[n]):
            d1, d2 = find_d1d2(z1, z2, z, R * lab[n])
            pot[:, n] = bessells(x, y, z1, z2, lab[n], order, d1, d2).real
    return pot


@numba.njit(nogil=True, cache=True)
def disbeslsv(x, y, z1, z2, lab, order, ilap, naq, R=8):
    z = x + y * 1j
    qxqy = np.zeros((2 * (order + 1), naq))
    if ilap:
        wdis = lapls_int_ho_wdis(x, y, z1, z2, order)
        qxqy[: order + 1, 0] = wdis.real
        qxqy[order + 1 :, 0] = -wdis.imag
    for n in range(ilap, len(lab)):
        if isinside(z1, z2, z, R * lab[n]):
            d1, d2 = find_d1d2(z1, z2, z, R * lab[n])
            qxqylab = bessellsqxqy(x, y, z1, z2, lab[n], order, d1, d2).real
            qxqy[: order + 1, n] = qxqylab[0 : order + 1]
            qxqy[order + 1 :, n] = qxqylab[order + 1 :]
    return qxqy

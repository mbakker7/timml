"""Control point utilities.

Helper functions to place control points.
"""

import numpy as np


def controlpoints(Ncp, z1, z2, eps=0, include_ends=False, dely=0):
    # thetacp = np.arange(np.pi, 0, -np.pi/self.Ncp) - 0.5 * np.pi/self.Ncp
    # The following works MUCH better for a uniform head along the line
    delY = dely / np.abs(z2 - z1) * 2
    thetacp = np.linspace(np.pi, 0, Ncp + 2)[1:-1]
    if include_ends:
        Zcp = np.zeros(Ncp + 2, "D")
        Zcp[0] = -1
        Zcp[-1] = 1
        Zcp[1:-1] = np.cos(thetacp)
    else:
        # thetacp = np.arange(np.pi, 0, -np.pi/Ncp) - 0.5 * np.pi/Ncp
        Zcp = np.zeros(Ncp, "D")
        Zcp.real = np.cos(thetacp)
    Zcp.imag = eps + delY  # control point just on positive side (this is handy later on)
    zcp = Zcp * (z2 - z1) / 2.0 + 0.5 * (z1 + z2)
    return zcp.real, zcp.imag


def strengthinf_controlpoints(Ncp, Nlayers):
    # include_ends is False in comparison to function above
    thetacp = np.linspace(np.pi, 0, Ncp + 2)[1:-1]
    Xcp = np.cos(thetacp)
    s = np.zeros((Ncp, Ncp))
    for i in range(Ncp):
        s[i] = Xcp[i] ** np.arange(Ncp)
    rv = np.zeros((Ncp * Nlayers, Ncp * Nlayers))
    for i in range(Ncp):
        for j in range(Nlayers):
            rv[i * Nlayers + j, j::Nlayers] = s[i]
    return rv

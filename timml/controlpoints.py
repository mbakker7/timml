import numpy as np

def controlpoints(Ncp, z1, z2, eps=0, include_ends=False):
    #thetacp = np.arange(np.pi, 0, -np.pi/self.Ncp) - 0.5 * np.pi/self.Ncp
    # The following works MUCH better for a uniform head along the line
    thetacp = np.linspace(np.pi, 0, Ncp+2)[1:-1]
    if include_ends:
        Zcp = np.zeros(Ncp+2, 'D')
        Zcp[0] = -1
        Zcp[-1] = 1
        Zcp[1:-1] = np.cos(thetacp)
    else:
        #thetacp = np.arange(np.pi, 0, -np.pi/Ncp) - 0.5 * np.pi/Ncp
        Zcp = np.zeros(Ncp, 'D')
        Zcp.real = np.cos(thetacp)
    Zcp.imag = eps  # control point just on positive side (this is handy later on)
    zcp = Zcp * (z2 - z1) / 2.0 + 0.5 * (z1 + z2)
    return zcp.real, zcp.imag

def strengthinf(Ncp):
    # include_ends is False in comparison to function above
    thetacp = np.linspace(np.pi, 0, Ncp+2)[1:-1]
    Xcp = np.cos(thetacp)
    rv = np.zeros((Ncp, Ncp))
    for i in range(Ncp):
        rv[i] = Xcp[i] ** np.arange(Ncp)
    return rv
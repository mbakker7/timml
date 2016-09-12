import numpy as np

def timtraceline(ml, xstart, ystart, zstart, step, nstep):
    v0 = ml.velocity(xstart, ystart, zstart)
    xtrace = [xstart]
    ystart = [ystart]
    zstart = [zstart]
    tstart = [0.0]
    #

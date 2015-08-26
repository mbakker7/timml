from numpy import *
from besselaes import disbeslsho, disbesonlylsho, omegalaplsho


def qxqybase(x,y):
    qx = zeros(1)
    qy = zeros(1)
    disbeslsho(x, y,-1.0, 0.0, 1.0, 0.0, 1, zeros(1), 0, qx, qy)
    return qx, qy    

def qxqylabbase(x,y):
    qx = zeros(1)
    qy = zeros(1)
    disbesonlylsho(x, y,-1.0, 0.0, 1.0, 0.0, 2, ones(1), 1, qx, qy)
    return qx, qy

def omega(x,y):
    pot = zeros(1)
    psi = zeros(1)
    omegalaplsho(x, y, -1.0, 0.0, 1.0, 0.0, 0, pot, psi)
    return pot, psi

qxqy = vectorize(qxqybase)
qxqylab = vectorize(qxqylabbase)

def qnbase(x,y,theta):
    qx, qy = qxqy(x,y)
    qn = qx * cos(theta) + qy * sin(theta)
    return qn


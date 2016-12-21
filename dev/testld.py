from timml import *
import numpy as np



def numlap(func, x, y, d):
    numlap = (func(x + d, y) - 2.0 * func(x, y) + func(x - d, y)) / d ** 2 + \
             (func(x, y + d) - 2.0 * func(x, y) + func(x, y - d)) / d ** 2
    return numlap

ml = ModelMaq()
ld = LineDoubletHoBase(ml, x1=-2, y1=0, x2=2, y2=0, delp=1.0, order=0)
ml.solve()
d = 1e-6
x = 2.0
y = 3.0 
num = numlap(ml.potential, x, y, d)
print 'Numerical Laplacian confined:', num
print 'pot above and below, diff:', ml.potential(1, d), ml.potential(1, -d), ml.potential(1, d) - ml.potential(1, -d)


ml = ModelMaq(z=[2, 1, 0], c = 100, top='semi', hstar=0)
ld = LineDoubletHoBase(ml, x1=-2, y1=0, x2=2, y2=0, delp=1.0, order=0)
ml.solve()
d = 1e-3
x = 2.0
y = 3.0 
num = numlap(ml.potential, x, y, d)
print 'Numerical Laplacian confined:', num
print 'Pot / lab ** 2:', ml.potential(x, y) / ml.aq.lab ** 2
print 'pot above and below, diff:', ml.potential(1, d), ml.potential(1, -d), ml.potential(1, d) - ml.potential(1, -d)


ml = ModelMaq(z=[2, 1, 0], c = 100, top='semi', hstar=0)
ld = LineDoubletHoBase(ml, x1=-2, y1=0, x2=2, y2=0, delp=np.array([0, 1.0]), order=1)
ml.solve()
print 'pot above and below, diff:', ml.potential(1, d), ml.potential(1, -d), ml.potential(1, d) - ml.potential(1, -d)

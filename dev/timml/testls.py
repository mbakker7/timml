from besselaesnew import *
from pylab import *

besselaesnew.initialize()


x = 2.0
y = 1.0
z1 = -1 + 0j
z2 = 1 + 0j
lab = np.array([0, 7.0])
naq = 2
order = 1
ilap = 1
potv = zeros((2,2), 'd')
print potv
potv[:,:] = besselaesnew.potbeslsv(x, y, z1, z2, lab, order, ilap)
print potv

pot = zeros(2, 'd')
for i in range(0,2):
    pot[:] = besselaesnew.potbeslsho(x, y, z1, z2, lab, i, ilap)
    print pot
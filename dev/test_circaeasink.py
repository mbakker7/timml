from timml import *
from pylab import *

ml = ModelMaq(kaq=[1, 2], z=[10, 5, 4, 0], c=20)

c = CircAreaSink(ml, 0, 0, 100, 0.01)

ml.solve()

print 'circle is large:', c.islarge

d = 1e-5

x = 20
y = 10

qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)

qx, qy = ml.disvec(x, y)

print 'x, y:', x, y
print 'qx   , qy   :', qx, qy
print 'qxnum, qynum:', qxnum, qynum

x = 200
y = -50

qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)

qx, qy = ml.disvec(x, y)

print 'x, y:', x, y
print 'qx   , qy   :', qx, qy
print 'qxnum, qynum:', qxnum, qynum

################

print 'small resistance'

ml = ModelMaq(kaq=[1, 2], z=[10, 5, 4, 0], c=1e-3)

c = CircAreaSink(ml, 0, 0, 100, 0.01)

ml.solve()

print 'circle is large:', c.islarge

d = 1e-5

x = 99.5 * cos(pi / 3)
y = 99.5 * sin(pi / 3)

qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)

qx, qy = ml.disvec(x, y)

print 'x, y:', x, y
print 'qx   , qy   :', qx, qy
print 'qxnum, qynum:', qxnum, qynum

x = 105.5 * cos(pi / 3)
y = 105.5 * sin(pi / 3)

qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)

qx, qy = ml.disvec(x, y)

print 'x, y:', x, y
print 'qx   , qy   :', qx, qy
print 'qxnum, qynum:', qxnum, qynum

################

print 'multiple resistances'

ml = ModelMaq(kaq=[1, 2, 4, 2], z=[20, 16, 14, 12, 10, 5, 4, 0], c=[1, 100, 4])

c = CircAreaSink(ml, 0, 0, 1000, 0.01)

ml.solve()

print 'circle is large:', c.islarge

d = 1e-5

x = 99.5 * cos(pi / 3)
y = 99.5 * sin(pi / 3)

qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)

qx, qy = ml.disvec(x, y)

print 'x, y:', x, y
print 'qx   , qy   :', qx, qy
print 'qxnum, qynum:', qxnum, qynum

x = 105.5 * cos(pi / 3)
y = 105.5 * sin(pi / 3)

qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)

qx, qy = ml.disvec(x, y)

print 'x, y:', x, y
print 'qx   , qy   :', qx, qy
print 'qxnum, qynum:', qxnum, qynum

x = 995 * cos(pi / 3)
y = -995 * sin(pi / 3)

qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)

qx, qy = ml.disvec(x, y)

print 'x, y:', x, y
print 'qx   , qy   :', qx, qy
print 'qxnum, qynum:', qxnum, qynum

x = 1005 * cos(pi / 3)
y = -1005 * sin(pi / 3)

qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)

qx, qy = ml.disvec(x, y)

print 'x, y:', x, y
print 'qx   , qy   :', qx, qy
print 'qxnum, qynum:', qxnum, qynum

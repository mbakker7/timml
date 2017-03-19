from timml import *

ml = ModelMaq(kaq=[4, 2], z=[10, 5, 4, 0], c=20)

rf = Constant(ml, 0, 100, 50)

uf = Uflow(ml, 0.001, 45)

#ld = ImpLineDoubletString(ml, [(10, -10), (10, 10), (0, 20)], layers=0, order=5)
#ld = LeakyLineDoublet(ml, 10, -10, 10, 10, res=5, layers=[1], order=5)
ld = LeakyLineDoubletString(ml, [(10, -10), (10, 10), (0, 20)], res=5, layers=[0, 1], order=5)
ml.solve()

x, y = 10, 3
qx, qy = ml.disvec(x, y)
d = 1e-6
h1 = ml.head(x - d, y)
h2 = ml.head(x + d, y)
print 'qx:        ', qx
print 'T(h1-h2)/c:', ld.resfac * (h1 - h2)

ml.solve()
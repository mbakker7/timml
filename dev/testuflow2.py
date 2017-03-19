from timml import *
from pylab import *

ml = ModelMaq(kaq=[1, 2], z=[10, 5, 4, 0], c=20)

xy = [(-5,-5), (5,-3), (3,5), (-4,4)]
p1 = PolygonInhomMaq(ml, xy=xy, kaq=[0.2, 8], z=[11,10,5,4,0], c=[10, 50], top='semi', hstar=1.0, order=5, ndeg=3)
#p1 = PolygonInhomMaq(ml, xy=xy, kaq=[0.2, 8], z=[10,5,4,0], c=[20], top='conf', order=3, ndeg=1)

rf = Constant(ml, xr=50, yr=0, hr=1)
uf = Uflow(ml, 0.01, 45)

ml.solve()
xg,yg = linspace(-10,10,50), linspace(-10,10,50)
h = ml.headgrid(xg, yg)

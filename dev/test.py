from timml import *
from pylab import *

ml = ModelMaq(kaq=[1, 2], z=[10, 5, 4, 0], c=20)
#ml = ModelMaq(kaq=[1, 2], z=[12, 10, 5, 4, 0], c=[20,20], top='semi', hstar=2.0)
xy = [(-5,0), (5,0), (5,8), (-5,8)]
p1 = PolygonInhomMaq(ml, xy=[(-5,0), (5,0), (5,8), (-5,8)], kaq=[0.2, 8], z=[11,10,5,4,0], c=[2, 20], top='semi', hstar=1.0, order=3, ndeg=1)
#p1 = PolygonInhomMaq(ml, xy=[(-5,0), (5,0), (5,8), (-5,8)], kaq=[0.2, 8], z=[10,5,4,0], c=[20], top='conf', order=3, ndeg=3)
#p2 = PolygonInhomMaq(ml, xy=[(-5,8),(5,8),(5,14),(-5,14)], kaq=[5, 2], z=[10,5,4,0], c=[20], top='conf', order=3, ndeg=3)
# no shared boundary
#p2 = PolygonInhomMaq(ml, xy=[(-5,10),(5,10),(5,14),(-5,14)], kaq=[5, 2], z=[10,5,4,0], c=[20], top='conf', order=3, ndeg=3)

w = Well(ml, xw=0, yw=-10, Qw=100, layers=1)
#w2 = HeadWell(ml, xw=1, yw=4, hw=-1, rw=0.1, layers=[0])
rf = Constant(ml, xr=0, yr=-100, hr=2)

#from timml.constant import ConstantStar
#c = ConstantStar(ml, hstar=1, aq=p1)

ml.solve()
xg,yg = linspace(-10,10,50), linspace(-5,15,50)

h = ml.headgrid(xg, yg)

timcontour(ml, -10, 10, 50, -10, 10, 50, [0, 1], 20)
#show()
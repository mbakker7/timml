from timml import *

ml = ModelMaq(kaq=[1, 2], z=[10, 5, 4, 0], c=20)
xy = [(-5,0), (5,0), (5,8), (-5,8)]
p1 = PolygonInhomMaq(ml, xy=xy, kaq=[0.2, 8], z=[10,5,4,0], c=[20], top='conf', order=3, ndeg=3)
p2 = PolygonInhomMaq(ml, xy=[(-5,8),(5,8),(5,14),(-5,14)], kaq=[5, 2], z=[10,5,4,0], c=[20], top='conf', order=3, ndeg=3)
# no shared boundary
#p2 = PolygonInhomMaq(ml, xy=[(-5,10),(5,10),(5,14),(-5,14)], kaq=[5, 2], z=[10,5,4,0], c=[20], top='conf', order=3, ndeg=3)

w = WellBase(ml, xw=0, yw=-10, Qw=100, layers=1)
rf = Constant(ml, xr=0, yr=-100, hr=2)
ml.solve()

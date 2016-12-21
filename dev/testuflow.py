from timml import *
from pylab import *

ml = ModelMaq(kaq=[1, 2], z=[10, 5, 4, 0], c=20)

w = HeadWell(ml, xw=10, yw=0, hw=-1, rw=0.1, layers=[0])
rf = Constant(ml, xr=0, yr=0, hr=1)
uf = Uflow(ml, 0.01, 45)

ml.solve()
xg,yg = linspace(-50,50,50), linspace(-50,50,50)
h = ml.headgrid(xg, yg)

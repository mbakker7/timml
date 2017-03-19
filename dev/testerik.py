from timml import *
from pylab import *

xyeast = np.loadtxt('rivereast.txt')
xywest = np.loadtxt('riverwest.txt')
xyriver = np.vstack((xywest, xyeast[-1::-1]))
xw = 734197.1
yw = 4072815.9
x1,x2,y1,y2 = (733800.0, 735000.0, 4071800.0, 4073600.0)
xg = linspace(x1, x2, 50)
yg = linspace(y1, y2, 50)
xs1,xs2,ys1,ys2 = (734034.375, 734423.4375, 4072615.625, 4073000.0)
xg2 = linspace(xs1, xs2, 50)
yg2 = linspace(ys1, ys2, 50)

c = 34

#print 'using smaller resistance of the river'
#c = 1.0

ml = ModelMaq(kaq=24, z=[0, -15])
w = WellBase(ml, xw=xw, yw=yw, Qw=3145, layers=0)
rf = Constant(ml, xr=xw+1e6, yr=yw, hr=0)
p = PolygonInhomMaq(ml, xy=xyriver, kaq=24, z=[1, 0, -15], c=c, top='semi', order=2, ndeg=3, hstar=0)
ml.solve()

#ml2 = ModelMaq(kaq=24, z=[0, -15])
#w = WellBase(ml2, xw=xw, yw=yw, Qw=3145, layers=0)
#rf = Constant(ml2, xr=xw+1e6, yr=yw, hr=0)
#for i in range(len(xyeast) - 1):
#    HeadLineSink(ml2, xyeast[i,0], xyeast[i,1], xyeast[i+1,0], xyeast[i+1,1], 0, res=c, wh=20)
#for i in range(len(xyeast) - 1):
#    HeadLineSink(ml2, xywest[i,0], xywest[i,1], xywest[i+1,0], xywest[i+1,1], 0, res=c, wh=20)
#ml2.solve()
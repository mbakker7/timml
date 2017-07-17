timmlpath = '/Users/mark/git/timml'
import sys
if timmlpath not in sys.path: sys.path.append(timmlpath)
from timml import *

ml = Model(k=[10], zb=[0.0], zt=[100.0])

xylist0=[(-30, -30), (30, -30), (30, 30), (-30, 30)]
PolygonInhom(ml, k=[10], c=[], zb=[0.0], zt=[100.0], xylist=xylist0)

xylist1=[(-70, -70), (70, -70), (70, 70), (-70, 70)]
PolygonInhom(ml, k=[1], c=[], zb=[0.0], zt=[100.0], xylist=xylist1)

MakeInhomPolySide(ml, xylist=xylist0, order=5, closed=True)
MakeInhomPolySide(ml, xylist=xylist1, order=5, closed=True)

Constant(ml, 500, 500, 100.0, [0], label="constant")
Well(ml, 0, 0, 1000, 0.5, 0, label="well")

ml.solve()
timcontour(ml, -100, 100, 100, -100, 100, 100, 1, 50)
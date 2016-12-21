from timml import *
import numpy as np

ml = ModelMaq(kaq=[2,6,4], z= [165, 140, 120, 80, 60, 0], c = [2000, 20000], npor=0.3)
rf = Constant(ml, xr=20000, yr=20000, hr=175, layer=0)
p = CircAreaSink(ml, xc=10000, yc=10000, R=15000, N=0.0002)
w1 = Well(ml, xw=10000, yw=8000, Qw=1000, rw=0.3, layers=[1], label='well 1')
w2 = Well(ml, xw=12100, yw=10700, Qw=5000, rw=0.3, layers=2, label='well 2')
w3 = Well(ml, xw=10000, yw=4600, Qw=5000, rw=0.3, layers=[1,2], label='maq well')
#
xyls = [(9510, 19466), (12620, 17376), (12753, 14976), (13020, 12176), (15066, 9466), (16443, 7910), (17510, 5286), (17600, 976)]
hls = [170, 168, 166, 164, 162, 160, 158]
ls1 = HeadLineSinkString(ml, xyls, hls=hls, layers=0)
#
xyls2 = [(356, 6976), (4043, 7153), (6176, 8400), (9286, 9820), (12266, 9686), (15066, 9466)]
hls2 = [174, 171, 168, 166, 164]
ls2 = HeadLineSinkString(ml, xyls2, hls=hls2, layers=0)
#
xyls3 = [(1376,  1910), (4176, 2043), (6800, 1553), (9953, 2086), (14043, 2043), (17600, 976)]
hls3 = [170, 166, 162, 160, 158]
ls3 = HeadLineSinkString(ml, xyls3, hls=hls3, layers=0)

##
#
#def trace_example1():
#    an = arange(0,2*pi,pi/5.0)
#    timtracelines(ml,w1.xw+cos(an),w1.yw+sin(an),100*ones(len(an)),-100,Nmax=500)
#    timtracelines(ml,w2.xw+cos(an),w2.yw+sin(an),30*ones(len(an)),-100,Nmax=500)
#    timtracelines(ml,w3.xw+cos(an),w3.yw+sin(an),100*ones(len(an)),-100,tmax=200*365,Nmax=200)
#
#try:
#    # equivalent to %matplotlib in IPython
#    get_ipython().magic('matplotlib')
#except:
#    pass
## Solve model
ml.solve()
#xg, yg = np.linspace(0, 20000, 50), np.linspace(0, 20000, 50)
#h = ml.headgrid(xg, yg)
## Contour results
#timcontour(ml, 0, 20000, 50, 0, 20000, 50, 3, 20)
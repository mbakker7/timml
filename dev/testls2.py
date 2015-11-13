from timml import *
import numpy as np

ml = ModelMaq(kaq=[2,6,4], z= [165, 140, 120, 80, 60, 0], c = [2000, 20000], npor=0.3)

rf = Constant(ml, xr=20000, yr=20000, hr=20, layer=0)

ls1 = HeadLineSinkHo(ml, x1=-500, y1=-100, x2=700, y2=200, hls = 10, order=5, layers=0)

ml.initialize()

xg, yg = np.linspace(-1000, 1000, 50), np.linspace(-500, 500, 100)

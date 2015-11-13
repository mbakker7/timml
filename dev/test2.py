from timml import *

ml = ModelMaq(kaq=[4, 2], z=[10, 5, 4, 0], c=20)

rf = Constant(ml, 0, 100, 50)

ls = HeadLineSink(ml, -10, 0, 10, 0, 40, res=2, wh = 10)
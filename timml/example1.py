from timml import *
ml = Model(k = [2,6,4], zb = [140,80,0], zt = [165,120,60], c = [2000,20000], n = [0.3,0.25,0.3], nll = [0.2,0.25])
rf = Constant(ml,20000,20000,175,[0])
p=CircAreaSink(ml,10000,10000,15000,0.0002,[0])
w1=Well(ml,10000,8000,3000,.3,[1],'well 1')
w2=Well(ml,12000,8000,5000,.3,[2],'well 2')
w3=Well(ml,10000,4600,5000,.3,[1,2],'maq well')
#
HeadLineSink(ml, 9510,  19466, 12620, 17376, 170,[0])
HeadLineSink(ml, 12620, 17376, 12753, 14976, 168,[0])
HeadLineSink(ml, 12753, 14976, 13020, 12176, 166,[0])
HeadLineSink(ml, 13020, 12176, 15066, 9466,  164,[0])
HeadLineSink(ml, 15066, 9466,  16443, 7910,  162,[0])
HeadLineSink(ml, 16443, 7910,  17510, 5286,  160,[0])
HeadLineSink(ml, 17510, 5286,  17600, 976,   158,[0])
#
HeadLineSink(ml, 356,   6976,  4043,  7153, 174,[0])
HeadLineSink(ml, 4043,  7153,  6176,  8400, 171,[0])
HeadLineSink(ml, 6176,  8400,  9286,  9820, 168,[0])
HeadLineSink(ml, 9286,  9820,  12266, 9686, 166,[0])
HeadLineSink(ml, 12266, 9686,  15066, 9466, 164,[0])
#
HeadLineSink(ml, 1376,  1910,  4176,  2043, 170,[0])
HeadLineSink(ml, 4176,  2043,  6800,  1553, 166,[0])
HeadLineSink(ml, 6800,  1553,  9953,  2086, 162,[0])
HeadLineSink(ml, 9953,  2086,  14043, 2043, 160,[0])
HeadLineSink(ml, 14043, 2043,  17600, 976 , 158,[0])
#

def trace_example1():
    an = arange(0,2*pi,pi/5.0)
    timtracelines(ml,w1.xw+cos(an),w1.yw+sin(an),100*ones(len(an)),-100,Nmax=500)
    timtracelines(ml,w2.xw+cos(an),w2.yw+sin(an),30*ones(len(an)),-100,Nmax=500)
    timtracelines(ml,w3.xw+cos(an),w3.yw+sin(an),100*ones(len(an)),-100,tmax=200*365,Nmax=200)


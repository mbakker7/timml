
# coding: utf-8

# # TimML Exercises

# ## Exercise 2: A system with wells, rivers, and recharge
# 
# Consider a system of three aquifers. The aquifer parameters are presented in Table 1; note that an average thickness is specified for the top unconfined aquifer. A river with three branches cuts through the upper aquifer. The river is modeled with 7 head-specified line-sinks and each branch is modeled with 5 line-sinks. The Northern branch is modeled with resistance line-sinks with a resistance $c$ of 5 days and a width $w$ of 10 meters. The heads are specified at the centers of the line-sinks and are shown in Figure 1. 
# 
# Three wells are present. Well 1 is screened in aquifer 0 and has a discharge of 1000 m$^3$/d, well 2 is screened in aquifer 2 and has a discharge of 5000 m$^3$/d, and well 3 is screened in aquifers 1 and 2 and has a total discharge of 5000 m$^3$/d. A constant recharge through the upper boundary of aquifer 0 is simulated by one large circular infiltration area that covers the entire model area; the recharge rate is 0.2 mm/day. A head of 175 m is specified in layer 0 at the upper righthand corner of the model domain. A layout of all analytic elements, except the boundary of the infiltration area, is shown in Figure 1. 
# 
# #### Table 1: Aquifer data for Exercise 2
# |             | $k$ (m/d) | $z_b$ (m) | $z_t$ | $c$ (days) | $n$ (-) | $n_{ll}$ (-) |
# |------------:|----------:|----------:|------:|-----------:|--------:|----------:|
# |Aquifer 0    |   2       |   140     | 165   |            |  0.3    |           | 
# |Leaky Layer 1|           |   120     | 140   |    2000    |         |   0.2     |    
# |Aquifer 1    |   6       |   80      | 120   |            |  0.25   |           |  
# |Leaky Layer 2|           |   60      | 80    |    20000   |         |   0.25    |  
# |Aquifer 2    |   4       |   0       | 60    |            |  0.3    |           ||
# 
# <img src="Layout_exercise2.png"> </img>
# 
# #### Figure 1: Layout of elements for Exercise 2. Heads at centers of line-sinks are indicated. 

# In[1]:

from timml import *
from pylab import *


# In[2]:

# Create basic model elements
ml = Model(k=[2, 6, 4],
           zb=[140, 80, 0],
           zt=[165, 120, 60],
           c=[2000, 20000],
           n=[0.3, 0.3, 0.3],
           nll=[0.3, 0.3])
rf = Constant(ml, xr=20000, yr=20000, head=175, layer=0)
p = CircAreaSink(ml, xp=10000, yp=10000, Rp=15000, infil=0.0002, layer=0)
w1 = Well(ml, xw=10000, yw=8000, Qw=1000, rw=0.3, layers=1, label='well 1')
w2 = Well(ml, xw=12100, yw=10700, Qw=5000, rw=0.3, layers=2, label='well 2')
w3 = Well(ml, xw=10000, yw=4600, Qw=5000, rw=0.3, layers=[1,2], label='maq well')
#
HeadLineSink(ml, x1=9510, y1=19466, x2=12620, y2=17376, head=170, layers=0)
HeadLineSink(ml, 12620, 17376, 12753, 14976, 168, [0])
HeadLineSink(ml, 12753, 14976, 13020, 12176, 166, [0])
HeadLineSink(ml, 13020, 12176, 15066, 9466,  164, [0])
HeadLineSink(ml, 15066, 9466,  16443, 7910,  162, [0])
HeadLineSink(ml, 16443, 7910,  17510, 5286,  160, [0])
HeadLineSink(ml, 17510, 5286,  17600, 976,   158, [0])
#
HeadLineSink(ml, 356,   6976,  4043,  7153, 174, [0])
HeadLineSink(ml, 4043,  7153,  6176,  8400, 171, [0])
HeadLineSink(ml, 6176,  8400,  9286,  9820, 168, [0])
HeadLineSink(ml, 9286,  9820,  12266, 9686, 166, [0])
HeadLineSink(ml, 12266, 9686,  15066, 9466, 164, [0])
#
HeadLineSink(ml, 1376,  1910,  4176,  2043, 170, [0])
HeadLineSink(ml, 4176,  2043,  6800,  1553, 166, [0])
HeadLineSink(ml, 6800,  1553,  9953,  2086, 162, [0])
HeadLineSink(ml, 9953,  2086,  14043, 2043, 160, [0])
HeadLineSink(ml, 14043, 2043,  17600, 976 , 158, [0])

ml.solve()
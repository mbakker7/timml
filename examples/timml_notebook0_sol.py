# ruff: noqa
# coding: utf-8

# # TimML Notebook 0
# ## Single layer flow


# Consider uniform flow from East to West. The gradient is 0.001.
# The hydraulic conductivity is $k=10$ m/d. The aquifer bottom and top are at 0 m
# and 10 m. The head at $x=-1000$ m and $y=0$ is fixed at 41 m.

# In[1]:

import numpy as np

import timml


# In[2]:


ml = timml.ModelMaq(kaq=10, z=[10, 0])


# In[3]:


rf = timml.Constant(ml, xr=-1000, yr=0, hr=41)


# In[4]:


uf = timml.Uflow(ml, slope=0.001, angle=0)


# In[5]:


ml.solve()


# In[6]:


ml.plots.contour(
    x1=-1200,
    x2=200,
    nx=50,
    y1=-500,
    y2=500,
    ny=50,
    levels=10,
    labels=True,
    decimals=2,
    legend=True,
)


# The default contour levels are not what we want for this example, so let's
# specify the levels to go from 39 to 42 with steps of 0.1 (not all those levels
# are present in the current window).

# In[7]:


ml.plots.contour(
    x1=-1200,
    x2=200,
    nx=50,
    y1=-500,
    y2=500,
    ny=50,
    levels=np.arange(39, 42, 0.1),
    labels=True,
    decimals=1,
)


# A well is located at $(x,y)=(-400,0)$ with a discharge $Q=50$ m$^3$ and a radius of 0.2 m.

# In[8]:


w = timml.Well(ml, xw=-400, yw=0, Qw=50.0, rw=0.2)

# After the well is added (or any other elements), the model needs to be solved again. A contour plot is created and a 10 strace line are added. The stepsize is given in meters and represents the largest space step that is taken, but it is reduced when certain accuracy constraints are not met. Note that, after running the code cell below, for each trace line it is printed to the screen what the reason was that the traceline was aborted. In this case it was either because the trace line reached a well or reached the maximum number of steps (the default is 100 steps, but this can be changed by specifying the `nstepmax` keyword).

# In[9]:


ml.solve()
ml.plots.contour(-1000, 100, 50, -500, 500, 50, levels=np.arange(39, 42, 0.1))
ml.plots.tracelines(-800 * np.ones(1), -200 * np.ones(1), np.zeros(1), hstepmax=20)


# ### Exercise a
# Draw 10 tracelines from $x=-800$ and different values of $y$.

# In[10]:


ml.plots.contour(-1000, 100, 50, -500, 500, 50, levels=np.arange(39, 42, 0.1))
ml.plots.tracelines(
    -800 * np.ones(10), np.linspace(-500, 500, 10), np.zeros(10), hstepmax=20
)


# ### Exercise b
# Quadruple the discharge of the well and reproduce the same figure

# In[11]:


ml = timml.ModelMaq(kaq=10, z=[10, 0])
rf = timml.Constant(ml, xr=-1000, yr=0, hr=41)
uf = timml.Uflow(ml, slope=0.001, angle=0)
w = timml.Well(ml, xw=-400, yw=0, Qw=200, rw=0.2)
ml.solve()
ml.plots.contour(-1000, 100, 50, -500, 500, 50, levels=np.arange(39, 42, 0.1))
ml.plots.tracelines(
    -800 * np.ones(10), np.linspace(-500, 500, 10), np.zeros(10), hstepmax=20
)

print(("head at well:", w.headinside()))


# ### Add a river
# A river runs along $x=0$. The water level in the river is at 40 m.

# In[12]:


ml = timml.ModelMaq(kaq=10, z=[10, 0])
rf = timml.Constant(ml, xr=-1000, yr=0, hr=41)
uf = timml.Uflow(ml, slope=0.001, angle=0)
w = timml.Well(ml, xw=-400, yw=0, Qw=200, rw=0.2)
ls1 = timml.HeadLineSink(ml, 0, -500, 0, 500, 40)
ml.solve()
ml.plots.contour(-1000, 100, 50, -500, 500, 50, levels=np.arange(39, 42, 0.1))
print(("head at well:", w.headinside()))


# ### Exercise c
# Simulate the river with 20 line-sinks from $y=-800$ to $y=800$.

# In[13]:


ml = timml.ModelMaq(kaq=10, z=[10, 0])
rf = timml.Constant(ml, xr=-1000, yr=0, hr=41)
uf = timml.Uflow(ml, slope=0.001, angle=0)
w = timml.Well(ml, xw=-400, yw=0, Qw=200, rw=0.2)
xls = np.zeros(21)
yls = np.linspace(-800, 800, 21)
ls = timml.HeadLineSinkString(ml, xy=list(zip(xls, yls)), hls=40, layers=0)
ml.solve()
ml.plots.contour(-1000, 100, 50, -500, 500, 50, levels=np.arange(39, 42, 0.1))
ml.plots.tracelines(
    -800 * np.ones(10), np.linspace(-500, 500, 10), np.zeros(10), hstepmax=20, color="b"
)
ml.plots.tracelines(
    -0.01 * np.ones(5), np.linspace(-150, 150, 5), np.zeros(5), hstepmax=20, color="r"
)


# ### Capture zone
# Create a five year capture zone. You may want to create a contour plot first.

# In[14]:


ml.plots.contour(-1000, 100, 50, -500, 500, 50, levels=np.arange(39, 42, 0.1), layers=0)
w.plotcapzone(hstepmax=20, nt=20, zstart=0, tmax=5 * 365.25, color="b")


# ### Exercise d
# Create a 20 year capture zone using 20 tracelines.

# In[15]:


ml.plots.contour(-1000, 100, 50, -500, 500, 50, levels=np.arange(39, 42, 0.1), layers=0)
w.plotcapzone(hstepmax=20, nt=20, zstart=0, tmax=20 * 365.25, color="b")


# In[ ]:

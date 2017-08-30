import numpy as np
import matplotlib.pyplot as plt
from .trace import timtracelines

class PlotTim:
    def plot(self):
        for e in self.elementlist:
            e.plot()
    
    def contour(self, x1, x2, nx, y1, y2, ny, layers, levels, layout=True, labels=False, decimals=0, color=None, newfig=True):
        layers = np.atleast_1d(layers)
        xg = np.linspace(x1, x2, nx)
        yg = np.linspace(y1, y2, ny)
        h = self.headgrid2(x1, x2, nx, y1, y2, ny, layers)
        if newfig: plt.figure()
        for ilayer in range(len(layers)):
            cs = plt.contour(xg, yg, h[ilayer], levels, colors=color)
            if labels:
                fmt = '%1.' + str(decimals) + 'f'
                plt.clabel(cs, fmt=fmt)
        plt.axis('scaled')
        if layout: self.plot()
        #plt.show()
        
    def vcontour(self, x1, x2, y1, y2, n, levels, labels=False, decimals=0, color=None, newfig=True):
        h = self.headalongline(np.linspace(x1, x2, n), np.linspace(y1, y2, n))
        L = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        xg = np.linspace(0, L, n)
        zg = 0.5 * (self.aq.zaqbot + self.aq.zaqtop)
        zg = np.hstack((self.aq.zaqtop[0], zg, self.aq.zaqbot[-1]))
        h = np.vstack((h[0], h, h[-1]))
        if newfig: plt.figure()
        cs = plt.contour(xg, zg, h, levels, colors=color)
        if labels:
            fmt = '%1.' + str(decimals) + 'f'
            plt.clabel(cs, fmt=fmt)
        
    def tracelines(self, xstart, ystart, zstart, hstepmax, vstepfrac=0.2, tmax=1e12, nstepmax=100, silent='.', color=None, horizontal=True, win=[-1e30, 1e30, -1e30, 1e30]):
        xyztlist = timtracelines(self, xstart, ystart, zstart, hstepmax=hstepmax, vstepfrac=vstepfrac, tmax=tmax, nstepmax=nstepmax, silent=silent, win=win)
        if horizontal:
            for xyzt in xyztlist:
                plt.plot(xyzt[:, 0], xyzt[:, 1], color=color)
        else:
            for xyzt in xyztlist:
                plt.plot(xyzt[:, 0], xyzt[:, 2], color=color)
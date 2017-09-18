import numpy as np
import matplotlib.pyplot as plt
from .trace import timtracelines

class PlotTim:
    def plot(self, win=None, newfig=True, figsize=None, orientation='hor'):
        if newfig:
            plt.figure(figsize=figsize)
            ax1 = None
            ax2 = None
            if orientation == 'both':
                ax1 = plt.subplot(211)
                ax2 = plt.subplot(212, sharex=ax1)
            elif orientation[:3] == 'hor':
                ax1 = plt.subplot()
            elif orientation[:3] == 'ver':
                ax2 = plt.subplot()
        else:
            if orientation == 'both' or orientation == 'hor':
                fig = plt.gcf()
                ax1 = fig.axes[0]
        if ax1 is not None:
            plt.axes(ax1)
            for e in self.elementlist:
                e.plot()
            if orientation == 'hor':
                plt.axis('scaled')
            elif orientation == 'both':
                plt.axis('equal')  # cannot be 'scaled' when sharing axes
            if win is not None:
                plt.axis(win)
            
    def contour(self, x1, x2, nx, y1, y2, ny, layers, levels, layout=True, labels=False, decimals=0, color=None, newfig=True, figsize=None):
        layers = np.atleast_1d(layers)
        xg = np.linspace(x1, x2, nx)
        yg = np.linspace(y1, y2, ny)
        h = self.headgrid2(x1, x2, nx, y1, y2, ny, layers)
        if newfig:
            plt.figure(figsize=figsize)
        # color
        if type(color) is list:
            if len(color) < len(layers):
                print('len(color) smaller than len(layers)')
        if type(color) is str:
            color = len(layers) * [color]
        if color is None:
            color = len(layers) * [color]
        # contour
        for i in range(len(layers)):
            cs = plt.contour(xg, yg, h[i], levels, colors=color[i])
            if labels:
                fmt = '%1.' + str(decimals) + 'f'
                plt.clabel(cs, fmt=fmt)
        plt.axis('scaled')
        if layout:
            self.plot(win=[x1, x2, y1, y2], newfig=False)
        #plt.show()
        
    def vcontour(self, x1, x2, y1, y2, n, levels, labels=False, decimals=0, color=None, newfig=True, figsize=None):
        h = self.headalongline(np.linspace(x1, x2, n), np.linspace(y1, y2, n))
        L = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        xg = np.linspace(0, L, n)
        zg = 0.5 * (self.aq.zaqbot + self.aq.zaqtop)
        zg = np.hstack((self.aq.zaqtop[0], zg, self.aq.zaqbot[-1]))
        h = np.vstack((h[0], h, h[-1]))
        if newfig:
            plt.figure(figsize=figsize)
        cs = plt.contour(xg, zg, h, levels, colors=color)
        if labels:
            fmt = '%1.' + str(decimals) + 'f'
            plt.clabel(cs, fmt=fmt)
        
    def tracelines(self, xstart, ystart, zstart, hstepmax, vstepfrac=0.2,
                   tmax=1e12, nstepmax=100, silent='.', color=None, orientation='hor',
                   win=[-1e30, 1e30, -1e30, 1e30], newfig=False, figsize=None):
        fig = plt.gcf()
        assert len(fig.axes) > 0, 'Error: Need to specify axes in figure before invoking tracelines'
        ax1 = None
        ax2 = None
        if orientation == 'both':
            ax1 = fig.axes[0]
            ax2 = fig.axes[1]
        elif orientation[:3] == 'hor':
            ax1 = fig.axes[0]
        elif orientation[:3] == 'ver':
            ax2 = fig.axes[1]
        xyztlist = timtracelines(self, xstart, ystart, zstart, hstepmax=hstepmax, vstepfrac=vstepfrac, tmax=tmax, nstepmax=nstepmax, silent=silent, win=win)
        if ax1 is not None:
            plt.axes(ax1)
            for xyzt in xyztlist:
                ax1.plot(xyzt[:, 0], xyzt[:, 1], color=color)
        if ax2 is not None:
            plt.axes(ax1)
            for xyzt in xyztlist:
                ax2.plot(xyzt[:, 0], xyzt[:, 2], color=color)
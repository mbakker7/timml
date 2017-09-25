import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from .trace import timtracelines, timtraceline

class PlotTim:
    def plot(self, win=None, newfig=True, figsize=None, orientation='hor', topfigfrac=0.8):
        if newfig:
            plt.figure(figsize=figsize)
            ax1 = None
            ax2 = None
            if orientation == 'both':
                ax1 = plt.axes([0.125, 0.18 + (1 - topfigfrac) * 0.7, (0.9 - 0.125), topfigfrac * 0.7])
                ax2 = plt.axes([0.125, 0.11, (0.9 - 0.125), (1 - topfigfrac) * 0.7], sharex=ax1)
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
            
    def contour(self, x1, x2, nx, y1, y2, ny, layers=0, levels=20, layout=True, labels=False,
                decimals=0, color=None, newfig=True, figsize=None, legend=True):
        layers = np.atleast_1d(layers)
        xg = np.linspace(x1, x2, nx)
        yg = np.linspace(y1, y2, ny)
        h = self.headgrid2(x1, x2, nx, y1, y2, ny, layers)
        if newfig:
            plt.figure(figsize=figsize)
        # color
        if color is None:
            c = plt.rcParams['axes.prop_cycle'].by_key()['color']
        elif type(color) is str:
            c = len(layers) * [color]
        elif type(color) is list:
            c = color
        if len(c) < len(layers):
            n = np.ceil(self.aq.naq / len(c))
            c = n * c 
        # contour
        cscollectionlist = []
        for i in range(len(layers)):
            cs = plt.contour(xg, yg, h[i], levels, colors=c[i])
            cscollectionlist.append(cs.collections[0])
            if labels:
                fmt = '%1.' + str(decimals) + 'f'
                plt.clabel(cs, fmt=fmt)
        if legend == True:
            legendlist = ['layer ' + str(i) for i in layers]
            plt.legend(cscollectionlist, legendlist)
        elif type(legend) is list:
            plt.legend(cscollectionlist, legend)
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
        if color is None:
            c = plt.rcParams['axes.prop_cycle'].by_key()['color']
        elif type(color) is str:
            c = self.aq.naq * [color]
        elif type(color) is list:
            c = color
        if len(c) < self.aq.naq:
            n = np.ceil(self.aq.naq / len(c))
            c = n * c 
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
            xyztlist = []
        for i in range(len(xstart)):
            xyzt, layerlist = timtraceline(self, xstart[i], ystart[i], zstart[i], hstepmax=hstepmax,
                                           vstepfrac=vstepfrac, tmax=tmax, nstepmax=nstepmax,
                                           silent=silent, win=win, returnlayers=True)
            if silent == '.':
                print('.', end='', flush=True)
            if ax1 is not None:
                #plt.axes(ax1)
                color = [c[self.aq.layernumber[i]] if self.aq.ltype[i] == 'a' else 'k' for i in layerlist]
                points = np.array([xyzt[:,0], xyzt[:,1]]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = LineCollection(segments, colors=color)
                ax1.add_collection(lc)
                #ax1.plot(xyzt[:, 0], xyzt[:, 1], color=color)
            if ax2 is not None:
                color = [c[self.aq.layernumber[i]] if self.aq.ltype[i] == 'a' else 'k' for i in layerlist]
                points = np.array([xyzt[:,0], xyzt[:,2]]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = LineCollection(segments, colors=color)
                ax2.add_collection(lc)
                ax2.set_ylim(self.aq.z[-1], self.aq.z[0])
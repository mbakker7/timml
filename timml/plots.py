from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from timml.aquifer import SimpleAquifer
from timml.trace import timtraceline

plt.rcParams["contour.negative_linestyle"] = "solid"


class PlotTim:
    def __init__(self, ml):
        self._ml = ml

    def __repr__(self):
        """Return string representation of Plots submodule."""
        methods = "".join(
            [f"\n - {meth}" for meth in dir(self) if not meth.startswith("_")]
        )
        return "TimML plots, available methods:" + methods

    def topview(
        self,
        win=None,
        newfig=True,
        figsize=None,
        orientation="hor",
        topfigfrac=0.8,
        layer=None,
    ):
        """Plot the model layout.

        Other features such as pathlines or capture zones may be added to the plot with
        separate commands.

        Parameters
        ----------
        win : list or tuple
            [xmin, xmax, ymin, ymax]
        newfig : boolean (default True)
            create new figure
        figsize : tuple of 2 values (default is mpl default)
            size of figure
        orientation : ('hor', 'ver', 'both')
            'hor' for horizontal, 'ver' for vertical
            'both' for horizontal above vertical
        topfigfrac : float
            relative size of top figure when orientation='both'
        layer : integer
            layer for which plot is created

        Returns
        -------
        None
        """
        if newfig:
            plt.figure(figsize=figsize)
            ax1 = None
            ax2 = None
            if orientation == "both":
                ax1 = plt.axes(
                    [
                        0.125,
                        0.18 + (1 - topfigfrac) * 0.7,
                        (0.9 - 0.125),
                        topfigfrac * 0.7,
                    ]
                )
                ax2 = plt.axes(
                    [0.125, 0.11, (0.9 - 0.125), (1 - topfigfrac) * 0.7], sharex=ax1
                )
            elif orientation[:3] == "hor":
                ax1 = plt.subplot()
            elif orientation[:3] == "ver":
                ax2 = plt.subplot()
        else:
            if orientation == "both":
                fig = plt.gcf()
                ax1 = fig.axes[0]
                ax2 = fig.axes[1]
            elif orientation[:3] == "hor":
                fig = plt.gcf()
                ax1 = fig.axes[0]
                ax2 = None
            elif orientation[:3] == "ver":
                fig = plt.gcf()
                ax1 = None
                ax2 = fig.axes[0]
        if ax1 is not None:
            plt.sca(ax1)
            for e in self._ml.elementlist:
                e.plot(layer=layer)
            if orientation[:3] == "hor":
                plt.axis("scaled")
            elif orientation == "both":
                plt.axis("equal")  # cannot be 'scaled' when sharing axes
            if win is not None:
                plt.axis(win)
        if ax2 is not None:
            plt.sca(ax2)
            for i in range(self._ml.aq.nlayers):
                if self._ml.aq.ltype[i] == "l":
                    plt.axhspan(
                        ymin=self._ml.aq.z[i + 1],
                        ymax=self._ml.aq.z[i],
                        color=[0.8, 0.8, 0.8],
                    )
            for i in range(1, self._ml.aq.nlayers):
                if self._ml.aq.ltype[i] == "a" and self._ml.aq.ltype[i - 1] == "a":
                    plt.axhspan(
                        ymin=self._ml.aq.z[i],
                        ymax=self._ml.aq.z[i],
                        color=[0.8, 0.8, 0.8],
                    )

    def xsection(
        self,
        xy: Optional[list[tuple[float]]] = None,
        labels=True,
        params=False,
        ax=None,
    ):
        """Plot cross-section of model.

        Note: this method does not plot elements at this time. It does plot
        cross-section inhoms if the model is a cross-section model (ModelXsection).

        Parameters
        ----------
        xy : list of tuples, optional
            list of tuples with coordinates of the form [(x0, y0), (x1, y1)]. If not
            provided, a cross section with length 1 is plotted.
        labels : bool, optional
            add layer numbering labels to plot
        params : bool, optional
            add parameter values to plot
        ax : matplotlib.Axes, optional
            axes to plot on, default is None which creates a new figure

        Returns
        -------
        ax : matplotlib.Axes
            axes with plot
        """
        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=(8, 4))

        # if SimpleAquifer, plot inhoms and return
        if isinstance(self._ml.aq, SimpleAquifer):
            for e in self._ml.elementlist:
                e.plot(ax=ax)
            if xy is not None:
                (x1, _), (x2, _) = xy
            else:
                x1, x2 = ax.get_xlim()
            for inhom in self._ml.aq.inhomlist:
                inhom.plot(ax=ax, labels=labels, params=params, x1=x1, x2=x2)
            ax.set_xlim(x1, x2)
            ax.set_ylabel("elevation")
            ax.set_xlabel("x")
            return ax

        # else get cross-section line
        if xy is not None:
            (x0, y0), (x1, y1) = xy
            r = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
            if y0 == 0 and y1 == 0:
                r0 = x0
                ax.set_xlim(x0, x1)
            elif x0 == 0 and x1 == 0:
                ax.set_ylim(y0, y1)
                r0 = y0
            else:
                ax.set_xlim(0, r)
                r0 = 0.0
        else:
            r0 = 0.0
            r = 1.0

        # get values for layer and aquifer numbering
        if labels:
            lli = 1 if self._ml.aq.ltype[0] == "a" else 0
            aqi = 0
        else:
            lli = None
            aqi = None

        # plot layers
        for i in range(self._ml.aq.nlayers):
            # leaky layers
            if self._ml.aq.ltype[i] == "l":
                ax.axhspan(
                    ymin=self._ml.aq.z[i + 1],
                    ymax=self._ml.aq.z[i],
                    color=[0.8, 0.8, 0.8],
                )
                if labels:
                    ax.text(
                        r0 + 0.5 * r if not params else r0 + 0.25 * r,
                        np.mean(self._ml.aq.z[i : i + 2]),
                        f"leaky layer {lli}",
                        ha="center",
                        va="center",
                    )
                if params:
                    ax.text(
                        r0 + 0.75 * r if labels else r0 + 0.5 * r,
                        np.mean(self._ml.aq.z[i : i + 2]),
                        (f"$c$ = {self._ml.aq.c[lli]}"),
                        ha="center",
                        va="center",
                    )
                if labels or params:
                    lli += 1

            # aquifers
            if labels and self._ml.aq.ltype[i] == "a":
                ax.text(
                    r0 + 0.5 * r if not params else r0 + 0.25 * r,
                    np.mean(self._ml.aq.z[i : i + 2]),
                    f"aquifer {aqi}",
                    ha="center",
                    va="center",
                )
            if params and self._ml.aq.ltype[i] == "a":
                if aqi == 0:
                    paramtxt = f"$k_h$ = {self._ml.aq.kaq[aqi]}"
                if self._ml.name == "Model3D":
                    paramtxt += f", $k_z/k_h$ = {self._ml.aq.kzoverkh[aqi]:.2f}"
                ax.text(
                    r0 + 0.75 * r if labels else r0 + 0.5 * r,
                    np.mean(self._ml.aq.z[i : i + 2]),
                    paramtxt,
                    ha="center",
                    va="center",
                )
            if (labels or params) and self._ml.aq.ltype[i] == "a":
                aqi += 1

        # aquifer-aquifer boundaries (for e.g. Model3D)
        for i in range(1, self._ml.aq.nlayers):
            if self._ml.aq.ltype[i] == "a" and self._ml.aq.ltype[i - 1] == "a":
                ax.axhspan(
                    ymin=self._ml.aq.z[i], ymax=self._ml.aq.z[i], color=[0.8, 0.8, 0.8]
                )
        # top and bottom
        ax.axhline(self._ml.aq.z[0], color="k", lw=0.75)
        ax.axhline(self._ml.aq.z[-1], color="k", lw=3.0)
        # add y-label
        ax.set_ylabel("elevation")
        return ax

    def contour(
        self,
        win,
        ngr=20,
        layers=0,
        levels=20,
        layout=True,
        labels=True,
        decimals=0,
        color=None,
        newfig=True,
        figsize=None,
        legend=True,
        **kwargs,
    ):
        """Head contour plot.

        Parameters
        ----------
        win : list or tuple
            [xmin, xmax, ymin, ymax]
        ngr : scalar, tuple or list
            if scalar: number of grid points in x and y direction
            if tuple or list: nx, ny, number of grid points in x and y
            directions
        layers : integer, list or array
            layers for which grid is returned
        levels : integer or array (default 20)
            levels that are contoured
        layout : boolean (default True)
            plot layout of elements
        labels : boolean (default True)
            print labels along contours
        decimals : integer (default 0)
            number of decimals of labels along contours
        color : str or list of strings
            color of contour lines
        newfig : boolean (default True)
            create new figure
        figsize : tuple of 2 values (default is mpl default)
            size of figure
        legend : list or boolean (default True)
            add legend to figure
            if list of strings: use strings as names in legend

        Returns
        -------
        cs : list of contour sets for each contoured layer
        """
        x1, x2, y1, y2 = win
        if np.isscalar(ngr):
            nx = ny = ngr
        else:
            nx, ny = ngr
        layers = np.atleast_1d(layers)
        xg = np.linspace(x1, x2, nx)
        yg = np.linspace(y1, y2, ny)
        h = self._ml.headgrid(xg, yg, layers)
        if newfig:
            plt.figure(figsize=figsize)
        # color
        if color is None:
            c = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        elif isinstance(color, str):
            c = len(layers) * [color]
        elif isinstance(color, list):
            c = color
        if len(c) < len(layers):
            n = np.ceil(self._ml.aq.naq / len(c))
            c = n * c
        # contour
        cslist = []
        cshandlelist = []
        for i in range(len(layers)):
            cs = plt.contour(xg, yg, h[i], levels, colors=c[i], **kwargs)
            cslist.append(cs)
            handles, _ = cs.legend_elements()
            cshandlelist.append(handles[0])
            if labels:
                fmt = "%1." + str(decimals) + "f"
                plt.clabel(cs, fmt=fmt)
        if isinstance(legend, list):
            plt.legend(cshandlelist, legend)
        elif legend:
            legendlist = ["layer " + str(i) for i in layers]
            plt.legend(cshandlelist, legendlist)
        plt.axis("scaled")
        if layout:
            self.topview(win=[x1, x2, y1, y2], newfig=False, layer=layers)
        return cslist

    def vcontour(
        self,
        win,
        n,
        levels=20,
        labels=True,
        decimals=0,
        color=None,
        vinterp=True,
        nudge=1e-6,
        newfig=True,
        figsize=None,
        layout=True,
    ):
        """Head contour plot in vertical cross-section.

        Parameters
        ----------
        win : list or tuple
            [xmin, xmax, ymin, ymax]
        n : integer
            number of grid points along cross-section
        levels : integer or array (default 20)
            levels that are contoured
        labels : boolean (default True)
            print labels along contours
        decimals : integer (default 0)
            number of decimals of labels along contours
        color : str or list of strings
            color of contour lines
        vinterp : boolean
            when True, interpolate between centers of layers
            when False, constant value vertically in each layer
        nudge : float
            first value is computed nudge from the specified window
        newfig : boolean (default True)
            create new figure
        figsize : tuple of 2 values (default is mpl default)
            size of figure
        layout : boolean
            plot layout if True

        Returns
        -------
        cs : contour set
        """
        x1, x2, y1, y2 = win
        h = self._ml.headalongline(
            np.linspace(x1 + nudge, x2 - nudge, n),
            np.linspace(y1 + nudge, y2 - nudge, n),
        )
        L = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        xg = np.linspace(0, L, n)
        if vinterp:
            zg = 0.5 * (self._ml.aq.zaqbot + self._ml.aq.zaqtop)
            zg = np.hstack((self._ml.aq.zaqtop[0], zg, self._ml.aq.zaqbot[-1]))
            h = np.vstack((h[0], h, h[-1]))
        else:
            zg = np.empty(2 * self._ml.aq.naq)
            for i in range(self._ml.aq.naq):
                zg[2 * i] = self._ml.aq.zaqtop[i]
                zg[2 * i + 1] = self._ml.aq.zaqbot[i]
            h = np.repeat(h, 2, 0)
        if newfig:
            plt.figure(figsize=figsize)
        cs = plt.contour(xg, zg, h, levels, colors=color)
        if labels:
            fmt = "%1." + str(decimals) + "f"
            plt.clabel(cs, fmt=fmt)
        if layout:
            self.topview(win=[x1, x2, y1, y2], orientation="ver", newfig=False)
        return cs

    def tracelines(
        self,
        xstart,
        ystart,
        zstart,
        hstepmax,
        vstepfrac=0.2,
        tmax=1e12,
        nstepmax=100,
        silent=".",
        color=None,
        orientation="hor",
        win=None,
        newfig=False,
        figsize=None,
        *,
        return_traces=False,
        metadata=False,
    ):
        """Function to trace multiple pathlines.

        Parameters
        ----------
        ml : Model object
            model to which the element is added
        xstart : array
            x-coordinates of starting locations
        ystart : array
            y-coordinates of starting locations
        zstart : array
            z-coordinates of starting locations
        hstepmax : scalar
            maximum horizontal step size [L]
        vstepfrac : scalar
            maximum vertical step as fraction of layer thickness
        tmax : scalar
            maximum travel time
        nstepmax : int
            maximum number of steps
        silent : string
            if '.', prints dot upon completion of each traceline
        color : string
            matplotlib color of traceline
        orientation : ('hor', 'ver', 'both')
            'hor' for horizontal, 'ver' for vertical
            'both' for horizontal above vertical
        win : list
            list with [xmin, xmax, ymin, ymax]
        figsize : 2tuple
            figure size
        return_traces : boolean
            return traces if True
        metadata: boolean
            if False, return list of xyzt arrays
            if True, return list of result dicionaries

        Returns
        -------
        traces : result
            only if return_traces = True
        """
        if win is None:
            win = [-1e30, 1e30, -1e30, 1e30]
        if color is None:
            c = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        elif isinstance(color, str):
            c = self._ml.aq.naq * [color]
        elif isinstance(color, list):
            c = color
        if len(c) < self._ml.aq.naq:
            n = int(np.ceil(self._ml.aq.naq / len(c)))
            c = n * c
        fig = plt.gcf()
        assert len(fig.axes) > 0, (
            "Error: Need to specify axes in figure before invoking tracelines"
        )

        axes = {}
        if orientation == "both":
            axes["hor"] = fig.axes[0]
            axes["ver"] = fig.axes[1]
        elif orientation[:3] == "hor":
            axes["hor"] = fig.axes[0]
        elif orientation[:3] == "ver":
            axes["ver"] = fig.axes[-1]

        if return_traces:
            traces = []
        else:
            metadata = True  # suppress future warning from timtraceline
        for i, _ in enumerate(xstart):
            trace = timtraceline(
                self._ml,
                xstart[i],
                ystart[i],
                zstart[i],
                hstepmax=hstepmax,
                vstepfrac=vstepfrac,
                tmax=tmax,
                nstepmax=nstepmax,
                silent=silent,
                win=win,
                returnlayers=True,
                metadata=metadata,
            )
            if return_traces:
                traces.append(trace)
            if metadata:
                xyzt, layerlist = trace["trace"], trace["layers"]
            else:
                xyzt, layerlist = trace
            if silent == ".":
                print(".", end="", flush=True)
            if "hor" in axes:
                color = []
                for ixyzt, ilayer in zip(xyzt, layerlist, strict=False):
                    aq = self._ml.aq.find_aquifer_data(ixyzt[0], ixyzt[1])
                    color.append(
                        c[aq.layernumber[ilayer]] if aq.ltype[ilayer] == "a" else "k"
                    )
                points = np.array([xyzt[:, 0], xyzt[:, 1]]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = LineCollection(segments, colors=color)
                axes["hor"].add_collection(lc)
                # ax1.plot(xyzt[:, 0], xyzt[:, 1], color=color)
            if "ver" in axes:
                color = []
                for ixyzt, ilayer in zip(xyzt, layerlist, strict=False):
                    aq = self._ml.aq.find_aquifer_data(ixyzt[0], ixyzt[1])
                    color.append(
                        c[aq.layernumber[ilayer]] if aq.ltype[ilayer] == "a" else "k"
                    )
                points = np.array([xyzt[:, 0], xyzt[:, 2]]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = LineCollection(segments, colors=color)
                axes["ver"].add_collection(lc)
                axes["ver"].set_ylim(aq.z[-1], aq.z[0])
        if silent == ".":
            print("")  # Print the final newline after the dots
        if return_traces:
            return traces

    def vcontoursf1D(
        self,
        x1,
        x2,
        nx,
        levels,
        labels=False,
        decimals=0,
        color=None,
        nudge=1e-6,
        newfig=True,
        figsize=None,
        layout=True,
        ax=None,
    ):
        """Contour plot in vertical cross-section of 1D model.

        Parameters
        ----------
        x1 : scalar
            left edge of contour domain
        x2 : scalar
            right edge of contour domain
        nx : integer
            number of grid points along cross-section
        levels : integer or array (default 20)
            levels that are contoured
        labels : boolean (default True)
            print labels along contours
        decimals : integer (default 0)
            number of decimals of labels along contours
        color : str or list of strings
            color of contour lines
        nudge : float
            first value is computed nudge from the specified x1 and x2
        newfig : boolean (default True)
            create new figure. if False, ax must be supplied
        figsize : tuple of 2 values (default is mpl default)
            size of figure
        layout : boolean
            plot layout if True
        ax : matplotlib axis
            add plot to specified axis

        Returns
        -------
        ax : axis
        """
        naq = self._ml.aq.naq
        xflow = np.linspace(x1 + nudge, x2 - nudge, nx)
        Qx = np.empty((naq, nx))
        for i in range(nx):
            Qx[:, i], _ = self._ml.disvec(xflow[i], 0)
        zflow = np.empty(2 * naq)
        for i in range(self._ml.aq.naq):
            aq = self._ml.aq.find_aquifer_data(xflow[0], 0)  # use first x as reference
            zflow[2 * i] = aq.zaqtop[i]
            zflow[2 * i + 1] = aq.zaqbot[i]
        Qx = Qx[::-1]  # set upside down
        Qxgrid = np.empty((2 * naq, nx))
        Qxgrid[0] = 0
        for i in range(naq - 1):
            Qxgrid[2 * i + 1] = Qxgrid[2 * i] - Qx[i]
            Qxgrid[2 * i + 2] = Qxgrid[2 * i + 1]
        Qxgrid[-1] = Qxgrid[-2] - Qx[-1]
        Qxgrid = Qxgrid[::-1]  # index 0 at top
        if figsize is not None:
            plt.figure(figsize=figsize)
        if ax:
            ax = ax
        else:
            ax = plt.subplot(111)
        cs = ax.contour(xflow, zflow, Qxgrid, levels, colors=color)
        if labels:
            fmt = "%1." + str(decimals) + "f"
            plt.clabel(cs, fmt=fmt)
        if layout:
            self.xsection(xy=[(x1, 0), (x2, 0)], labels=False, ax=ax)
        return ax

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from .trace import timtraceline, timtracelines

plt.rcParams["contour.negative_linestyle"] = "solid"


class PlotTim:
    def plot(
        self,
        win=None,
        newfig=True,
        figsize=None,
        orientation="hor",
        topfigfrac=0.8,
        layer=None,
    ):
        """Plot layout

        Parameters
        ----------

        win : list or tuple
            [x1, x2, y1, y2]

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
            for e in self.elementlist:
                e.plot(layer=layer)
            if orientation[:3] == "hor":
                plt.axis("scaled")
            elif orientation == "both":
                plt.axis("equal")  # cannot be 'scaled' when sharing axes
            if win is not None:
                plt.axis(win)
        if ax2 is not None:
            plt.sca(ax2)
            for i in range(self.aq.nlayers):
                if self.aq.ltype[i] == "l":
                    plt.axhspan(
                        ymin=self.aq.z[i + 1], ymax=self.aq.z[i], color=[0.8, 0.8, 0.8]
                    )
            for i in range(1, self.aq.nlayers):
                if self.aq.ltype[i] == "a" and self.aq.ltype[i - 1] == "a":
                    plt.axhspan(
                        ymin=self.aq.z[i], ymax=self.aq.z[i], color=[0.8, 0.8, 0.8]
                    )

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
        """Contour plot

        Parameters
        ----------

        win : list or tuple
            [x1, x2, y1, y2]
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
        h = self.headgrid(xg, yg, layers)
        if newfig:
            plt.figure(figsize=figsize)
        # color
        if color is None:
            c = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        elif type(color) is str:
            c = len(layers) * [color]
        elif type(color) is list:
            c = color
        if len(c) < len(layers):
            n = np.ceil(self.aq.naq / len(c))
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
        if type(legend) is list:
            plt.legend(cshandlelist, legend)
        elif legend:
            legendlist = ["layer " + str(i) for i in layers]
            plt.legend(cshandlelist, legendlist)
        plt.axis("scaled")
        if layout:
            self.plot(win=[x1, x2, y1, y2], newfig=False, layer=layers)
        return cslist

    def vcontour(
        self,
        win,
        n,
        levels,
        labels=False,
        decimals=0,
        color=None,
        vinterp=True,
        nudge=1e-6,
        newfig=True,
        figsize=None,
        layout=True,
    ):
        """
        Vertical contour

        Returns
        -------

        cs : contour set
        """

        x1, x2, y1, y2 = win
        h = self.headalongline(
            np.linspace(x1 + nudge, x2 - nudge, n),
            np.linspace(y1 + nudge, y2 - nudge, n),
        )
        L = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        xg = np.linspace(0, L, n)
        if vinterp:
            zg = 0.5 * (self.aq.zaqbot + self.aq.zaqtop)
            zg = np.hstack((self.aq.zaqtop[0], zg, self.aq.zaqbot[-1]))
            h = np.vstack((h[0], h, h[-1]))
        else:
            zg = np.empty(2 * self.aq.naq)
            for i in range(self.aq.naq):
                zg[2 * i] = self.aq.zaqtop[i]
                zg[2 * i + 1] = self.aq.zaqbot[i]
            h = np.repeat(h, 2, 0)
        if newfig:
            plt.figure(figsize=figsize)
        cs = plt.contour(xg, zg, h, levels, colors=color)
        if labels:
            fmt = "%1." + str(decimals) + "f"
            plt.clabel(cs, fmt=fmt)
        if layout:
            self.plot(win=[x1, x2, y1, y2], orientation="ver", newfig=False)
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
        win=[-1e30, 1e30, -1e30, 1e30],
        newfig=False,
        figsize=None,
        *,
        return_traces=False,
        metadata=False,
    ):
        """Draw trace lines"""
        if color is None:
            c = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        elif type(color) is str:
            c = self.aq.naq * [color]
        elif type(color) is list:
            c = color
        if len(c) < self.aq.naq:
            n = int(np.ceil(self.aq.naq / len(c)))
            c = n * c
        fig = plt.gcf()
        assert (
            len(fig.axes) > 0
        ), "Error: Need to specify axes in figure before invoking tracelines"
        ax1 = None
        ax2 = None
        if orientation == "both":
            ax1 = fig.axes[0]
            ax2 = fig.axes[1]
        elif orientation[:3] == "hor":
            ax1 = fig.axes[0]
        elif orientation[:3] == "ver":
            ax2 = fig.axes[1]
        if return_traces:
            traces = []
        else:
            metadata = True  # suppress future warning from timtraceline
        for i in range(len(xstart)):
            trace = timtraceline(
                self,
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
            if ax1 is not None:
                # plt.axes(ax1)
                color = [
                    c[self.aq.layernumber[i]] if self.aq.ltype[i] == "a" else "k"
                    for i in layerlist
                ]
                points = np.array([xyzt[:, 0], xyzt[:, 1]]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = LineCollection(segments, colors=color)
                ax1.add_collection(lc)
                # ax1.plot(xyzt[:, 0], xyzt[:, 1], color=color)
            if ax2 is not None:
                color = [
                    c[self.aq.layernumber[i]] if self.aq.ltype[i] == "a" else "k"
                    for i in layerlist
                ]
                points = np.array([xyzt[:, 0], xyzt[:, 2]]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = LineCollection(segments, colors=color)
                ax2.add_collection(lc)
                ax2.set_ylim(self.aq.z[-1], self.aq.z[0])
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
        """
        Vertical contour for 1D model

        """
        naq = self.aq.naq
        xflow = np.linspace(x1 + nudge, x2 - nudge, nx)
        Qx = np.empty((naq, nx))
        for i in range(nx):
            Qx[:, i], _ = self.disvec(xflow[i], 0)
        zflow = np.empty(2 * naq)
        for i in range(self.aq.naq):
            zflow[2 * i] = self.aq.zaqtop[i]
            zflow[2 * i + 1] = self.aq.zaqbot[i]
        Qx = Qx[::-1]  # set upside down
        Qxgrid = np.empty((2 * naq, nx))
        Qxgrid[0] = 0
        for i in range(naq - 1):
            Qxgrid[2 * i + 1] = Qxgrid[2 * i] - Qx[i]
            Qxgrid[2 * i + 2] = Qxgrid[2 * i + 1]
        Qxgrid[-1] = Qxgrid[-2] - Qx[-1]
        Qxgrid = Qxgrid[::-1]  # index 0 at top
        if newfig:
            _, ax = plt.subplots(1, 1, figsize=figsize)
        else:
            ax = ax
        cs = ax.contour(xflow, zflow, Qxgrid, levels, colors=color)
        if labels:
            fmt = "%1." + str(decimals) + "f"
            plt.clabel(cs, fmt=fmt)
        return ax
        # if layout:
        #    self.plot(win=[x1, x2, y1, y2], orientation='ver', newfig=False)


def compute_z1z2(xy):
    # Returns z1 and z2 of polygon, in clockwise order
    x, y = list(zip(*xy))
    if x[0] == x[-1] and y[0] == y[-1]:  # In case last point is repeated
        x = x[:-1]
        y = y[:-1]
    z1 = np.array(x) + np.array(y) * 1j
    index = list(range(1, len(z1))) + [0]
    z2 = z1[index]
    Z = 1e-6j
    z = Z * (z2[0] - z1[0]) / 2.0 + 0.5 * (z1[0] + z2[0])
    bigZ = (2.0 * z - (z1 + z2)) / (z2 - z1)
    bigZmin1 = bigZ - 1.0
    bigZplus1 = bigZ + 1.0
    angle = np.sum(np.log(bigZmin1 / bigZplus1).imag)
    if angle < np.pi:  # reverse order
        z1 = z1[::-1]
        z2 = z1[index]
    return z1, z2


def refine_n_segments(xy, shape_type, n_segments):
    """Refine line segments into n_segments each.

    Use controlpoints half-circle approach to determine new segment lengths.

    Parameters
    ----------
    xy : list of tuple or np.array
        list of coordinates or 2d-array containing x-coordinates in the first column
        and y-coordinates in the second column
    shape_type : str
        shape type, either "line" or "polygon".
    n_segments : int
        number of segments to split each line segment into.

    Returns
    -------
    xynew : np.array
        array containing refined coordinates
    reindexer : np.array
        array containing index to original line segment, useful for obtaining element
        parameters from original input.
    """

    if shape_type == "polygon":
        z1, z2 = compute_z1z2(xy)
    elif shape_type == "line":
        z = xy[:, 0] + 1j * xy[:, 1]
        z1 = z[:-1]
        z2 = z[1:]
    else:
        raise ValueError("shptype must be one of 'polygon' or 'line'.")
    xpts = []
    ypts = []
    reindexer = []
    for i in range(len(z1)):
        xcpi, ycpi = controlpoints(n_segments - 1, z1[i], z2[i], include_ends=True)
        if i < len(z1) - 1:
            xpts.append(xcpi[:-1])
            ypts.append(ycpi[:-1])
        else:
            xpts.append(xcpi)
            ypts.append(ycpi)
        reindexer.append(i * np.ones(len(xcpi[:-1]), dtype=int))
    return (
        np.vstack([np.concatenate(xpts), np.concatenate(ypts)]).T,
        np.concatenate(reindexer),
    )

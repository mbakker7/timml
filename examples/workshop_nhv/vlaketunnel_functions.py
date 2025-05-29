# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

import timml as tml


def create_model(
    kaq=[0.1, 5.0, 15.0, 5.0],
    c=[1000.0, 2.0, 2.0, 2.0],
    hstar=0,
    c_channel_bot=30,
    do_plot=True,
    df_dh=None,
):
    """Create a TimML model for Vlaketunnel case.

    Parameters
    ----------
    kaq : list, optional
        Kh of aquifers. The default is [0.1, 5.0, 15.0, 5.0].
    c : list, optional
        c of aquitards. The default is [1000.0, 2.0, 2.0, 2.0].
    hstar : float, optional
        Top boundary condition of semi-confining toplayer. The default is 0.
    c_channel_bot : float, optional
        resistance of Kanaal door Zuid-Beveland. The default is 30.
    do_plot : boolean, optional
        Plot results? The default is True.
    df_dh: pd.DataFrame, optional
        Information about observed drawdowns, required for plotting.
        The default is None.

    Returns
    -------
    ml : timml model
        The model

    """
    # create model
    if c is None:
        c = [1000.0, 2.0, 2.0, 2.0]
    if kaq is None:
        kaq = [0.1, 5.0, 15.0, 5.0]
    ml = tml.ModelMaq(
        kaq=kaq,
        z=[1.0, -3.0, -7.0, -7.0, -14.0, -14.0, -30.0, -30.0, -40.0],
        c=c,
        topboundary="semi",
        npor=[None, None, None, None, None, None, None, None],
        hstar=hstar,
    )

    # add dewatering
    dewatering_east_xys = [
        [59224, 387382],
        [59359, 387375],
        [59360, 387311],
        [59234, 387298],
    ]
    q_east_total = 325 * 24

    q_west_total = 75 * 24
    dewatering_west_xys = [
        [58781, 387375],
        [58785, 387307],
    ]

    for dewatering_xys, q_total in zip(
        [dewatering_east_xys, dewatering_west_xys],
        [q_east_total, q_west_total],
        strict=False,
    ):
        # loop over both dewatering locations
        for dewatering_xy in dewatering_xys:
            # loop over the modelled wells, in pratice a lot of more wells are used.
            # Current model has focus on regional effect, therefore limited number
            # of wells are considered sufficient

            # dewatering_east:
            tml.Well(
                xw=dewatering_xy[0],
                yw=dewatering_xy[1],
                Qw=q_total / len(dewatering_xys),
                rw=0.5,
                res=1.0,
                layers=1,
                label=None,
                model=ml,
            )

    c_channel = ml.aq.c.copy()
    c_channel[0] = c_channel_bot

    # channel_0:
    tml.PolygonInhomMaq(
        kaq=ml.aq.kaq,
        z=ml.aq.z,
        c=c_channel,
        topboundary="semi",
        npor=[None, None, None, None, None, None, None, None],
        hstar=0.0,
        # compared to QGIS-Tim export the channel is extended to the north in order to
        # cover the northern observation wells better
        xy=[
            [58921, 390500],
            [59065, 390500],
            [59110, 387996],
            [59146, 387447],
            [59263, 386809],
            [59317, 386260],
            [59110, 386251],
            [58966, 386863],
            [58921, 388617],
        ],
        order=4,
        ndeg=6,
        model=ml,
    )
    ml.solve()

    if do_plot and (df_dh is not None):
        plot_model_results(ml, df_dh)

    return ml


def plot_model_input(ml):
    """Plot model input in schematic section.

    Parameters
    ----------
    ml : timml Model
        The model

    Returns
    -------
    None.

    """
    # some plotting constants
    xmin = -1
    xchannel = -0.25
    xhinter = -0.2
    xmax = 1
    zaqmid = np.mean([ml.aq.zaqtop, ml.aq.zaqbot], axis=0)

    # plot layers
    plt.hlines(y=ml.aq.zlltop, xmin=xmin, xmax=xmax, color="darkgray")
    plt.hlines(y=ml.aq.zaqbot, xmin=xmin, xmax=xmax, color="darkgray")

    # plot kh
    for kh, z in zip(ml.aq.kaq, zaqmid, strict=False):
        plt.annotate(f"kh={kh:0.1f}m/d", (0, z), ha="center")
    # plot c
    for c, z in zip(ml.aq.c, ml.aq.zaqtop, strict=False):
        plt.annotate(f"c={c:0.1f}d", (0.5, z), ha="center", va="center")
    # plot channel
    plt.plot([xmin, xchannel], [ml.aq.inhomlist[0].hstar] * 2, color="blue")
    plt.annotate(
        f"h_ch={ml.aq.inhomlist[0].hstar:0.1f}",
        (xchannel, ml.aq.inhomlist[0].hstar),
        ha="right",
        va="bottom",
    )
    plt.annotate(
        f"c_ch={ml.aq.inhomlist[0].c[0]:0.1f}",
        (xchannel, ml.aq.zaqtop[0]),
        ha="right",
        va="bottom",
    )

    # plot hinterland
    plt.plot([xhinter, xmax], [ml.aq.hstar] * 2, color="darkblue")
    plt.annotate(
        f"h_polder={ml.aq.hstar:0.1f}", (xhinter, ml.aq.hstar), ha="left", va="bottom"
    )

    plt.xlim([xmin, xmax])


def plot_model_results(ml, df_dh):
    """Plot results of TimML model of Vlaketunnel case.

    Parameters
    ----------
    ml : timml Model,
        The model.
    df_dh : pd.DataFrame
        Observed drawdowns

    Returns
    -------
    None.

    """
    # contour plot
    plt.subplot(221)
    ml.plots.contour(
        win=[57000, 60000, 386900, 389100],
        ngr=50,
        layers=1,
        levels=[-5, -2, -1, -0.5, -0.1],
        labels=True,
        decimals=2,
        legend=False,
        newfig=False,
    )
    plt.scatter(df_dh.x, df_dh.y, 20, c=df_dh.color)
    for _, row in df_dh.iterrows():
        plt.annotate(f"{row.dh_obs:0.2f}", (row.x, row.y), ha=row.ha, va=row.va)
    plt.title("contours in layer 1")
    # plot model input
    plt.subplot(222)
    plot_model_input(ml)

    for plotid in (223, 224):
        plt.subplot(plotid)
        if plotid == 223:
            # first plot, get model results
            df_dh["ml_layer"] = None
            df_dh["dh_calc"] = None
            for index, row in df_dh.iterrows():
                df_dh.loc[index, "ml_layer"] = np.where(ml.aq.zaqtop > row.screen_top)[
                    0
                ][-1]
                df_dh.loc[index, "dh_calc"] = ml.headalongline(
                    row.x, row.y, row.ml_layer
                )[0][0]
            # plot all model results
            plot_df = df_dh
        else:
            # second plot, only plot outside dewatering area
            plot_df = df_dh.loc[df_dh.r > 100]

        plt.scatter(
            plot_df.r, plot_df.dh_obs, 50, c=plot_df.color, alpha=0.3, label="observed"
        )
        plt.scatter(plot_df.r, plot_df.dh_calc, 40, marker="+", label="modelled")
        plt.legend()
        plt.title("heads from screened modellayer")
        plt.grid()

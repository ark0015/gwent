import astropy.units as u
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import z_at_value


def Plot_SNR(
    var_x,
    sample_x,
    var_y,
    sample_y,
    SNRMatrix,
    fig=None,
    ax=None,
    display=True,
    return_plt=False,
    dl_axis=False,
    lb_axis=False,
    smooth_contours=True,
    cfill=True,
    display_cbar=True,
    x_axis_label=True,
    y_axis_label=True,
    x_axis_line=None,
    y_axis_line=None,
    logLevels_min=-1.0,
    logLevels_max=0.0,
    hspace=0.15,
    wspace=0.1,
    contour_kwargs={},
    contourf_kwargs={},
    xticklabels_kwargs={},
    xlabels_kwargs={},
    xline_kwargs={},
    yticklabels_kwargs={},
    ylabels_kwargs={},
    yline_kwargs={},
):
    """Plots the SNR contours from calcSNR

    Parameters
    ----------
    var_x: str
        x-axis variable
    sample_x: array
        samples at which ``SNRMatrix`` was calculated corresponding to the x-axis variable
    var_y: str
        y-axis variable
    sample_y: array
        samples at which ``SNRMatrix`` was calculated corresponding to the y-axis variable
    SNRMatrix: array-like
        the matrix at which the SNR was calculated corresponding to the particular x and y-axis variable choices

    fig: object, optional
        matplotlib figure object on which to collate the individual plots
    ax: object, optional
        matplotlib axes object on which to plot the individual plot
    display: bool, optional
        Option to turn off display if saving multiple plots to a file
    return_plt: bool, optional
        Option to return ``fig`` and ``ax``
    dl_axis: bool, optional
        Option to turn on the right hand side labels of luminosity distance
    lb_axis: bool, optional
        Option to turn on the right hand side labels of lookback time
    smooth_contours: bool, optional
        Option to have contours appear smooth instead of tiered (depending on sample size the edges appear boxey).
    cfill: bool, optional
        Option to use filled contours or not, default is ``True``
    display_cbar: bool, optional
        Option to display the colorbar on the axes object
    x_axis_label: bool, optional
        Option to display the x axis label
    y_axis_label: bool, optional
        Option to display the y axis label
    x_axis_line: int,float, optional
        Option to display a line on the x axis if not None
    y_axis_line: int,float, optional
        Option to display a line on the y axis if not None
    logLevels_min: float, optional
        Sets the minimum log level of the colorbar, default is -1.0 which set the minimum to the log minimum of the given ``SNRMatrix``
    logLevels_max: float, optional
        Sets the maximum log level of the colorbar, default is 0.0, which sets the maximum to the log maximum value of the given ``SNRMatrix``
    hspace: float, optional
        Sets the vertical space between axes objects, default is 0.15
    wspace: float, optional
        Sets the horizontal space between axes objects, default is 0.1
    contour_kwargs: dict, optional
        Sets additional kwargs taken by contour in matplotlib
    contourf_kwargs: dict, optional
        Sets additional kwargs taken by contourf in matplotlib
    xticklabels_kwargs: dict, optional
        Sets additional kwargs taken by xticklabel in matplotlib
    xlabels_kwargs=: dict, optional
        Sets additional kwargs taken by xlabel in matplotlib
    xline_kwargs: dict, optional
        Sets additional kwargs taken by ax.axvline in matplotlib
    yticklabels_kwargs: dict, optional
        Sets additional kwargs taken by yticklabel in matplotlib
    ylabels_kwargs: dict, optional
        Sets additional kwargs taken by ylabel in matplotlib
    yline_kwargs: dict, optional
        Sets additional kwargs taken by ax.axhline in matplotlib

    """
    if fig is not None:
        if ax is not None:
            pass
        else:
            fig, ax = plt.subplots()
    else:
        fig, ax = plt.subplots()

    if "colors" not in contour_kwargs.keys() and "cmap" not in contour_kwargs.keys():
        contour_kwargs["colors"] = "k"
    if "linewidths" not in contour_kwargs.keys():
        contour_kwargs["linewidths"] = 2.0

    if "cmap" not in contourf_kwargs.keys():
        contourf_kwargs["cmap"] = "viridis"

    logSNR = np.log10(SNRMatrix)
    if logLevels_min == -1.0:
        logLevels_min = np.log10(np.array([1.0]))
    if logLevels_max == 0.0:
        logLevels_max = np.ceil(np.amax(logSNR))
    if logLevels_max < logLevels_min:
        raise ValueError("All SNRs are lower than 5.")

    logLevels_add = np.log10(np.array([3.0, 10.0, 31.0]))
    print_logLevels = np.concatenate(
        (logLevels_min, logLevels_add, np.arange(2.0, logLevels_max + 1.0))
    )

    logLevels = print_logLevels

    ylabel_min = min(sample_y)
    ylabel_max = max(sample_y)
    xlabel_min = min(sample_x)
    xlabel_max = max(sample_x)

    # Set whether log or linearly spaced axes
    if xlabel_max < 0.0 or xlabel_min < 0.0 or var_x in ["n_p", "T_obs"]:
        xaxis_type = "lin"
        step_size = int(xlabel_max - xlabel_min + 1)
        x_labels = np.linspace(xlabel_min, xlabel_max, step_size)
    else:
        x_log_range = np.log10(xlabel_max) - np.log10(xlabel_min)
        if x_log_range >= 2.0:
            xaxis_type = "log"
            step_size = int(np.log10(xlabel_max) - np.log10(xlabel_min) + 1)
            x_labels = np.logspace(
                np.log10(xlabel_min), np.log10(xlabel_max), step_size
            )
        else:
            xaxis_type = "lin"
            x_scale = 10 ** round(np.log10(xlabel_min))
            x_labels = (
                np.arange(
                    round(xlabel_min / x_scale), round(xlabel_max / x_scale) + 1, 1
                )
                * x_scale
            )
            if x_labels[0] < xlabel_min:
                x_labels[0] = xlabel_min
            if x_labels[-1] > xlabel_max:
                x_labels[-1] = xlabel_max

    if ylabel_max < 0.0 or ylabel_min < 0.0 or var_y in ["n_p", "T_obs"]:
        yaxis_type = "lin"
        step_size = int(ylabel_max - ylabel_min + 1)
        y_labels = np.linspace(ylabel_min, ylabel_max, step_size)
    else:
        y_log_range = np.log10(ylabel_max) - np.log10(ylabel_min)
        if y_log_range >= 2.0:
            yaxis_type = "log"
            step_size = int(np.log10(ylabel_max) - np.log10(ylabel_min) + 1)
            y_labels = np.logspace(
                np.log10(ylabel_min), np.log10(ylabel_max), step_size
            )
        else:
            yaxis_type = "lin"
            y_scale = 10 ** round(np.log10(ylabel_min))
            y_labels = (
                np.arange(
                    round(ylabel_min / y_scale), round(ylabel_max / y_scale) + 1, 1
                )
                * y_scale
            )
            if y_labels[0] < ylabel_min:
                y_labels[0] = ylabel_min
            if y_labels[-1] > ylabel_max:
                y_labels[-1] = ylabel_max

    # Set axis scales based on what data sampling we used
    if yaxis_type == "lin" and xaxis_type == "log":
        if not cfill:
            CS1 = ax.contour(
                np.log10(sample_x), sample_y, logSNR, print_logLevels, **contour_kwargs
            )
        else:
            if smooth_contours:
                # cmap = mpl.cm.get_cmap(name=contourf_kwargs["cmap"])
                cmap = mpl.colormaps[contourf_kwargs["cmap"]]
                cmap.set_under(color="white")
                CS1 = ax.imshow(
                    logSNR,
                    extent=[
                        np.log10(xlabel_min),
                        np.log10(xlabel_max),
                        ylabel_min,
                        ylabel_max,
                    ],
                    vmin=logLevels_min,
                    vmax=logLevels_max,
                    origin="lower",
                    aspect="auto",
                    cmap=cmap,
                )
            else:
                CS1 = ax.contourf(
                    np.log10(sample_x), sample_y, logSNR, logLevels, **contourf_kwargs
                )
            ax.contour(
                np.log10(sample_x), sample_y, logSNR, print_logLevels, **contour_kwargs
            )
        ax.set_xlim(np.log10(xlabel_min), np.log10(xlabel_max))
        ax.set_ylim(ylabel_min, ylabel_max)

    elif yaxis_type == "log" and xaxis_type == "lin":
        if not cfill:
            CS1 = ax.contour(
                sample_x, np.log10(sample_y), logSNR, print_logLevels, **contour_kwargs
            )
        else:
            if smooth_contours:
                # cmap = mpl.cm.get_cmap(name=contourf_kwargs["cmap"])
                cmap = mpl.colormaps[contourf_kwargs["cmap"]]
                cmap.set_under(color="white")
                CS1 = ax.imshow(
                    logSNR,
                    extent=[
                        xlabel_min,
                        xlabel_max,
                        np.log10(ylabel_min),
                        np.log10(ylabel_max),
                    ],
                    vmin=logLevels_min,
                    vmax=logLevels_max,
                    origin="lower",
                    aspect="auto",
                    cmap=cmap,
                )
            else:
                CS1 = ax.contourf(
                    sample_x, np.log10(sample_y), logSNR, logLevels, **contourf_kwargs
                )
            ax.contour(
                sample_x, np.log10(sample_y), logSNR, print_logLevels, **contour_kwargs
            )
        ax.set_xlim(xlabel_min, xlabel_max)
        ax.set_ylim(np.log10(ylabel_min), np.log10(ylabel_max))
    elif yaxis_type == "lin" and xaxis_type == "lin":
        if not cfill:
            CS1 = ax.contour(
                sample_x, sample_y, logSNR, print_logLevels, **contour_kwargs
            )
        else:
            if smooth_contours:
                # cmap = mpl.cm.get_cmap(name=contourf_kwargs["cmap"])
                cmap = mpl.colormaps[contourf_kwargs["cmap"]]
                cmap.set_under(color="white")
                CS1 = ax.imshow(
                    logSNR,
                    extent=[xlabel_min, xlabel_max, ylabel_min, ylabel_max],
                    vmin=logLevels_min,
                    vmax=logLevels_max,
                    origin="lower",
                    aspect="auto",
                    cmap=cmap,
                )
            else:
                CS1 = ax.contourf(
                    sample_x, sample_y, logSNR, logLevels, **contourf_kwargs
                )
            ax.contour(sample_x, sample_y, logSNR, print_logLevels, **contour_kwargs)
        ax.set_xlim(xlabel_min, xlabel_max)
        ax.set_ylim(ylabel_min, ylabel_max)
    else:
        if not cfill:
            CS1 = ax.contour(
                np.log10(sample_x),
                np.log10(sample_y),
                logSNR,
                print_logLevels,
                **contour_kwargs,
            )
        else:
            if smooth_contours:
                # cmap = mpl.cm.get_cmap(name=contourf_kwargs["cmap"])
                cmap = mpl.colormaps[contourf_kwargs["cmap"]]
                cmap.set_under(color="white")
                CS1 = ax.imshow(
                    logSNR,
                    extent=[
                        np.log10(xlabel_min),
                        np.log10(xlabel_max),
                        np.log10(ylabel_min),
                        np.log10(ylabel_max),
                    ],
                    vmin=logLevels_min,
                    vmax=logLevels_max,
                    origin="lower",
                    aspect="auto",
                    cmap=cmap,
                    interpolation="None",
                )
            else:
                CS1 = ax.contourf(
                    np.log10(sample_x),
                    np.log10(sample_y),
                    logSNR,
                    logLevels,
                    **contourf_kwargs,
                )
            ax.contour(
                np.log10(sample_x),
                np.log10(sample_y),
                logSNR,
                print_logLevels,
                **contour_kwargs,
            )
        ax.set_xlim(np.log10(xlabel_min), np.log10(xlabel_max))
        ax.set_ylim(np.log10(ylabel_min), np.log10(ylabel_max))

    Get_Axes_Labels(
        ax,
        "x",
        var_x,
        xaxis_type,
        x_labels,
        x_axis_line,
        xlabels_kwargs,
        xticklabels_kwargs,
        xline_kwargs,
    )
    Get_Axes_Labels(
        ax,
        "y",
        var_y,
        yaxis_type,
        y_labels,
        y_axis_line,
        ylabels_kwargs,
        yticklabels_kwargs,
        yline_kwargs,
    )

    if not x_axis_label:
        ax.set_xticklabels("")
        ax.set_xlabel("")
    if not y_axis_label:
        ax.set_yticklabels("")
        ax.set_ylabel("")

    # If true, display luminosity distance on right side of plot
    if dl_axis:
        if var_y != "z":
            raise ValueError(
                "Sorry, we can only plot luminosity distance when redshift is on the y axis."
            )

        # Set other side y-axis for luminosity distance scalings
        ax2 = ax.twinx()
        # Set axis scales based on what data sampling we used
        if yaxis_type == "lin" and xaxis_type == "log":
            ax2.contour(
                np.log10(sample_x), sample_y, logSNR, print_logLevels, **contour_kwargs
            )
        elif yaxis_type == "log" and xaxis_type == "lin":
            ax2.contour(
                sample_x, np.log10(sample_y), logSNR, print_logLevels, **contour_kwargs
            )
        else:
            ax2.contour(
                np.log10(sample_x),
                np.log10(sample_y),
                logSNR,
                print_logLevels,
                **contour_kwargs,
            )

        dists_min = cosmo.luminosity_distance(ylabel_min).to("Gpc")
        dists_min = np.ceil(np.log10(dists_min.value))
        dists_max = cosmo.luminosity_distance(ylabel_max).to("Gpc")
        dists_max = np.ceil(np.log10(dists_max.value))
        dists = np.arange(dists_min, dists_max)
        dists = 10**dists * u.Gpc

        distticks = [
            z_at_value(cosmo.luminosity_distance, dist, method="Bounded")
            for dist in dists
        ]
        # Set other side y-axis for lookback time scalings
        ax2.set_yticks(np.log10(distticks))
        ax2.set_yticklabels(
            [
                r"$10^{%i}$" % np.log10(dist)
                if np.abs(int(np.log10(dist))) > 1
                else "{:g}".format(dist)
                for dist in dists.value
            ]
        )
        ax2.set_ylabel(r"$D_{L}$ [Gpc]")
    elif lb_axis:
        if var_y != "z":
            raise ValueError(
                "Sorry, we can only plot lookback time when redshift is on the y axis."
            )
        # Set other side y-axis for lookback time scalings
        ax2 = ax.twinx()
        # Set axis scales based on what data sampling we used
        if yaxis_type == "lin" and xaxis_type == "log":
            ax2.contour(
                np.log10(sample_x), sample_y, logSNR, print_logLevels, **contour_kwargs
            )
        elif yaxis_type == "log" and xaxis_type == "lin":
            ax2.contour(
                sample_x, np.log10(sample_y), logSNR, print_logLevels, **contour_kwargs
            )
        else:
            ax2.contour(
                np.log10(sample_x),
                np.log10(sample_y),
                logSNR,
                print_logLevels,
                **contour_kwargs,
            )

        ages1 = np.array([13.5, 13, 10, 5, 1]) * u.Gyr
        ages2 = np.array([500, 100, 10, 1]) * u.Myr
        ages2 = ages2.to("Gyr")
        ages = np.hstack((ages1.value, ages2.value))
        ages = ages * u.Gyr
        ageticks = [z_at_value(cosmo.age, age, method="Bounded") for age in ages]

        # Set axes limits
        ax2.set_yticks(np.log10(ageticks))
        ax2.set_yticklabels(["{:g}".format(age) for age in ages.value])
        ax2.set_ylabel(r"$t_{\rm cosmic}$ [Gyr]")
        ax2.yaxis.set_label_coords(1.2, 0.5)

    if display_cbar:
        if lb_axis or dl_axis:
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.9, 0.15, 0.025, 0.7])
            # Make colorbar
            if not cfill:
                # Make colorbar
                norm = colors.Normalize(vmin=logLevels_min, vmax=logLevels_max)
                tick_levels = np.linspace(
                    float(logLevels_min), logLevels_max, len(print_logLevels)
                )
                cbar = mpl.colorbar.ColorbarBase(
                    cbar_ax,
                    ax=(ax, ax2),
                    pad=0.01,
                    cmap=CS1.cmap,
                    norm=norm,
                    boundaries=tick_levels,
                    ticks=tick_levels,
                    spacing="proportional",
                )
            else:
                cbar = fig.colorbar(CS1, cax=cbar_ax, ax=(ax, ax2), pad=0.01)
        else:
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.82, 0.15, 0.025, 0.7])
            if not cfill:
                # Make colorbar
                norm = colors.Normalize(vmin=logLevels_min, vmax=logLevels_max)
                tick_levels = np.linspace(
                    float(logLevels_min), logLevels_max, len(print_logLevels)
                )
                cbar = mpl.colorbar.ColorbarBase(
                    cbar_ax,
                    cmap=CS1.cmap,
                    norm=norm,
                    boundaries=tick_levels,
                    ticks=tick_levels,
                    spacing="proportional",
                )
            else:
                # Make colorbar
                cbar = fig.colorbar(CS1, cax=cbar_ax, ticks=print_logLevels)

        cbar.set_label(r"SNR")
        cbar.ax.set_yticklabels(
            [
                r"$10^{%i}$" % x if int(x) > 1 else r"$%i$" % (10**x)
                for x in print_logLevels
            ],
            **yticklabels_kwargs,
        )

    if display:
        # fig.tight_layout()
        fig.subplots_adjust(hspace=hspace, wspace=wspace)
        plt.show()

    if return_plt:
        return fig, ax


def Get_Axes_Labels(
    ax,
    var_axis,
    var,
    var_scale,
    orig_labels,
    line_val,
    label_kwargs,
    tick_label_kwargs,
    line_kwargs,
):
    """Gives paper plot labels for given axis

    Parameters
    ----------
    ax: object
        The current axes object
    var_axis: str
        The axis to change labels and ticks, can either be ``'y'`` or ``'x'``
    var: str
        The variable to label
    orig_labels: list,np.ndarray
        The original labels for the particular axis, may be updated depending on parameter
    line_val: int,float
        Value of line plotted on ``var_axis`` if not None. Assumed to be non-log10 value
    label_kwargs: dict
        The dictionary adjusting the particular axis' label kwargs
    tick_label_kwargs: dict
        The dictionary adjusting the particular axis' tick label kwargs
    line_kwargs: dict
        The dictionary associated with the line displayed on ``var_axis``

    """

    # Set axes labels and whether log or linearly spaced
    if var_axis not in ["y", "x"]:
        raise ValueError("var_axis can only by x or y")

    ax_dict = {}

    if var == "M":
        ax_dict[var_axis + "ticks"] = np.log10(orig_labels)
        ax_dict[var_axis + "label"] = r"$M_{\mathrm{tot}}~[M_{\odot}]$"
        ax_dict[var_axis + "ticklabels"] = [
            r"$10^{%i}$" % x if int(x) > 1 else r"$%i$" % (10**x)
            for x in np.log10(orig_labels)
        ]
    elif var == "q":
        new_labels = orig_labels[::2]
        ax_dict[var_axis + "ticks"] = new_labels
        ax_dict[var_axis + "label"] = r"$\mathrm{Mass~Ratio}~q$"
        ax_dict[var_axis + "ticklabels"] = [r"$%i$" % int(x) for x in new_labels]
    elif var == "z":
        ax_dict[var_axis + "ticks"] = np.log10(orig_labels)
        ax_dict[var_axis + "label"] = r"$\mathrm{Redshift}~z$"
        ax_dict[var_axis + "ticklabels"] = [
            x if int(x) < 1 else int(x) for x in orig_labels
        ]
    elif var in ["chi1", "chi2", "chii"]:
        new_labels = (
            np.arange(round(min(orig_labels) * 10), round(max(orig_labels) * 10) + 1, 1)
            / 10
        )
        new_labels = new_labels[::2]
        ax_dict[var_axis + "ticks"] = new_labels
        ax_dict[var_axis + "label"] = r"$\mathrm{Spin}~\chi_{i}$"
        ax_dict[var_axis + "ticklabels"] = [r"$%.1f$" % x for x in new_labels]
    elif var == "L":
        # Proposed Value = 2.5Gm: L3 LISA
        ax_dict[var_axis + "ticks"] = np.log10(orig_labels)
        ax_dict[var_axis + "label"] = r"Arm Length [m]"
        ax_dict[var_axis + "ticklabels"] = [
            r"$10^{%i}$" % x if int(x) > 1 else r"$%i$" % (10**x)
            for x in np.log10(orig_labels)
        ]
    elif var == "A_acc":
        # Proposed Value = 3x10^{-15}: L3 LISA
        ax_dict[var_axis + "ticks"] = np.log10(orig_labels)
        ax_dict[var_axis + "label"] = r"$A_{\mathrm{acc}} [\mathrm{m~s^{-2}}]$"
        ax_dict[var_axis + "ticklabels"] = [
            r"$10^{%.0f}$" % x for x in np.log10(orig_labels)
        ]
    elif var == "A_IFO":
        # Proposed Value = 10^{-12}: L3 LISA
        if var_scale == "log":
            ax_dict[var_axis + "ticks"] = np.log10(orig_labels)
            ax_dict[var_axis + "label"] = r"$A_{\mathrm{IFO}}$ [m]"
            ax_dict[var_axis + "ticklabels"] = [
                r"$10^{%.0f}$" % x for x in np.log10(orig_labels)
            ]
        elif var_scale == "lin":
            scale = 10 ** round(np.log10(min(orig_labels)))
            new_labels = (
                np.arange(
                    round(min(orig_labels) / scale),
                    round(max(orig_labels) / scale) + 1,
                    3,
                )
                * scale
            )
            ax_dict[var_axis + "ticks"] = new_labels
            ax_dict[var_axis + "label"] = r"$A_{\mathrm{IFO}}$ [pm]"
            ax_dict[var_axis + "ticklabels"] = [
                r"$%.1f$" % (x) for x in new_labels / scale
            ]
    elif var == "f_acc_break_low":
        # Proposed Value = 0.4mHz: L3 LISA
        scale = 10 ** round(np.log10(min(orig_labels)))
        new_labels = (
            np.arange(
                round(min(orig_labels) / scale), round(max(orig_labels) / scale) + 1, 1
            )
            * scale
        )
        ax_dict[var_axis + "ticks"] = new_labels
        ax_dict[var_axis + "label"] = r"$f_{\mathrm{acc,low}}$ [mHz]"
        ax_dict[var_axis + "ticklabels"] = [r"$%.1f$" % x for x in new_labels * 1e3]
    elif var == "f_acc_break_high":
        # Proposed Value = 8mHz: L3 LISA
        scale = 10 ** round(np.log10(min(orig_labels)))
        new_labels = (
            np.arange(
                round(min(orig_labels) / scale), round(max(orig_labels) / scale) + 1, 1
            )
            * scale
        )
        ax_dict[var_axis + "ticks"] = new_labels
        ax_dict[var_axis + "label"] = r"$f_{\mathrm{acc,high}}$ [mHz]"
        ax_dict[var_axis + "ticklabels"] = [r"$%.1f$" % x for x in new_labels * 1e3]
    elif var == "f_IFO_break":
        # Proposed Value = 2mHz: L3 LISA
        scale = 10 ** round(np.log10(min(orig_labels)))
        new_labels = (
            np.arange(
                round(min(orig_labels) / scale), round(max(orig_labels) / scale) + 1, 1
            )
            * scale
        )
        ax_dict[var_axis + "ticks"] = new_labels
        ax_dict[var_axis + "label"] = r"$f_{\mathrm{IFO,break}}$ [mHz]"
        ax_dict[var_axis + "ticklabels"] = [r"$%.1f$" % x for x in new_labels * 1e3]
    elif var == "n_p":
        sample_range = max(orig_labels) - min(orig_labels)
        sample_rate = max(2, int(sample_range / 10))
        new_labels = orig_labels[::sample_rate]
        ax_dict[var_axis + "ticks"] = new_labels
        ax_dict[var_axis + "label"] = r"$\mathrm{Number~of~Pulsars}$"
        ax_dict[var_axis + "ticklabels"] = [r"$%i$" % int(x) for x in new_labels]
    elif var == "cadence":
        new_labels = np.arange(round(min(orig_labels)), round(max(orig_labels)) + 1, 5)
        ax_dict[var_axis + "ticks"] = new_labels
        ax_dict[
            var_axis + "label"
        ] = r"$\mathrm{Observation~Cadence}$ $[\mathrm{yr}^{-1}]$"
        ax_dict[var_axis + "ticklabels"] = [r"$%i$" % int(x) for x in new_labels]
    elif var == "sigma":
        scale = 10 ** round(np.log10(min(orig_labels)))
        new_labels = (
            np.arange(
                round(min(orig_labels) / scale), round(max(orig_labels) / scale) + 1, 1
            )
            * scale
        )
        ax_dict[var_axis + "ticks"] = new_labels
        ax_dict[var_axis + "label"] = r"TOA Error RMS [ns]"
        ax_dict[var_axis + "ticklabels"] = [r"$%.0f$" % x for x in new_labels * 1e9]
    elif var == "T_obs":
        new_labels = orig_labels[::2]
        ax_dict[var_axis + "ticks"] = new_labels
        ax_dict[var_axis + "label"] = r"${\rm T_{obs}}$ [yr]"
        ax_dict[var_axis + "ticklabels"] = [r"$%i$" % int(x) for x in new_labels]
    elif var == "Infrastructure Length":
        # Proposed Value = 3995: aLIGO, Voyager
        # Proposed Value = 40000: CE1
        if var_scale == "log":
            ax_dict[var_axis + "ticks"] = np.log10(orig_labels)
            ax_dict[var_axis + "label"] = r"Infrastructure Length [m]"
            ax_dict[var_axis + "ticklabels"] = [
                r"$10^{%.0f}$" % y if abs(int(y)) > 1 else r"$%.1f$" % (10**y)
                for y in np.log10(orig_labels)
            ]
        elif var_scale == "lin":
            ax_dict[var_axis + "ticks"] = orig_labels
            ax_dict[var_axis + "label"] = r"Infrastructure Length [km]"
            ax_dict[var_axis + "ticklabels"] = [
                r"$%.1f$" % (x / 1e3) for x in orig_labels
            ]
    elif var == "Laser Power":
        # Proposed Value = 125: aLIGO
        # Proposed Value = 145: Voyager
        # Proposed Value = 150: CE1
        if var_scale == "log":
            ax_dict[var_axis + "ticks"] = np.log10(orig_labels)
            ax_dict[var_axis + "label"] = r"Laser Power [W]"
            ax_dict[var_axis + "ticklabels"] = [
                r"$10^{%.0f}$" % x if abs(int(x)) > 1 else r"$%.1f$" % (10**x)
                for x in np.log10(orig_labels)
            ]
        elif var_scale == "lin":
            ax_dict[var_axis + "ticks"] = orig_labels
            ax_dict[var_axis + "label"] = r"Laser Power [W]"
            ax_dict[var_axis + "ticklabels"] = [r"$%.1f$" % x for x in orig_labels]
    elif var == "Seismic Gamma":
        # Proposed Value = 0.8: aLIGO, Voyager, CE1
        if var_scale == "log":
            ax_dict[var_axis + "ticks"] = np.log10(orig_labels)
            ax_dict[var_axis + "label"] = r"Seismic Gamma"
            ax_dict[var_axis + "ticklabels"] = [
                r"$10^{%.0f}$" % y if abs(int(y)) > 1 else r"$%.1f$" % (10**y)
                for y in np.log10(orig_labels)
            ]
        elif var_scale == "lin":
            ax_dict[var_axis + "ticks"] = orig_labels
            ax_dict[var_axis + "label"] = r"Seismic Gamma"
            ax_dict[var_axis + "ticklabels"] = [r"$%.1f$" % y for y in orig_labels]
    elif var == "Materials Substrate Temp":
        # Proposed Value = 295: aLIGO, CE1
        # Proposed Value = 123: Voyager
        if var_scale == "lin":
            ax_dict[var_axis + "ticks"] = orig_labels
            ax_dict[var_axis + "label"] = r"Mirror Substrate Temp [K]"
            ax_dict[var_axis + "ticklabels"] = [
                r"$%.1f \times 10^{%i}$" % (x / 10 ** int(np.log10(x)), np.log10(x))
                if np.abs(int(np.log10(x))) > 1
                else "{:g}".format(x)
                for x in orig_labels
            ]
        elif var_scale == "log":
            ax_dict[var_axis + "ticks"] = np.log10(orig_labels)
            ax_dict[var_axis + "label"] = r"Mirror Substrate Temp [K]"
            ax_dict[var_axis + "ticklabels"] = [
                r"$10^{%.0f}$" % y if abs(int(y)) > 1 else r"$%.1f$" % (10**y)
                for y in np.log10(orig_labels)
            ]
    else:
        if var_scale == "lin":
            ax_dict[var_axis + "ticks"] = orig_labels
            ax_dict[var_axis + "label"] = str(var)
            ax_dict[var_axis + "ticklabels"] = [
                r"$%.1f \times 10^{%i}$" % (x / 10 ** int(np.log10(x)), np.log10(x))
                if np.abs(int(np.log10(x))) > 1
                else "{:g}".format(x)
                for x in orig_labels
            ]
        elif var_scale == "log":
            ax_dict[var_axis + "ticks"] = np.log10(orig_labels)
            ax_dict[var_axis + "label"] = str(var)
            ax_dict[var_axis + "ticklabels"] = [
                r"$10^{%.0f}$" % y if abs(int(y)) > 1 else r"$%.1f$" % (10**y)
                for y in np.log10(orig_labels)
            ]
    if line_val is not None:
        if "linestyle" not in line_kwargs.keys():
            line_kwargs["linestyle"] = "--"
        if "color" not in line_kwargs.keys():
            line_kwargs["color"] = "k"
        if "label" not in line_kwargs.keys():
            line_kwargs["label"] = "Proposed Value"

        if var_scale == "log":
            if var_axis == "y":
                ax.axhline(y=np.log10(line_val), **line_kwargs)
            elif var_axis == "x":
                ax.axvline(x=np.log10(line_val), **line_kwargs)
        elif var_scale == "lin":
            if var_axis == "y":
                ax.axhline(y=line_val, **line_kwargs)
            elif var_axis == "x":
                ax.axvline(x=line_val, **line_kwargs)

    ax.update(ax_dict)
    if label_kwargs:
        if var_axis == "y":
            ax.set_ylabel(ax.get_ylabel(), **label_kwargs)
        elif var_axis == "x":
            ax.set_xlabel(ax.get_xlabel(), **label_kwargs)

    if tick_label_kwargs:
        if var_axis == "y":
            ax.set_yticklabels(ax.get_yticklabels(), **tick_label_kwargs)
        elif var_axis == "x":
            ax.set_xticklabels(ax.get_xticklabels(), **tick_label_kwargs)

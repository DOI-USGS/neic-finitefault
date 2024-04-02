import pathlib
from typing import List, Literal, Optional, Union

import numpy as np
from matplotlib import pyplot as plt  # type: ignore
from matplotlib import ticker  # type: ignore
from obspy import read  # type: ignore
from scipy.signal import butter, filtfilt  # type: ignore

plt.rc("axes", titlesize=14)
plt.rc("axes", labelsize=12)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("font", size=12)


def plot_waveforms(
    axes: List[plt.Axes],
    times: List[Union[list, np.ndarray]],
    waveforms: List[Union[list, np.ndarray]],
    weights: List[float],
    type_str: Optional[str] = None,
    comp: Optional[str] = None,
    color: str = "blue",
    custom: Optional[Literal["fill", "syn"]] = None,
) -> List[plt.Axes]:
    """Plot timeseries waveform data

    :param axes: The axes to plot the data on
    :type axes: List[plt.Axes]
    :param times: The time
    :type times: List[np.ndarray]
    :param waveforms: The waveform values
    :type waveforms: List[np.ndarray]
    :param weights: The waveform weight
    :type weights: List[float]
    :param type_str: The data type, defaults to None
    :type type_str: Optional[str], optional
    :param comp: The component, defaults to None
    :type comp: Optional[str], optional
    :param color: The color of the plot line, defaults to "blue"
    :type color: str, optional
    :param custom: The custom option (fill), defaults to None
    :type custom: Optional[Literal[&quot;fill&quot;, &quot;syn&quot;]], optional
    :return: The updated axes
    :rtype: List[plt.Axes]
    """
    nPlot = 0
    for ax, time, waveform, weight in zip(axes, times, waveforms, weights):
        nPlot += 1
        if weight == 0.0:
            ax.plot(time, waveform, color=color, linewidth=2, linestyle="dashed")
        else:
            ax.plot(time, waveform, color=color, linewidth=2 * weight)
        min_time, max_time = ax.get_xlim()
        min_time = np.minimum(np.min(time), min_time)
        max_time = np.maximum(np.max(time), max_time)
        if custom == "fill":
            min_val, max_val = ax.get_ylim()
            min_val = np.minimum(np.min(waveform), min_val)
            max_val = np.maximum(np.max(waveform), max_val)
            ax.set_ylim((-(max(abs(min_val), max_val)), max(abs(min_val), max_val)))
            min_val, max_val = ax.get_ylim()
            ax.vlines(0, min_val, max_val, "k", lw=1)
            if type_str == "body":
                ax.text(
                    np.max(time),
                    0.6 * max_val,
                    "{:0.1f}".format(max(abs(min_val), max_val)),
                    ha="right",
                    va="center",
                )
                ax.hlines(0, -20, np.max(time), "k", lw=1)
                ax.set_xlim((-20, np.max(time)))
            elif type_str == "surf":
                ax.text(
                    np.max(time),
                    0.6 * max_val,
                    "{:0.3f}".format(max(abs(min_val), max_val)),
                    ha="right",
                    va="center",
                )
                ax.hlines(0, -350, np.max(time), "k", lw=1)
                ax.set_xlim((-350, np.max(time)))
            elif type_str == "cgps" or type_str == "strong":
                min_wval = np.min(waveform)
                max_wval = np.max(waveform)
                if max_wval > abs(min_wval):
                    ax.text(
                        np.max(time),
                        0.6 * max_val,
                        "{:0.2f}".format(max_wval),
                        ha="right",
                        va="center",
                    )
                else:
                    ax.text(
                        np.max(time),
                        0.6 * max_val,
                        "{:0.2f}".format(min_wval),
                        ha="right",
                        va="center",
                    )
                ax.hlines(0, -15, np.max(time), "k", lw=1)
                ax.set_xlim((-15, np.max(time)))
            min_time, max_time = ax.get_xlim()
            if type_str == "body" and comp == "BHZ":
                ax.text(
                    1.1 * min_time,
                    0.2 * max(abs(min_val), max_val),
                    "P",
                    ha="right",
                    va="bottom",
                )
            if type_str == "body" and comp == "SH":
                ax.text(
                    1.1 * min_time,
                    0.2 * max(abs(min_val), max_val),
                    "SH",
                    ha="right",
                    va="bottom",
                )
            if type_str == "surf" and comp == "BHZ":
                ax.text(
                    1.2 * min_time,
                    0.2 * max(abs(min_val), max_val),
                    "Z",
                    ha="right",
                    va="bottom",
                )
            if type_str == "surf" and comp == "SH":
                ax.text(
                    1.2 * min_time,
                    0.2 * max(abs(min_val), max_val),
                    "T",
                    ha="right",
                    va="bottom",
                )
        if custom == "syn":
            max_val = np.maximum(abs(min(waveform)), max(waveform))
            tmin, tmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            if type_str == "body":
                ax.text(
                    tmax,
                    0.6 * ymin,
                    "{:0.1f}".format(max_val),
                    ha="right",
                    va="center",
                    color="red",
                )
            elif type_str == "surf":
                ax.text(
                    tmax,
                    0.6 * ymin,
                    "{:0.3f}".format(max_val),
                    ha="right",
                    va="center",
                    color="red",
                )
            elif type_str == "cgps" or type_str == "strong":
                ax.text(
                    tmax,
                    0.6 * ymin,
                    "{:0.2f}".format(max_val),
                    ha="right",
                    va="center",
                    color="red",
                )
        if type_str == "body":
            ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5, min_n_ticks=5))
            ax.yaxis.set_major_locator(ticker.NullLocator())
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
        elif type_str == "surf":
            ax.xaxis.set_major_locator(ticker.MultipleLocator(500))
            ax.yaxis.set_major_locator(ticker.NullLocator())
        elif type_str == "cgps" or type_str == "strong":
            ax.xaxis.set_major_locator(ticker.MultipleLocator(40))
            ax.yaxis.get_major_locator().set_params(integer=True)  # type:ignore
        if nPlot > len(weights) - 3:
            ax.set_xlabel("Time After OT (s)")
        ax.grid(axis="x", which="both", linestyle="dotted", color="0.5")
    return axes


def add_metadata(
    axes: List[plt.Axes],
    azimuths: Optional[List[float]] = None,
    comps: Optional[List[str]] = None,
    distances: Optional[List[float]] = None,
    names: Optional[List[str]] = None,
    type_str: Optional[List[float]] = None,
    weights: Optional[List[float]] = None,
) -> List[plt.Axes]:
    """Add metadata to axes

    :param axes: The axes to add the metadata to
    :type axes: List[plt.Axes]
    :param azimuths: A list of azimuth values, defaults to None
    :type azimuths: Optional[List[float]], optional
    :param comps: A list of components, defaults to None
    :type comps: Optional[List[float]], optional
    :param distances: A list of distance values, defaults to None
    :type distances: Optional[List[float]], optional
    :param names: A list of names, defaults to None
    :type names: Optional[List[str]], optional
    :param type_str: A list of data types, defaults to None
    :type type_str: Optional[List[float]], optional
    :param weights: A list of weight values, defaults to None
    :type weights: Optional[List[float]], optional
    :return: The updated axes
    :rtype: List[plt.Axes]
    """
    if type_str is not None:
        if type_str == "cgps" or type_str == "strong":
            if distances is not None:
                for ax, dist in zip(axes, distances):
                    ax.text(
                        0.01,
                        0.46,
                        "{:0.2f}".format(dist),
                        ha="left",
                        va="top",
                        transform=ax.transAxes,
                    )
            if comps is not None:
                if names is not None:
                    for ax, comp, name in zip(axes, comps, names):
                        if comp[-1] == "E":
                            ax.text(
                                -0.18,
                                0.5,
                                name,
                                ha="right",
                                va="center",
                                transform=ax.transAxes,
                                rotation=90,
                                fontweight="bold",
                            )
                for ax, comp in zip(axes, comps):
                    ax.text(
                        -0.13,
                        0.5,
                        comp,
                        ha="right",
                        va="center",
                        transform=ax.transAxes,
                        rotation=90,
                    )
        else:
            if names is not None:
                for ax, name in zip(axes, names):
                    ax.text(
                        -0.02,
                        0.50,
                        name,
                        ha="right",
                        va="center",
                        transform=ax.transAxes,
                    )
            if distances is not None:
                for ax, dist in zip(axes, distances):
                    ax.text(
                        0.01,
                        0.46,
                        "{:0.0f}".format(dist),
                        ha="left",
                        va="top",
                        transform=ax.transAxes,
                    )
    if azimuths is not None:
        for ax, az in zip(axes, azimuths):
            ax.text(
                0.01,
                0.5,
                "{:0.0f}".format(az),
                ha="left",
                va="bottom",
                transform=ax.transAxes,
            )
    if weights is not None:
        for ax, weight in zip(axes, weights):
            alpha = 1 if weight > 0 else 0.25
            lines = ax.get_lines()
            for line in lines:
                line.set_alpha(alpha)
    return axes


def plot_waveform_fits(
    files: List[dict],
    components: list,
    type_str: str,
    start_margin: int = 10,
    plot_directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot fit of observed to synthetic data for selected channels

    :param files: Waveform files to plot
    :type files: List[dict]
    :param components: Components (channels) of data selected for plotting
    :type components: list
    :param type_str: Data type of given waveform files
    :type type_str: str
    :param start_margin: Start margin of data for plotting, defaults to 10
    :type start_margin: int, optional
    :param plot_directory: Directory where plots should be saved, defaults to pathlib.Path()
    :type plot_directory: Union[pathlib.Path, str], optional
    """
    plot_directory = pathlib.Path(plot_directory)
    if type_str == "body" or type_str == "surf":
        files = [file for file in files if file["component"] in components]
        print("Creating Waveform Fit Plot: " + str(type_str) + " " + str(components[0]))
    if type_str == "cgps" or type_str == "strong":
        files = [file for file in files]
        print("Creating Waveform Fit Plot: " + str(type_str))
    if not len(files):
        print(f"No files found matching components ({components}). "
              "Not plotting waveform fits.")
        return
    files = sorted(files, key=lambda k: (k["azimuth"], k["component"]))
    sampling = [file["dt"] for file in files]
    names = [file["name"] for file in files]
    azimuths = [file["azimuth"] for file in files]
    distances = [file["distance"] for file in files]
    weights = [file["trace_weight"] for file in files]
    obs_waveforms = [file["observed"] for file in files]
    syn_waveforms = [file["synthetic"] for file in files]
    comp = [file["component"] for file in files]
    zipped = zip(sampling, syn_waveforms)
    syn_times = [dt * np.arange(0, len(synthetic)) for dt, synthetic in zipped]
    start_waveform: List[float] = []
    for file in files:
        dt = file["dt"]
        nstart = file["start_signal"]
        margin = int(start_margin / dt)
        margin = min(nstart, margin)
        start_waveform = start_waveform + [margin]
    obs_times = [
        dt * np.arange(-start, len(observed) - start)
        for dt, start, observed in zip(sampling, start_waveform, obs_waveforms)
    ]
    numrows_phase = len(files) // 3 + 1
    fig, axes = plt.subplots(
        max(3, numrows_phase), 3, figsize=(20, int(2.2 * numrows_phase))
    )
    axes2 = axes.ravel()
    for ax in axes2[len(files) :]:
        ax.axis("off")

    if type_str == "body" or type_str == "surf":
        axes2 = plot_waveforms(
            axes2,
            obs_times,
            obs_waveforms,
            weights,
            type_str=type_str,
            comp=comp[0],
            color="black",
            custom="fill",
        )
        axes2 = plot_waveforms(
            axes2,
            syn_times,
            syn_waveforms,
            weights,
            type_str=type_str,
            comp=comp[0],
            color="red",
            custom="syn",
        )
    if type_str == "cgps" or type_str == "strong":
        axes2 = plot_waveforms(
            axes2,
            obs_times,
            obs_waveforms,
            weights,
            type_str=type_str,
            comp=comp[0],
            color="black",
            custom="fill",
        )
        axes2 = plot_waveforms(
            axes2,
            syn_times,
            syn_waveforms,
            weights,
            type_str=type_str,
            comp=comp[0],
            color="red",
            custom="syn",
        )

    dict = {
        "weights": weights,
        "azimuths": azimuths,
        "names": names,
        "distances": distances,
        "type_str": type_str,
        "comps": comp[0],
    }
    axes2 = add_metadata(axes2, **dict)

    if type_str == "body":
        if "BHZ" in components:
            plot_name = "P_body_waves"
        if "SH" in components:
            plot_name = "SH_body_waves"

    if type_str == "surf":
        if "BHZ" in components:
            plot_name = "Rayleigh_surf_waves"
        if "SH" in components:
            plot_name = "Love_surf_waves"

    if type_str == "cgps":
        plot_name = "cGPS_waves"

    if type_str == "strong":
        plot_name = "strong_motion_waves"

    plt.savefig(plot_directory / (plot_name + ".png"), dpi=300)  # bbox_inches='tight')
    plt.savefig(plot_directory / (plot_name + ".ps"))
    plt.close()
    return


def filt_waveform(file: dict, high_freq: float) -> dict:
    """Filter the waveform

    :param file: The file information
    :type file: dict
    :param high_freq: The high frequency to filter
    :type high_freq: float
    :return: The updated file information
    :rtype: dict
    """
    stream = read(file["file"])
    data = stream[0][1000:]
    high_freq = high_freq / 25
    b, a = butter(2, high_freq, btype="lowpass")
    filt_data = filtfilt(b, a, data)
    file["synthetic"] = filt_data
    return file


def plot_spectra(
    files: List[dict],
    dt: float,
    plot_directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot the spectra of the waveforms

    :param files: The file information dictionaries
    :type files: List[dict]
    :param dt: The dt
    :type dt: float
    :param plot_directory: The directory where plots are saved, defaults to pathlib.Path()
    :type plot_directory: Union[pathlib.Path, str], optional
    """
    plot_directory = pathlib.Path(plot_directory)
    for file in files:
        waveform = file["synthetic"]
        name = file["name"]
        comp = file["component"]
        fft = np.fft.fft(waveform)
        n = len(waveform)
        freq = np.fft.fftfreq(n, d=dt)
        plt.loglog(freq[: n // 2], np.abs(fft[: n // 2]))
        plt.title("{} {}".format(name, comp))
        plt.savefig(plot_directory / "spectra_{}_{}".format(name, comp))
        plt.close()

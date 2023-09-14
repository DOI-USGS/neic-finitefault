import pathlib
from typing import List, Literal, Optional, Union

import numpy as np
from matplotlib import pyplot as plt  # type: ignore
from matplotlib import ticker

from wasp.waveform_plots_NEIC import filt_waveform, plot_spectra  # type: ignore


def plot_waveforms(
    axes: List[plt.Axes],
    times: List[Union[list, np.ndarray]],
    waveforms: List[Union[list, np.ndarray]],
    color: str = "blue",
    custom: Optional[Literal["fill"]] = None,
) -> List[plt.Axes]:
    """Plot timeseries waveform data

    :param axes: The axes to plot the data on
    :type axes: List[plt.Axes]
    :param times: The time
    :type times: List[np.ndarray]
    :param waveforms: The waveform values
    :type waveforms: List[np.ndarray]
    :param color: The color of the plot line, defaults to "blue"
    :type color: str, optional
    :param custom: The custom option (fill), defaults to None
    :type custom: Optional[Literal[&quot;fill&quot;]], optional
    :return: The updated axes
    :rtype: List[plt.Axes]
    """
    for ax, time, waveform in zip(axes, times, waveforms):
        if waveform is None:
            continue
        if len(waveform) == 0:
            continue
        ax.plot(time, waveform, color=color, linewidth=0.6)
        min_time, max_time = ax.get_xlim()
        min_time = np.minimum(np.min(time), min_time)
        max_time = np.maximum(np.max(time), max_time)
        ax.set_xlim([min_time, max_time])
        if custom == "fill":
            min_val, max_val = ax.get_ylim()
            min_val = np.minimum(np.min(waveform), min_val)
            max_val = np.maximum(np.max(waveform), max_val)
            ax.set_ylim([min_val, max_val])
            ax.vlines(0, min_val, max_val)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3, min_n_ticks=3))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, min_n_ticks=3))
    return axes


def add_metadata(
    axes: List[plt.Axes],
    azimuths: Optional[List[float]] = None,
    distances: Optional[List[float]] = None,
    names: Optional[List[str]] = None,
    weights: Optional[List[float]] = None,
) -> plt.Axes:
    """Add metadata to axes

    :param axes: The axes to add the metadata to
    :type axes: List[plt.Axes]
    :param azimuths: A list of azimuth values, defaults to None
    :type azimuths: Optional[List[float]], optional
    :param distances: A list of distance values, defaults to None
    :type distances: Optional[List[float]], optional
    :param names: A list of names, defaults to None
    :type names: Optional[List[str]], optional
    :param weights: A list of weight values, defaults to None
    :type weights: Optional[List[float]], optional
    :return: The updated axes
    :rtype: plt.Axes
    """
    if names is not None:
        for ax, name in zip(axes, names):
            if name is None:
                continue
            ax.text(0.9, 0.9, name, ha="center", va="center", transform=ax.transAxes)
    if distances is not None:
        for ax, dist in zip(axes, distances):
            if dist is None:
                continue
            ax.text(
                0.1,
                0.1,
                "{:0.1f}".format(dist),
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
    if azimuths is not None:
        for ax, az in zip(axes, azimuths):
            if az is None:
                continue
            ax.text(
                0.1,
                0.9,
                "{:0.1f}".format(az),
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
    if weights is not None:
        for ax, weight in zip(axes, weights):
            if weight is None:
                continue
            alpha = 1 if weight > 0 else 0.1
            lines = ax.get_lines()
            for line in lines:
                line.set_alpha(alpha)
    return axes


def plot_waveform_fits(
    files: List[dict],
    components: List[str],
    type_str: str,
    start_margin: int = 10,
    event: Optional[int] = None,
    plot_directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot fit of observed to synthetic data for selected channels.

    :param files: waveform files to plot
    :param components: components (channels) of data selected for plotting
    :param type_str: data type of given waveform files
    :param start_margin: start margin of data for plotting
    :type files: list
    :type components: list
    :type type_str: string
    :type start_margin: float, optinoal
    """
    plot_directory = pathlib.Path(plot_directory)
    files = [file for file in files if file["component"] in components]
    files = sorted(files, key=lambda k: k["azimuth"])
    sampling = [file["dt"] for file in files]
    names = [file["name"] for file in files]
    azimuths = [file["azimuth"] for file in files]
    distances = [file["distance"] for file in files]
    weights = [file["trace_weight"] for file in files]
    obs_waveforms = [file["observed"] for file in files]
    syn_waveforms = [file["synthetic"] for file in files]
    zipped = zip(obs_waveforms, syn_waveforms)
    syn_waveforms = [
        syn_waveform[: len(obs_waveform)] for obs_waveform, syn_waveform in zipped
    ]
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
    ]  # type: ignore
    numrows_phase = len(files) // 4 + 1
    fig, axes = plt.subplots(max(4, numrows_phase), 4, figsize=(13, 9))
    axes2 = axes.ravel()
    for ax in axes2[len(files) :]:
        ax.axis("off")
    obs_waveforms = [waveform for waveform in obs_waveforms]
    syn_waveforms = [waveform for waveform in syn_waveforms]

    axes2 = plot_waveforms(axes2, obs_times, obs_waveforms, color="black")
    axes2 = plot_waveforms(axes2, syn_times, syn_waveforms, color="red", custom="fill")
    dict = {
        "weights": weights,
        "azimuths": azimuths,
        "names": names,
        "distances": distances,
    }
    axes2 = add_metadata(axes2, **dict)

    if type_str == "cgps":
        if "LXZ" in components:
            plot_name = "LXZ_cgps_waves"
        if "LXN" in components:
            plot_name = "LXN_cgps_waves"
        if "LXE" in components:
            plot_name = "LXE_cgps_waves"

    if type_str == "strong_motion":
        if "HNZ" in components:
            plot_name = "HNZ_strong_motion_waves"
        if "HNN" in components:
            plot_name = "HNN_strong_motion_waves"
        if "HNE" in components:
            plot_name = "HNE_strong_motion_waves"

    if type_str == "tele_body":
        if "BHZ" in components:
            plot_name = "P_body_waves"
        if "SH" in components:
            plot_name = "SH_body_waves"

    if type_str == "surf_tele":
        if "BHZ" in components:
            plot_name = "Rayleigh_surf_waves"
        if "SH" in components:
            plot_name = "Love_surf_waves"

    if event is not None:
        plot_name = "{}_event{}".format(plot_name, event)
    plt.savefig(plot_directory / plot_name, bbox_inches="tight")
    plt.close()
    return


if __name__ == "__main__":
    files = [
        {
            "file": "/home/pkoch/folder_plot16/STR.MT07.HNE.C1.ACC",
            "name": "MT07",
            "component": "HNE",
            "synthetic": [],
        },
        {
            "file": "/home/pkoch/folder_plot16/STR.MT07.HNN.C1.ACC",
            "name": "MT07",
            "component": "HNN",
            "synthetic": [],
        },
        {
            "file": "/home/pkoch/folder_plot16/STR.MT07.HNZ.C1.ACC",
            "name": "MT07",
            "component": "HNZ",
            "synthetic": [],
        },
    ]
    high_freq = 15
    for file in files:
        file = filt_waveform(file, high_freq)
    plot_spectra(files, 0.01)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    waveforms = [file["synthetic"] for file in files]
    times = [np.arange(len(waveform)) * 0.01 for waveform in waveforms]
    axes = plot_waveforms(axes, times, waveforms)  # type: ignore
    axes[0].set_title("N")
    axes[1].set_title("E")
    axes[2].set_title("Z")
    plot_name = "MT07_lowpass_{}".format(high_freq)
    plt.savefig(plot_name, bbox_inches="tight")
    print(1)

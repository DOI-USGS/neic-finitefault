# -*- coding: utf-8 -*-
"""
Script for implementing the baseline removal scheme proposed by
Wang et al. [2011]_

.. [2011] Wang, R., B. Schurr, C. Milkereit, Z. Shao, and M. Jin (2011),
     An improved automatic scheme for empirical baseline correction of digital
     strong motion records, Bull. Seism. Soc. Am.,101(5), 2029–2044,
     doi:10.1785/0120110039.
"""


import argparse
import glob
import os
import pathlib
import time
from typing import Optional, Tuple, Union

import matplotlib.pyplot as plt  # type:ignore
import numpy as np
import scipy.integrate as integrate  # type:ignore
from obspy import read  # type:ignore
from obspy.core.trace import Trace  # type:ignore
from scipy.optimize import shgo  # type:ignore


def wang_process(
    file: Union[pathlib.Path, str],
    plot: bool = False,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> Trace:
    """Remove baselines from strong motion data and optionally plot the results

    :param file: File with strong motion data in sac format
    :type file: Union[pathlib.Path, str]
    :param plot: Whether to plot the results of the baseline removal, defaults to False
    :type plot: bool, optional
    :param directory: Where to save plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The updated strong motion data
    :rtype: Trace

    .. note::

       This has only been tested with waveform files in SAC format. We can't
       ensure this works with other formats.

    .. note::

       The waveform should have some data prior to the start of shaking.
       We have found it is best if there is data up to at least 20 seconds
       prior to origin time
    """
    directory = pathlib.Path(directory)
    trace, vel_trace, disp_trace, constants = _wang_baselines(file)
    if not constants:
        return None
    a0, a1, a2, a3, t1, t2, tp, t_d0, t_pgd, t_pga, t_end = constants
    disp_data = disp_trace.data
    delta = disp_trace.stats.delta
    disp_data2, gps, t_jump = _gps_jump(disp_data, t_end, t1, t2, a1, a2, a3, tp, delta)
    disp_trace2 = disp_trace.copy()
    disp_trace2.data = disp_data2
    vel_trace2 = disp_trace2.copy()
    vel_trace2.differentiate()
    if plot:
        constants = [a0, a1, a2, a3, t1, t2, tp, t_d0, t_pgd, t_pga, t_end, gps, t_jump]
        _optional_plots(vel_trace, disp_trace, constants, directory=directory)
        _opt_plots2(trace, vel_trace, disp_trace, tp, directory=directory)
    return vel_trace2


def _wang_baselines(
    file: Union[pathlib.Path, str]
) -> Tuple[Trace, Trace, Trace, Optional[list]]:
    """Perform baseline removal

    :param file: File with strong motion data in sac format
    :type file: Union[pathlib.Path, str]
    :return: _description_
    :rtype: Tuple[Trace, Trace, Trace, Optional[list]]
    """
    st = read(file)
    st1 = st.copy()
    trace = st1[0]
    if np.max(np.abs(trace.data)) < 10**-8:
        return trace, trace, trace, None
    delta = trace.stats.delta
    #
    # we define some relevant times
    #
    time_p = _first_arrival(trace)
    time_end, time_pga = _time_end_pga(trace)
    if time_p >= min(time_end, time_pga):
        return trace, trace, trace, None
    preevent_baseline = np.mean(trace.data[:time_p])  # -  int(1 / delta)])
    trace.data = trace.data - preevent_baseline
    trace.data = trace.data[: time_p + 4 * (time_end - time_p)]
    #
    # we integrate to velocity and displacement and remove quadratic trend at the
    # end of the signal.
    #
    new_trace = trace.copy()
    vel_trace = new_trace.integrate(method="cumtrapz")
    vel_trace.data = vel_trace.data - vel_trace.data[time_p]  # - 1 * int(1/delta)]
    new_trace = vel_trace.copy()
    disp_trace = new_trace.integrate(method="cumtrapz")
    disp_trace.data = disp_trace.data - disp_trace.data[time_p]  # - 1 * int(1 / delta)]
    time_d0, time_pgd = _time_d0_pgd(disp_trace, time_p, time_pga, time_end)
    disp_tail = disp_trace.data[time_end:]
    print("TAIL", disp_tail)
    if len(disp_tail) < 5:
        return trace, trace, trace, None
    a0, a1, a2, aic3 = _polinomial_fit(disp_tail, delta, 2)
    b0, b1, b2, b3, aic4 = _polinomial_fit(disp_tail, delta, 3)
    a3 = 0
    if aic4 < aic3:
        a0, a1, a2, a3 = [b0, b1, b2, b3]
    #
    # we search for a baseline correction such that corrected displacement best
    # resembles a heaviside function.
    #
    disp_data = disp_trace.data
    max_val = max(time_pga, time_d0)
    constraint0 = {"type": "ineq", "fun": _constraint}
    bounds = [(time_pgd * delta, time_end * delta), (max_val * delta, time_end * delta)]
    if time_pgd > time_end or max_val > time_end:
        return trace, trace, trace, None
    args = (disp_data, a1, a2, a3, time_p, time_end, delta)

    res = shgo(_function, bounds, iters=6, constraints=constraint0, args=args)
    time1, time2 = res.x
    success = res.success
    if not success:
        print(file)
        print(res)
    time1 = int(time1 / delta)
    time2 = int(time2 / delta)
    constants = [
        a0,
        a1,
        a2,
        a3,
        time1,
        time2,
        time_p,
        time_d0,
        time_pgd,
        time_pga,
        time_end,
    ]
    return trace, vel_trace, disp_trace, constants


def _function(
    time: list,
    disp_data: np.ndarray,
    a1: float,
    a2: float,
    a3: float,
    time_p: int,
    time_end: int,
    delta: float,
) -> float:
    """Compute misfit of corrected displacement to a heaviside

    :param time: List of times
    :type time: list
    :param disp_data: The displacement data
    :type disp_data: np.ndarray
    :param a1: The first polynomial of the fit
    :type a1: float
    :param a2: The second polynomial of the fit
    :type a2: float
    :param a3: The third polynomial of the fit
    :type a3: float
    :param time_p: The arrival time
    :type time_p: int
    :param time_end: The end time
    :type time_end: int
    :param delta: The data delta
    :type delta: float
    :return: The minimum value of the heaviside fit
    :rtype: float
    """
    time1, time2 = time
    time1 = int(time1 / delta)
    time2 = int(time2 / delta)
    length = len(disp_data)
    disp_corr = _correct_disp(length, delta, a1, a2, a3, time1, time2, time_end)
    disp_data2 = disp_data - disp_corr
    gps = np.mean(disp_data2[time_end:])
    comparison = disp_data2 - gps / 2
    after = comparison[1:]
    before = comparison[:-1]
    crossings = np.nonzero(comparison == 0)[0]
    crossings = list(crossings) + list(np.nonzero(before * after < 0)[0])  # type:ignore
    crossings = [index for index in crossings if index > time_p]  # type:ignore
    if not crossings:
        return np.sum(disp_data2**2)
    results = [_heaviside_fit(int(index), disp_data2, gps) for index in crossings]
    return min(results)


def _gps_jump(
    disp_data: np.ndarray,
    time_end: int,
    time1: int,
    time2: int,
    a1: float,
    a2: float,
    a3: float,
    time_p: int,
    delta: float,
) -> Tuple[np.ndarray, float, int]:
    """Get time and amplitude of heaviside jump

    :param disp_data: The displacement data
    :type disp_data: np.ndarray
    :param time_end: The end time
    :type time_end: int
    :param time1: Time 1 of the optimized result
    :type time1: int
    :param time2: Time 2 of the optimized result
    :type time2: int
    :param a1: The first polynomial of the fit
    :type a1: float
    :param a2: The second polynomial of the fit
    :type a2: float
    :param a3: The third polynomial of the fit
    :type a3: float
    :param time_p: The arrival time
    :type time_p: int
    :param delta: The data's delta
    :type delta: float
    :return: The updated displacement data, the mean of the updated displacement data,
            and the time of the jump
    :rtype: Tuple[np.ndarray, float, int]
    """
    length = len(disp_data)
    disp_corr = _correct_disp(length, delta, a1, a2, a3, time1, time2, time_end)
    disp_data2 = disp_data - disp_corr
    gps = np.mean(disp_data2[time_end:])
    comparison = disp_data2 - gps / 2
    after = comparison[1:]
    before = comparison[:-1]
    crossings = np.nonzero(comparison == 0)[0]
    crossings = list(crossings) + list(np.nonzero(before * after < 0)[0])  # type:ignore
    crossings = [index for index in crossings if index > time_p]  # type:ignore
    if not crossings:
        return np.sum(disp_data**2), gps, 0
    results = [_heaviside_fit(int(index), disp_data2, gps) for index in crossings]
    best_result = min(results)
    zipped = zip(crossings, results)
    time_jump = next(index for index, result in zipped if result == best_result)
    return disp_data2, gps, time_jump


def _constraint(x: np.ndarray) -> float:
    """Simple subtraction method

    :param x: The array
    :type x: np.ndarray
    :return: The difference of the first and second value
    :rtype: float
    """
    return x[1] - x[0]


def _integral_arias(trace: Trace, delta: float) -> np.ndarray:
    """Compute `\int_0^T \ddot{x}^2(t)dt`

    :param trace: The strong motion trace
    :type trace: Trace
    :param delta: The data's delta
    :type delta: float
    :return: The arias integral
    :rtype: np.ndarray
    """
    trace = trace**2
    integral = integrate.cumtrapz(trace, dx=delta, initial=0)
    return integral


def _first_arrival(trace: Trace) -> int:
    """A binary search hand-made picker

    :param trace: The strong motion trace
    :type trace: Trace
    :return: The arrival time
    :rtype: int
    """
    accel = trace.data
    delta = trace.stats.delta
    trace0 = accel[0 : int(8 / delta)]
    val0 = max(_integral_arias(trace0, delta))
    for k in range(8, int(len(accel) / delta), 8):
        trace1 = accel[int(k / delta) : int((k + 8) / delta)]
        if len(trace1) == 0:
            continue
        val = max(_integral_arias(trace1, delta))
        if val > 10 * val0:
            break
    start0 = k
    trace0 = accel[int((start0 - 12) / delta) : int((start0 - 8) / delta)]
    if len(trace0) == 0:
        return int(start0 / delta)
    val0 = max(_integral_arias(trace0, delta))
    for k in range(start0 - 8, start0 + 8, 4):
        trace1 = accel[int(k / delta) : int((k + 4) / delta)]
        val = max(_integral_arias(trace1, delta))
        if val > 10 * val0:
            break
    start0 = k
    trace0 = accel[int((start0 - 6) / delta) : int((start0 - 4) / delta)]
    val0 = max(_integral_arias(trace0, delta))
    for k in range(start0 - 4, start0 + 4, 2):
        trace1 = accel[int(k / delta) : int((k + 2) / delta)]
        val = max(_integral_arias(trace1, delta))
        if val > 10 * val0:
            break
    start0 = k
    trace0 = accel[int((start0 - 3) / delta) : int((start0 - 2) / delta)]
    val0 = max(_integral_arias(trace0, delta))
    for k in range(start0 - 2, start0 + 2, 1):
        trace1 = accel[int(k / delta) : int((k + 1) / delta)]
        val = max(_integral_arias(trace1, delta))
        if val > 10 * val0:
            break
    start0 = k  # - 1
    return int(start0 / delta)


def _time_end_pga(trace: Trace) -> Tuple[int, int]:
    """Get the end and pga times

    :param trace: The strong motion trace
    :type trace: Trace
    :return: The end time and pga time
    :rtype: Tuple[int,int]
    """
    accel = trace.data
    delta = trace.stats.delta
    pga = np.max(np.abs(accel))
    time_pga = next(index for index, v in enumerate(accel) if abs(v) == pga)
    integral = _integral_arias(accel, delta)
    integral = integral / integral[-1]
    time_end = next(index for index, v in enumerate(integral) if v > 0.95)
    return time_end, time_pga


def _time_d0_pgd(
    trace: Trace, time_p: int, time_pga: int, time_end: int
) -> Tuple[int, int]:
    """Find the last zero crossing of displacement prior to the estimated trace
        end, for uncorrected displacement. Then Find the last local optimum of
        displacement prior to this instant.

    :param trace: The strong motion trace
    :type trace: Trace
    :param time_p: The arrival time
    :type time_p: int
    :param time_pga: The pga time
    :type time_pga: int
    :param time_end: The end time
    :type time_end: int
    :return: The zero crossing and local optimum times
    :rtype: Tuple[int,int]
    """
    zero_crossings = [
        index
        for index in range(time_p, time_end)
        if _iz_zero_crossing(trace.data, index)
    ]
    if zero_crossings:
        time_d0 = zero_crossings[-1]
        pgd = np.max(np.abs(trace.data[time_p:time_d0]))
        time_pgd = next(
            index for index in range(time_d0) if abs(trace.data[index]) == pgd
        )
    else:
        #
        # In case we can't find these values, we search them in some value previous to
        # the PGA time
        #
        vel_trace = np.gradient(trace, trace.stats.delta)
        pgv = np.max(np.abs(vel_trace.data[time_p:time_pga]))
        time_pgv = next(
            index for index in range(time_pga) if abs(vel_trace.data[index]) == pgv
        )
        time_d0 = time_pgv
        zero_crossings = [
            index
            for index in range(time_p, time_pgv)
            if _iz_zero_crossing(vel_trace.data, index)
        ]
        #
        # In case we still can't find a t_pgd estimation, we use a gross estimation.
        #
        time_pgd = (
            zero_crossings[-1] if zero_crossings else int((time_p + time_pga) / 2)
        )
    return time_d0, time_pgd


def _iz_zero_crossing(disp: np.ndarray, index: int) -> bool:
    """Find the zero crossing of signals

    :param disp: The displacement data
    :type disp: np.ndarray
    :param index: The index to evaluate
    :type index: int
    :return: Whether the index is at a zero crossing
    :rtype: bool
    """
    if disp[index] * disp[index - 1] > 0:
        return False
    elif disp[index] * disp[index - 1] < 0:
        return True
    elif disp[index] == 0:
        if np.max(np.abs(disp[:index])) == 0:
            return False
        if np.max(np.abs(disp[index:])) == 0:
            return False
        nonzero0 = next(i for i in range(index + 1) if abs(disp[index - i]) > 0)
        val0 = disp[nonzero0]
        nonzero1 = next(i for i in range(len(disp) - index) if abs(disp[index + i]) > 0)
        val1 = disp[nonzero1]
        if val1 * val0 < 0:
            return True
    return False


def _polinomial_fit(disp_tail: np.ndarray, delta: int, n: int) -> np.ndarray:
    """Find best fitting polinomial of degree n, and Akaike criteria of
    polinomial model

    :param disp_tail: The tail of displacement
    :type disp_tail: np.ndarray
    :param delta: The data's delta
    :type delta: int
    :param n: Power
    :type n: int
    :return: The polynomial fit values
    :rtype: np.ndarray
    """
    time = np.arange(len(disp_tail)) * delta
    matrix = np.fliplr(np.vander(time, n + 1))
    sym_mat = np.dot(np.transpose(matrix), matrix)
    inv_mat = np.linalg.inv(sym_mat)
    full_mat = np.dot(inv_mat, np.transpose(matrix))
    coeffs = np.dot(full_mat, disp_tail)
    polin = np.zeros(len(disp_tail))
    for i, coeff in enumerate(coeffs):
        polin = polin + coeff * time**i
    log_likelihood = -np.sum((disp_tail - polin) ** 2) * 5000
    aic = _akaike_criteria(len(disp_tail), n + 1, log_likelihood)
    return np.concatenate((np.dot(full_mat, disp_tail), [aic]))


def _akaike_criteria(nsize: int, nparam: int, log_likelihood: float) -> float:
    """The Akaike criteria

    :param nsize: The length of the data
    :type nsize: int
    :param nparam: The power
    :type nparam: int
    :param log_likelihood: The log likelihood
    :type log_likelihood: float
    :return: The Akaika criteria value
    :rtype: float
    """
    if nsize > nparam + 1:
        return (
            2 * nparam
            - 2 * log_likelihood
            + 2 * nparam * (1 + nparam) / (nsize - nparam - 1)
        )
    else:
        return 2 * nparam - 2 * log_likelihood


def _heaviside_fit(t_jump: int, disp: np.ndarray, static: float) -> float:
    """Fit to heaviside

    :param t_jump: The time of the jump
    :type t_jump: int
    :param disp: The displacement data
    :type disp: np.ndarray
    :param static: The mean of the tail of the displacement data
    :type static: float
    :return: The sum of the error
    :rtype: float
    """
    length = len(disp)
    step = static * np.ones(length)
    step[:t_jump] = np.zeros(t_jump)
    error = (disp - step) ** 2
    return np.sum(error)


def _correct_disp(
    length: int,
    delta: float,
    a1: float,
    a2: float,
    a3: float,
    t1: int,
    t2: int,
    t_end: int,
) -> np.ndarray:
    """Compute estimated baselines of displacement signal

    :param length: The length of the data
    :type length: int
    :param delta: The data's delta
    :type delta: float
    :param a1: The first polynomial of the fit
    :type a1: float
    :param a2: The second polynomial of the fit
    :type a2: float
    :param a3: The third polynomial of the fit
    :type a3: float
    :param t1: The first time
    :type t1: int
    :param t2: The second time
    :type t2: int
    :param t_end: The end time
    :type t_end: int
    :return: Estimated baselines
    :rtype: np.ndarray
    """
    new_time = np.arange(length) * delta
    disp_corr = np.zeros(length)
    p0 = a1 + 2 * a2 * (t2 - t_end) * delta + 3 * a3 * (t2 - t_end) ** 2 * delta**2
    if t1 < t2:
        disp_corr[t1:t2] = (p0 / (t2 - t1) / delta / 2) * (
            new_time[t1:t2] - t1 * delta
        ) ** 2
    disp_corr[t2:] = (
        p0 * (t2 - t1) * delta / 2
        + p0 * (new_time[t2:] - t2 * delta)
        + (a2 + 3 * a3 * (t2 - t_end) * delta) * (new_time[t2:] - t2 * delta) ** 2
        + a3 * (new_time[t2:] - t2 * delta) ** 3
    )
    return disp_corr


def _opt_plots2(
    trace: Trace,
    vel_trace: Trace,
    disp_trace: Trace,
    tp: int,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plots to show effect of this baseline removal procedure on selected
    strong motion data

    :param trace: The trace data
    :type trace: Trace
    :param vel_trace: The velocity trace
    :type vel_trace: Trace
    :param disp_trace: The displacement trace
    :type disp_trace: Trace
    :param tp: The arrival time
    :type tp: int
    :param directory: Where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    station = disp_trace.stats.station
    channel = disp_trace.stats.channel
    delta = disp_trace.stats.delta
    time = np.arange(len(disp_trace.data[:tp])) * delta
    fig, axes = plt.subplots(3, 1, figsize=(6, 8))
    ax1, ax2, ax3 = axes
    ax1.set_title("Trend of accelerometer before signal")
    ax1.plot(time, trace.data[:tp])
    ax2.set_title("Trend of velocity before signal")
    ax2.plot(time, vel_trace.data[:tp])
    ax3.set_title("Trend of displacement before signal")
    ax3.plot(time, disp_trace.data[:tp])
    plt.savefig(directory / "plots" / f"{station}_{channel}_plot4")
    plt.close(fig)


def _optional_plots(
    vel: Trace,
    disp: Trace,
    constants: list,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plots to show effect of this baseline removal procedure on selected
    strong motion data

    :param vel: The velocity trace
    :type vel: Trace
    :param disp: The displacement trace
    :type disp: Trace
    :param constants: The polynomial fits and relevant times
    :type constants: list
    :param directory: Where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    station = disp.stats.station
    channel = disp.stats.channel
    delta = disp.stats.delta
    a0, a1, a2, a3, t1, t2, tp, t_d0, t_pgd, t_pga, t_end, gps, t_jump = constants

    length = len(disp.data)
    disp_corr = _correct_disp(length, delta, a1, a2, a3, t1, t2, t_end)
    fig2, ax2 = plt.subplots(1, 1, figsize=(8, 5))
    time = np.arange(len(disp.data)) * delta
    plt.title("Fit of displacement trend to uncorrected displacement")
    ax2.set_xlabel("time[s]")
    ax2.set_ylabel("disp[mt]")
    ax2.plot(time, disp.data)
    ax2.plot(time, disp_corr, "k", linewidth=2)
    ax2.axvline(x=tp * delta, color="k", label="$t_p$")
    ax2.axvline(x=t1 * delta, color="r", label="$t_1$")
    ax2.axvline(x=t2 * delta, color="g", label="$t_2$")
    ax2.axvline(x=t_end * delta, label="$t_f$")
    ax2.legend()
    plt.savefig(directory / "plots" / f"{station}_{channel}_plot2")
    plt.close(fig2)

    fig3, ax3 = plt.subplots(1, 1, figsize=(8, 5))
    plt.title("Fit of corrected displacement to heaviside")
    ax3.set_xlabel("time[s]")
    ax3.set_ylabel("disp[mt]")
    ax3.plot(time, disp.data - disp_corr)
    gps_trace = np.zeros(len(disp.data))
    gps_trace[t_jump:] = gps
    ax3.plot(time, gps_trace, "k", linewidth=2)
    ax3.axvline(x=tp * delta, color="k", label="$t_p$")
    ax3.axvline(x=t1 * delta, color="r", label="$t_1$")
    ax3.axvline(x=t2 * delta, color="g", label="$t_2$")
    ax3.axvline(x=t_end * delta, label="$t_f$")
    ax3.legend()
    plt.savefig(directory / "plots" / f"{station}_{channel}_plot3")
    plt.close(fig3)

    vel = disp.copy()
    vel.differentiate()
    delta2 = int(15 / delta)
    fig5, ax5 = plt.subplots(1, 1, figsize=(5, 3))
    plt.title("Detect arrival of P-wave to accelerometer")
    ax5.set_xlabel("time[s]")
    ax5.set_ylabel("vel[mt/s]")
    ax5.plot(time[: tp + delta2], vel.data[: tp + delta2])
    ax5.axvline(x=tp * delta, color="k", label="$t_p$")
    ax5.legend()
    plt.savefig(directory / "plots" / f"{station}_{channel}_plot5")
    plt.close(fig5)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(), help="folder where there are input files"
    )
    parser.add_argument(
        "-p",
        "--plot_results",
        action="store_true",
        help="plot results for baseline removal procedure",
    )
    args = parser.parse_args()
    os.chdir(args.folder)
    if not os.path.isdir("plots"):
        os.mkdir("plots")
    if not os.path.isdir("../int_STR"):
        os.mkdir("../int_STR")
    time0 = time.time()
    files = glob.glob("acc*")
    results: list = []
    for file in files:
        wang_process(file, plot=True)
    print("Time spent: ", time.time() - time0)

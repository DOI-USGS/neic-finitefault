# -*- coding: utf-8 -*-
"""
Script for implementing the baseline removal scheme proposed by
Wang et al. [2011]_

.. [2011] Wang, R., B. Schurr, C. Milkereit, Z. Shao, and M. Jin (2011),
     An improved automatic scheme for empirical baseline correction of digital
     strong motion records, Bull. Seism. Soc. Am.,101(5), 2029–2044,
     doi:10.1785/0120110039.
"""


from obspy import read
import numpy as np
from scipy.optimize import shgo
import os
import argparse
import matplotlib.pyplot as plt
import time
import glob
import scipy.integrate as integrate


def wang_process(file, plot=False):
    """Scheme to remove baselines from strong motion data,
    and eventually plot the results.
    
    :param file: file with strong motion data, in sac format
    :param plot: whether to plot results of baseline removal procedure
    :type file: string
    :type plot: bool, optional

    .. note::

       This has only been tested with waveform files in SAC format. We can't
       ensure this works with other formats.

    .. note::

       The waveform should have some data prior to the start of shaking.
       We have found it is best if there is data up to at least 20 seconds
       prior to origin time
    """
    trace, vel_trace, disp_trace, constants = _wang_baselines(file)
    if not constants:
        return None
    a0, a1, a2, a3, t1, t2, tp, t_d0, t_pgd, t_pga, t_end = constants
    disp_data = disp_trace.data
    delta = disp_trace.stats.delta
    disp_data2, gps, t_jump = _gps_jump(
            disp_data, t_end, t1, t2, a1, a2, a3, tp, delta)
    disp_trace2 = disp_trace.copy()
    disp_trace2.data = disp_data2
    vel_trace2 = disp_trace2.copy()
    vel_trace2.differentiate()
    if plot:
#        t1 = int(t1 / delta)
#        t2 = int(t2 / delta)
        constants = [a0, a1, a2, a3, t1, t2, tp, t_d0, t_pgd, t_pga, t_end,
                     gps, t_jump]                
        _optional_plots(vel_trace, disp_trace, disp_trace2, constants)
        _opt_plots2(trace, vel_trace, disp_trace, tp)
    return vel_trace2


def _wang_baselines(file):
    """
    """
    st = read(file)
    st1 = st.copy()   
    trace = st1[0]
    if np.max(np.abs(trace.data)) < 10 ** -8:
        return trace, trace, trace, None
    delta = trace.stats.delta
#
# we define some relevant times
#
    time_p = _first_arrival(trace)
    time_end, time_pga = _time_end_pga(trace)
    if time_p >= min(time_end, time_pga):
        return trace, trace, trace, None
    preevent_baseline = np.mean(trace.data[:time_p])# -  int(1 / delta)])
    trace.data = trace.data - preevent_baseline  
    trace.data = trace.data[:time_p + 4*(time_end - time_p)]
#
# we integrate to velocity and displacement and remove cuadratic trend at the
# end of the signal.
#
    new_trace = trace.copy()
    vel_trace = new_trace.integrate(method='cumtrapz')
    vel_trace.data = vel_trace.data - vel_trace.data[time_p]# - 1 * int(1/delta)]
    new_trace = vel_trace.copy()
    disp_trace = new_trace.integrate(method='cumtrapz')
    disp_trace.data = disp_trace.data - disp_trace.data[time_p]# - 1 * int(1 / delta)]
    time_d0, time_pgd = _time_d0_pgd(disp_trace, time_p, time_pga, time_end)
    disp_tail = disp_trace.data[time_end:]
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
    constraint0 = ({'type': 'ineq', 'fun': _constraint})
    bounds = [(time_pgd * delta, time_end * delta),
              (max_val * delta, time_end * delta)]
    if time_pgd > time_end or max_val > time_end:
        return trace, trace, trace, None
    args = (disp_data, a1, a2, a3, time_p, time_end, time_pgd, delta)

    res = shgo(_function, bounds, iters=6, constraints=constraint0, args=args)
    time1, time2 = res.x
    success = res.success
    if not success:
        print(file)
        print(res)
    time1 = int(time1 / delta)
    time2 = int(time2 / delta)
    constants = [a0, a1, a2, a3, time1, time2, time_p, time_d0, time_pgd,
                 time_pga, time_end]     
    return trace, vel_trace, disp_trace, constants


def _function(time, disp_data, a1, a2, a3, time_p, time_end, time_pgd, delta):
    """Compute misfit of corrected displacement to a heaviside.
    """
    time1, time2 = time
    time1 = int(time1 / delta)
    time2 = int(time2 / delta)
    length = len(disp_data)
    disp_corr = _correct_disp(
        length, delta, a1, a2, a3, time1, time2, time_end)
    disp_data2 = disp_data - disp_corr
    gps = np.mean(disp_data2[time_end:])
    comparison = disp_data2 - gps / 2
    after = comparison[1:]
    before = comparison[:-1]
    crossings = np.nonzero(comparison == 0)[0]
    crossings = list(crossings) + list(np.nonzero(before * after < 0)[0])
    crossings = [index for index in crossings if index > time_p]
    if not crossings:
        return np.sum(disp_data2 ** 2)
    results = [_heaviside_fit(int(index), disp_data2, gps)\
        for index in crossings]
    return min(results)


def _gps_jump(disp_data, time_end, time1, time2, a1, a2, a3, time_p, delta):
    """Get time and amplitude of heaviside jump
    """
    length = len(disp_data)
    disp_corr = _correct_disp(
            length, delta, a1, a2, a3, time1, time2, time_end)
    disp_data2 = disp_data - disp_corr
    gps = np.mean(disp_data2[time_end:])
    comparison = disp_data2 - gps / 2
    after = comparison[1:]
    before = comparison[:-1]
    crossings = np.nonzero(comparison == 0)[0]
    crossings = list(crossings) + list(np.nonzero(before * after < 0)[0])
    crossings = [index for index in crossings if index > time_p]
    if not crossings:
        return np.sum(disp_data ** 2), gps, 0
    results = [_heaviside_fit(int(index), disp_data2, gps)\
        for index in crossings]
    best_result = min(results)
    zipped = zip(crossings, results)
    time_jump = next(index for index, result in zipped if result==best_result)
    return disp_data2, gps, time_jump


def _constraint(x):
    """
    """
    return x[1] - x[0]


def _integral_arias(trace, delta):
    r"""We compute :math:`\int_0^T \ddot{x}^2(t)dt`. 
    
    :param trace: input waveform
    :param delta: sampling period of waveform.
    """
    trace = trace ** 2
    integral = integrate.cumtrapz(trace, dx=delta, initial=0)
    return integral
    
    
def _first_arrival(trace):
    """A binary search hand-made picker.
    """
    accel = trace.data
    delta = trace.stats.delta
    trace0 = accel[0:int(8 / delta)]
    val0 = max(_integral_arias(trace0, delta))
    for k in range(8, int(len(accel) / delta), 8):
        trace1 = accel[int(k / delta):int((k + 8) / delta)]
        if len(trace1) == 0:
            continue
        val = max(_integral_arias(trace1, delta))
        if val > 10 * val0: break
    start0 = k
    trace0 = accel[int((start0 - 12) / delta):int((start0 - 8) / delta)]
    if len(trace0) == 0:
        return int(start0 / delta)
    val0 = max(_integral_arias(trace0, delta))
    for k in range(start0 - 8, start0 + 8, 4):
        trace1 = accel[int(k / delta):int((k + 4) / delta)]
        val = max(_integral_arias(trace1, delta))
        if val > 10 * val0: break
    start0 = k
    trace0 = accel[int((start0 - 6) / delta):int((start0 - 4) / delta)]
    val0 = max(_integral_arias(trace0, delta))
    for k in range(start0 - 4, start0 + 4, 2):
        trace1 = accel[int(k / delta):int((k + 2) / delta)]
        val = max(_integral_arias(trace1, delta))
        if val > 10 * val0: break
    start0 = k
    trace0 = accel[int((start0 - 3) / delta):int((start0 - 2) / delta)]
    val0 = max(_integral_arias(trace0, delta))
    for k in range(start0 - 2, start0 + 2, 1):
        trace1 = accel[int(k / delta):int((k + 1) / delta)]
        val = max(_integral_arias(trace1, delta))
        if val > 10 * val0: break
    start0 = k# - 1
    return int(start0 / delta)
    
    
def _time_end_pga(trace):
    """
    """
    accel = trace.data
    delta = trace.stats.delta
    pga = np.max(np.abs(accel))
    time_pga = next(index for index, v in enumerate(accel) if abs(v) == pga)
    integral = _integral_arias(accel, delta)
    integral = integral / integral[-1]
    time_end = next(index for index, v in enumerate(integral) if v > 0.95)
    return time_end, time_pga
    
    
def _time_d0_pgd(trace, time_p, time_pga, time_end):
    """
    """
    zero_crossings = [index for index in range(time_p, time_end)\
        if _iz_zero_crossing(trace.data, index)]
#
# we find the last zero crossing of displacement prior to the estimated trace
# end, for uncorrected displacement. Then we find the last local optimum of
# displacement prior to this instant.
#
    if zero_crossings:
        time_d0 = zero_crossings[-1]
        pgd = np.max(np.abs(trace.data[time_p:time_d0]))
        time_pgd = next(index for index in range(time_d0)\
                        if abs(trace.data[index])==pgd)
    else:
#
# In case we can't find these values, we search them in some value previous to
# the PGA time
#
        vel_trace = np.gradient(trace, trace.stats.delta)
#        print(time_p, time_pga)
#        print(trace.stats.station, trace.stats.channel)
        pgv = np.max(np.abs(vel_trace.data[time_p:time_pga]))
        time_pgv = next(index for index in range(time_pga)\
                        if abs(vel_trace.data[index])==pgv)
        time_d0 = time_pgv
        zero_crossings = [index for index in range(time_p, time_pgv)\
                          if _iz_zero_crossing(vel_trace.data, index)]
#
# In case we still can't find a t_pgd estimation, we use a gross estimation.
#
        time_pgd = zero_crossings[-1] if zero_crossings\
        else int((time_p + time_pga) / 2)
    return time_d0, time_pgd
    
    
def _iz_zero_crossing(disp, index):
    """find zero crossings of signals
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
        nonzero0 = next(
            i for i in range(index + 1) if abs(disp[index - i]) > 0)
        val0 = disp[nonzero0]
        nonzero1 = next(
            i for i in range(len(disp) - index) if abs(disp[index + i]) > 0)
        val1 = disp[nonzero1]
        if val1 * val0 < 0:
            return True
    else:
        return False


def _polinomial_fit(disp_tail, delta, n):
    """Find best fitting polinomial of degree n, and Akaike criteria of
    polinomial model
    """
    time = np.arange(len(disp_tail)) * delta
    matrix = np.fliplr(np.vander(time, n + 1))
    sym_mat = np.dot(np.transpose(matrix), matrix)
    inv_mat = np.linalg.inv(sym_mat)
    full_mat = np.dot(inv_mat, np.transpose(matrix))
    coeffs = np.dot(full_mat, disp_tail)
    polin = np.zeros(len(disp_tail))
    for i, coeff in enumerate(coeffs):
        polin = polin + coeff*time**i
    log_likelihood = - np.sum((disp_tail - polin)**2)*5000
    aic = _akaike_criteria(len(disp_tail), n + 1, log_likelihood)
    return np.concatenate((np.dot(full_mat, disp_tail), [aic]))
    
    
def _akaike_criteria(nsize, nparam, log_likelihood):
    """
    """
    if nsize > nparam + 1:
        return 2 * nparam - 2 * log_likelihood\
            + 2 * nparam * (1 + nparam) / (nsize - nparam - 1)
    else:
        return 2 * nparam - 2 * log_likelihood
    
    
def _heaviside_fit(t_jump, disp, static):
    """Fit to heaviside
    """
    length = len(disp)
    step = static * np.ones(length)
    step[:t_jump] = np.zeros(t_jump)
    error = (disp - step) ** 2
    return np.sum(error)
        
        
def _correct_disp(length, delta, a1, a2, a3, t1, t2, t_end):
    """Compute estimated baselines of displacement signal
    """
    new_time = np.arange(length) * delta
    disp_corr = np.zeros(length)
    p0 = a1 + 2*a2*(t2 - t_end)*delta + 3*a3*(t2 - t_end)**2 * delta**2
    if t1 < t2:
        disp_corr[t1:t2]\
            = (p0/(t2 - t1)/delta/2) * (new_time[t1:t2] - t1*delta)**2
    disp_corr[t2:]\
        = p0*(t2 - t1)*delta/2 + p0*(new_time[t2:] - t2*delta)\
        + (a2 + 3*a3*(t2 - t_end)*delta)*(new_time[t2:] - t2*delta)**2\
        + a3*(new_time[t2:] - t2*delta)**3
    return disp_corr
    
    
def _opt_plots2(trace, vel_trace, disp_trace, tp):
    """Plots to show effect of this baseline removal procedure on selected
    strong motion data.
    """
    station = disp_trace.stats.station
    channel = disp_trace.stats.channel
    delta = disp_trace.stats.delta
    time = np.arange(len(disp_trace.data[:tp])) * delta
    fig, axes = plt.subplots(3, 1, figsize=(6, 8))
    ax1, ax2, ax3 = axes
    ax1.set_title('Trend of accelerometer before signal')
    ax1.plot(time, trace.data[:tp])
    ax2.set_title('Trend of velocity before signal')
    ax2.plot(time, vel_trace.data[:tp])
    ax3.set_title('Trend of displacement before signal')
    ax3.plot(time, disp_trace.data[:tp])
    plt.savefig('plots/{}_{}_plot4'.format(station, channel))
    plt.close(fig)
    
    
def _optional_plots(vel, disp, corr_disp, constants):
    """Plots to show effect of this baseline removal procedure on selected
    strong motion data.
    """
    station = disp.stats.station
    channel = disp.stats.channel
    delta = disp.stats.delta
    a0, a1, a2, a3, t1, t2, tp, t_d0, t_pgd, t_pga, t_end, gps, t_jump = constants
    
    length = len(disp.data)
    disp_corr = _correct_disp(length, delta, a1, a2, a3, t1, t2, t_end)
    fig2, ax2 = plt.subplots(1, 1, figsize=(8, 5))
    time = np.arange(len(disp.data)) * delta
    plt.title('Fit of displacement trend to uncorrected displacement')
    ax2.set_xlabel('time[s]')
    ax2.set_ylabel('disp[mt]')
    ax2.plot(time, disp.data)
    ax2.plot(time, disp_corr, 'k', linewidth = 2)
    ax2.axvline(x=tp * delta, color='k', label='$t_p$')
    ax2.axvline(x=t1 * delta, color='r', label='$t_1$')
    ax2.axvline(x=t2 * delta, color='g', label='$t_2$')
    ax2.axvline(x=t_end * delta, label='$t_f$')
    ax2.legend()
    plt.savefig('plots/{}_{}_plot2'.format(station, channel))
    plt.close(fig2)
    
    fig3, ax3 = plt.subplots(1, 1, figsize=(8, 5))
    plt.title('Fit of corrected displacement to heaviside')
    ax3.set_xlabel('time[s]')
    ax3.set_ylabel('disp[mt]')
    ax3.plot(time, disp.data - disp_corr)
    gps_trace = np.zeros(len(disp.data))
    gps_trace[t_jump:] = gps
    ax3.plot(time, gps_trace, 'k', linewidth = 2)
    ax3.axvline(x=tp * delta, color='k', label='$t_p$')
    ax3.axvline(x=t1 * delta, color='r', label='$t_1$')
    ax3.axvline(x=t2 * delta, color='g', label='$t_2$')
    ax3.axvline(x=t_end * delta, label='$t_f$')
    ax3.legend()
    plt.savefig('plots/{}_{}_plot3'.format(station, channel))
    plt.close(fig3)

    vel = disp.copy()
    vel.differentiate()
    delta2 = int(15 / delta)
    fig5, ax5 = plt.subplots(1, 1, figsize=(5, 3))
    plt.title('Detect arrival of P-wave to accelerometer')
    ax5.set_xlabel('time[s]')
    ax5.set_ylabel('vel[mt/s]')
    ax5.plot(time[:tp + delta2], vel.data[:tp + delta2])
    ax5.axvline(x=tp * delta, color='k', label='$t_p$')
    ax5.legend()
    plt.savefig('plots/{}_{}_plot5'.format(station, channel))
    plt.close(fig5)
    return
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-p", "--plot_results", action="store_true",
                        help="plot results for baseline removal procedure")
    args = parser.parse_args()
    os.chdir(args.folder)
    if not os.path.isdir('plots'):
        os.mkdir('plots')
    if not os.path.isdir('../int_STR'):
        os.mkdir('../int_STR')
    time0 = time.time()
    files = glob.glob('acc*')
    results = []
    for file in files:
        wang_process(file, plot=True)
    print('Time spent: ', time.time() - time0)

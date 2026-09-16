import json
import os
import pathlib
import subprocess
from io import StringIO
from os import chdir
from shlex import split
from shutil import copy
from typing import Union
from urllib.request import Request, urlopen

import cartopy.crs as ccrs  # type: ignore
import cartopy.feature as cfeature  # type: ignore
import imageio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  # type: ignore
import requests
from alive_progress import alive_bar
from clawpack.geoclaw import dtopotools  # type: ignore
from clawpack.geoclaw import topotools  # type: ignore
from obspy import UTCDateTime  # type: ignore
from obspy.core.stream import Stream  # type: ignore
from obspy.core.trace import Trace  # type: ignore
from obspy.signal import filter  # type: ignore
from pyproj import Geod
from scipy.interpolate import interp1d  # type: ignore
from scipy.signal import butter, find_peaks, sosfilt, sosfiltfilt  # type: ignore

import ffm.fault_plane as pf
import ffm.management as mng
import ffm.plane_management as pl_mng
import ffm.plot_graphic_NEIC as pg
from ffm.data_management import _dict_trace
from ffm.management import default_dirs
from ffm.read_config import PROJECT_DIRECTORY


class EmptyDataError(Exception):
    """Raised when API request is successful but data returned is empty."""

    pass


class RequestError(Exception):
    """Raised when API request fails."""

    pass


def make_DART_GFs_DS_SS(
    directory: Union[pathlib.Path, str],
    tensor_info: dict,
    segments_data: dict,
    dart_data: dict,
    dt: float = 60,
    tfinal: int = 10800,
    subfault_start: int = 0,
    time_domain: bool = False,
    plots: bool = False,
):
    """
    Use subfault locations and DART location to make tsunami Green's functions

    :param directory: Path to working directory, with dart_waves.json file.
    :type directory: Union[pathlib.Path, str]
    :param tensor_info: Event origin information
    :type tensor_info: dict
    :param segments_data: Fault plane segment information
    :type segments_data: dict
    :param dart_data: Dart station information
    :type dart_data: dict
    :param dt: Desired sampling rate of the Green's functions
    :type dt: float
    :param tfinal: Length of Green's functions, in seconds
    :type tfinal: int
    :param subfault_start: Subfault to start calculation (use only if previous run failed midway)
    :type subfault_start: int
    :param time_domain: Whether to save Green's functions in time domain, default is False
    :type time_domain: bool
    :param plots: Whether to save plot of subfault slip/rake and resulting seafloor deformation, default is False
    :type plots: bool
    """

    print("Computing DART GFs. This will take a while.")
    #####################
    ### Get Organized ###
    #####################
    directory = pathlib.Path(directory).resolve()
    run_directory = os.getcwd()  # What directory is this being run from
    GF_directory = directory / "_GF"
    if not os.path.exists(GF_directory):
        os.makedirs(GF_directory)
    geoclaw_parameter_file = directory / "geoclaw_parameter_file.json"
    ##############################
    ### NUMBER OF OUTPUT TIMES ###
    ##############################
    num_output_times = 0  ### Just making GFs at DART stations, so don't need any general tsunami output
    clines_zeta = [1, 3, 5, 7, 10, 15, 20, 25, 30, 50]  # max amplitude for plot, in m

    ### Find DART locations and create "gauge" and "high_res_region" dicts
    high_res_regions, gauges, lon_dart, lat_dart = dart_high_res_locations(
        dart_data, segments_data, tensor_info
    )
    ### Find locations of subfaults in model space ###
    lon_sub, lat_sub, dep_sub, sf_len, sf_wid, strike, dip, rake = (
        get_subfault_locations(segments_data, tensor_info)
    )

    ### Find model region, with 0.5 degrees of space beyond your fault and stations ###
    min_lon = min(min(lon_sub), min(lon_dart)) - 0.5
    max_lon = max(max(lon_sub), max(lon_dart)) + 0.5
    min_lat = min(min(lat_sub), min(lat_dart)) - 0.5
    max_lat = max(max(lat_sub), max(lat_dart)) + 0.5
    model_region = {
        "west": min_lon,
        "east": max_lon,
        "south": min_lat,
        "north": max_lat,
    }
    ### get topography for this region of interest ###
    topo_path = make_topo(min_lon, max_lon, min_lat, max_lat, directory)
    copy(topo_path, GF_directory)

    ##########################
    ### Initialize geoclaw ###
    ##########################
    geoclaw_init = "make .exe"
    with open(directory / "output.log", "a") as outfile:
        subprocess.call(
            split(geoclaw_init), stdout=outfile
        )  ### Pipe all the unneccesary output to a file
    outfile.close()

    source_dx = 2 / 60.0  # dtopo spatial sampling (2 arcmin)
    fcorner2 = [
        1.0 / (120 * 60),
        1.0 / (5 * 60),
    ]  # we will bandpass filter the GFs from 5 to 120 mins for consistency with data

    #######################################
    ###       Loop over subfaults       ###
    ### First dip-slip then strike-slip ###
    #######################################
    gf_rake = [90, 0]  # pure dip-slip, pure strike-slip
    gf_type = ["dip", "stk"]  # dip-slip, strike-slip, for naming

    for rake_direction in range(len(gf_type)):
        print(f"Computing {gf_type[rake_direction]} GFs")
        ### Loop over subfaults ###
        with alive_bar(int(len(lon_sub) - int(subfault_start))) as bar:
            for subflt in range(len(lon_sub) - int(subfault_start)):
                subf = int(subflt + subfault_start)

                ##############################################
                ### Create dtopo for slip on each subfault ###
                ##############################################
                fault = dtopotools.Fault()
                fault.subfaults = []
                fault.rupture_type = "static"

                subfault = dtopotools.SubFault()
                subfault.length = sf_len[subf] * 1000.0  # convert from km to m
                subfault.width = sf_wid[subf] * 1000.0  # convert from km to m
                subfault.depth = dep_sub[subf] * 1000.0  # convert from km to m
                subfault.strike = strike[subf]
                subfault.slip = 1.0  # GF slip in m
                subfault.rake = gf_rake[rake_direction]
                subfault.dip = dip[subf]
                subfault.longitude = lon_sub[subf]
                subfault.latitude = lat_sub[subf]
                subfault.coordinate_specification = "centroid"
                fault.subfaults = [subfault]
                if plots == True:
                    # Quick GeoClaw plot of slip and rake on subfault #
                    fault.plot_subfaults(plot_rake=True, slip_color=True)
                    plt.savefig(f"{GF_directory}/{str(subf).rjust(4, '0')}_GF_slip.png")

                x, y = fault.create_dtopo_xy(dx=source_dx, buffer_size=1.5)
                gf_dtopo = fault.create_dtopography(x, y, times=[1.0])
                dtopo_fname = f"{GF_directory}/{str(subf).rjust(4, '0')}_{gf_type[rake_direction]}.tt1"
                gf_dtopo.write(dtopo_fname, dtopo_type=1)

                if plots == True:
                    # Quick GeoClaw plot of seafloor deformation
                    fig = plt.figure(figsize=(12, 5))
                    gf_dtopo.plot_dZ_colors(t=1)  # , cmax_dZ=0.01)
                    plt.savefig(
                        f"{GF_directory}/{str(subf).rjust(4, '0')}_{gf_type[rake_direction]}_SeafloorDeformation.png"
                    )
                    plt.close("all")
                ##############################
                ### Fill setrun parameters ###
                ##############################
                setrun_parameter = {
                    "model_region": model_region,
                    "high_res_regions": high_res_regions,
                    "gauges": gauges,
                    "topo_path": f"{GF_directory}/{os.path.basename(topo_path)}",
                    "dtopo_path": f"{GF_directory}/{os.path.basename(dtopo_fname)}",
                    "dtopo_type": 1,
                    "dt_max_dtopo": 0.2,
                    "tfinal": tfinal,
                    "num_output_times": num_output_times,
                    "wave_tolerance": 1.0e-5,
                    "clines_zeta": clines_zeta,
                }
                with open(geoclaw_parameter_file, "w") as fout:
                    json_dumps_str = json.dumps(setrun_parameter, indent=4)
                    print(json_dumps_str, file=fout)
                fout.close()

                ##################################
                ### Move to _GF folder and run ###
                ##################################
                copy(geoclaw_parameter_file, GF_directory)
                copy(PROJECT_DIRECTORY / "src/ffm/geoclaw_setrun.py", GF_directory)
                copy(PROJECT_DIRECTORY / "src/ffm/Makefile", GF_directory)
                copy(PROJECT_DIRECTORY / "src/ffm/tsunami_plot.py", GF_directory)
                chdir(GF_directory)

                geoclaw_run = "make .output"
                with open(directory / "output.log", "a") as outfile:
                    subprocess.call(
                        split(geoclaw_run), stdout=outfile
                    )  ### Pipe all the unneccesary output to a file
                outfile.close()
                chdir(run_directory)
                ### Make geoclaw output into usable format ###
                tlims = [0, tfinal]
                num_dart = len(dart_data)
                for kdart in range(num_dart):
                    # load in dart output
                    id_dart = dart_data[kdart]["name"]
                    gf_file_freq = f"{directory}/{id_dart}_{gf_type[rake_direction]}_gf.txt"  # save GF in frequency domain for inversion
                    if (
                        id_dart.isalpha()
                    ):  # GeoClaw can't handle letters as station IDs, so make a dummy name with numbers
                        alpha_ind = ord(id_dart[-1]) - 64  # Convert NZ[A-L] to number
                        id_dart = str(alpha_ind).rjust(5, "0")
                    gf = f"{GF_directory}/_output/gauge{id_dart}.txt"
                    greens_func = np.genfromtxt(gf)
                    # Get time vector
                    time_gf = greens_func[:, 1]
                    # Get amplitude of GF
                    eta_gf = greens_func[:, 5]
                    # Trim length and resample time series to regular interval
                    itrim = np.intersect1d(
                        np.where(time_gf >= tlims[0])[0],
                        np.where(time_gf <= tlims[1])[0],
                    )
                    time_gf = time_gf[itrim]
                    eta_gf = eta_gf[itrim]
                    time_0 = np.arange(tlims[0], tlims[1] + dt, dt)
                    eta_gf = np.interp(time_0, time_gf, eta_gf)
                    time_out = np.arange(0, tlims[1] + dt, dt)
                    eta_gf = np.interp(time_out, np.r_[0, time_0], np.r_[0, eta_gf])

                    gf_fil = filter.bandpass(
                        eta_gf,
                        fcorner2[0],
                        fcorner2[1],
                        1 / dt,
                        corners=2,
                        zerophase=True,
                    )
                    # Divide by 100 to get from 1 m back to 1 cm
                    gf_fil = gf_fil / 100.0

                    if time_domain == True:
                        # Write out the filtered time domain GFs to text file
                        with open(
                            f"{GF_directory}/dart.{id_dart}.{str(subf).rjust(4, '0')}.GF.{gf_type[rake_direction]}.txt",
                            "w",
                        ) as gf_out:
                            for kint in range(len(time_out)):
                                gf_out.write(
                                    "%10.4f \t %10.20e \t %10.20e\n"
                                    % (time_out[kint], gf_fil[kint], eta_gf[kint])
                                )
                        gf_out.close()
                        if plots == True:
                            plt.figure(figsize=(8, 6))
                            plt.plot(time_out / 3600.0, gf_fil, "k", label="bandpass")
                            plt.plot(
                                time_out / 3600.0, eta_gf / 100.0, "b", label="raw"
                            )
                            plt.legend()
                            plt.xlabel("Time (hrs)")
                            plt.ylabel("Amplitude (cm)")
                            plt.title(
                                f"{gf_type[rake_direction]} GF, subfault {str(subf).rjust(4, '0')} to DART {id_dart}"
                            )
                            plt.savefig(
                                f"{GF_directory}/dart.{id_dart}.{str(subf).rjust(4, '0')}.GF.{gf_type[rake_direction]}.png"
                            )
                    # Convert GF from time domain to frequency domain and add to running txt file
                    dart_gf_time_freq_domain(
                        id_dart,
                        int(dt),
                        subf,
                        time_out,
                        gf_fil,  # send bandpass filtered GF in cm to convert to frequency domain
                        gf_file_freq,
                        GF_directory,
                    )
                bar()  # update progress bar, one subfault complete


def get_subfault_locations(segments_data, tensor_info):
    """
    Get subfault location and orientation information from segments_data.json
    :param segments_data: Dictionary of fault segment parameters
    :type segments_data: dict
    """
    subfault_lons = []
    subfault_lats = []
    subfault_depths = []
    subfault_length = []
    subfault_width = []
    subfault_strike = []
    subfault_dip = []
    subfault_rake = []
    segments = segments_data["segments"]
    rupture_vel = segments[0]["rupture_vel"]
    delta_strike = segments[0]["delta_strike"]
    delta_dip = segments[0]["delta_dip"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )
    subfaults = {"delta_strike": delta_strike, "delta_dip": delta_dip}
    subfaults2 = pf._point_sources_def(rise_time, rupture_vel, subfaults)
    strike_ps = int(subfaults2["strike_ps"] / 2)
    dip_ps = int(subfaults2["dip_ps"] / 2)
    latitudes = [ps_segment[:, :, dip_ps, strike_ps, 0] for ps_segment in point_sources]
    longitudes = [
        ps_segment[:, :, dip_ps, strike_ps, 1] for ps_segment in point_sources
    ]
    depths = [ps_segment[:, :, dip_ps, strike_ps, 2] for ps_segment in point_sources]
    if len(segments) == 1:
        subfault_lats = latitudes[0].flatten()
        subfault_lons = longitudes[0].flatten()
        subfault_depths = depths[0].flatten()
        subfault_length = np.ones(len(subfault_lats)) * delta_strike
        subfault_width = np.ones(len(subfault_lats)) * delta_dip
        subfault_strike = np.ones(len(subfault_lats)) * segments[0]["strike"]
        subfault_dip = np.ones(len(subfault_lats)) * segments[0]["dip"]
        subfault_rake = np.ones(len(subfault_lats)) * segments[0]["rake"]
    else:  # multi-segment
        zipped2 = zip(
            segments,
            point_sources,
            latitudes,
            longitudes,
            depths,
        )
        for i_segment, fault_segment_data in enumerate(zipped2):
            plane_info = segments[i_segment]
            strike = plane_info["strike"]
            dip = plane_info["dip"]
            rake = plane_info["rake"]
            lat_fault = fault_segment_data[2].flatten()
            lon_fault = fault_segment_data[3].flatten()
            depth_fault = fault_segment_data[4].flatten()
            zipped3 = zip(
                lon_fault,
                lat_fault,
                depth_fault,
            )
            for line in zipped3:
                lon, lat, dep = line
                subfault_lons.append(lon)
                subfault_lats.append(lat)
                subfault_depths.append(dep)
                subfault_length.append(delta_strike)
                subfault_width.append(delta_dip)
                subfault_strike.append(strike)
                subfault_dip.append(dip)
                subfault_rake.append(rake)

    return (
        subfault_lons,
        subfault_lats,
        subfault_depths,
        subfault_length,
        subfault_width,
        subfault_strike,
        subfault_dip,
        subfault_rake,
    )


def make_dart_json_no_observations(
    directory: Union[pathlib.Path, str],
    dart_station_file: Union[pathlib.Path, str],
):
    """
    Create a dart_waves.json file without fetching data (for use for fwd modeling without data going into inversion)

    :param directory: Working directory (where tensor_info.json is located)
    :type directory: Union[pathlib.Path, str]
    :param dart_station_file: Path to file with list of DART locations
    :type dart_station_file: Union[pathlib.Path, str]
    """
    directory = pathlib.Path(directory)
    tensor_info = json.load(open(directory / "tensor_info.json"))
    hypo_lat = tensor_info["lat"]
    hypo_lon = tensor_info["lon"]
    df_dart = pd.read_csv(dart_station_file, sep=r"\s+", index_col=False)
    dart_ID = df_dart["Station_ID"]
    longitude = df_dart["Longitude"]
    latitude = df_dart["Latitude"]
    info_traces = []
    for idx, station in enumerate(dart_ID):
        dist_km, azimuth, baz = mng._distazbaz(
            latitude[idx], longitude[idx], hypo_lat, hypo_lon
        )
        dist_degrees = dist_km / 111.11
        if dist_degrees < 20:
            file: str = "None"
            info = _dict_trace(
                file=file,
                name_station=str(dart_ID[idx]),
                azimuth=azimuth,
                channel="dart",
                distance=dist_degrees,
                dt=60.0,
                duration=10800,  # TO DO: CHECK THIS
                n_start_obs=0,
                trace_weight=1.0,
                wavelet_weight="0 0 2 2 2 2 2 2\n",
                observed_trace=None,
                synthetic_trace=None,
                location=[latitude[idx], longitude[idx]],
            )
            info_traces.append(info)
    with open(directory / "dart_waves.json", "w") as f:
        json.dump(
            info_traces,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )


def query_dart_api(origin_time, station, dt=60, fcorner=[5, 120]):
    """
    Makes query to NOAA or New Zealand API for DART observations

    :param origin_time: Earthquake origin time
    :type origin_time:
    :param station: DART station from which to query data
    :type station: str
    :param dt: Data sampling rate, in seconds, for interpolation. The default is 60.
    :type dt: float, optional
    :param fcorner: bandpass filter corner frequencies in minutes
    :type fcorner: List[float]. The default is [5,120]
    :raises EmptyDataError: No data available for query.
    :raises RequestError: Issue with query request.
    """
    # Use origin time to determine query range
    startmonth = origin_time.month
    startday = origin_time.day
    startyear = origin_time.year
    if origin_time.hour > 16:
        end_time = (
            origin_time + 86400
        )  # If event was toward end of day, grab next day of observations as well
    else:
        end_time = origin_time
    endmonth = end_time.month
    endday = end_time.day
    endyear = end_time.year
    fcorner2 = [1.0 / (fcorner[1] * 60), 1.0 / (fcorner[0] * 60)]
    if station.isnumeric():  # Stations with numeric IDs are NOAA stations
        source = "http://www.ndbc.noaa.gov/dart_data.php"
        url = f"{source}?station={station}&startmonth={startmonth}&startday={startday}&startyear={startyear}&endmonth={endmonth}&endday={endday}&endyear={endyear}"
        r = requests.get(url, headers={"Accept": "text/csv"})
        try:
            if r.status_code == 200:
                req = Request(url)
                csv_len = len(urlopen(req).readlines())
                if csv_len <= 10:
                    raise EmptyDataError
                else:
                    df = pd.read_csv(StringIO(r.text), skiprows=8, sep=r"\s+")
                    df = df.drop(0)
                    df = df[::-1]
                    df["Date"] = df[["#YY", "MM", "DD"]].agg("-".join, axis=1)
                    df["Time"] = df[["hh", "mm", "ss"]].agg(":".join, axis=1)
                    df["Datetime"] = df[["Date", "Time"]].agg("T".join, axis=1)
                    dates = pd.to_datetime(df["Datetime"], utc=True)
                    # calculate unix datetime
                    df["Datetime"] = (
                        dates - pd.Timestamp(str(origin_time))
                    ) // pd.Timedelta("1s")
                    df["HEIGHT"] = df["HEIGHT"].astype(float)
                    df = df.drop(df[df.HEIGHT == 9999.000].index)
                    df = df.reset_index(drop=True)
                    if len(df) < 10:
                        raise EmptyDataError
                    else:
                        print(f"{station}: Data acquired")
                        time = df["Datetime"] - df["Datetime"][0]
                        data = np.array(df["HEIGHT"])
                        peak_idx = find_peaks(data, prominence=1)[0]
                        for idx in range(len(peak_idx)):
                            data[peak_idx[idx]] = data[peak_idx[idx] - 1]
                        t_in = np.arange(time[0], time[len(time) - 1], dt)
                        # Bandpass filter from 5-120 mins #
                        f4 = interp1d(time, data, fill_value=0)
                        model_in = f4(t_in)

                        model_in = model_in - np.mean(model_in[2:])
                        data_fil = filter.bandpass(
                            model_in,
                            fcorner2[0],
                            fcorner2[1],
                            1 / dt,
                            corners=2,
                            zerophase=True,
                        )
                        return (t_in, data_fil)
            else:
                raise RequestError("Request Error encountered")
        except RequestError:
            print(f"Got non-200 status code")
            return (np.array([0]), np.array([0]))
        except EmptyDataError:
            print(f"{station}: Returned empty data")
            return (np.array([0]), np.array([0]))
    elif station.isalpha():  # Stations with alpha IDs are New Zealand stations
        # TILDE requires 2-digit month and day
        if int(startmonth) < 10:
            startmonth = f"0{startmonth}"
        if int(startday) < 10:
            startday = f"0{startday}"
        if int(endmonth) < 10:
            endmonth = f"0{endmonth}"
        if int(endday) < 10:
            endday = f"0{endday}"
        source = "https://tilde.geonet.org.nz/v3/data/"
        domain = "dart/"
        key = station + "/water-height/-/15s/-"
        startdate = f"/{startyear}-{startmonth}-{startday}/"
        enddate = f"{endyear}-{endmonth}-{endday}"
        url = source + domain + key + startdate + enddate
        # print(url)
        r = requests.get(url, headers={"Accept": "text/csv"})
        try:
            if r.status_code == 200:
                df = pd.read_csv(
                    StringIO(r.text), parse_dates=["timestamp"]
                )  # , index_col='timestamp')
                df["Datetime"] = df["timestamp"]
                dates = pd.to_datetime(df["Datetime"], utc=True)
                # calculate unix datetime
                df["Datetime"] = (
                    dates - pd.Timestamp(str(origin_time))
                ) // pd.Timedelta("1s")
                df["HEIGHT"] = df["value"].astype(float)
                if len(df) < 10:
                    raise EmptyDataError
                else:
                    print(f"{station}: Data acquired")
                    time = df["Datetime"] - df["Datetime"][0]
                    t_in = np.arange(time[0], time[len(time) - 1], dt)

                    f4 = interp1d(time, df["HEIGHT"], "cubic")
                    model_in = f4(t_in)

                    model_in = model_in - np.mean(model_in[2:])
                    # data_fil = bandpass(model_in, fcorner2, 1 / dt)
                    data_fil = filter.bandpass(
                        model_in,
                        fcorner2[0],
                        fcorner2[1],
                        1 / dt,
                        corners=2,
                        zerophase=True,
                    )
                    t_in = t_in + df["Datetime"][0]
                    return (t_in, data_fil)
            else:
                raise RequestError("Request Error")
        except RequestError:
            print("Got non-200 status code")
            return (np.array([0]), np.array([0]))
        except EmptyDataError:
            print("Empty data returned.")
            return (np.array([0]), np.array([0]))


def get_dart_observations(
    tensor_info,
    data_directory,
    dt=60.0,
    max_distance=20,
    plot=False,
    fcorner=[5, 120],
    make_json=False,
):
    """
    Gets DART observations from station file with DART Station ID, Longitude, and Latitude.

    :param tensor_info: Moment tensor information
    :type tensor_info: dict
    :param data_directory: Directory where DART observations should be written.
    :type data_directory: Union[pathlib.Path, str]
    :param dt: Desired data sampling rate in seconds. The default is 60.0.
    :type dt: float, optional
    :param max_distance: Maximum hypocenteral distance in degrees. Default is 20.
    :type max_distance: float, optional
    :param plot: Whether to plot observed data. Default is False.
    :type plot: Bool, optional
    :param fcorner: Bandpass corner frequencies, in seconds. Default is 5 to 120 minutes.
    :type fcorner: array, optional
    :param make_json: Whether to make dart_waves.json from all available stations. Default is False.
    :type make_json: Bool, optional
    """

    hypo_lat = tensor_info["lat"]
    hypo_lon = tensor_info["lon"]
    origin_time = UTCDateTime(tensor_info["datetime"])
    dart_station_file = PROJECT_DIRECTORY / "DART_station_locations.txt"

    data_directory = pathlib.Path(data_directory).resolve()
    if not os.path.exists(f"{data_directory}"):
        os.makedirs(f"{data_directory}")
    df_dart = pd.read_csv(dart_station_file, sep=r"\s+", index_col=False)
    dart_ID = df_dart["Station_ID"]
    longitude = df_dart["Longitude"]
    latitude = df_dart["Latitude"]

    if make_json == True:
        info_traces = []
    for idx, station in enumerate(dart_ID):
        dist_km, azimuth, baz = mng._distazbaz(
            latitude[idx], longitude[idx], hypo_lat, hypo_lon
        )
        dist_degrees = dist_km / 111.11
        if dist_degrees < max_distance:

            # Reading from API
            t_in, data_fil = query_dart_api(
                origin_time,
                str(station),
                dt=dt,
                fcorner=fcorner,
            )
            if len(t_in) > 1:
                ### WRITE TO SAC FILE ###
                new_name = "processed_{}_dart.sac".format(station)
                # print(f"...Data retrieved, writing to {data_directory}/{new_name}")
                tr = Trace(data=data_fil)  # leave as m
                tr.stats.delta = dt
                tr.stats.station = str(station)
                tr.stats.channel = "dart"
                utcdt_origin_time = UTCDateTime(
                    origin_time.year,
                    origin_time.month,
                    origin_time.day,
                    origin_time.hour,
                    origin_time.minute,
                    origin_time.second,
                )
                tr.stats.starttime = (
                    UTCDateTime(origin_time.year, origin_time.month, origin_time.day)
                    + t_in[0]
                )
                sac_dict = {"stla": latitude[idx], "stlo": longitude[idx]}
                tr.stats.sac = sac_dict
                stream = Stream(
                    traces=[
                        tr,
                    ]
                )
                # Trim from one hour before origin time to 7 hours later #
                stream.trim(origin_time - 1 * 3600, origin_time + 7 * 3600)
                tr.data = tr.data - np.mean(tr.data[0:60])
                sac_start = tr.stats.starttime
                stream.write(f"{data_directory}/{new_name}", format="SAC", byteorder=0)
                if make_json == True:
                    n_start_obs = int(np.floor((origin_time - sac_start) / dt))
                    duration = len(tr.data)
                    info = _dict_trace(
                        file=f"{data_directory}/DART/{new_name}",
                        name_station=str(dart_ID[idx]),
                        azimuth=azimuth,
                        channel="dart",
                        distance=dist_degrees,
                        dt=dt,
                        duration=duration,
                        n_start_obs=n_start_obs,
                        trace_weight=1.0,
                        wavelet_weight="0 0 2 2 2 2 2 2\n",
                        observed_trace=None,
                        synthetic_trace=None,
                        location=[latitude[idx], longitude[idx]],
                    )
                    info_traces.append(info)
                if plot == True:
                    ### MAKE PLOT OF OBSERVATIONS ###
                    plt.figure(figsize=(15, 5))
                    plt.plot(
                        tr.times() / 3600 - (origin_time - sac_start) / 3600,
                        tr.data * 100.0,
                    )
                    plt.xlim([-1, 7])
                    plt.xlabel("Time After Origin (hrs)")
                    plt.ylabel("Tsunami Amplitude (cm)")
                    plt.savefig(f"{data_directory}/{station}_dart.png")
                    plt.close()
            else:
                continue

    if make_json == True:
        print(f"Writing to {data_directory}/dart_waves.json")
        with open(data_directory / "dart_waves.json", "w") as f:
            json.dump(
                info_traces,
                f,
                sort_keys=True,
                indent=4,
                separators=(",", ": "),
                ensure_ascii=False,
            )


def make_ffm_kinematic_dtopo(
    ffm_solution, directory: Union[pathlib.Path, str] = pathlib.Path(), makeplots=True
):
    """
    Creates kinematic dtopo file from Solution.txt output

    :param ffm_solution: Solution.txt, finite fault solution file
    :type ffm_solution: str
    :param directory: Working directory
    :type directory: Union[pathlib.Path, str]
    :param makeplots: Whether to create dtopo animation, default is True.
    :type makeplots: Bool, optional
    :return: the filename for the kinamtic dtopo
    :rtype: str
    """
    directory = pathlib.Path(directory)
    fault = dtopotools.Fault()
    fault.rupture_type = "kinematic"
    dtopo_fname = directory / "kinematic_dtopo.tt1"
    plots_out = directory / "kinematic_dtopo_plots"

    source_dx = 2 / 60.0  # dtopo sampling (2 arcmins)
    num_subfaults = 0
    max_stf = 0
    LAT = []
    LON = []
    T_RUP = []
    with open(ffm_solution, "r") as sol:
        for line in sol:
            if line.startswith("#Fault_segment"):
                sf_length = float(line.split()[7].split("km")[0])
                sf_width = float(line.split()[12].split("km")[0])
            elif "#" in line:  # HEADER LINES
                continue
            elif len(np.array(line.split())) < 4:  # FAULT BOUNDARY LINES
                continue
            else:  # ACTUAL SUBFAULT DETAILS
                (
                    lat,
                    lon,
                    dep,
                    slip,
                    rake,
                    strike,
                    dip,
                    t_rup,
                    t_ris,
                    t_fal,
                    mo,
                ) = [float(v) for v in line.split()]
                # Make rake be -180 to 180 (not 0-360)
                if float(rake) > 180:
                    rake = float(rake) - 360

                # save info to find epicenter location
                LAT.append(lat)
                LON.append(lon)
                T_RUP.append(t_rup)

                # make into dtopo format
                subfault = dtopotools.SubFault()
                subfault.length = sf_length * 1000.0  # units in m
                subfault.width = sf_width * 1000.0  # units in m
                subfault.depth = dep * 1000.0  # units in m
                subfault.strike = strike
                subfault.slip = slip / 100.0  # units in m
                subfault.rake = rake
                subfault.dip = dip
                subfault.longitude = lon
                subfault.latitude = lat
                subfault.coordinate_specification = "centroid"
                subfault.rupture_time = t_rup
                subfault.rise_time = t_ris + t_fal

                if num_subfaults == 0:
                    fault.subfaults = [subfault]
                else:
                    fault.subfaults.append(subfault)

                num_subfaults += 1
                max_stf = np.max([max_stf, t_rup + t_ris + t_fal])
    ind_hyp = np.where(T_RUP == np.min(T_RUP))[0][0]
    print(f"{num_subfaults} total subfaults")
    max_stf = np.ceil(max_stf / 10) * 10  # round max_stf up to the nearest 10
    print(f"Max STF: {max_stf}")
    output_times = np.arange(0, int(max_stf) + 10, 10)

    x, y = fault.create_dtopo_xy(dx=source_dx, buffer_size=3)
    fault.create_dtopography(x, y, times=output_times)
    dtopo = fault.dtopo
    dtopo.write(dtopo_fname, dtopo_type=1, dZ_format="%.9f")
    if makeplots == True:
        # Plot dtopo GIF
        plot_dtopo_snapshots(
            dtopo_fname, LON[ind_hyp], LAT[ind_hyp], output_dir=plots_out
        )

    return dtopo_fname


def plot_dtopo_snapshots(
    dtopo_file,
    hyp_lon,
    hyp_lat,
    output_dir="dtopo_pngs",
    gif_name="dtopo_animation.gif",
):
    """
    Reads a dtopo file and saves a PNG for each time step, showing deformation as a function of lon/lat.
    Also creates a GIF of all frames.
    Assumes columns: time, lon, lat, deformation.

    :param dtopo_file: filename of dtopo (seafloor deformation)
    :type dtopo_file: str
    :param hyp_lon: hypocenter longitude
    :type hyp_lon: float
    :param hyp_lat: hypocenter latitude
    :type hyp_lat: float
    :param output_dir: directory to save snapshots. Default is 'dtopo_pngs'.
    :type output_dir: str
    :param gif_name: filename of GIF of kinematic seafloor displacement. Default is 'dtopo_animation.gif'.
    :type gif_name: str
    """

    # Load data
    data = np.loadtxt(dtopo_file)
    times = np.unique(data[:, 0])

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    frame_files = []

    for t in times:
        mask = data[:, 0] == t
        lons = data[mask, 1]
        lats = data[mask, 2]
        deformation = data[mask, 3]

        plt.figure(figsize=(8, 6))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        sc = plt.scatter(
            lons, lats, c=deformation, cmap="seismic", vmin=-4, vmax=4, s=5
        )
        epi = plt.plot(hyp_lon, hyp_lat, "*", c="k")
        plt.colorbar(sc, label="Deformation")
        gl = ax.gridlines(draw_labels=True, linestyle="--", color="gray")

        # Optional: Customize label appearance and visibility
        gl.top_labels = False  # Disable labels at the top
        gl.right_labels = False  # Disable labels on the right
        plt.title(f"Deformation at t={t:.2f} s")
        plt.grid()
        plt.tight_layout()
        frame_path = os.path.join(output_dir, f"dtopo_t{t:.2f}.png")
        plt.savefig(frame_path)
        plt.close()
        frame_files.append(frame_path)

    # Create GIF
    images = [imageio.imread(f) for f in frame_files]
    gif_path = os.path.join(output_dir, gif_name)
    imageio.mimsave(gif_path, images, duration=0.5)


def run_ffm_tsunami(
    directory: Union[pathlib.Path, str],
    ffm_solution: Union[pathlib.Path, str] = "Solution.txt",
    tfinal: int = 10800,
    model_type: str = "kinematic",
    run_geoclaw: bool = True,
):
    """
    Runs tsunami forward model calculation from FFM solution assuming static slip

    :param directory: Directory that holds the FFM solution param file. The default is pathlib.Path().
    :type directory: Union[pathlib.Path, str]
    :param ffm_solution: Path to FFM solution in PARAM format (e.g., Solution.txt)
    :type ffm_solution: Union[pathlib.Path, str]
    :param tfinal: End time for calculating tsunami propagation, in seconds. The default is 10800.
    :type tfinal: int, optional
    :param model_type: Whether to use static or kinematic seafloor deformation, default is kinematic
    :type model_type: str, optional
    :param run_geoclaw: Whether to run the model (otherwise, just make the input files), default is True
    :type run_geoclaw: bool, optional
    """
    # Set up directories to stay organized
    directory = pathlib.Path(directory).resolve()
    working_dir = os.getcwd()
    dart_station_file = PROJECT_DIRECTORY / "DART_station_locations.txt"
    tsunami_output_directory = directory / "_TsunamiOutput"
    if not os.path.exists(tsunami_output_directory):
        os.makedirs(tsunami_output_directory)
    geoclaw_parameter_file = directory / "geoclaw_parameter_file.json"
    with open(directory / "segments_data.json") as s:
        segments_data = json.load(s)
    dart_stations = directory / "dart_waves.json"
    with open(dart_stations) as d:
        dart_data = json.load(d)
    with open(directory / "tensor_info.json") as f:
        tensor_info = json.load(open(directory / "tensor_info.json"))
    #########################################################
    ### CREATE DICT WITH PARAMETERS FOR GEOCLAW_SETRUN.PY ###
    #########################################################
    setrun_parameter = {}
    ### Find DART locations and create "gauge" and "high_res_region" dicts
    if os.path.isfile(directory / "dart_waves.json") == False:
        make_dart_json_no_observations(directory, dart_station_file)
    high_res_regions, gauges, lon_dart, lat_dart = dart_high_res_locations(
        dart_data, segments_data, tensor_info
    )
    ### Find locations of subfaults in model space ###
    lon_sub, lat_sub, dep_sub, _, _, _, _, _ = get_subfault_locations(
        segments_data, tensor_info
    )
    ###################################################################################
    ### Find model region, with 0.5 degrees of space beyond your fault and stations ###
    ###################################################################################
    min_lon = min(min(lon_sub), min(lon_dart)) - 25
    max_lon = max(max(lon_sub), max(lon_dart)) + 25
    min_lat = min(min(lat_sub), min(lat_dart)) - 25
    max_lat = max(max(lat_sub), max(lat_dart)) + 25
    model_region = {
        "west": min_lon,
        "east": max_lon,
        "south": min_lat,
        "north": max_lat,
    }
    clines_zeta = [1, 3, 5, 7, 10, 15, 20, 25, 30, 50]  # max amplitude for plot, in m
    ##################
    ### TOPO/DTOPO ###
    ##################
    if model_type == "kinematic":
        dtopo_fname = make_ffm_kinematic_dtopo(ffm_solution, directory, makeplots=True)
    elif model_type == "static":
        dtopo_fname = make_ffm_static_dtopo(ffm_solution, directory, makeplots=True)
    else:
        raise ValueError("model_type must be either static or kinematic")
    dtopo_type = 1
    topo_path = make_topo(min_lon, max_lon, min_lat, max_lat, directory)
    copy(dtopo_fname, tsunami_output_directory)
    copy(topo_path, tsunami_output_directory)
    ##############################
    ### NUMBER OF OUTPUT TIMES ###
    ##############################
    num_output_times = int(
        np.ceil(tfinal / 60.0 / 5.0)
    )  ### Enough to have one output every 5.0 minutes
    ###################################################
    ### WRITE PARAMETERS FOR SETRUN.PY TO JSON FILE ###
    ###################################################
    setrun_parameter = {
        "model_region": model_region,
        "high_res_regions": high_res_regions,
        "gauges": gauges,
        "topo_path": f"{tsunami_output_directory}/{os.path.basename(topo_path)}",
        "dtopo_path": f"{tsunami_output_directory}/{os.path.basename(dtopo_fname)}",
        "dtopo_type": dtopo_type,
        "dt_max_dtopo": 5,
        "tfinal": tfinal,
        "num_output_times": num_output_times,
        "wave_tolerance": 1.0e-4,
        "clines_zeta": clines_zeta,
    }
    with open(geoclaw_parameter_file, "w") as fout:
        json_dumps_str = json.dumps(setrun_parameter, indent=4)
        print(json_dumps_str, file=fout)
    fout.close()
    ###################################################################
    ### Copy everything to the TsunamiOutput folder and run GeoClaw ###
    ###################################################################
    if run_geoclaw == True:
        copy(geoclaw_parameter_file, tsunami_output_directory)
        copy(PROJECT_DIRECTORY / "src/ffm/geoclaw_setrun.py", tsunami_output_directory)
        copy(PROJECT_DIRECTORY / "src/ffm/Makefile", tsunami_output_directory)
        copy(PROJECT_DIRECTORY / "src/ffm/tsunami_plot.py", tsunami_output_directory)
        chdir(tsunami_output_directory)

        geoclaw_init = "make .exe"
        with open(directory / "output.log", "a") as outfile:
            subprocess.call(
                split(geoclaw_init), stdout=outfile
            )  ### Pipe all the unneccesary output to a file
        outfile.close()

        subprocess.call(split(geoclaw_init))
        print("RUN GEOCLAW")
        geoclaw_plot = "make .plots"
        with open(directory / "output.log", "a") as outfile:
            subprocess.call(split(geoclaw_plot), stdout=outfile)
        outfile.close()
        chdir(working_dir)
        ### Make output plots comparing synthetic waveforms to data ###
        num_dart = len(dart_data)
        fig = plt.figure(figsize=(10, 10))
        for kdart in range(num_dart):
            # load in dart output
            id_dart = dart_data[kdart]["name"]
            dart_out = np.genfromtxt(
                f"{tsunami_output_directory}/_output/gauge{id_dart}.txt"
            )
            time_out = dart_out[:, 1]
            eta_out = dart_out[:, 5]
            ax = fig.add_subplot(num_dart, 1, kdart + 1)
            ax.plot(time_out / 3600, eta_out)
            ax.set_ylabel(id_dart)
        plt.savefig(f"{directory}/dart_fwdmodel.png")


def make_ffm_static_dtopo(ffm_solution, directory, makeplots=True):
    """
    Create dtopo data file for deformation of seafloor due to earthquake.
    Uses the Okada model with fault parameters and mesh specified below.

    :param ffm_solution: FFM solution file in PARAM format.
    :type ffm_solution: str
    :param makeplots: Whether to make plots of seafloor deformation. The default is True.
    :type makeplots: bool, optional
    :return: Filename where dtopo is saved.
    :rtype: str
    """

    directory = pathlib.Path(directory)
    dtopo_fname = directory / "static_dtopo.tt1"
    plots_out = directory / "static_dtopo_plots"
    fault = dtopotools.Fault()
    fault.subfaults = []
    fault.rupture_type = "static"
    source_dx = 2 / 60.0  # dtopo spatial sampling (2 arcmin)
    event_m0 = 0
    num_subfaults = 0
    with open(ffm_solution, "r") as sol:
        for line in sol:
            if line.startswith("#Fault_segment"):
                sf_length = float(line.split()[7].split("km")[0])
                sf_width = float(line.split()[12].split("km")[0])
            elif "#" in line:  # HEADER LINES
                continue
            elif len(np.array(line.split())) < 4:  # FAULT BOUNDARY LINES
                continue
            else:  # ACTUAL SUBFAULT DETAILS
                (
                    lat,
                    lon,
                    dep,
                    slip,
                    rake,
                    strike,
                    dip,
                    t_rup,
                    t_ris,
                    t_fal,
                    mo,
                ) = [float(v) for v in line.split()]
                # Make rake be -180 to 180 (not 0-360)
                if float(rake) > 180:
                    rake = float(rake) - 360
                event_m0 += mo

                subfault = dtopotools.SubFault()
                subfault.length = sf_length * 1000.0  # units in m
                subfault.width = sf_width * 1000.0  # units in m
                subfault.depth = dep * 1000.0  # units in m
                subfault.strike = strike
                subfault.slip = slip / 100.0  # units in m
                subfault.rake = rake
                subfault.dip = dip
                subfault.longitude = lon
                subfault.latitude = lat
                subfault.coordinate_specification = "centroid"
                subfault.rupture_time = t_rup
                subfault.rise_time = t_ris + t_fal

                if num_subfaults == 0:
                    fault.subfaults = [subfault]
                else:
                    fault.subfaults.append(subfault)
                num_subfaults += 1

    x, y = fault.create_dtopo_xy(dx=source_dx, buffer_size=3)
    fault.create_dtopography(x, y, times=[1.0])
    dtopo = fault.dtopo
    dtopo.write(dtopo_fname, dtopo_type=1)

    event_Mw = (2.0 / 3) * (np.log10(event_m0 * 1e-7) - 9.1)
    print(f"Event Mw: {event_Mw}")

    # Make into dtopo format
    if makeplots == True:
        os.makedirs(plots_out, exist_ok=True)
        print("Plotting Seafloor Deformation")
        fig = plt.figure(figsize=(12, 5))
        fault.plot_subfaults(plot_rake=False, slip_color=True)
        plt.savefig(plots_out / "FFM.png")
        print("The seismic moment is %g N-m" % fault.Mo())
        print("The Moment magnitude is %g" % fault.Mw())

        fig = plt.figure(figsize=(12, 5))
        dtopo.plot_dZ_colors(t=1)
        plt.savefig(plots_out / "SeafloorDeformation.png")

    return dtopo_fname


def dart_high_res_locations(dart_data, segments_data, tensor_info):
    """
    Determine regions around each DART station for higher resolution in GeoClaw

    :param dart_data: Dictionary information from dart_waves.json
    :type dart_data: dict
    :param segments_data: Dictionary information from segments_data.json
    :type segments_data: dict
    :return: high-resolution areas around each DART, DART station info, DART longitudes, DART latitudes
    :rtype: Tuple[dict, dict, list, list]
    """

    ############################################################################
    ### Set high-resolution spots around the locations of your dart stations ###
    ############################################################################
    high_res_regions = {}
    gauges = {}
    ### READ IN DART STATION INFO ###
    #### Add dart info to json ###
    num_dart = len(dart_data)
    lon_dart = []
    lat_dart = []
    for kdart in range(num_dart):
        lat_dart.append(float(dart_data[kdart]["location"][0]))
        lon_dart.append(float(dart_data[kdart]["location"][1]))
        id_dart = dart_data[kdart]["name"]
        if id_dart.isalpha():
            print(f"Alpha DART: {id_dart}")
            alpha_ind = ord(id_dart[-1]) - 64  # Convert NZ[A-L] to number
            id_dart = str(alpha_ind).rjust(5, "0")
            print(f"... becomes Number DART: {id_dart}")
        region = {
            "minlevel": 3,
            "maxlevel": 3,
            "t1": 0.0,
            "t2": 10000.0,
            "x1": lon_dart[kdart] - 0.5,
            "x2": lon_dart[kdart] + 0.5,
            "y1": lat_dart[kdart] - 0.5,
            "y2": lat_dart[kdart] + 0.5,
        }
        dart_info = {
            "gauge_id": id_dart,
            "lon": lon_dart[kdart],
            "lat": lat_dart[kdart],
            "t1": 0,
            "t2": 1e10,
        }
        high_res_regions.update({str(kdart): region})
        gauges.update({str(kdart): dart_info})

    ### Add high-res region around fault itself ###
    subfault_lons, subfault_lats, _, _, _, _, _, _ = get_subfault_locations(
        segments_data, tensor_info
    )
    fault_region = {
        "minlevel": 3,
        "maxlevel": 3,
        "t1": 0.0,
        "t2": 300.0,
        "x1": np.min(subfault_lons) - 0.2,
        "x2": np.max(subfault_lons) + 0.2,
        "y1": np.min(subfault_lats) - 0.2,
        "y2": np.max(subfault_lats) + 0.2,
    }
    high_res_regions.update({str(num_dart): fault_region})
    return high_res_regions, gauges, lon_dart, lat_dart


def make_topo(
    west: float,
    east: float,
    south: float,
    north: float,
    directory: Union[pathlib.Path, str],
):
    """
    Makes topo file for GeoClaw

    :param west: western limit of the study region, degrees longitude
    :type west: float
    :param east: eastern limit of the study region, degrees longitude
    :type east: float
    :param south: southern limit of the study region, degrees latitude
    :type south: float
    :param north: northern limit of the study region, degrees latitude
    :type north: float
    :param directory: directory where modeling is being performed
    :type directory: Union[pathlib.Path, str]
    :return: GeoClaw topo file in ESRI Arc/Info ASCII format
    :rtype: string
    """
    directory = pathlib.Path(directory)
    topogrd = pathlib.Path(os.environ["CLAW"]) / "topo30.grd"
    crop_grd = directory / "topo_crop.nc"
    ei_topo_out = directory / "topo.txt"
    ### Round boundaries to the nearest 10th of a degree
    west = np.floor(west * 10) / 10.0
    east = np.ceil(east * 10) / 10.0
    south = np.floor(south * 10) / 10.0
    north = np.ceil(north * 10) / 10.0
    # print(west, east, south, north)
    ### Use gmt grdcut to crop topo30.grd to region of interest
    grdcut_command = (
        f"gmt grdcut {topogrd} -G{crop_grd} -R{west}/{east}/{south}/{north}"
    )
    os.system(grdcut_command)
    ### Use grdconvert to convert format from grd to ESRI Arc/Info ASCII Grid Interchange format
    grdconvert_command = f"gmt grdconvert {crop_grd}=nf -G{ei_topo_out}=ei"
    os.system(grdconvert_command)
    ### Invert the columns of the header to have number before description
    topo = topotools.Topography(ei_topo_out, topo_type=2)
    topo_crop_out = directory / "topo.topo"
    topo.write(topo_crop_out, no_data_value=-999999)
    return topo_crop_out


def dart_gf_time_freq_domain(
    station,
    dt,
    subfault,
    time_out,
    eta,
    gf_file_freq,
    GF_directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """
    Computes text files of Green's function in both time and frequency domain

    :param station: DART station ID.
    :type station: str
    :param dt: Sampling interval of Green's function.
    :type dt: float
    :param subfault: Subfault index.
    :type subfault: int
    :param time_out: Time vector for Green's functions.
    :type time_out: array
    :param eta: Amplitude vector for Green's functions.
    :type eta: array
    :param gf_file_time: Filename for time domain Green's functions to be written to.
    :type gf_file_freq: str
    :param GF_directory: Directory where Green's functions should be written. The default is pathlib.Path().
    :type GF_directory: Union[pathlib.Path, str], optional
    """

    # If this is the first subfault and a frequency GF file already exists (i.e., from a previous run) delete it.
    if subfault == 0:
        if os.path.exists(gf_file_freq):
            os.remove(gf_file_freq)
    # Write out to file
    with open(gf_file_freq, "a") as out_freq:
        waveform = np.array(eta, dtype=complex)

        fft_waveform = np.fft.fft(waveform, n=1024)
        out_freq.write("{} 256\n".format(int(subfault + 1)))
        for value in fft_waveform[:256]:
            real = np.real(value)
            imag = np.imag(value)
            out_freq.write("{} {}\n".format(real, imag))
    out_freq.close()

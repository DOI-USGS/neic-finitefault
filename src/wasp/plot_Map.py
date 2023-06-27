# -*- coding: utf-8 -*-
"""
Map plot with PyGMT
"""
import argparse
from matplotlib import pyplot as plt
from matplotlib import gridspec, ticker, patches, colors
import matplotlib
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cf
from obspy.imaging.beachball import beach
import numpy as np
import os
import get_outputs
import load_ffm_model
import json
import errno
from datetime import datetime
from shutil import copy2, move
from glob import glob
from cartopy.io.img_tiles import Stamen
from matplotlib.patches import Rectangle
import pandas as pd
#from clawpack.geoclaw import dtopotools
#
# local modules
#
import fault_plane as pf
import velocity_models as mv
import plane_management as pl_mng
import seismic_tensor as tensor
import shakemap_tools as shakemap
from waveform_plots_NEIC import plot_waveform_fits
from plot_maps_NEIC import plot_map, set_map_cartopy, plot_borders
from matplotlib import cm
from matplotlib.colors import ListedColormap
from static2fsp import static_to_fsp
from plot_graphic_NEIC import __redefine_lat_lon


def _PlotMap(tensor_info, segments, point_sources, solution, default_dirs, convex_hulls=[], 
            files_str=None, stations_gps=None, stations_cgps=None, option='Solucion.txt',
            max_slip=None, legend_len=None, scale=None, limits=[None,None,None,None]):

    import pygmt
    import numpy as np
    import collections

    ################################
    ### GET DESIRED PLOT REGION ####
    ################################

    lat0 = tensor_info['lat']
    lon0 = tensor_info['lon']

    segments_lats, segments_lons, segments_deps = __redefine_lat_lon(segments, point_sources)
    min_lats = [min(segment_lat.flatten()) for segment_lat in segments_lats]
    max_lats = [max(segment_lat.flatten()) for segment_lat in segments_lats]
    min_lons = [min(segment_lon.flatten()) for segment_lon in segments_lons]
    max_lons = [max(segment_lon.flatten()) for segment_lon in segments_lons]
    min_lat = np.min(min_lats)# - 0.5
    max_lat = np.max(max_lats)# + 0.5
    min_lon = np.min(min_lons)# - 0.5
    max_lon = np.max(max_lons)# + 0.5

    if files_str is not None:
        for file in files_str:
            name = file['name']
            latp, lonp = file['location']
            min_lat = min(min_lat, latp)
            max_lat = max(max_lat, latp)
            min_lon = min(min_lon, lonp)
            max_lon = max(max_lon, lonp)
    if stations_cgps is not None:
        for file in stations_cgps:
            name = file['name']
            latp, lonp = file['location']
            min_lat = min(min_lat, latp)
            max_lat = max(max_lat, latp)
            min_lon = min(min_lon, lonp)
            max_lon = max(max_lon, lonp)
    if stations_gps is not None:
        max_obs = np.zeros(3)
        stations_gps2 = []
        for name, sta_lat, sta_lon, obs, syn, error in stations_gps:
            min_lat = min(min_lat, sta_lat)
            max_lat = max(max_lat, sta_lat)
            min_lon = min(min_lon, sta_lon)
            max_lon = max(max_lon, sta_lon)
            stations_gps2 = stations_gps2\
                + [[name, sta_lat, sta_lon, obs, syn, error]]
            max_obs = np.maximum([abs(float(v)) for v in obs], max_obs)
        max_obs = np.max(max_obs)
        if legend_len==None:
            if max_obs < 5:
                legend_len = 1
            elif max_obs < 10:
                legend_len = 5
            elif max_obs < 20:
                legend_len = 10
            else:
                legend_len = 20
        if scale==None:
            scale=2
        max_obs = max_obs/scale


    if limits == [None,None,None,None]:
        region = [min_lon - 0.5, max_lon + 0.5, min_lat - 0.5, max_lat + 0.5]
    else:
        region = [limits[0], limits[1], limits[2], limits[3]]

    ################################
    ### PLOT BASEMAP/ COASTLINES ###
    ################################

    fig = pygmt.Figure()
    resolution='03s'
    map_scale_len='50k'
    if region[1]-region[0] > 5 or region[3]-region[2] > 5:
        resolution='30s'
        map_scale_len='100k'
    if region[1]-region[0] > 10 or region[3]-region[2] > 10:
        resolution='01m'
        map_scale_len='200k'
    grid = pygmt.datasets.load_earth_relief(resolution=resolution, region=region)
    projection = "M15c"
    map_scale = "g"+str(region[0])+"/"+str(region[2])+"+c17.40+w"+map_scale_len+"+ar+l+jBL+o0.5/0.5+f"

    pygmt.config(MAP_FRAME_TYPE="plain")
    fig.basemap(region=region, projection=projection, frame=["WSne", "afg"])
    fig.grdimage(grid=grid, cmap="oleron", shading=True, transparency=20)
    fig.colorbar(yshift='a-10p',
                 xshift='a-5p',
                 position="n0.05/-0.1+jTL+w100p/8%+h",
                 frame="x+lTopography (km)",
                 #box="+p2p,black+ggray80",
                 scale=0.001
                 )
    fig.plot(str(default_dirs['root_dir'])+'/pb2002_boundaries.gmt', style="f10/3p", region=region, pen="2p,black")
    fig.plot(str(default_dirs['root_dir'])+'/pb2002_boundaries.gmt', style="f10/3p", region=region, pen="1p,white")
    fig.basemap(region=region, projection=projection, frame="ag1", map_scale=map_scale)
    fig.coast(resolution="h", shorelines=True)

    faults = glob('*.fault')
    if len(faults) > 0:
        for kfault in range(len(faults)):
            print('...Adding fault trace from: '+str(faults[kfault]))
            fault_trace = np.genfromtxt(faults[kfault])
            fig.plot(x=fault_trace[:,0],y=fault_trace[:,1],pen='2p,black')
            fig.plot(x=fault_trace[:,0],y=fault_trace[:,1],pen='1p,white')

    ###############################
    ### PLOT FINITE FAULT MODEL ###
    ###############################
    segments_data = json.load(open('segments_data.json'))
    segments = segments_data['segments']
    solution = get_outputs.read_solution_static_format(segments)
    plane_info = segments[0]
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
            = pl_mng.__unpack_plane_data(plane_info)

    segments_lats, segments_lons, segments_deps = __redefine_lat_lon(segments, point_sources)
    slip=solution['slip']


    ### MAKE CPT OF SLIP ###
    if max_slip==None:
        maxslip=0
        for segment in range(len(segments_lats)):
            slips = slip[segment].flatten()
            maxslip = np.max([maxslip,np.array(slips).max()/100])
    else:
        maxslip=max_slip
    pygmt.makecpt(cmap=str(default_dirs['root_dir'])+'/python_code/fault2.cpt', series=[0, maxslip])

    plane_info = segments[0]
    strike = plane_info['strike']
    dip = plane_info['dip']

    rupture_vel = segments[0]['rupture_vel']
    subfaults = {'delta_strike': delta_strike, 'delta_dip': delta_dip}
    subfaults2 = pf._point_sources_def(rise_time, rupture_vel, subfaults)
    strike_ps = int(subfaults2['strike_ps'] / 2)
    dip_ps = int(subfaults2['dip_ps'] / 2)
    latitudes = [ps_segment[:, :, dip_ps, strike_ps, 0] for ps_segment in point_sources]
    longitudes = [ps_segment[:, :, dip_ps, strike_ps, 1] for ps_segment in point_sources]
    
    for segment in range(len(segments_lats)):
        plane_info = segments[segment]
        strike = plane_info['strike']
        dip = plane_info['dip']
        lons = longitudes[segment].flatten()
        lats = latitudes[segment].flatten()
        slips = slip[segment].flatten()
        fig.plot(x=lons,
             y=lats,
             style="J"+str(strike)+"/"+str(delta_strike)+"/"+str(delta_dip*np.cos(np.radians(dip))),
             cmap=True,
             color=slips/100.
             )

    fig.colorbar(yshift='a-10p',
                 xshift='a135p',
                 position="n0.05/-0.1+jTL+w100p/8%+h",
                 frame="x+lSlip (m)",
                 #box="+p2p,black+ggray80"
                 )

    ### PLOT AFTERSHOCKS OVER TOP ###
    aftershocks = glob('*aftershock*')
    if len(aftershocks) > 0:
        for kafter in range(len(aftershocks)):
            print('...Adding aftershocks from: '+str(aftershocks[kafter]))
            aftershock = np.genfromtxt(aftershocks[kafter], delimiter='\t', skip_header=1)
            aftershock_lat = aftershock[:,3]
            aftershock_lon = aftershock[:,4]
            aftershock_mag = aftershock[:,6]
            fig.plot(x=aftershock_lon, y=aftershock_lat,style="cp", size=aftershock_mag, color='grey', pen='black', transparency=60)
    
    #################################
    # outline the segments in black #
    #################################
    for segment in range(len(segments_lats)):
        depths = segments_deps[segment].flatten()
        min_lats_idx = np.where(segments_lats[segment].flatten()==min_lats[segment])[0][0]
        cornerA = [segments_lons[segment].flatten()[min_lats_idx], min_lats[segment], depths[min_lats_idx], 0]
        min_lons_idx = np.where(segments_lons[segment].flatten()==min_lons[segment])[0][0]
        cornerB = [min_lons[segment], segments_lats[segment].flatten()[min_lons_idx], depths[min_lons_idx], 1]
        max_lats_idx = np.where(segments_lats[segment].flatten()==max_lats[segment])[0][-1]
        cornerC = [segments_lons[segment].flatten()[max_lats_idx], max_lats[segment], depths[max_lats_idx], 2]
        max_lons_idx = np.where(segments_lons[segment].flatten()==max_lons[segment])[0][-1]
        cornerD = [max_lons[segment], segments_lats[segment].flatten()[max_lons_idx], depths[max_lons_idx], 3]
        corners = np.c_[cornerA, cornerB, cornerC, cornerD]
        max_dep = max(corners[2,:])
        idx_max = np.where(np.array(corners[2,:]) == max_dep)[0][0]
        min1 = min2 = corners[:,idx_max]
        for i in range(len(corners)):
            if corners[2,i] < min1[2]:
                min2 = min1
                min1 = corners[:,i]
            elif corners[2,i] < min2[2] and collections.Counter(corners[:,i]) != collections.Counter(min1):
                min2 = corners[:,i]

        updip = np.c_[min1,min2]
        corners = np.c_[corners,cornerA]

        fig.plot(x=corners[0,:], y=corners[1,:], pen="1p,black")
        fig.plot(x=updip[0,:], y=updip[1,:], pen="1p,red")

        #plot multiple segment hypocenters if any
        fig.plot(x=lon0, y=lat0, style='a7p', color="white", pen="black")
        if 'hypocenter' in segments[segment]:
            hyp = segments[segment]['hypocenter']
            fig.plot(x=hyp['lon'], y=hyp['lat'], style='a7p', color="white", pen="black")

    ###########################
    ### PLOT LOCAL STATIONS ###
    ###########################

    ### STRONG-MOTION ACCELEROMETER ###
    if files_str is not None:
        print('...Plotting Strong Motion Stations')
        for file in files_str:
            name = file['name']
            latp, lonp = file['location']
            fig.plot(x=lonp, y=latp, style="i10p", color="white", pen="black")
            fig.text(x=lonp, y=latp, text=name)
        ### ADD TO LEGEND ###
        fig.plot(x=region[1], y=region[2],
            xshift="a-130p", yshift="a-54p",
            no_clip=True, style="i10p", color="white", pen="black")
        fig.text(x=region[1], y=region[2], text="Accelerometer",
            xshift="a-123p", yshift="a-55p", no_clip=True, justify="ML")

    ### HIGH-RATE GNSS ###
    if stations_cgps is not None:
        print('...Plotting cGPS Stations')
        for file in stations_cgps:
            name = file['name']
            latp, lonp = file['location']
            fig.plot(x=lonp, y=latp, style="t10p", color="navy", pen="black")
        ### ADD TO LEGEND ###
        fig.plot(x=region[1], y=region[2],
            xshift="a-55p", yshift="a-56p",
            no_clip=True, style="t10p", color="navy", pen="black")
        fig.text(x=region[1], y=region[2], text="HR GNSS",
            xshift="a-48p", yshift="a-55p", no_clip=True, justify="ML", offset=0/10)

    ### STATIC GNSS ####
    if stations_gps is not None:
        print('...Plotting Static GPS Stations')
        for name, sta_lat, sta_lon, obs, syn, error in stations_gps2:
            gps_z_obs, gps_n_obs, gps_e_obs = obs
            gps_z_syn, gps_n_syn, gps_e_syn = syn
            err_z, err_n, err_e = error

            staticv_obs = pd.DataFrame(
                 data={
                     "x": [sta_lon],
                     "y": [sta_lat],
                     "east_velocity": [float(gps_e_obs)/max_obs],
                     "north_velocity": [float(gps_n_obs)/max_obs],
                     "east_sigma": [float(err_e)/max_obs],
                     "north_sigma": [float(err_n)/max_obs],
                     "correlation_EN": [0.0],
                 }
             )
            v_obs = "0.3c+a45+p0.5p,grey+e+h0.5+gblack"
            # Plot thick white arrow behind, to get white outline on black arrow
            fig.velo(
                data=staticv_obs,
                pen="0.07c,grey",
                line="grey",
                color="grey",
                spec="e1/0",
                vector=v_obs
            )
            fig.velo(
                data=staticv_obs,
                pen="0.05c,black",
                line="BLACK",
                color="BLACK",
                spec="e1/0.34",
                vector=v_obs
            )

            staticv_syn = pd.DataFrame(
                  data={
                      "x": [sta_lon],
                      "y": [sta_lat],
                      "east_velocity": [float(gps_e_syn)/max_obs],
                      "north_velocity": [float(gps_n_syn)/max_obs],
                      "east_sigma": [0.0],
                      "north_sigma": [0.0],
                      "correlation_EN": [0.0],
                  }
              )
            v_syn = "0.3c+p0.5p,black+e+h0.5+gred"
            # Plot thick black arrow behind, to get black outline on red arrow
            fig.velo(
                data=staticv_syn,
                pen=".05c,black",
                line="black",
                color="black",
                spec="e1/0",
                vector=v_syn
            )
            #overlay thin red arrow
            fig.velo(
                data=staticv_syn,
                pen=".03c,red",
                line="red",
                color="red",
                spec="e1/0",
                vector=v_syn
            ) 
        ### ADD TO LEGEND ###
        fig.text(x=region[1], y=region[2], text="Observed GNSS",
            xshift="a-133p", yshift="a-30p", no_clip=True, justify="ML")
        fig.text(x=region[1], y=region[2], text="Synthetic  GNSS",
            xshift="a-133p", yshift="a-40p", no_clip=True, justify="ML")
        static_legend = pd.DataFrame(
             data={
                 "x": [region[1]],
                 "y": [region[2]],
                 "east_velocity": [legend_len/max_obs],
                 "north_velocity": [0],
                 "east_sigma": [legend_len/max_obs/10],
                 "north_sigma": [legend_len/max_obs/10],
                 "correlation_EN": [0],
                 }
             )
        # Plot thick white arrow behind, to get white outline on black arrow

        fig.velo(
                data=static_legend,
                pen="0.07c,grey",
                line="grey",
                color="grey",
                spec="e1/0",
                vector=v_obs,
                xshift="a-45p", yshift="a-30p",
                no_clip=True
            )
        fig.velo(
                data=static_legend,
                pen="0.05c,black",
                line="BLACK",
                color="BLACK",
                spec="e1/0.34",
                vector=v_obs,
                xshift="a-45p", yshift="a-30p",
                no_clip=True
            )
        # Plot thick black arrow behind, to get black outline on red arrow
        fig.velo(
                data=static_legend,
                pen=".05c,black",
                line="black",
                color="black",
                spec="e1/0",
                vector=v_syn,
                xshift="a-45p", yshift="a-40p",
                no_clip=True
            )
        #overlay thin red arrow
        fig.velo(
                data=static_legend,
                pen=".03c,red",
                line="red",
                color="red",
                spec="e1/0",
                vector=v_syn,
                xshift="a-45p", yshift="a-40p",
                no_clip=True
            )
        fig.text(x=region[1], y=region[2], text=str(legend_len*10)+"+/-"+str(legend_len)+" mm",
            xshift="a-60p", yshift="a-20p", no_clip=True, justify="ML")


    fig.savefig("PyGMT_Map.png")


if __name__ == '__main__':
    """
    """
    import management as mng
    import eventpage_downloads as dwnlds
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
    parser.add_argument("-ffms", "--ffm_solution", action="store_true",
                        help="plot FFM solution slip maps, rise time")
    parser.add_argument("-t", "--tele", action="store_true",
                        help="plot misfit of teleseismic data")
    parser.add_argument("-su", "--surface", action="store_true",
                        help="plot misfit of surface waves data")
    parser.add_argument("-st", "--strong", action="store_true",
                        help="plot strong motion stations and strong motion misfit")
    parser.add_argument("--cgps", action="store_true",
                        help="plot misfit of cGPS data")
    parser.add_argument("--gps", action="store_true", help="plot GPS data")
    parser.add_argument(
        "-in", "--insar", action="store_true", help="plot InSar data")
    parser.add_argument("-o","--option",choices=['autoscale','noscale'],
                        help="choose whether Rupture plot needs to be scaled")
    parser.add_argument("-mr","--mrtime", default='0', type=float,
                        help="choose cutoff time for Moment Rate plot")
    parser.add_argument("-shakemap","--shakemappolygon",action="store_true",
                        help="create shakemap_polygon.txt")
    parser.add_argument("-ev","--EventID", nargs='?', const='Not Provided', type=str,
                        help="Provide event ID")
    parser.add_argument("-d","--downloads",action="store_true",
                        help="create event page download files")
    parser.add_argument("-pub","--publish",action="store_true",
                        help="rename files for use with sendproduct")
    parser.add_argument("-check","--checkerboard",action="store_true",
                        help="plot comparison for checkerboard test")
    parser.add_argument("-max","--maxvalue",default=None, type=float,
                        help="Choose maximum slip value for plot")
    parser.add_argument("-legend","--legend_length",default=None, type=float,
                        help="Length of static GNSS vector for legend (in cm)")
    parser.add_argument("-scale","--scale_factor",default=None, type=float,
                        help="Scale factor for static GNSS vector lengths (larger means shorter vectors)")
    parser.add_argument("-limits","--map_limits", type=float, nargs=4,
                        help="Specify map limits [W,E,N,S] from edges of plotted features. eg: 0.5 0.5 0.5 0.5\
                        gives a 0.5 degree buffer on each side. Negative numbers will cut off plotted features.")
    args = parser.parse_args()
    os.chdir(args.folder)
    used_data = []
    used_data = used_data + ['strong_motion'] if args.strong else used_data
    used_data = used_data + ['cgps'] if args.cgps else used_data
    used_data = used_data + ['tele_body'] if args.tele else used_data
    used_data = used_data + ['surf_tele'] if args.surface else used_data
    default_dirs = mng.default_dirs()
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    segments_data = json.load(open('segments_data.json'))
    segments = segments_data['segments']
    rise_time = segments_data['rise_time']
    connections = None
    if 'connections' in segments_data:
        connections = segments_data['connections']
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections)
    solution = get_outputs.read_solution_static_format(segments)

    traces_info, traces_info_cgps, stations_gps = [None, None, None]
    if args.gps:
        names, lats, lons, observed, synthetic, error\
                = get_outputs.retrieve_gps()
        stations_gps = zip(names, lats, lons, observed, synthetic, error)
    if args.cgps:
        traces_info_cgps = json.load(open('cgps_waves.json'))
    if args.strong:
        traces_info = json.load(open('strong_motion_waves.json'))
    if args.maxvalue != None:
        maxval=args.maxvalue
    else:
        maxval=None
    if args.legend_length != None:
        legend_len=args.legend_length
    else:
        legend_len=None
    if args.scale_factor != None:
        scale=args.scale_factor
    else:
        scale=None
    if args.map_limits:
        limits=[args.map_limits[0],args.map_limits[1],args.map_limits[2],args.map_limits[3]]
        print(f'Axes limits: {limits}')
    else:
        limits=[None,None,None,None]
    if args.ffm_solution:
        if not os.path.isfile('velmodel_data.json'):
            vel_model = mv.select_velmodel(tensor_info, default_dirs)
        else:
            vel_model = json.load(open('velmodel_data.json'))
        shear = pf.shear_modulous(point_sources, velmodel=vel_model)
        if args.option == 'autoscale':
            autosize = True
        else:
            autosize = False
        if args.mrtime > 0:
            mr_time = args.mrtime
        else:
            mr_time = None
        if args.EventID:
            evID = args.EventID
        else:
            evID = None

        _PlotMap(tensor_info, segments, point_sources, solution, default_dirs, convex_hulls=[],
            files_str=traces_info, stations_gps=stations_gps, stations_cgps=traces_info_cgps, option='Solucion.txt',
            max_slip=maxval, legend_len=legend_len, scale=scale, limits=limits)

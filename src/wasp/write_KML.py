#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 11:21:45 2022

@author: degoldberg
"""
import numpy as np
import collections
from plot_graphic_NEIC import __redefine_lat_lon
import plane_management as pl_mng
import argparse
import os
import get_outputs
import json
import fault_plane as pf
import management as mng
import seismic_tensor as tensor
import cartopy.crs as ccrs
from matplotlib import pyplot as plt
import cartopy.feature as cf
import cartopy.io.shapereader as shpreader
from matplotlib.patches import Rectangle
from plot_maps_NEIC import plot_map, plot_borders
from matplotlib import cm 
from matplotlib.colors import ListedColormap
import cartopy
from glob import glob

"""
Set colorbar for slip
"""
rm = 100 # amount of lines to remove on black end of magma_r
ad = 50 # how much at the zero end should be *just* white before transitioning to meet colors
magma_cpt = cm.get_cmap('magma_r', 512) #start with magma_r
white_bit = np.array([255/256, 250/256, 250/256, 1]) #create array of white
slip_cpt = magma_cpt(np.linspace(0,1,512)) #initialize slip_cpt
slip_cpt[rm:,:] = slip_cpt[0:-rm,:] #move beginning up to remove black end
r_s = np.linspace(white_bit[0],slip_cpt[rm][0],rm-ad) # gradient from white to beginning of new magma
g_s = np.linspace(white_bit[1],slip_cpt[rm][1],rm-ad)
b_s = np.linspace(white_bit[2],slip_cpt[rm][2],rm-ad)
slip_cpt[ad:rm,:][:,0] = r_s
slip_cpt[ad:rm,:][:,1] = g_s
slip_cpt[ad:rm,:][:,2] = b_s
slip_cpt[:ad,:] = white_bit
slipcpt = ListedColormap(slip_cpt)
#####

def _write_KML(segments, point_sources, evID, margins, directory=None):
    print('Writing KML file...')
    if directory == None:
        directory = os.getcwd()

    kmlfile = str(evID)+'.kml'    
    kml = open(directory + '/' + kmlfile, 'w')
    kml.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    kml.write("<kml xmlns=\"http://earth.google.com/kml/2.2\">\n")
    kml.write("<Document>\n")
    kml.write("	<name>"+str(evID)+"</name>\n")
    kml.write("	<StyleMap id=\"msn_ylw-pushpin\">\n")
    kml.write("	<Pair>\n")
    kml.write("	<key>normal</key>\n")
    kml.write("	<styleUrl>#sn_ylw-pushpin</styleUrl>\n")
    kml.write("		</Pair>\n")
    kml.write("		<Pair>\n")
    kml.write("		<key>highlight</key>\n")
    kml.write("	<styleUrl>#sh_ylw-pushpin</styleUrl>\n")
    kml.write("		</Pair>\n")
    kml.write("	</StyleMap>\n")
    kml.write("	<Style id=\"sn_ylw-pushpin\">\n")
    kml.write("		<IconStyle>\n")
    kml.write("			<scale>1.1</scale>\n")
    kml.write("			<Icon>\n")
    kml.write("				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>\n")
    kml.write("			</Icon>\n")
    kml.write("			<hotSpot x=\"20\" y=\"2\" xunits=\"pixels\" yunits=\"pixels\"/>\n")
    kml.write("		</IconStyle>\n")
    kml.write("		<LineStyle>\n")
    kml.write("		<color>ff00ffff</color>\n")
    kml.write("	<width>5</width>\n")
    kml.write("		</LineStyle>\n")
    kml.write("	</Style>\n")
    kml.write("	<Style id=\"sh_ylw-pushpin\">\n")
    kml.write("		<IconStyle>\n")
    kml.write("		<scale>1.3</scale>\n")
    kml.write("	<Icon>\n")
    kml.write("				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>\n")
    kml.write("		</Icon>\n")
    kml.write("			<hotSpot x=\"20\" y=\"2\" xunits=\"pixels\" yunits=\"pixels\"/>\n")
    kml.write("		</IconStyle>\n")
    kml.write("	<LineStyle>\n")
    kml.write("	<color>ff00ffff</color>\n")
    kml.write("			<width>5</width>\n")
    kml.write("	</LineStyle>\n")
    kml.write("	</Style>\n")

    segments_lats, segments_lons, segments_deps = __redefine_lat_lon(segments, point_sources)
    min_lats = [min(segment_lat.flatten()) for segment_lat in segments_lats]
    max_lats = [max(segment_lat.flatten()) for segment_lat in segments_lats]
    min_lons = [min(segment_lon.flatten()) for segment_lon in segments_lons]
    max_lons = [max(segment_lon.flatten()) for segment_lon in segments_lons]

    # find outline of segments #
    segment_number=0
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

        #updip = np.c_[min1,min2]
        corners = np.c_[corners,cornerA]
        # write outline of segments to kml file: #
        kml.write(" <Placemark>\n")
        kml.write("		<name>Fault Segment "+str(segment_number)+"</name>\n")
        kml.write("		<styleUrl>#msn_ylw-pushpin</styleUrl>\n")
        kml.write("		<LineString>\n")
        kml.write("			<tessellate>1</tessellate>\n")
        kml.write("			<coordinates>\n")
        kml.write(" "+str(corners[0,0])+","+str(corners[1,0])+","+str(corners[2,0])+"   "+str(corners[0,1])+","+str(corners[1,1])+","+str(corners[2,1])+"   "+str(corners[0,2])+","+str(corners[1,2])+","+str(corners[2,2])+"   "+str(corners[0,3])+","+str(corners[1,3])+","+str(corners[2,3])+"   "+str(corners[0,4])+","+str(corners[1,4])+","+str(corners[2,4]) + "   </coordinates>\n")
        kml.write("		</LineString>\n")
        kml.write("	</Placemark>\n")
        segment_number+=1

    kml.write("     <GroundOverlay>\n")
    kml.write("             <name>"+str(evID)+" FFM Image Overlay</name>\n")
    kml.write("             <Icon>\n")
    kml.write("                     <href>Map_kml.png</href>\n")
    kml.write("                     <viewBoundScale>0.75</viewBoundScale>\n")
    kml.write("             </Icon>\n")
    kml.write("		<color>bfffffff</color>\n")
    kml.write("             <LatLonBox>\n")
    kml.write("                     <north>"+str(margins[3])+"</north>\n")
    kml.write("                     <south>"+str(margins[2])+"</south>\n")
    kml.write("                     <east>"+str(margins[1])+"</east>\n")
    kml.write("                     <west>"+str(margins[0])+"</west>\n")
    kml.write("             </LatLonBox>\n")
    kml.write("     </GroundOverlay>\n")




    kml.write(" </Document>\n")
    kml.write(" </kml>\n")        
    kml.close()
    return

def _PlotMap_KML(tensor_info, segments, point_sources, solution, default_dirs,
             convex_hulls=[], stations_str=None, stations_gps=None, stations_cgps=None,
             option='Solucion.txt', max_slip=None, legend_len=None, scale=None, limits=[None,None,None,None], evID=None):
    """We plot slip map.
    """
    import collections
    print('Creating Slip Map...')
    plane_info = segments[0]
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
            = pl_mng.__unpack_plane_data(plane_info)
    x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
    y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
    slip = solution['slip']
    #
    # accurate plot coordinates
    #
    segments_lats, segments_lons, segments_deps = __redefine_lat_lon(segments, point_sources)
    min_lats = [min(segment_lat.flatten()) for segment_lat in segments_lats]
    max_lats = [max(segment_lat.flatten()) for segment_lat in segments_lats]
    min_lons = [min(segment_lon.flatten()) for segment_lon in segments_lons]
    max_lons = [max(segment_lon.flatten()) for segment_lon in segments_lons]
    min_lat = np.min(min_lats)# - 0.5
    max_lat = np.max(max_lats)# + 0.5
    min_lon = np.min(min_lons)# - 0.5
    max_lon = np.max(max_lons)# + 0.5

    margin = 1.3 * (stk_subfaults * delta_strike) / 111.19#min(3 * (stk_subfaults * delta_strike) / 111.19, 10)
    lat0 = tensor_info['lat']
    lon0 = tensor_info['lon']
    tectonic = '{}.shp'.format(default_dirs['trench_graphics'])
    dictn = {
            'projection':ccrs.AlbersEqualArea(lon0,lat0,standard_parallels=(lat0-3., lat0+3)),
            'transform': ccrs.PlateCarree(),
            #'facecolor': '#eafff5'
            'facecolor': 'None'
    }

    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111,projection=dictn['projection'],facecolor=dictn['facecolor'], frameon=False)
    ax.spines['geo'].set_linewidth(2)
    gl = ax.gridlines()
    gl.top_labels=False
    gl.bottom_labels=False
    gl.left_labels=False
    gl.right_labels=False
    fig.subplots_adjust(hspace=0, wspace=0, top=0.9, bottom=0.1, right=0.8)
    tectonic=None
    shpfilename=None
    countries=None
    
    if stations_str is not None:
        for file in stations_str:
            name = file['name']
            latp, lonp = file['location']
            min_lat = min(min_lat, latp)
            max_lat = max(max_lat, latp)
            min_lon = min(min_lon, lonp)
            max_lon = max(max_lon, lonp)
            distance = max(np.abs(latp - lat0), np.abs(lonp - lon0))
            margin = max(margin, 1.2 * distance)
            ax.plot(lonp, latp, 'w', marker='v', markersize=15, lw=0.3, markeredgecolor='k',
                    transform=dictn['transform'], zorder=4)
            #ax.text(lonp + 0.02, latp + 0.02, '{}'.format(name),
            #        transform=dictn['projection'], zorder=4)
    if stations_cgps is not None:
        for file in stations_cgps:
            name = file['name']
            latp, lonp = file['location']
            min_lat = min(min_lat, latp)
            max_lat = max(max_lat, latp)
            min_lon = min(min_lon, lonp)
            max_lon = max(max_lon, lonp)
            distance = max(np.abs(latp - lat0), np.abs(lonp - lon0))
            margin = max(margin, 1.2 * distance)
            ax.plot(lonp, latp, 'b', marker='^', markersize=10, lw=0.3, markeredgecolor='k',
                    transform=dictn['transform'], zorder=4)
            #ax.text(lonp + 0.02, latp + 0.02, '{}'.format(name),
            #        transform=dictn['projection'], zorder=4)

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
            distance = max(np.abs(sta_lat - lat0), np.abs(sta_lon - lon0))
            margin = max(margin, 1.2 * distance)
        max_obs = np.max(max_obs)
        for name, sta_lat, sta_lon, obs, syn, error in stations_gps2:
            if scale == None:
                scale=2  # bigger number here makes arrows look longer
            scale2=2
            gps_z, gps_n, gps_e = syn
            east_west = float(gps_e) / max_obs / (1./scale)
            north_south = float(gps_n) / max_obs / (1./scale)
            plt.arrow(sta_lon, sta_lat, east_west, north_south, facecolor='r',
                      zorder=7, linewidth=0.5, width=0.02*scale2, head_width=0.05*scale2, head_length=0.05*scale2,
                      transform=dictn['transform'], edgecolor='k')
            up_down = float(gps_z) / max_obs#/ 100
            gps_z, gps_n, gps_e = obs
            east_west = float(gps_e) / max_obs / (1./scale)
            north_south = float(gps_n) / max_obs / (1./scale)
            plt.arrow(sta_lon, sta_lat, east_west, north_south, facecolor='k', zorder=5,
                      linewidth=0.5, width=0.02*scale2, head_width=0.05*scale2, head_length=0.05*scale2,
                      transform=dictn['transform'])
            up_down = float(gps_z) / max_obs / (1./scale) #/ 100
            ### INCLUDE GNSS STATION NAMES ON PLOT? ###
    #        plt.text(sta_lon + 0.02, sta_lat + 0.02, '{}'.format(name),
    #                 transform=ccrs.PlateCarree())
            err_z, err_n, err_e = error
            width = float(err_e) / max_obs/ (1./scale)
            height = float(err_n) / max_obs/ (1./scale)
        if legend_len == None:
            legend_len=20
        plt.arrow(0.83*(max_lon-min_lon)+(min_lon+0.5), 0.08*(max_lat-min_lat)+(min_lat-0.5), legend_len / max_obs / (1./scale), 0, facecolor='k', edgecolor='k', zorder=10, linewidth=0.5, width=0.02*scale2, head_width=0.05*scale2, head_length=0.05*scale2,transform=dictn['transform'])
        plt.arrow(0.83*(max_lon-min_lon)+(min_lon+0.5), 0.04*(max_lat-min_lat)+(min_lat-0.5), legend_len / max_obs / (1./scale), 0, facecolor='r', edgecolor='k', zorder=10, linewidth=0.5, width=0.02*scale2, head_width=0.05*scale2, head_length=0.05*scale2,transform=dictn['transform'])
        plt.text(0.9, 0.09, f'{int(legend_len)} cm', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
        plt.text(0.82, 0.06, 'Observed GNSS', horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
        plt.text(0.82, 0.03, 'Synthetic GNSS', horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    max_slip = max([np.amax(slip_fault) for slip_fault in slip])\
        if not max_slip else max_slip
    if limits == [None,None,None,None]:
        margins = [min_lon - 0.5, max_lon + 0.5, min_lat - 0.5, max_lat + 0.5]
    else:
        margins = [min_lon - limits[0], max_lon + limits[1], min_lat - limits[2], max_lat + limits[3]]

    ax.add_patch(Rectangle((margins[0]+.05, margins[2]+.05), 0.975*(margins[1]-margins[0]), 0.1*(margins[3]-margins[2]), edgecolor='k', facecolor='0.5', alpha=0.5, transform=dictn['transform'],zorder=3))
    ax.set_extent(margins)
    ax = set_KML_map_cartopy(ax, margins, tectonic=tectonic, countries=countries, bathymetry=None, faults=True, aftershocks=True, transform=dictn['transform'])
    ax.plot(lon0, lat0, '*', markersize=15, transform=dictn['transform'],
            zorder=5, markerfacecolor="None", markeredgecolor='k', markeredgewidth=1.5)

    # outline the segments in black #
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

        ax.plot(corners[0,:],corners[1,:], 'k',lw=2, zorder=4, transform=dictn['transform'])
        ax.plot(updip[0,:],updip[1,:], 'r', lw=1.5, zorder=4, transform=dictn['transform'])

        #plot multiple segment hypocenters if any
        if 'hypocenter' in segments[segment]:
            hyp = segments[segment]['hypocenter']
            ax.plot(hyp['lon'], hyp['lat'], '*', markersize=15, transform=dictn['transform'], zorder=5, markerfacecolor="None", markeredgecolor='k', markeredgewidth=1.5)
    #
    # plot slip map
    #
    ax, cs = plot_map(ax, segments_lats, segments_lons, slip, max_val=max_slip,
                      transform = dictn['transform'])
    sm = plt.cm.ScalarMappable(cmap=slipcpt,norm=plt.Normalize(vmin=0.,vmax=max_slip/100.))

    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#    cb_ax = inset_axes(ax, width = "30%", height = "5%", loc = 'lower left')
    cb_ax = ax.inset_axes([0.06, 0.03, 0.3, 0.02])
    cbar = plt.colorbar(sm, cax=cb_ax, orientation='horizontal')
    cbar.outline.set_linewidth(3)
    cbar.set_label('Slip (m)')
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    
    plt.savefig('Map_kml.png', dpi=300, bbox_inches='tight')
    plt.savefig('Map_kml.ps')
    plt.close()

    _write_KML(segments, point_sources, evID, margins, directory=None)
    return

def set_KML_map_cartopy(ax, margins, tectonic=None, countries=None, bathymetry=None, faults=True, aftershocks=True, transform=None):
    """
    """
    ax.set_extent(margins)
#    ax.coastlines(resolution='10m', zorder=3)
    ax.spines["bottom"].set_linewidth(10)
    ax.spines["top"].set_linewidth(10)
    ax.spines["left"].set_linewidth(10)
    ax.spines["right"].set_linewidth(10)
    gl = ax.gridlines(linewidth=1, color='black', alpha=0.3, draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels=False
    gl.bottom_labels=False
    #ax.add_feature(cartopy.feature.OCEAN)
    ocean = cartopy.feature.NaturalEarthFeature('physical', 'ocean', \
        scale='110m', edgecolor='none', facecolor=cartopy.feature.COLORS['water'])
#    ax.add_feature(ocean)
    if tectonic:
        ax.add_feature(tectonic)
    if countries:
        ax.add_feature(countries)
    if bathymetry:
        ax.add_feature(bathymetry)
    if faults==True:
        faults = glob('*.fault')
        if len(faults) > 0:
            for kfault in range(len(faults)):
                print('...Adding fault trace from: '+str(faults[kfault]))
                fault_trace = np.genfromtxt(faults[kfault])
                ax.plot(fault_trace[:,0],fault_trace[:,1],'k',zorder=100,transform=transform)
    if aftershocks==True:
        aftershocks = glob('*aftershock*')
        if len(aftershocks) > 0:
            for kafter in range(len(aftershocks)):
                print('...Adding aftershocks from: '+str(aftershocks[kafter]))
                aftershock = np.genfromtxt(aftershocks[kafter], delimiter="\t")
                aftershock_lat = aftershock[:,3]
                aftershock_lon = aftershock[:,4]
                aftershock_mag = aftershock[:,6]
                ax.scatter(aftershock_lon, aftershock_lat, s=aftershock_mag*10, c='0.5', zorder = 200,transform=transform)

    min_lon, max_lon, min_lat, max_lat = margins
    return ax


if __name__ == '__main__':
    # directory = '/Users/degoldberg/Desktop/'
    # eventID='testevent'

    parser = argparse.ArgumentParser()
    parser.add_argument("-ev","--EventID", nargs='?', const='Not Provided', type=str,
                        help="Provide event ID")
    parser.add_argument("-st", "--strong", action="store_true",
                        help="plot strong motion stations and strong motion misfit")
    parser.add_argument("--cgps", action="store_true",
                        help="plot misfit of cGPS data")
    parser.add_argument("--gps", action="store_true", help="plot GPS data")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
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

    default_dirs = mng.default_dirs()

    stations_str, stations_cgps, stations_gps = [None, None, None]
    if args.gps:
        names, lats, lons, observed, synthetic, error\
                = get_outputs.retrieve_gps()
        stations_gps = zip(names, lats, lons, observed, synthetic, error)
    if args.cgps:
        traces_info_cgps = json.load(open('cgps_waves.json'))
    if args.strong:
        traces_info = json.load(open('strong_motion_waves.json'))

    if args.EventID:
        evID = args.EventID
    else:
        evID = None

    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
#    segments, rise_time, point_sources = pl_mng.__read_planes_info() # Loads point sources and segments information
    segments_data = json.load(open('segments_data.json'))
    segments = segments_data['segments']
    rise_time = segments_data['rise_time']
    connections = None
    if 'connections' in segments_data:
        connections = segments_data['connections']
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections)
    solution = get_outputs.read_solution_static_format(segments)

    if args.maxvalue != None:
        max_val=args.maxvalue
    else:
        max_val=None
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

    _PlotMap_KML(tensor_info, segments, point_sources, solution, default_dirs, stations_str=stations_str, stations_gps=stations_gps, stations_cgps=stations_cgps, max_slip=max_val, legend_len=legend_len, scale=scale, limits=limits, evID=evID)


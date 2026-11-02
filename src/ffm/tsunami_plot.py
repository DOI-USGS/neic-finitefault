"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""
import json
import os

import matplotlib.pyplot as plt
import numpy as np
from clawpack.geoclaw import fgmax_tools  # type: ignore
from clawpack.geoclaw import topotools  # type: ignore
from clawpack.visclaw import (  # type: ignore [import-untyped, import-not-found, attr-defined]
    gaugetools,
    geoplot,
)
from clawpack.visclaw.data import ClawPlotData  # type: ignore
from pylab import axis, floor, legend, mod, plot, title, xlabel, xticks  # type: ignore


# --------------------------
def setplot(plotdata=None, parameter_file="geoclaw_parameter_file.json"):
    # --------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    """

    # ------------------------------------------------------------------
    # Parameters from json file:
    # ------------------------------------------------------------------
    f = open(parameter_file)
    parameter = json.load(f)

    # ------------------------------------------------------------------

    if plotdata is None:
        plotdata = ClawPlotData()

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = "ascii"  # 'ascii', 'binary', 'netcdf'

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):

        gaugetools.plot_gauge_locations(
            current_data.plotdata, gaugenos="all", format_string="ko", add_labels=True
        )

    # To add title with time format hours:minutes:seconds...

    def title_hours(current_data):

        t = current_data.t
        hours = int(t / 3600.0)
        tmin = mod(t, 3600.0)
        min = int(tmin / 60.0)
        tsec = mod(tmin, 60.0)
        sec = int(mod(tmin, 60.0))
        timestr = "%s:%s:%s" % (hours, str(min).zfill(2), str(sec).zfill(2))
        title("%s after earthquake" % timestr)

    # -----------------------------------------
    # Figure for surface
    # -----------------------------------------
    plotfigure = plotdata.new_plotfigure(name="Surface", figno=0)
    plotfigure.kwargs = {"figsize": (10, 10)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes("pcolor")
    plotaxes.axescmd = "subplot(111)"
    plotaxes.title = "Surface"
    plotaxes.scaled = True

    def fixup(current_data):
        addgauges(current_data)
        title_hours(current_data)

    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.2
    plotitem.pcolor_cmax = 0.2
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0, 0, 0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0, 0, 0]
    plotitem.patchedges_show = 0
    plotaxes.xlimits = [
        parameter["model_region"]["west"],
        parameter["model_region"]["east"],
    ]
    plotaxes.ylimits = [
        parameter["model_region"]["south"],
        parameter["model_region"]["north"],
    ]

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type="2d_contour")
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = np.linspace(-3000, -3000, 1)
    plotitem.amr_contour_colors = ["y"]  # color on each level
    plotitem.kwargs = {"linestyles": "solid", "linewidths": 2}
    plotitem.amr_contour_show = [1, 0, 0]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    # -----------------------------------------
    # Figures for gauges
    # -----------------------------------------
    plotfigure = plotdata.new_plotfigure(
        name="Surface at gauges", figno=300, type="each_gauge"
    )
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = "auto"
    plotaxes.ylimits = "auto"
    plotaxes.title = "Surface"

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type="1d_plot")
    plotitem.plot_var = 3
    plotitem.plotstyle = "b-"

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type="1d_plot")
    plotitem.show = False

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0, :]
        eta = q[3, :]
        topo = eta - h
        return topo

    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = "g-"

    def add_zeroline(current_data):

        t = current_data.t
        gaugeno = current_data.gaugeno

        if gaugeno == 32412:
            try:
                plot(TG32412[:, 0], TG32412[:, 1], "r")
                legend(["GeoClaw", "Obs"], loc="lower right")
            except:
                pass
            axis((0, t.max(), -0.3, 0.3))

        plot(t, 0 * t, "k")
        n = int(floor(t.max() / 3600.0) + 2)
        xticks([3600 * i for i in range(n)], ["%i" % i for i in range(n)])
        xlabel("time (hours)")

    plotaxes.afteraxes = add_zeroline

    # -----------------------------------------
    # Figures for fgmax - max values on fixed grids
    # -----------------------------------------
    otherfigure = plotdata.new_otherfigure(
        name="max amplitude and arrival times", fname="amplitude_times.png"
    )

    def makefig(plotdata):

        fg = fgmax_tools.FGmaxGrid()
        fg.outdir = plotdata.outdir
        data_file = os.path.join(plotdata.outdir, "fgmax_grids.data")
        fg.read_fgmax_grids_data(fgno=1, data_file=data_file)
        fg.read_output()

        clines_zeta = [1, 3, 5, 7, 10, 15, 20, 25, 30, 50]
        colors = geoplot.discrete_cmap_1(clines_zeta)
        plt.figure(1)
        plt.clf()
        zeta = np.where(fg.B > 0, fg.h, fg.h + fg.B)  # surface elevation in ocean
        plt.contourf(fg.X, fg.Y, zeta * 100, clines_zeta, colors=colors)
        cbar = plt.colorbar()
        cbar.set_label("Maximum Amplitude (cm)")
        plt.contour(fg.X, fg.Y, fg.B, [0.0], colors="k")  # coastline

        # plot arrival time contours and label:
        arrival_t = fg.arrival_time / 3600.0  # arrival time in hours
        clines_t = np.linspace(0, 8, 17)  # hours
        clines_t_label = clines_t[2::2]  # which ones to label
        clines_t_colors = ([0.5, 0.5, 0.5],)
        con_t = plt.contour(fg.X, fg.Y, arrival_t, clines_t, colors=clines_t_colors)
        plt.clabel(con_t, clines_t_label)

        # fix axes:
        plt.ticklabel_format(style="plain", useOffset=False)
        plt.xticks(rotation=20)
        plt.gca().set_aspect(1.0 / np.cos(fg.Y.mean() * np.pi / 180.0))
        plt.title("Maximum amplitude / arrival times")

        if not os.path.isdir(plotdata.plotdir):
            os.mkdir(plotdata.plotdir)
        fname = os.path.join(plotdata.plotdir, "amplitude_times.pdf")
        plt.savefig(fname)
        print("Created ", fname)

    otherfigure.makefig = makefig

    # -----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True  # print figures
    plotdata.print_format = "png"  # file format
    plotdata.print_framenos = "all"  # list of frames to print
    plotdata.print_gaugenos = "all"  # list of gauges to print
    plotdata.print_fignos = "all"  # list of figures to print
    plotdata.html = True  # create html files of plots?
    plotdata.html_homelink = "../README.html"  # pointer for top of index
    plotdata.latex = True  # create latex file of plots?
    plotdata.latex_figsperline = 2  # layout of plots
    plotdata.latex_framesperline = 1  # layout of plots
    plotdata.latex_makepdf = False  # also run pdflatex?
    plotdata.parallel = True  # make multiple frame png's at once

    return plotdata

"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

from __future__ import absolute_import, print_function

from dotenv import load_dotenv
load_dotenv()

import json
import os

import numpy as np

try:
    CLAW = os.environ["CLAW"]
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(CLAW, "geoclaw", "scratch")


# ---------------------------------------------------------------
def setrun(claw_pkg="geoclaw", parameter_file="geoclaw_parameter_file.json"):
    # ---------------------------------------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data  # type: ignore

    assert claw_pkg.lower() == "geoclaw", "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    # ------------------------------------------------------------------
    # Parameters from json file:
    # ------------------------------------------------------------------
    f = open(parameter_file)
    parameter = json.load(f)

    # ------------------------------------------------------------------
    # GeoClaw specific parameters:
    # ------------------------------------------------------------------
    rundata = setgeo(rundata, parameter_file)

    # ------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    # ------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # Set single grid parameters first.
    # See below for AMR parameters.

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = parameter["model_region"]["west"]  # west longitude
    clawdata.upper[0] = parameter["model_region"]["east"]  # east longitude

    clawdata.lower[1] = parameter["model_region"]["south"]  # south latitude
    clawdata.upper[1] = parameter["model_region"]["north"]  # north latitude

    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] = 30
    clawdata.num_cells[1] = 30

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0

    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False  # True to restart from prior results
    clawdata.restart_file = "foo.bar"  # File to use for restart data

    # -------------
    # Output times:
    # --------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style == 1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = parameter["num_output_times"]
        clawdata.tfinal = parameter["tfinal"]
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0.5, 1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 3
        clawdata.output_t0 = True

    clawdata.output_format = "binary"  # 'ascii' or 'binary'

    clawdata.output_q_components = "all"  # need all
    clawdata.output_aux_components = "none"  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False  # output aux arrays each frame

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.2

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2

    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = "unsplit"

    # For unsplit method, transverse_waves can be
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ["mc", "mc", "mc"]

    clawdata.use_fwaves = True  # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = "godunov"

    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = "extrap"
    clawdata.bc_upper[0] = "extrap"

    clawdata.bc_lower[1] = "extrap"
    clawdata.bc_upper[1] = "extrap"

    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif np.abs(clawdata.checkpt_style) == 1:
        # Checkpoint only at tfinal.
        pass

    elif np.abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.
        clawdata.checkpt_times = [0.1, 0.15]

    elif np.abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5

    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # maximum size of patches in each direction (matters in parallel):
    amrdata.max1d = 60

    # max number of refinement levels:
    amrdata.amr_levels_max = 3

    # List of refinement ratios at each level (length at least mxnest-1)
    amrdata.refinement_ratios_x = [2, 6]
    amrdata.refinement_ratios_y = [2, 6]
    amrdata.refinement_ratios_t = [2, 6]

    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ["center", "capacity", "yleft"]

    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False  # use Richardson?
    amrdata.flag_richardson_tol = 0.002  # Richardson tolerance
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0

    #  ----- For developers -----
    # Toggle debugging print statements:
    amrdata.dprint = False  # print domain flags
    amrdata.eprint = False  # print err est flags
    amrdata.edebug = False  # even more err est flags
    amrdata.gprint = False  # grid bisection/clustering
    amrdata.nprint = False  # proper nesting output
    amrdata.pprint = False  # proj. of tagged points
    amrdata.rprint = False  # print regridding summary
    amrdata.sprint = False  # space/memory output
    amrdata.tprint = True  # time step reporting each level
    amrdata.uprint = False  # update/upbnd reporting

    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    # ---------------
    # Regions:
    # ---------------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    for kregion in range(len(parameter["high_res_regions"])):
        kstr = str(kregion)
        minlevel = parameter["high_res_regions"][kstr]["minlevel"]
        maxlevel = parameter["high_res_regions"][kstr]["maxlevel"]
        t1 = parameter["high_res_regions"][kstr]["t1"]
        t2 = parameter["high_res_regions"][kstr]["t2"]
        x1 = parameter["high_res_regions"][kstr]["x1"]
        x2 = parameter["high_res_regions"][kstr]["x2"]
        y1 = parameter["high_res_regions"][kstr]["y1"]
        y2 = parameter["high_res_regions"][kstr]["y2"]
        rundata.regiondata.regions.append([minlevel, maxlevel, t1, t2, x1, x2, y1, y2])

    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    for kgauge in range(len(parameter["gauges"])):
        kstr = str(kgauge)
        gauge_id = int(parameter["gauges"][kstr]["gauge_id"])
        lon = parameter["gauges"][kstr]["lon"]
        lat = parameter["gauges"][kstr]["lat"]
        t1 = parameter["gauges"][kstr]["t1"]
        t2 = parameter["gauges"][kstr]["t2"]
        rundata.gaugedata.gauges.append([gauge_id, lon, lat, t1, t2])

    return rundata
    # end of function setrun
    # ----------------------


# -------------------
def setgeo(rundata, parameter_file):
    # -------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    # ------------------------------------------------------------------
    # Parameters from json file:
    # ------------------------------------------------------------------
    f = open(parameter_file)
    parameter = json.load(f)

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.0e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 1e6

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = parameter["wave_tolerance"]

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, fname]
    topo_path = os.path.join(parameter["topo_path"])
    topo_data.topofiles.append([2, topo_path])

    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, fname]
    dtopo_path = os.path.join(parameter["dtopo_path"])
    dtopo_type = int(parameter["dtopo_type"])
    dtopo_data.dtopofiles.append([dtopo_type, dtopo_path])
    dtopo_data.dt_max_dtopo = parameter["dt_max_dtopo"]

    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [fname]

    # == fgout grids ==
    # new style as of v5.9.0 (old rundata.fixed_grid_data is deprecated)
    # set rundata.fgout_data.fgout_grids to be a
    # list of objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:
    # rundata.fgout_data.fgout_grids = []

    # == fgmax_grids.data values ==
    # NEW STYLE STARTING IN v5.7.0
    from clawpack.geoclaw import fgmax_tools  # type: ignore
    from clawpack.geoclaw import fgout_tools  # type: ignore

    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 2  # Save depth and speed

    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list to start

    # Now append to this list objects of class fgmax_tools.FGmaxGrid
    # specifying any fgmax grids.

    # Points on a uniform 2d grid:
    dx_fine = 2.0 / (3.0 * 4.0)  # grid resolution at finest level
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2  # uniform rectangular x-y grid
    fg.x1 = (
        parameter["model_region"]["west"] + dx_fine / 2.0
    )  # specify pts to align with FV cell centers
    fg.x2 = parameter["model_region"]["east"] - dx_fine / 2.0
    fg.y1 = parameter["model_region"]["south"] + dx_fine / 2.0
    fg.y2 = parameter["model_region"]["north"] - dx_fine / 2.0
    fg.dx = dx_fine
    fg.tstart_max = 10.0  # when to start monitoring max values
    fg.tend_max = parameter["tfinal"]  # when to stop monitoring max values
    fg.dt_check = 60.0  # target time (sec) increment between updating
    # max values
    fg.min_level_check = 3  # which levels to monitor max on
    fg.arrival_tol = 1.0e-5  # tolerance for flagging arrival

    fg.interp_method = 0  # 0 ==> pw const in cells, recommended
    fgmax_grids.append(fg)  # written to fgmax_grids.data

    # == fgout_grids.data values ==
    # NEW IN v5.9.0
    # Set rundata.fgout_data.fgout_grids to be a list of
    # objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:
    fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially

    fgout = fgout_tools.FGoutGrid()
    fgout.fgno = 1
    fgout.point_style = 2  # will specify a 2d grid of points
    fgout.output_format = "binary32"  # 4-byte, float32
    fgout.nx = 200
    fgout.ny = 250
    fgout.x1 = parameter["model_region"][
        "west"
    ]  # specify edges (fgout pts will be cell centers)
    fgout.x2 = parameter["model_region"]["east"]
    fgout.y1 = parameter["model_region"]["south"]
    fgout.y2 = parameter["model_region"]["north"]
    fgout.tstart = 0.0
    fgout.tend = parameter["tfinal"]
    fgout.nout = 25
    fgout_grids.append(fgout)  # written to fgout_grids.data

    return rundata
    # end of function setgeo
    # ----------------------


if __name__ == "__main__":
    # Set up run-time parameters and write all data files.
    import sys

    from clawpack.geoclaw import kmltools  # type: ignore

    rundata = setrun(*sys.argv[1:])
    rundata.write()

    kmltools.make_input_data_kmls(rundata)

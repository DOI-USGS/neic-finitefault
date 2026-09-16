import json
import os
import pathlib
import shutil
import tempfile

from ffm.geoclaw_setrun import setgeo, setrun

from .testutils import RESULTS_DIR


def test_setrun():
    tempdir = tempfile.mkdtemp()
    try:
        geoclaw_parameter_file = RESULTS_DIR / "NP1" / "geoclaw_parameter_file.json"
        f = open(geoclaw_parameter_file)
        parameter = json.load(f)
        rundata = setrun(parameter_file=geoclaw_parameter_file)
        assert rundata.clawdata.lower == [
            parameter["model_region"]["west"],
            parameter["model_region"]["south"],
        ]
        assert rundata.clawdata.upper == [
            parameter["model_region"]["east"],
            parameter["model_region"]["north"],
        ]
        assert rundata.clawdata.num_cells == [30, 30]
        assert rundata.clawdata.num_eqn == 3
        assert rundata.clawdata.num_aux == 3
        assert rundata.clawdata.capa_index == 2
        assert rundata.clawdata.t0 == 0.0
        assert rundata.clawdata.output_style == 1
        assert rundata.clawdata.num_output_times == parameter["num_output_times"]
        assert rundata.clawdata.tfinal == parameter["tfinal"]
        assert rundata.clawdata.output_t0 == True
        assert rundata.clawdata.output_format == "binary"
        assert rundata.clawdata.output_q_components == "all"
        assert rundata.clawdata.output_aux_components == "none"
        assert rundata.clawdata.output_aux_onlyonce == False
        assert rundata.clawdata.dt_variable == True
        assert rundata.clawdata.dt_initial == 0.2
        assert rundata.clawdata.dt_max == 1e99
        assert rundata.clawdata.cfl_desired == 0.75
        assert rundata.clawdata.cfl_max == 1.0
        assert rundata.clawdata.steps_max == 5000
        assert rundata.clawdata.order == 2
        assert rundata.clawdata.dimensional_split == "unsplit"
        assert rundata.clawdata.transverse_waves == 2
        assert rundata.clawdata.num_waves == 3
        assert rundata.clawdata.limiter == ["mc", "mc", "mc"]
        assert rundata.clawdata.use_fwaves == True
        assert rundata.clawdata.source_split == "godunov"
        assert rundata.clawdata.num_ghost == 2
        assert rundata.clawdata.checkpt_style == 0
        assert rundata.amrdata.max1d == 60
        assert rundata.amrdata.amr_levels_max == 3
        assert rundata.amrdata.refinement_ratios_x == [2, 6]
        assert rundata.amrdata.refinement_ratios_y == [2, 6]
        assert rundata.amrdata.refinement_ratios_t == [2, 6]
        assert rundata.amrdata.aux_type == ["center", "capacity", "yleft"]
        assert rundata.amrdata.flag_richardson == False
        assert rundata.amrdata.flag2refine == True
        assert rundata.amrdata.regrid_interval == 3
        assert rundata.amrdata.regrid_buffer_width == 2
        assert rundata.amrdata.clustering_cutoff == 0.700000
        for kregion in range(len(parameter["high_res_regions"])):
            kstr = str(kregion)
            assert rundata.regiondata.regions[kregion] == [
                parameter["high_res_regions"][kstr]["minlevel"],
                parameter["high_res_regions"][kstr]["maxlevel"],
                parameter["high_res_regions"][kstr]["t1"],
                parameter["high_res_regions"][kstr]["t2"],
                parameter["high_res_regions"][kstr]["x1"],
                parameter["high_res_regions"][kstr]["x2"],
                parameter["high_res_regions"][kstr]["y1"],
                parameter["high_res_regions"][kstr]["y2"],
            ]
        for kgauge in range(len(parameter["gauges"])):
            kstr = str(kgauge)
            assert rundata.gaugedata.gauges[kgauge] == [
                int(parameter["gauges"][kstr]["gauge_id"]),
                parameter["gauges"][kstr]["lon"],
                parameter["gauges"][kstr]["lat"],
                parameter["gauges"][kstr]["t1"],
                parameter["gauges"][kstr]["t2"],
            ]
    finally:
        shutil.rmtree(tempdir)


def test_setgeo():
    tempdir = tempfile.mkdtemp()
    try:
        geoclaw_parameter_file = RESULTS_DIR / "NP1" / "geoclaw_parameter_file.json"
        f = open(geoclaw_parameter_file)
        parameter = json.load(f)
        rundata = setrun(parameter_file=geoclaw_parameter_file)
        rundata = setgeo(rundata, geoclaw_parameter_file)
        assert rundata.geo_data.gravity == 9.81
        assert rundata.geo_data.coordinate_system == 2
        assert rundata.geo_data.earth_radius == 6367.5e3
        assert rundata.geo_data.coriolis_forcing == False
        assert rundata.geo_data.sea_level == 0.0
        assert rundata.geo_data.dry_tolerance == 1.0e-3
        assert rundata.geo_data.friction_forcing == True
        assert rundata.geo_data.manning_coefficient == 0.025
        assert rundata.geo_data.friction_depth == 1e6
        assert rundata.refinement_data.variable_dt_refinement_ratios == True
        assert rundata.refinement_data.wave_tolerance == parameter["wave_tolerance"]
        assert rundata.topo_data.topofiles[0] == [2, "topo.topo"]
        assert rundata.dtopo_data.dtopofiles[0] == [
            int(parameter["dtopo_type"]),
            os.path.join(parameter["dtopo_path"]),
        ]
        assert rundata.dtopo_data.dt_max_dtopo == parameter["dt_max_dtopo"]
        assert rundata.qinit_data.qinit_type == 0
        assert rundata.qinit_data.qinitfiles == []
        assert rundata.fgmax_data.num_fgmax_val == 2
    finally:
        shutil.rmtree(tempdir)

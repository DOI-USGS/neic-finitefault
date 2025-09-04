import json
import os
import pathlib
import shutil
import tempfile

from wasp.fault_plane import point_sources_param
from wasp.get_outputs import read_solution_static_format
from wasp.plot_Map import PlotMap

from .testutils import (
    RESULTS_DIR,
    get_cgnss_json,
    get_segments_data,
    get_tele_waves_json,
    get_tensor_info,
    update_manager_file_locations,
)


def testPlotMap():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solution.txt", tempdir / "Solution.txt")
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "segments_data.json", tempdir / "segments_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "segments_data.json", tempdir / "segments_data.json"
        )
        cgnss_waves = get_cgnss_json()
        new_cgnss_waves = update_manager_file_locations(
            cgnss_waves, tempdir, replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(tempdir / "data")
        os.mkdir(tempdir / "cGNSS")
        for o, n in zip(cgnss_waves, new_cgnss_waves):
            shutil.copyfile(o["file"], n["file"])
        tele_waves = get_tele_waves_json()
        new_tele_waves = update_manager_file_locations(
            tele_waves, tempdir / "data", replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(tempdir / "data" / "P")
        for o, n in zip(tele_waves, new_tele_waves):
            shutil.copyfile(o["file"], n["file"])
        segments = get_segments_data()
        tensor = get_tensor_info()
        point_sources = point_sources_param(
            segments["segments"],
            tensor,
            segments["rise_time"],
        )
        solution = read_solution_static_format(
            segments["segments"],
            tempdir,
        )
        PlotMap(
            tensor,
            segments["segments"],
            point_sources,
            solution,
            {"root_dir": pathlib.Path(__file__).parent.parent.parent},
            files_str=new_tele_waves,
            stations_cgnss=new_cgnss_waves,
            directory=tempdir,
        )
        assert (tempdir / "PyGMT_Map.png").exists()
    finally:
        shutil.rmtree(tempdir)

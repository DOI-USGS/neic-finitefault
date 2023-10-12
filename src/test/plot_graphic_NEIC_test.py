import json
import os
import pathlib
import shutil
import tempfile
from test.testutils import (
    RESULTS_DIR,
    get_cgps_json,
    get_insar_json,
    get_segments_data,
    get_strong_motion_json,
    get_surf_waves_json,
    get_tele_waves_json,
    get_tensor_info,
    get_velmodel_data,
    update_manager_file_locations,
)

import pytest

from wasp.fault_plane import point_sources_param, shear_modulous
from wasp.get_outputs import get_insar, read_solution_static_format, retrieve_gps
from wasp.plot_graphic_NEIC import (
    PlotComparisonMap,
    PlotInsar,
    PlotSlipDist_Compare,
    _plot_vel_model,
    _PlotCumulativeSlip,
    _PlotMultiSlipDist,
    _PlotRiseTime,
    _PlotRuptTime,
    calculate_cumulative_moment_tensor,
    plot_beachball,
    plot_ffm_sol,
    plot_misfit,
    shakemap_polygon,
)

SEGMENTS = get_segments_data()
TENSOR = get_tensor_info()


def __get_solution():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solucion.txt", tempdir / "Solucion.txt")
        point_sources = point_sources_param(
            SEGMENTS["segments"],
            TENSOR,
            SEGMENTS["rise_time"],
        )
        solution = read_solution_static_format(
            SEGMENTS["segments"],
            tempdir,
        )
        return point_sources, solution
    finally:
        shutil.rmtree(tempdir)


POINT_SOURCES, SOLUTION = __get_solution()


def test_calculate_cumulative_moment_tensor():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        with open(tempdir / "segments_data.json", "w") as f:
            json.dump(SEGMENTS, f)
        calculate_cumulative_moment_tensor(SOLUTION, directory=tempdir)
        assert (tempdir / "Cumulative_Moment_Tensor.png").exists()
    finally:
        shutil.rmtree(tempdir)


def test_plot():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solucion.txt", tempdir / "Solucion.txt")
        cgps_waves = get_cgps_json()
        new_cgps_waves = update_manager_file_locations(
            cgps_waves,
            tempdir,
            replace_dir=str(RESULTS_DIR / "data"),
        )
        os.mkdir(tempdir / "data")
        os.mkdir(tempdir / "cGPS")
        for o, n in zip(cgps_waves, new_cgps_waves):
            shutil.copyfile(o["file"], n["file"])
        tele_waves = get_tele_waves_json(all=True)
        new_tele_waves = update_manager_file_locations(
            tele_waves, tempdir / "data", replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(tempdir / "data" / "P")
        os.mkdir(tempdir / "data" / "SH")
        os.mkdir(tempdir / "data" / "LONG")
        os.mkdir(tempdir / "data" / "STR")
        for o, n in zip(tele_waves, new_tele_waves):
            shutil.copyfile(o["file"], n["file"])
        surf_waves = get_surf_waves_json(all=True)
        new_surf_waves = update_manager_file_locations(
            surf_waves, tempdir / "data", replace_dir=str(RESULTS_DIR / "data")
        )
        for o, n in zip(surf_waves, new_surf_waves):
            shutil.copyfile(o["file"], n["file"])
        strong_waves = get_strong_motion_json(all=True)
        new_strong_waves = update_manager_file_locations(
            strong_waves, tempdir / "data", replace_dir=str(RESULTS_DIR / "data")
        )
        for o, n in zip(strong_waves, new_strong_waves):
            shutil.copyfile(o["file"], n["file"])
        new_tele_waves = update_manager_file_locations(
            tele_waves, tempdir / "data", replace_dir=str(RESULTS_DIR / "data")
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "segments_data.json", tempdir / "segments_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "static_synthetics.txt",
            tempdir / "static_synthetics.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "static_data.json", tempdir / "static_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "fault&rise_time.txt", tempdir / "fault&rise_time.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "synthetics_body.txt", tempdir / "synthetics_body.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "synthetics_surf.txt", tempdir / "synthetics_surf.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "synthetics_strong.txt",
            tempdir / "synthetics_strong.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "synthetics_cgps.txt", tempdir / "synthetics_cgps.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "velmodel_data.json", tempdir / "velmodel_data.json"
        )
        with open(tempdir / "tele_waves.json", "w") as f:
            json.dump(new_tele_waves, f)
        with open(tempdir / "surf_waves.json", "w") as f:
            json.dump(new_surf_waves, f)
        with open(tempdir / "strong_motion_waves.json", "w") as f:
            json.dump(new_strong_waves, f)
        with open(tempdir / "cgps_waves.json", "w") as f:
            json.dump(new_cgps_waves, f)
        names, lats, lons, observed, synthetic, error = retrieve_gps(directory=tempdir)
        stations_gps = zip(names, lats, lons, observed, synthetic, error)
        shear = shear_modulous(POINT_SOURCES, velmodel=get_velmodel_data())
        plot_ffm_sol(
            TENSOR,
            SEGMENTS,
            POINT_SOURCES,
            shear,
            SOLUTION,
            {"root_dir": pathlib.Path(__file__).parent.parent.parent},
            autosize=False,
            files_str=new_strong_waves,
            stations_gps=stations_gps,
            stations_cgps=new_cgps_waves,
            max_val=None,
            legend_len=None,
            scale=None,
            separate_planes=False,
            label_stations=False,
            directory=tempdir,
        )
        assert (tempdir / "MomentRate.png").exists()
        assert (tempdir / "Map.png").exists()
        assert (tempdir / "Map.pdf").exists()
        assert (tempdir / "SlipDist_plane0.png").exists()
        assert (tempdir / "Map.eps").exists()
        assert (tempdir / "SlipDist_plane0.ps").exists()
        with open(tempdir / "SlipDist_plane0.ps", "rb") as d:
            data = b"".join(d.readlines()[10:])
        with open(RESULTS_DIR / "NP1" / "SlipDist_plane0.ps", "rb") as t:
            target = b"".join(t.readlines()[10:])
        assert data == target
        plot_misfit(
            used_data_type=[
                "cgps",
                "strong",
                "surf",
                "body",
            ],
            directory=tempdir,
        )
        assert (tempdir / "Rayleigh_surf_waves.png").exists()
        assert (tempdir / "cGPS_waves.png").exists()
        assert (tempdir / "strong_motion_waves.png").exists()
        assert (tempdir / "SH_body_waves.png").exists()
        assert (tempdir / "Love_surf_waves.png").exists()
        assert (tempdir / "P_body_waves.png").exists()
    finally:
        shutil.rmtree(tempdir)


def test_plot_beachball():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        tele_waves = get_tele_waves_json(all=True)
        new_tele_waves = update_manager_file_locations(
            tele_waves, tempdir / "data", replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(tempdir / "data")
        os.mkdir(tempdir / "data" / "P")
        os.mkdir(tempdir / "data" / "SH")
        for o, n in zip(tele_waves, new_tele_waves):
            shutil.copyfile(o["file"], n["file"])
        plot_beachball(SEGMENTS["segments"], new_tele_waves, directory=tempdir)
        assert (tempdir / "Tensor.png").exists()
    finally:
        shutil.rmtree(tempdir)


def test__plot_vel_model():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        _plot_vel_model(get_velmodel_data(), tempdir)
        assert (tempdir / "crust_body_wave_vel_model.png").exists()
    finally:
        shutil.rmtree(tempdir)


@pytest.mark.skipif(
    os.getenv("CI_REGISTRY") is not None or not os.getenv("RUN_ALL", False),
    reason="Build runner does not have the resources to run",
)
def test_PlotComparisonMap():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        PlotComparisonMap(
            TENSOR,
            SEGMENTS["segments"],
            POINT_SOURCES,
            SOLUTION,
            SOLUTION,
            directory=tempdir,
        )
        assert (tempdir / "Checkerboard_Map_Comparison.png").exists()
    finally:
        shutil.rmtree(tempdir)


def test___PlotCumulativeSlip():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        c1, c2, c3, c4 = _PlotCumulativeSlip(
            SEGMENTS["segments"], POINT_SOURCES, SOLUTION, TENSOR, directory=tempdir
        )
        assert (tempdir / "CumulativeSlip_plane0.png").exists()
        assert (tempdir / "CumulativeSlip_plane0.ps").exists()
        assert c1 == "-72.38 -31.74 -0.16"
        assert c2 == "-72.06 -29.34 -0.16"
        assert c3 == "-71.15 -29.43 30.18"
        assert c4 == "-71.47 -31.83 30.18"
    finally:
        shutil.rmtree(tempdir)


def test_PlotInsar():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_synthetics.txt",
            tempdir / "insar_synthetics.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_ascending.txt", tempdir / "insar_ascending.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_descending.txt",
            tempdir / "insar_descending.txt",
        )
        new_insar = update_manager_file_locations(
            get_insar_json(),
            tempdir,
            replace_dir=str(RESULTS_DIR / "NP1"),
            file_key="name",
        )
        with open(tempdir / "insar_data.json", "w") as f:
            json.dump(new_insar, f)

        insar_data = get_insar(data_dir=tempdir)
        PlotInsar(
            TENSOR,
            SEGMENTS["segments"],
            POINT_SOURCES,
            SOLUTION,
            insar_data["ascending"][0]["points"],
            "0",
            "ascending",
            directory=tempdir,
        )
        PlotInsar(
            TENSOR,
            SEGMENTS["segments"],
            POINT_SOURCES,
            SOLUTION,
            insar_data["descending"][0]["points"],
            "0",
            "descending",
            directory=tempdir,
        )
        assert (tempdir / "InSAR_ascending_fit_0.png").exists()
        assert (tempdir / "InSAR_ascending_fit_0.ps").exists()
        assert (tempdir / "InSAR_descending_fit_0.png").exists()
        assert (tempdir / "InSAR_descending_fit_0.ps").exists()
    finally:
        shutil.rmtree(tempdir)


def test__PlotMultiSlipDist():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        _PlotMultiSlipDist(
            SEGMENTS["segments"], POINT_SOURCES, SOLUTION, directory=tempdir
        )
        assert (tempdir / "SlipDist_plane0.png").exists()
        assert (tempdir / "SlipDist_plane0.ps").exists()
    finally:
        shutil.rmtree(tempdir)


def test__PlotRiseTime():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        _PlotRiseTime(SEGMENTS["segments"], POINT_SOURCES, SOLUTION, tempdir)
        assert (tempdir / "RiseTime_plane0.png").exists()

    finally:
        shutil.rmtree(tempdir)


def test__PlotRuptTime():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        _PlotRuptTime(SEGMENTS["segments"], POINT_SOURCES, SOLUTION, tempdir)
        assert (tempdir / "RuptTime_plane0.png").exists()

    finally:
        shutil.rmtree(tempdir)


def test_PlotSlipDist_Compare():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        PlotSlipDist_Compare(
            SEGMENTS["segments"], POINT_SOURCES, SOLUTION, SOLUTION, directory=tempdir
        )
        assert (tempdir / "Checkerboard_SlipDist_plane0.png").exists()
    finally:
        shutil.rmtree(tempdir)


def test_shakemap_polygon():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shakemap_polygon(
            SEGMENTS["segments"],
            POINT_SOURCES,
            SOLUTION,
            TENSOR,
            "us20003k7a",
            directory=tempdir,
        )
        with open(tempdir / "shakemap_polygon.txt", "rb") as d:
            data = b"".join(d.readlines()[10:])
        with open(RESULTS_DIR / "NP1" / "shakemap_polygon.txt", "rb") as t:
            target = b"".join(t.readlines()[10:])
        assert data == target
    finally:
        shutil.rmtree(tempdir)

import os
import pathlib
import shutil
import tempfile
from glob import glob

import pytest

from wasp.green_functions import fk_green_fun1, gf_retrieve

from .testutils import (
    RESULTS_DIR,
    get_sampling_filter,
    get_tele_waves_json,
    get_tensor_info,
)


def test_gf_retrieve_gps():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "fault&rise_time.txt", tempdir / "fault&rise_time.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "point_sources.txt", tempdir / "point_sources.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "vel_model.txt", tempdir / "vel_model.txt"
        )
        shutil.copyfile(RESULTS_DIR / "static_data.txt", tempdir / "static_data.txt")
        os.mkdir(tempdir / "logs")
        with open(tempdir / "logs" / "green_gps_log", "w"):
            pass
        fdir = pathlib.Path(__file__).parent.parent.parent / "fortran_code"
        gf_retrieve(
            ["gps"],
            {
                "tele_gf": fdir / "bin_inversion_gfortran_f95" / "green_tele",
                "strong_motion_gf": fdir / "bin_str_f95" / "get_strong_motion",
                "gps_gf": fdir / "src_dc_f95" / "gf_static_f95",
            },
            tempdir,
        )
        files = glob(str(tempdir) + "/*")
        ## TODO DETERMINE if test data should be added for the GF files
        # in the meantime just check that the files were created
        assert str(tempdir / "Green_static_subfault.txt") in files
    finally:
        shutil.rmtree(tempdir)


def test_gf_retrieve_tele():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "filtro_tele.txt", tempdir / "filtro_tele.txt")
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "fault&rise_time.txt", tempdir / "fault&rise_time.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "point_sources.txt", tempdir / "point_sources.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "vel_model.txt", tempdir / "vel_model.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "channels_body.txt", tempdir / "channels_body.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "waveforms_body.txt", tempdir / "waveforms_body.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "instrumental_response.txt",
            tempdir / "instrumental_response.txt",
        )
        os.mkdir(tempdir / "logs")
        with open(tempdir / "logs" / "green_tele_log", "w"):
            pass
        fdir = pathlib.Path(__file__).parent.parent.parent / "fortran_code"
        gf_retrieve(
            ["body"],
            {
                "tele_gf": fdir / "bin_inversion_gfortran_f95" / "green_tele",
                "strong_motion_gf": fdir / "bin_str_f95" / "get_strong_motion",
                "gps_gf": fdir / "src_dc_f95" / "gf_static_f95",
            },
            tempdir,
        )
        files = glob(str(tempdir) + "/*")
        ## TODO DETERMINE if test data should be added for the GF files
        # in the meantime just check that the files were created
        for f in [f["name"] for f in get_tele_waves_json(all=True)]:
            assert f"{tempdir}/{f}.GRE" in files
            assert f"{tempdir}/{f}.TDE" in files
    finally:
        shutil.rmtree(tempdir)


def test_fk_green_fun0():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        cgps = fk_green_fun1(
            get_sampling_filter(), get_tensor_info(), "test_location", cgps=True
        )
        gps = fk_green_fun1(
            get_sampling_filter(), get_tensor_info(), "test_location", cgps=False
        )
        assert cgps == {
            "location": "test_location",
            "min_depth": 1,
            "max_depth": 49.8,
            "min_dist": 0,
            "max_dist": 1000,
            "dt": 0.4,
            "time_corr": 25,
        }
        assert gps == {
            "location": "test_location",
            "min_depth": 1,
            "max_depth": 49.8,
            "min_dist": 0,
            "max_dist": 1000,
            "dt": 0.4,
            "time_corr": 10,
        }
    finally:
        shutil.rmtree(tempdir)

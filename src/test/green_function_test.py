import os
import pathlib
import shutil
import tempfile
from glob import glob

import pytest

from wasp.green_functions import fk_green_fun1, gf_retrieve

from .testutils import (
    RESULTS_DIR,
    get_cgps_json,
    get_sampling_filter,
    get_strong_motion_json,
    get_tele_waves_json,
    get_tensor_info,
)


@pytest.mark.skipif(
    os.getenv("CI_REGISTRY") is not None or os.getenv("RUN_ALL", False) == False,
    reason="Requires a 2+GB file to run",
)
def test_gf_retrieve_cgps():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "filtro_cgps.txt", tempdir / "filtro_cgps.txt")
        shutil.copyfile(
            RESULTS_DIR / "filtro_strong.txt", tempdir / "filtro_strong.txt"
        )
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
            RESULTS_DIR / "channels_cgps.txt", tempdir / "channels_cgps.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "wavelets_cgps.txt", tempdir / "wavelets_cgps.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "instrumental_response.txt",
            tempdir / "instrumental_response.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "Green_cgps.txt",
            tempdir / "Green_cgps.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "vel_model.txt", tempdir / "vel_model.txt"
        )
        with open(tempdir / "Green_cgps.txt", "a") as gf:
            gf.write(f"\n{str(tempdir/'GF_cgps')}\n")
        os.mkdir(tempdir / "logs")
        with open(tempdir / "logs" / "green_cgps_log", "w"):
            pass
        fdir = pathlib.Path(__file__).parent.parent.parent / "fortran_code"
        gf_retrieve(
            ["cgps"],
            {
                "tele_gf": fdir / "bin_inversion_gfortran_f95" / "green_tele",
                "strong_motion_gf": fdir / "bin_str_f95" / "get_strong_motion",
                "gps_gf": fdir / "src_dc_f95" / "gf_static_f95",
            },
            tempdir,
        )
        files = glob(str(tempdir) + "/*")
        ## TODO DETERMINE if test data should be added for the GRE/TDE files
        # in the meantime just check that the files were created
        for name, component in [
            (f["name"], f["component"]) for f in get_cgps_json(all=True)
        ]:
            assert f"{tempdir}/{name}.{component}" in files
    finally:
        shutil.rmtree(tempdir)


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


@pytest.mark.skipif(
    os.getenv("CI_REGISTRY") is not None or os.getenv("RUN_ALL", False) == False,
    reason="Build runner does not have the resources to run",
)
def test_gf_retrieve_insar():
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
        shutil.copyfile(RESULTS_DIR / "insar_data.txt", tempdir / "insar_data.txt")
        os.mkdir(tempdir / "logs")
        with open(tempdir / "logs" / "green_insar_log", "w"):
            pass
        fdir = pathlib.Path(__file__).parent.parent.parent / "fortran_code"
        gf_retrieve(
            ["insar"],
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
        assert str(tempdir / "Green_insar_subfault.txt") in files
    finally:
        shutil.rmtree(tempdir)


@pytest.mark.if(
    os.getenv("CI_REGISTRY") is not None or os.getenv("RUN_ALL", False) == False,
    reason="Requires a 2+GB file to run",
)
def test_gf_retrieve_strong_motion():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            RESULTS_DIR / "filtro_strong.txt", tempdir / "filtro_strong.txt"
        )
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
            RESULTS_DIR / "channels_strong.txt",
            tempdir / "channels_strong.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "wavelets_strong.txt",
            tempdir / "wavelets_strong.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "instrumental_response.txt",
            tempdir / "instrumental_response.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "Green_strong.txt",
            tempdir / "Green_strong.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "vel_model.txt", tempdir / "vel_model.txt"
        )
        with open(tempdir / "Green_strong.txt", "a") as gf:
            gf.write(f"\n{str(tempdir/'GF_strong')}\n")
        os.mkdir(tempdir / "logs")
        with open(tempdir / "logs" / "green_str_log", "w"):
            pass
        fdir = pathlib.Path(__file__).parent.parent.parent / "fortran_code"
        gf_retrieve(
            ["strong_motion"],
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
        for name, component in [
            (f["name"], f["component"]) for f in get_strong_motion_json(all=True)
        ]:
            assert f"{tempdir}/{name}.{component}" in files
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
            ["tele_body"],
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

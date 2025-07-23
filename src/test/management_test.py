import pathlib
import shutil
from tempfile import mkdtemp

import numpy as np
from obspy.taup import TauPyModel  # type:ignore

from wasp.management import (
    _distazbaz,
    coords2utm,
    correct_response_file,
    default_dirs,
    start_time_id,
    theoretic_arrivals,
    use_waveforms,
)


def test_correct_response_file():
    tempdir = mkdtemp()
    tensor_info = {"date_origin": "2023-07-08T12:10:09"}
    response_str = """START 2023-07-08T12:10:09Z
    START 2023-07-09T12:10:09Z
    START 2023-07-10T12:10:09Z
    END 2023-07-11T12:10:09Z
    END 2023-07-12T12:10:09Z
    END 2023-07-13T12:10:09Z
    NETWORK aa
    NETWORK bb
    NETWORK cc
    CONSTANT 11
    CONSTANT 22
    CONSTANT 33
    """
    try:
        response_file = pathlib.Path(tempdir) / "response.txt"
        with open(response_file, "w") as f:
            f.write(response_str)
        correct_response_file(tensor_info=tensor_info, pzfile=response_file)
        with open(response_file, "r") as f:
            updated_response = f.read()
        assert (
            updated_response.strip().replace("\t", "").replace(" ", "")
            == (
                """END 2023-07-13T12:10:09Z
        NETWORK aa
        NETWORK bb
        NETWORK cc
        CONSTANT 11"""
            )
            .strip()
            .replace("\t", "")
            .replace(" ", "")
        )
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_coofds2utm():
    coords2utm(32, 55, 40) == (1923.9229318287914, 3640.62059307422)


def test_distazbaz():
    assert _distazbaz(32, 55, 33, 56) == (
        145.35453183162141,
        220.41350123893108,
        39.876631986557655,
    )


def test_default_dirs():
    examples = pathlib.Path(__file__).parent / "data"
    config_path = examples / "config.ini"
    target_code_path = "/home/user/neic-finitefault"
    target_compute_near_gf = f"{target_code_path}/fortran_code/src_dc_f95"
    target_gf_bank = f"{target_code_path}/fortran_code/gfs_nm/long/low.in"
    target_info = f"{target_code_path}/fortran_code/info"
    target_modelling = f"{target_code_path}/fortran_code/bin_inversion_gfortran_f95"
    target_near_gf = f"{target_code_path}/fortran_code/bin_str_f95"
    target_tectonicplates = f"{target_code_path}/fortran_code/tectonicplates"
    assert default_dirs(config_path) == {
        "root_dir": target_code_path,
        "long_gf_bank": target_gf_bank,
        "crust_codes": f"{target_info}/CNtype2.txt",
        "models_codes": f"{target_info}/CNtype2_key.txt",
        "litho_model": f"{target_info}/LITHO1.0.nc",
        "gf_bank": target_gf_bank,
        "strong_motion_gf_bank2": f"{target_compute_near_gf}/green_bank_openmp_f95",
        "strong_motion_gf": f"{target_near_gf}/get_strong_motion",
        "cgps_gf_bank": f"{target_near_gf}/cgps",
        "gps_gf": f"{target_compute_near_gf}/gf_static_f95",
        "tele_gf": f"{target_modelling}/green_tele",
        "finite_fault": f"{target_modelling}/run_modelling",
        "forward": f"{target_modelling}/run_forward",
        "trench_graphics": f"{target_tectonicplates}/PB2002_plates",
    }


def test_start_time_id():
    assert start_time_id({"datetime": "2023-07-08T12:10:09"}) == "20230708121009"


def test_theoretical_arrivals():
    model = TauPyModel(model="ak135f_no_mud")
    arrivals = theoretic_arrivals(model, 20, 35)
    np.testing.assert_almost_equal(
        arrivals["p_arrival"][0].time, 269.4806309224751, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["p_arrival"][1].time, 271.3779514441007, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["p_arrival"][2].time, 271.5371397639714, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["p_arrival"][3].time, 274.5748317621238, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["p_arrival"][4].time, 274.9422784875743, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["pp_arrival"][0].time, 286.0118331807067, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["s_arrival"][0].time, 492.4317556718514, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["s_arrival"][1].time, 495.0629873974531, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["s_arrival"][2].time, 495.27848920935264, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["s_arrival"][3].time, 496.01914493591386, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["s_arrival"][4].time, 496.49098362769223, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["s_arrival"][5].time, 500.20455604428753, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["s_arrival"][6].time, 501.3508832952396, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["ss_arrival"][0].time, 509.84935358894506, decimal=5
    )
    np.testing.assert_almost_equal(
        arrivals["p_slowness"], 10.875067057735919, decimal=5
    )
    np.testing.assert_almost_equal(arrivals["s_slowness"], 19.92570580559402, decimal=5)


def test_use_waveforms():
    assert use_waveforms("strong") == True
    assert use_waveforms("cgps") == True
    assert use_waveforms("body") == True
    assert use_waveforms("surf") == True
    assert use_waveforms("other") == False

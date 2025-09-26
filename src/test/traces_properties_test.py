import json
import pathlib
import shutil
import tempfile

from wasp.traces_properties import (
    filter_strong,
    filter_surf,
    filter_tele,
    nyquist_frequency,
    properties_json,
    sampling,
    wavelet_scales,
)


def test_filter_tele():
    assert filter_tele({"time_shift": 5, "depth": 22.4}) == {
        "freq0": 0.003,
        "freq3": 1.2,
        "high_freq": 1.0,
        "low_freq": 0.01,
    }
    assert filter_tele({"time_shift": 15, "depth": 22.4}) == {
        "freq0": 0.003,
        "freq3": 1.2,
        "high_freq": 1.0,
        "low_freq": 0.006,
    }
    assert filter_tele({"time_shift": 40, "depth": 22.4}) == {
        "freq0": 0.002,
        "freq3": 1.2,
        "high_freq": 1.0,
        "low_freq": 0.004,
    }
    assert filter_tele({"time_shift": 90, "depth": 210}) == {
        "freq0": 0.001,
        "freq3": 0.9,
        "high_freq": 0.8,
        "low_freq": 0.002,
    }


def test_filter_surf():
    assert filter_surf() == {
        "freq1": 0.003,
        "freq2": 0.004,
        "freq3": 0.006,
        "freq4": 0.007,
    }


def test_filter_strong():
    assert filter_strong(
        {"time_shift": 5, "depth": 22.4},
    ) == {"high_freq": 0.125, "low_freq": 0.02}
    assert filter_strong(
        {"time_shift": 15, "depth": 22.4},
    ) == {"high_freq": 0.125, "low_freq": 0.01}
    assert filter_strong(
        {"time_shift": 35, "depth": 22.4},
    ) == {"high_freq": 0.125, "low_freq": 0.01}
    assert filter_strong(
        {"time_shift": 65, "depth": 22.4},
    ) == {"high_freq": 0.125, "low_freq": 0.01}
    assert filter_strong({"time_shift": 105, "depth": 22.4}, cgnss=True) == {
        "high_freq": 0.3,
        "low_freq": 0.01,
    }


def test_nyquist_frequency():
    assert nyquist_frequency(2000) == 0.00025
    assert nyquist_frequency(10) == 0.05
    assert nyquist_frequency(0.02) == 25


def test_properties_json():
    tempdir = tempfile.mkdtemp()
    try:
        d1 = properties_json(
            {"time_shift": 5, "depth": 22.4, "moment_mag": 2.9e30},
            0.2,
            data_directory=tempdir,
        )
        target = {
            "sampling": {"dt_tele": 0.4, "dt_strong": 0.2, "dt_cgnss": 0.2},
            "tele_filter": {
                "freq0": 0.003,
                "low_freq": 0.01,
                "high_freq": 1.0,
                "freq3": 1.2,
            },
            "surf_filter": {
                "freq1": 0.003,
                "freq2": 0.004,
                "freq3": 0.006,
                "freq4": 0.007,
            },
            "strong_filter": {"low_freq": 0.02, "high_freq": 0.125},
            "cgnss_filter": {"low_freq": 0.02, "high_freq": 0.3},
            "wavelet_scales": [1, 8],
        }
        with open(pathlib.Path(tempdir) / "sampling_filter.json", "r") as f:
            d2 = json.load(f)
        assert d1 == d2 == target
    finally:
        shutil.rmtree(tempdir)


def test_sampling():
    assert sampling({"moment_mag": 2.9e21, "depth": 22.4, "time_shift": 10}) == {
        "dt_cgnss": 0.2,
        "dt_strong": 0.2,
        "dt_tele": 0.2,
    }
    assert sampling({"moment_mag": 2.9e21, "depth": 300, "time_shift": 10}) == {
        "dt_cgnss": 0.2,
        "dt_strong": 0.2,
        "dt_tele": 0.4,
    }
    assert sampling({"moment_mag": 2.9e21, "depth": 500, "time_shift": 10}) == {
        "dt_cgnss": 0.8,
        "dt_strong": 0.8,
        "dt_tele": 0.5,
    }
    assert sampling({"moment_mag": 2.9e21, "depth": 22.4, "time_shift": 10}, 20) == {
        "dt_cgnss": 20,
        "dt_strong": 0.2,
        "dt_tele": 0.2,
    }
    assert sampling({"moment_mag": 2.9e21, "depth": 22.4, "time_shift": 50}) == {
        "dt_cgnss": 0.4,
        "dt_strong": 0.4,
        "dt_tele": 0.2,
    }
    assert sampling({"moment_mag": 2.9e21, "depth": 22.4, "time_shift": 100}) == {
        "dt_cgnss": 0.8,
        "dt_strong": 0.8,
        "dt_tele": 0.2,
    }


def test_wavelet_scales():
    assert wavelet_scales() == [1, 8]

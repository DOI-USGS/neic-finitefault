import json
import math
import pathlib
from copy import deepcopy
from typing import List, Optional

import numpy as np

from ffm.seismic_tensor import get_tensor

HOME = pathlib.Path(__file__).parent.parent.parent
DATA_DIR = pathlib.Path(__file__).parent / "data"
END_TO_END_DIR = DATA_DIR / "end_to_end"
RESULTS_DIR = END_TO_END_DIR / "results"


def assert_values_close(new, target, rtol: float = 1e-8):
    """Recursively assert equality of (possibly nested) values, comparing
    floats with a relative tolerance to absorb tiny platform-dependent
    differences (e.g. x86 vs arm64 libm/LAPACK). Everything else is exact."""
    if isinstance(new, dict) and isinstance(target, dict):
        assert new.keys() == target.keys(), (new.keys(), target.keys())
        for key in new:
            assert_values_close(new[key], target[key], rtol=rtol)
    elif isinstance(new, (list, tuple)) and isinstance(target, (list, tuple)):
        assert len(new) == len(target), (new, target)
        for new_item, target_item in zip(new, target):
            assert_values_close(new_item, target_item, rtol=rtol)
    elif isinstance(new, float) and isinstance(target, float):
        np.testing.assert_allclose(new, target, rtol=rtol)
    else:
        assert new == target


def assert_text_close(new_text: str, target_text: str, rtol: float = 1e-8):
    """Assert two text blobs are equal, allowing tiny platform-dependent
    floating point differences (e.g. x86 vs arm64 libm/LAPACK) in numeric
    tokens. Non-numeric tokens must match exactly."""
    new_lines = new_text.splitlines()
    target_lines = target_text.splitlines()
    assert len(new_lines) == len(target_lines)
    for new_line, target_line in zip(new_lines, target_lines):
        new_tokens = new_line.split()
        target_tokens = target_line.split()
        assert len(new_tokens) == len(target_tokens), (new_line, target_line)
        for new_tok, target_tok in zip(new_tokens, target_tokens):
            if new_tok == target_tok:
                continue
            np.testing.assert_allclose(
                float(new_tok),
                float(target_tok),
                rtol=rtol,
                err_msg=f"'{new_line}' vs '{target_line}'",
            )


def assert_solution_equivalent(
    new_text: str, target_text: str, mw_tol: float = 0.1, geometry_rtol: float = 1e-4
):
    """Assert two Solution.txt files describe equivalent inversions.

    Simulated annealing trajectories diverge across platforms/BLAS backends,
    so per-subfault slip, rake, rupture time and moment cannot be compared.
    Instead require: identical structure, matching fault geometry (lat, lon,
    depth, strike, dip), and total moment magnitude Mw within mw_tol."""
    new_lines = new_text.splitlines()
    target_lines = target_text.splitlines()
    assert len(new_lines) == len(target_lines)
    new_moment = 0.0
    target_moment = 0.0
    # columns of subfault lines: lat lon depth slip rake strike dip t_rup t_ris t_fal mo
    geometry_columns = [0, 1, 2, 5, 6]
    for new_line, target_line in zip(new_lines, target_lines):
        new_tokens = new_line.split()
        target_tokens = target_line.split()
        assert len(new_tokens) == len(target_tokens), (new_line, target_line)
        try:
            new_values = [float(t) for t in new_tokens]
            target_values = [float(t) for t in target_tokens]
        except ValueError:
            # header/label line: numeric tokens compared loosely, text exactly
            assert_text_close(new_line, target_line, rtol=geometry_rtol)
            continue
        if len(new_values) == 11:
            for col in geometry_columns:
                np.testing.assert_allclose(
                    new_values[col],
                    target_values[col],
                    rtol=geometry_rtol,
                    err_msg=f"'{new_line}' vs '{target_line}'",
                )
            new_moment += new_values[10]
            target_moment += target_values[10]
        else:
            # other numeric lines (e.g. segment corner coordinates) are geometry
            np.testing.assert_allclose(
                new_values, target_values, rtol=geometry_rtol, err_msg=new_line
            )
    new_mw = (2.0 / 3) * (math.log10(new_moment) - 16.1)
    target_mw = (2.0 / 3) * (math.log10(target_moment) - 16.1)
    assert abs(new_mw - target_mw) < mw_tol, (new_mw, target_mw)


class MockResponse:
    def __init__(self, text: str, status: int):
        self.text = text
        self.status_code = status

    def json(self):
        return json.loads(self.text)

    def raise_for_status(self):
        if self.status_code != 200:
            raise Exception()


def get_cgnss_json(min: int = 0, max: int = 3, all: bool = False):
    """Get cgnss_waves.json"""
    with open(RESULTS_DIR / "NP1" / "cgnss_waves.json", "r") as f:
        d = json.load(f)
        if not all:
            d = d[min:max]
        cgnss_waves = cgnss_waves = update_manager_file_locations(
            d, RESULTS_DIR / "data"
        )
    return cgnss_waves


def get_imagery_json():
    """Get imagery_data.json"""
    with open(RESULTS_DIR / "NP1" / "imagery_data.json", "r") as f:
        imagery_data = cgnss_waves = update_manager_file_locations(
            json.load(f), RESULTS_DIR / "NP1", file_key="name"
        )
    return imagery_data


def get_sampling_filter():
    """Get sampling_filter.json"""
    with open(RESULTS_DIR / "NP1" / "sampling_filter.json", "r") as f:
        sampling_filter = json.load(f)
    return sampling_filter


def get_segments_data():
    """Get segments_data.json"""
    with open(RESULTS_DIR / "NP1" / "segments_data.json", "r") as f:
        segments = json.load(f)
    return segments


def get_multisegments_data():
    """Get multisegments_data.json"""
    with open(RESULTS_DIR / "NP1" / "multisegments_data.json", "r") as f:
        segments = json.load(f)
    return segments


def get_static_json(min: int = 0, max: int = 3, all: bool = False):
    """Get static_data.json"""
    with open(RESULTS_DIR / "NP1" / "static_data.json", "r") as f:
        d = json.load(f)
        if not all:
            d = d[min:max]
        static_data = d
    return static_data


def get_strong_motion_json(min: int = 0, max: int = 3, all: bool = False):
    """Get strong_motion_waves.json"""
    with open(RESULTS_DIR / "NP1" / "strong_motion_waves.json", "r") as f:
        d = json.load(f)
        if not all:
            d = d[min:max]
        sm_waves = update_manager_file_locations(d, RESULTS_DIR / "data")
    return sm_waves


def get_surf_waves_json(min: int = 0, max: int = 3, all: bool = False):
    """Get surf_waves.json"""
    with open(RESULTS_DIR / "NP1" / "surf_waves.json", "r") as f:
        d = json.load(f)
        if not all:
            d = d[min:max]
        surf_waves = update_manager_file_locations(d, RESULTS_DIR / "data")
    return surf_waves


def get_tele_waves_json(min: int = 0, max: int = 3, all: bool = False):
    """Get tele_waves.json"""
    with open(RESULTS_DIR / "NP1" / "tele_waves.json", "r") as f:
        d = json.load(f)
        if not all:
            d = d[min:max]
        tele_waves = update_manager_file_locations(d, RESULTS_DIR / "data")
    return tele_waves


def get_tensor_info():
    """Get the tensor information"""
    return get_tensor(cmt_file=END_TO_END_DIR / "info" / "20003k7a_cmt_CMT")


def update_manager_file_locations(
    manage_dicts: List[dict],
    new_directory: pathlib.Path,
    replace_dir: Optional[str] = None,
    file_key: str = "file",
) -> List[dict]:
    """Update file locations for testing"""
    updated_managed_dicts = deepcopy(manage_dicts)
    if isinstance(updated_managed_dicts, list):
        for idx, d in enumerate(manage_dicts):
            if replace_dir is not None:
                updated_managed_dicts[idx][file_key] = manage_dicts[idx][
                    file_key
                ].replace(replace_dir, str(new_directory))
            else:
                updated_managed_dicts[idx][file_key] = str(
                    new_directory / manage_dicts[idx][file_key]
                )
    elif isinstance(updated_managed_dicts, dict):
        for key in updated_managed_dicts:
            for idx, d in enumerate(manage_dicts[key]):
                if replace_dir is not None:
                    updated_managed_dicts[key][idx][file_key] = manage_dicts[key][idx][
                        file_key
                    ].replace(replace_dir, str(new_directory))
                else:
                    updated_managed_dicts[key][idx][file_key] = str(
                        new_directory / manage_dicts[key][idx][file_key]
                    )
    return updated_managed_dicts


def get_velmodel_data():
    """Get velmodel_data.json"""
    with open(RESULTS_DIR / "NP1" / "velmodel_data.json", "r") as f:
        velmodel = json.load(f)
    return velmodel

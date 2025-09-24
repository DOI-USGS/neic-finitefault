import json
import pathlib
from copy import deepcopy
from typing import List, Optional

from wasp.seismic_tensor import get_tensor

HOME = pathlib.Path(__file__).parent.parent.parent
DATA_DIR = pathlib.Path(__file__).parent / "data"
END_TO_END_DIR = DATA_DIR / "end_to_end"
RESULTS_DIR = END_TO_END_DIR / "results"


def get_cgps_json(min: int = 0, max: int = 3, all: bool = False):
    """Get cgps_waves.json"""
    with open(RESULTS_DIR / "NP1" / "cgps_waves.json", "r") as f:
        d = json.load(f)
        if not all:
            d = d[min:max]
        cgps_waves = cgps_waves = update_manager_file_locations(d, RESULTS_DIR / "data")
    return cgps_waves


def get_insar_json():
    """Get insar_data.json"""
    with open(RESULTS_DIR / "NP1" / "insar_data.json", "r") as f:
        insar_data = cgps_waves = update_manager_file_locations(
            json.load(f), RESULTS_DIR / "NP1", file_key="name"
        )
    return insar_data


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

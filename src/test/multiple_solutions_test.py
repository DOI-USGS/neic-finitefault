import json
import os
import pathlib
import shutil
import tempfile
from copy import deepcopy
from unittest import mock

from wasp.multiple_solutions import (
    __worker,
    force_plane_above_ground,
    get_summary,
    get_summary_all_models,
    multiple_solutions,
)

from .testutils import HOME, get_segments_data, get_tensor_info

SEGMENTS = get_segments_data()["segments"]
TENSOR = get_tensor_info()


@mock.patch(target="wasp.inversion_chen_new.manual_modelling", return_value=None)
def test___worker(mock_inversion):
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        folder1 = tempdir / "folder1"
        os.mkdir(folder1)
        with open(folder1 / "segments_data.json", "w") as f:
            json.dump(
                {
                    "segments": [
                        {
                            "strike": 1,
                            "dip": 2,
                            "rupture_vel": 3,
                        },
                    ]
                },
                f,
            )
        __worker(TENSOR, "cgnss", {}, folder1)
    finally:
        shutil.rmtree(tempdir)


def test_force_plane_above_ground():
    updated_segments1 = force_plane_above_ground(TENSOR, SEGMENTS)
    assert updated_segments1 == SEGMENTS

    # test where plane intersects with the surface
    bad_segments = deepcopy(SEGMENTS)
    bad_segments[0]["dip_subfaults"] = 100
    updated_segments2 = force_plane_above_ground(TENSOR, bad_segments)
    assert updated_segments2 == [
        {
            "delay_segment": 0,
            "delta_dip": 13.431917887412927,
            "delta_strike": 17.92760869565217,
            "dip": 19.280827965117993,
            "dip_subfaults": 100,
            "hyp_dip": 5,
            "hyp_stk": 9,
            "neighbours": [],
            "rake": 109.27817171619564,
            "rupture_vel": 2.5,
            "min_vel": 1,
            "max_vel": 3.125,
            "stk_subfaults": 23,
            "strike": 6.613912311529926,
        }
    ]


def test_get_summary():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        with open(tempdir / "modelling_summary.txt", "w") as f:
            f.write(
                """Modelling Report

                averaged misfit error  0.1234
                objective function value  0.5678   
                """
            )
        assert get_summary(directory=tempdir) == {
            "misfit_error": 0.1234,
            "objective_error": 0.5678,
        }
    finally:
        shutil.rmtree(tempdir)
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        folder1 = tempdir / "folder1"
        folder2 = tempdir / "folder2"
        os.mkdir(folder1)
        os.mkdir(folder2)

    finally:
        shutil.rmtree(tempdir)


def test_get_summary_all_models():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        folder1 = tempdir / "folder1"
        folder2 = tempdir / "folder2"
        os.mkdir(folder1)
        os.mkdir(folder2)
        with open(folder1 / "modelling_summary.txt", "w") as f:
            f.write(
                """Modelling Report

                averaged misfit error  0.1234
                objective function value  0.5678   
                """
            )
        with open(folder2 / "modelling_summary.txt", "w") as f:
            f.write(
                """Modelling Report

                averaged misfit error  0.01911
                objective function value  0.121314   
                """
            )
        with open(folder1 / "segments_data.json", "w") as f:
            json.dump(
                {
                    "segments": [
                        {
                            "strike": 1,
                            "dip": 2,
                            "rupture_vel": 3,
                        },
                    ]
                },
                f,
            )
        with open(folder2 / "segments_data.json", "w") as f:
            json.dump(
                {
                    "segments": [
                        {
                            "strike": 4,
                            "dip": 5,
                            "rupture_vel": 6,
                        },
                    ]
                },
                f,
            )
        summary = get_summary_all_models(
            subfolders=[str(folder1), str(folder2)]
        ).to_dict()
        del summary["subfolders"]
        assert summary == {
            "strike": {0: 1, 1: 4},
            "dip": {0: 2, 1: 5},
            "rupt_vel": {0: 3, 1: 6},
            "objective_error": {0: 0.5678, 1: 0.121314},
            "misfit_error": {0: 0.1234, 1: 0.01911},
        }
    finally:
        shutil.rmtree(tempdir)


def _mock_worker():
    pass


def test_multiple_solutions():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        os.chdir(tempdir)
        folder1 = tempdir / "folder1"
        os.mkdir(folder1)
        with open(folder1 / "segments_data.json", "w") as f:
            json.dump(get_segments_data(), f)
        multiple_solutions(
            TENSOR,
            ["cgnss"],
            {},
            [folder1],
            strike=[1],
            dip=[1],
            rupt_vel=[1],
            directory=folder1,
            test=True,
        )
    finally:
        os.chdir(HOME)
        shutil.rmtree(tempdir)

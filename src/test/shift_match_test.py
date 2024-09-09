import json
import os
import pathlib
import shutil
import tempfile

from wasp.shift_match import print_arrival, shift_match2, shift_match_regional

from .testutils import (
    RESULTS_DIR,
    get_cgps_json,
    get_strong_motion_json,
    get_surf_waves_json,
    get_tele_waves_json,
    get_tensor_info,
    update_manager_file_locations,
)


def test_shift_match2():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_surf_waves = update_manager_file_locations(
            get_surf_waves_json(),
            tempdir / "data",
            replace_dir=str(RESULTS_DIR / "data"),
        )
        new_tele_waves = update_manager_file_locations(
            get_tele_waves_json(),
            tempdir / "data",
            replace_dir=str(RESULTS_DIR / "data"),
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "synthetics_body.txt", tempdir / "synthetics_body.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "synthetics_surf.txt", tempdir / "synthetics_surf.txt"
        )
        for f, d in zip(
            [
                "tele_waves.json",
                "surf_waves.json",
            ],
            [new_tele_waves, new_surf_waves],
        ):
            with open(tempdir / f, "w") as w:
                json.dump(d, w)
        os.mkdir(tempdir / "data")
        os.mkdir(tempdir / "data" / "SH")
        os.mkdir(tempdir / "data" / "P")
        os.mkdir(tempdir / "data" / "LONG")
        for (
            a,
            b,
            c,
            d,
        ) in zip(
            get_surf_waves_json(),
            new_surf_waves,
            get_tele_waves_json(),
            new_tele_waves,
        ):
            shutil.copyfile(a["file"], b["file"])
            shutil.copyfile(c["file"], d["file"])
        # TODO add example baseline data, for now make sure runs without error
        shift_match2("body", plot=True, directory=tempdir)
        for d in new_tele_waves:
            assert (tempdir / "tele_shift" / f"{d['name']}_{d['component']}.png").exists
        shift_match2("surf", plot=True, directory=tempdir)
        for d in new_surf_waves:
            assert (tempdir / "surf_shift" / f"{d['name']}_{d['component']}.png").exists
    finally:
        shutil.rmtree(tempdir)


def test_shift_match_regional():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_strong = update_manager_file_locations(
            get_strong_motion_json(),
            tempdir / "data",
            replace_dir=str(RESULTS_DIR / "data"),
        )
        new_cgps = update_manager_file_locations(
            get_cgps_json(),
            tempdir / "data",
            replace_dir=str(RESULTS_DIR / "data"),
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "synthetics_strong.txt",
            tempdir / "synthetics_strong.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "synthetics_cgps.txt", tempdir / "synthetics_cgps.txt"
        )
        for f, d in zip(
            [
                "strong_motion_waves.json",
                "cgps_waves.json",
            ],
            [new_strong, new_cgps],
        ):
            with open(tempdir / f, "w") as w:
                json.dump(d, w)
        os.mkdir(tempdir / "data")
        os.mkdir(tempdir / "data" / "STR")
        os.mkdir(tempdir / "data" / "cGPS")
        for (
            a,
            b,
            c,
            d,
        ) in zip(
            get_strong_motion_json(),
            new_strong,
            get_cgps_json(),
            new_cgps,
        ):
            shutil.copyfile(a["file"], b["file"])
            shutil.copyfile(c["file"], d["file"])
        # TODO add example baseline data, for now make sure runs without error
        shift_match_regional("strong", plot=True, directory=tempdir)
        for d in new_strong:
            assert (
                tempdir / "strong_shift" / f"{d['name']}_{d['component']}.png"
            ).exists
        shift_match_regional("cgps", plot=True, directory=tempdir)
        for d in new_cgps:
            assert (tempdir / "cgps_shift" / f"{d['name']}_{d['component']}.png").exists
    finally:
        shutil.rmtree(tempdir)


def test__print_arrival():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_tele_waves = update_manager_file_locations(
            get_tele_waves_json(),
            tempdir / "data",
            replace_dir=str(RESULTS_DIR / "data" / "SH"),
        )
        new_tele_waves = update_manager_file_locations(
            new_tele_waves,
            tempdir / "data",
            replace_dir=str(RESULTS_DIR / "data" / "P"),
        )

        for f, d in zip(
            [
                "tele_waves.json",
            ],
            [new_tele_waves],
        ):
            with open(tempdir / f, "w") as w:
                json.dump(d, w)
        print(new_tele_waves)
        os.mkdir(tempdir / "data")
        for (
            a,
            b,
        ) in zip(
            get_tele_waves_json(),
            new_tele_waves,
        ):
            shutil.copyfile(a["file"], b["file"])
        print_arrival(get_tensor_info(), directory=tempdir)
        for d in new_tele_waves:
            assert (tempdir / "tele_arrival" / f"{d['name']}_pick.png").exists
    finally:
        shutil.rmtree(tempdir)

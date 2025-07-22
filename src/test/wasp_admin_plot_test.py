import json
import os
import pathlib
import shutil
import tempfile
from glob import glob
from test.testutils import (
    DATA_DIR,
    END_TO_END_DIR,
    HOME,
    RESULTS_DIR,
    get_cgps_json,
    get_strong_motion_json,
    get_surf_waves_json,
    get_tele_waves_json,
    update_manager_file_locations,
)

from typer.testing import CliRunner

runner = CliRunner()


def test_kml():
    from wasp.wasp_admin.plot import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        with open(DATA_DIR / "config.ini") as f:
            config = f.read().replace("/home/user/neic-finitefault", str(HOME))
        with open(tempdir / "config.ini", "w") as wf:
            wf.write(config)
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solution.txt", tempdir / "Solution.txt")
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
        os.mkdir(tempdir / "data" / "RAYLEIGH")
        os.mkdir(tempdir / "data" / "LOVE")
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

        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )

        # test kml
        result = runner.invoke(
            app,
            [
                "kml",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "cgps",
                "-t",
                "surf",
                "-t",
                "strong",
                "-t",
                "body",
                "-c",
                str(tempdir / "config.ini"),
                "-e",
                "2003k7a",
            ],
        )
        assert result.exit_code == 0
        # validate results
        print(glob(str(tempdir) + "/*"))
        assert (tempdir / "2003k7a.kml").exists()
        assert (tempdir / "Map_kml.png").exists()
        assert (tempdir / "Map_kml.ps").exists()
    finally:
        shutil.rmtree(tempdir)


def test_map():
    from wasp.wasp_admin.plot import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        with open(DATA_DIR / "config.ini") as f:
            config = f.read().replace("/home/user/neic-finitefault", str(HOME))
        with open(tempdir / "config.ini", "w") as wf:
            wf.write(config)
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solution.txt", tempdir / "Solution.txt")
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
        os.mkdir(tempdir / "data" / "RAYLEIGH")
        os.mkdir(tempdir / "data" / "LOVE")
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

        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )

        # test map
        result = runner.invoke(
            app,
            [
                "map",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "cgps",
                "-t",
                "surf",
                "-t",
                "strong",
                "-t",
                "body",
                "-c",
                str(tempdir / "config.ini"),
            ],
        )
        print(result.exception)
        assert result.exit_code == 0
        # validate results
        assert (tempdir / "PyGMT_Map.png").exists()
    finally:
        shutil.rmtree(tempdir)


def test_neic():
    from wasp.wasp_admin.plot import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        with open(DATA_DIR / "config.ini") as f:
            config = f.read().replace("/home/user/neic-finitefault", str(HOME))
        with open(tempdir / "config.ini", "w") as wf:
            wf.write(config)
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solution.txt", tempdir / "Solution.txt")
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
        os.mkdir(tempdir / "data" / "RAYLEIGH")
        os.mkdir(tempdir / "data" / "LOVE")
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
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tensor_info.json", tempdir / "tensor_info.json"
        )

        # test neic
        result = runner.invoke(
            app,
            [
                "neic",
                str(tempdir),
                "-t",
                "cgps",
                "-t",
                "surf",
                "-t",
                "strong",
                "-t",
                "body",
                "-ffms",
                "-c",
                str(tempdir / "config.ini"),
            ],
        )
        print(result.exception)
        assert result.exit_code == 0
        # validate results
        plot_dir = tempdir / "plots"
        assert (plot_dir / "MomentRate.png").exists()
        assert (plot_dir / "Map.png").exists()
        assert (tempdir / "Map.pdf").exists()
        assert (plot_dir / "SlipDist_plane0.png").exists()
        assert (tempdir / "Map.eps").exists()
        assert (tempdir / "SlipDist_plane0.ps").exists()
        assert (plot_dir / "Rayleigh_surf_waves.png").exists()
        assert (plot_dir / "cGPS_waves.png").exists()
        assert (plot_dir / "strong_motion_waves.png").exists()
        assert (plot_dir / "SH_body_waves.png").exists()
        assert (plot_dir / "Love_surf_waves.png").exists()
        assert (plot_dir / "P_body_waves.png").exists()
    finally:
        shutil.rmtree(tempdir)

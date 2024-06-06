import json
import os
import pathlib
import shutil
import tempfile
from os import mkdir

import pytest

from .testutils import (
    END_TO_END_DIR,
    RESULTS_DIR,
    get_surf_waves_json,
    get_tele_waves_json,
)

RESULTS_DIR = pathlib.Path(__file__).parent / "data" / "end_to_end" / "results"
from wasp.eventpage_downloads import (
    make_waveproperties_json,
    temporary_file_reorganization_for_publishing,
    write_CMTSOLUTION_file,
    write_Coulomb_file,
    write_Okada_displacements,
)


def test_make_waveproperties_json():
    tempdir = tempfile.mkdtemp()
    try:
        tempdir = pathlib.Path(tempdir)
        with open(tempdir / "tele_waves.json", "w") as f:
            json.dump(get_tele_waves_json(all=True), f)
        with open(tempdir / "surf_waves.json", "w") as f:
            json.dump(get_surf_waves_json(all=True), f)
        make_waveproperties_json(tempdir)
        with open(tempdir / "wave_properties.json", "r") as f:
            wp = json.load(f)
        with open(RESULTS_DIR / "NP1" / "wave_properties.json") as f:
            targetwp = json.load(f)
        assert wp == targetwp
    finally:
        shutil.rmtree(tempdir)


def test_temporary_file_reorganization_for_publishing():
    files = [
        "CMTSOLUTION",
        "Coulomb.inp",
        "STF.txt",
        "Solucion.txt",
        ["plots", "MomentRate.png"],
        "insar_ascending.txt",
        "insar_descending.txt",
        "fsp_sol_file.txt",
        "shakemap_polygon.txt",
        "surface_deformation.disp",
        "wave_properties.json",
        ["plots", "Map.png"],
        ["plots", "SlipDist_plane0.png"],
        ["plots", "test_waves.png"],
        ["plots", "InSAR_test_fit.png"],
    ]
    tempdir = tempfile.mkdtemp()
    try:
        tempdir = pathlib.Path(tempdir)
        for f in files:
            if isinstance(f, list):
                dir = tempdir / f[0]
                if not dir.exists():
                    mkdir(tempdir / f[0])
                path = dir / f[1]
            else:
                path = tempdir / f
            with open(path, "w") as f:
                pass
        temporary_file_reorganization_for_publishing("eventid", tempdir)
        pub_dir = tempdir / "PublicationFiles"
        assert (pub_dir / "eventid.param").exists()
        assert (pub_dir / "eventid.fsp").exists()
        assert (pub_dir / "shakemap_polygon.txt").exists()
        assert (pub_dir / "surface_deformation.disp").exists()
        assert (pub_dir / "eventid_coulomb.inp").exists()
        assert (pub_dir / "CMTSOLUTION").exists()
        assert (pub_dir / "wave_properties.json").exists()
        assert (pub_dir / "eventid.mr").exists()
        assert (pub_dir / "mr.png").exists()
        assert (pub_dir / "eventid_basemap.png").exists()
        assert (pub_dir / "eventid_slip2.png").exists()
        assert (pub_dir / "fits").exists()
        assert (pub_dir / "fits.zip").exists()

    finally:
        shutil.rmtree(tempdir)


def test_write_CMTSOLUTION_file():
    tempdir = tempfile.mkdtemp()
    try:
        tempdir = pathlib.Path(tempdir)
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solucion.txt", tempdir / "Solucion.txt")
        write_CMTSOLUTION_file(END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir)
        with open(tempdir / "CMTSOLUTION", "r") as f:
            cmt = f.read()
        with open(RESULTS_DIR / "NP1" / "CMTSOLUTION") as f:
            cmt_target = f.read()
        assert cmt == cmt_target
    finally:
        shutil.rmtree(tempdir)


def test_write_Coulomb_file():
    tempdir = tempfile.mkdtemp()
    try:
        tempdir = pathlib.Path(tempdir)
        for f in ["fsp_sol_file.txt", "Solucion.txt"]:
            shutil.copyfile(RESULTS_DIR / "NP1" / f, tempdir / f)
        write_Coulomb_file(tempdir, "us20003k7a")
        with open(tempdir / "Coulomb.inp", "r") as f:
            cmt = f.read()
        with open(RESULTS_DIR / "NP1" / "Coulomb.inp") as f:
            cmt_target = f.read()
        assert cmt == cmt_target
    finally:
        shutil.rmtree(tempdir)


@pytest.mark.skipif(
    os.getenv("CI_REGISTRY") is not None,
    reason="Build runner does not have the resources to run",
)
def test_write_Okada_displacements():
    tempdir = tempfile.mkdtemp()
    try:
        tempdir = pathlib.Path(tempdir)
        for f in ["fsp_sol_file.txt", "Solucion.txt"]:
            shutil.copyfile(RESULTS_DIR / "NP1" / f, tempdir / f)
        write_Okada_displacements(END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir)
        with open(tempdir / "surface_deformation.disp", "r") as f:
            disp = f.read()
        with open(RESULTS_DIR / "NP1" / "surface_deformation.disp") as f:
            disp_target = f.read()
        assert disp == disp_target
        assert (tempdir / "Okada_Displacement.png").exists()
        # assert (tempdir / "Vertical_Surface_Displacement.png").exists()
    finally:
        shutil.rmtree(tempdir)

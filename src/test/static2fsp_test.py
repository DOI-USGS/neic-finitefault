import pathlib
import shutil
import tempfile

from wasp import get_outputs
from wasp.static2fsp import static_to_fsp

from .testutils import (
    RESULTS_DIR,
    get_segments_data,
    get_tensor_info,
    get_velmodel_data,
)


def test_static_to_fsp():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solucion.txt", tempdir / "Solucion.txt")
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "strong_motion_waves.json",
            tempdir / "strong_motion_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "cgps_waves.json",
            tempdir / "cgps_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "static_data.json",
            tempdir / "static_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tele_waves.json",
            tempdir / "tele_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "surf_waves.json",
            tempdir / "surf_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.json",
            tempdir / "insar_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.txt",
            tempdir / "insar_data.txt",
        )
        segments = get_segments_data()
        solution = get_outputs.read_solution_static_format(
            segments["segments"], tempdir
        )
        static_to_fsp(
            get_tensor_info(),
            segments,
            ["cgps", "gps", "insar", "strong", "surf", "body"],
            get_velmodel_data(),
            solution,
            tempdir,
        )
        with open(tempdir / "fsp_sol_file.txt") as f:
            data = f.readlines()[41:]
            dstr = f.read()
        with open(RESULTS_DIR / "NP1" / "fsp_sol_file.txt") as f:
            target = f.readlines()[41:]

        for l, t in zip(data, target):
            assert l == t
    finally:
        shutil.rmtree(tempdir)

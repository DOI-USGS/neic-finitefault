import pathlib
import shutil
import tempfile

from wasp import get_outputs
from wasp.fault_plane import point_sources_param
from wasp.write_KML import _PlotMap_KML

from .testutils import END_TO_END_DIR, RESULTS_DIR, get_segments_data, get_tensor_info

SEGMENTS = get_segments_data()
TENSOR = get_tensor_info()


def test_write_KML():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solucion.txt", tempdir / "Solucion.txt")
        point_sources = point_sources_param(
            SEGMENTS["segments"],
            TENSOR,
            SEGMENTS["rise_time"],
        )
        solution = get_outputs.read_solution_static_format(
            SEGMENTS["segments"], tempdir
        )
        _PlotMap_KML(
            TENSOR,
            SEGMENTS["segments"],
            point_sources,
            solution,
            {
                "trench_graphics": pathlib.Path(__file__).parent.parent.parent
                / "fortran_code"
                / "tectonicplates"
                / "PB2002_plates"
            },
            evID="20003k7a",
            directory=tempdir,
        )
        # TODO: add test data to check against; check exists for now
        assert (tempdir / "20003k7a.kml").exists()
        assert (tempdir / "Map_kml.png").exists()
        assert (tempdir / "Map_kml.ps").exists()
    finally:
        shutil.rmtree(tempdir)

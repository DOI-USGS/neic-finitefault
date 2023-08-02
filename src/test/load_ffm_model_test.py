import pathlib
import shutil
import tempfile

import numpy as np

from wasp.fault_plane import point_sources_param
from wasp.get_outputs import read_solution_static_format
from wasp.load_ffm_model import load_ffm_model

from .testutils import RESULTS_DIR, get_segments_data, get_tensor_info


def test_load_ffm_model():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        segments = get_segments_data()
        tensor = get_tensor_info()
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "Solucion.txt",
            tempdir / "Solucion.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "fault&rise_time.txt",
            tempdir / "fault&rise_time.txt",
        )
        point_sources = point_sources_param(
            segments["segments"],
            tensor,
            segments["rise_time"],
        )

        # test solucion option
        solucion_model = load_ffm_model(
            segments_data=segments, point_sources=point_sources, directory=tempdir
        )
        target_static = read_solution_static_format(
            segments["segments"],
            tempdir,
        )
        for key in target_static:
            if key in solucion_model:
                np.testing.assert_array_equal(solucion_model[key], target_static[key])

        # TODO: add test data for the following. Run through in the meantime
        # test fault&rise_time
        frt_model = load_ffm_model(
            segments_data=segments,
            point_sources=point_sources,
            option="fault&rise_time.txt",
            directory=tempdir,
        )
        # test Checkerboard
        checker_model = load_ffm_model(
            segments_data=segments,
            point_sources=point_sources,
            option="Checkerboard",
            directory=tempdir,
        )
        # test Patches
        patches_model = load_ffm_model(
            segments_data=segments,
            point_sources=point_sources,
            option="Patches",
            directory=tempdir,
        )
        # test point_source
        patches_model = load_ffm_model(
            segments_data=segments,
            point_sources=point_sources,
            option="point_source",
            directory=tempdir,
        )

    finally:
        shutil.rmtree(tempdir)

import json
import pathlib
import shutil
import tempfile

from wasp.modelling_parameters import modelling_prop

from .testutils import RESULTS_DIR, get_segments_data, get_tensor_info


def test_modelling_prop():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        dict, model_space = modelling_prop(
            get_tensor_info(),
            get_segments_data(),
            ["cgnss", "insar", "surf_waves", "strong_waves", "tele_waves"],
            directory=tempdir,
        )
        # compare annealing_prop
        with open(tempdir / "annealing_prop.json") as ad:
            annealing_data = json.load(ad)
        with open(RESULTS_DIR / "NP1" / "annealing_prop.json") as ad:
            target_annealing = json.load(ad)
        assert dict == annealing_data == target_annealing
        # compare model_space
        with open(tempdir / "model_space.json") as md:
            model_data = json.load(md)
        with open(RESULTS_DIR / "NP1" / "model_space.json") as md:
            target_model = json.load(md)
        assert model_space == model_data == target_model
    finally:
        shutil.rmtree(tempdir)

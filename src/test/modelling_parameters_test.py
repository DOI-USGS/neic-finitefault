import json
import pathlib
import shutil
import tempfile

from ffm.modelling_parameters import create_dataset_weights, modelling_prop

from .testutils import (
    RESULTS_DIR,
    assert_values_close,
    get_segments_data,
    get_tensor_info,
)


def test_modelling_prop():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        dict, model_space = modelling_prop(
            get_tensor_info(),
            get_segments_data(),
            ["cgnss", "imagery", "surf_waves", "strong_waves", "tele_waves"],
            directory=tempdir,
        )
        # compare annealing_prop
        with open(tempdir / "annealing_prop.json") as ad:
            annealing_data = json.load(ad)
        with open(RESULTS_DIR / "NP1" / "annealing_prop.json") as ad:
            target_annealing = json.load(ad)
        # seismic_moment derives from LAPACK eigenvalues, which differ in the
        # last ulps across platforms/BLAS backends, so compare with tolerance
        assert_values_close(dict, annealing_data)
        assert_values_close(annealing_data, target_annealing)
        # compare model_space
        with open(tempdir / "model_space.json") as md:
            model_data = json.load(md)
        with open(RESULTS_DIR / "NP1" / "model_space.json") as md:
            target_model = json.load(md)
        assert_values_close(model_space, model_data)
        assert_values_close(model_data, target_model)
        # check dataset_weights.json creation
        assert (tempdir / "dataset_weights.json").is_file()
        with open(tempdir / "dataset_weights.json") as dw:
            dw_data = json.load(dw)
        assert "weights" in dw_data
        assert dw_data["weights"]["body"] == 1.0
        assert dw_data["weights"]["cgnss"] == 1.0
        assert dw_data["weights"]["imagery"] == 1.0
        assert dw_data["weights"]["strong"] == 1.0
        assert dw_data["weights"]["surf"] == 1.0
    finally:
        shutil.rmtree(tempdir)


def test_create_dataset_weights():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        weights_data = create_dataset_weights(
            data_type=["tele_waves", "gnss_data"],
            directory=tempdir,
        )
        assert weights_data["weights"]["body"] == 1.0
        assert weights_data["weights"]["static"] == 1.0
        assert "strong" not in weights_data["weights"]
        with open(tempdir / "dataset_weights.json") as f:
            saved = json.load(f)
        assert saved == weights_data
    finally:
        shutil.rmtree(tempdir)


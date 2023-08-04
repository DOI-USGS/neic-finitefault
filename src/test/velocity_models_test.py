import pathlib
import shutil
import tempfile

from wasp.velocity_models import (
    __process_line,
    model2dict,
    model_depth2thick,
    select_velmodel,
)

from .testutils import get_tensor_info, get_velmodel_data

TENSOR = get_tensor_info()


def test_model2dict():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        with open(tempdir / "test_model.txt", "w") as f:
            f.write("0 1 2 3 4 5\n")
            f.write("6 7 8 9 10 11\n")
            f.write("12 13 14 15 16 17\n")
        assert model2dict(tempdir / "test_model.txt") == {
            "p_vel": ["0", "6", "12"],
            "s_vel": ["1", "7", "13"],
            "dens": ["2", "8", "14"],
            "thick": ["3", "9", "15"],
            "qa": ["4", "10", "16"],
            "qb": ["5", "11", "17"],
        }
    finally:
        shutil.rmtree(tempdir)


def test_model_depth2thick():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        with open(tempdir / "test_model.txt", "w") as f:
            f.write("0 1 2 3 4 5\n")
            f.write("6 7 8 9 10 11\n")
            f.write("12 13 14 15 16 17\n")
        assert model_depth2thick(tempdir / "test_model.txt") == {
            "p_vel": ["0.0", "8.0"],
            "s_vel": ["1.75", "9.1"],
            "dens": ["2", "8"],
            "thick": ["6", "6"],
            "qa": ["4", "10"],
            "qb": ["5", "11"],
        }
    finally:
        shutil.rmtree(tempdir)


def test_process_line():
    assert __process_line("1 2 3 4 5 10.23") == [1, 2, 3, 4, 5, 10.23]


def test_select_velmodel():
    velmodel = get_velmodel_data()
    assert (
        select_velmodel(
            TENSOR,
            {
                "litho_model": pathlib.Path(__file__).parent.parent.parent
                / "fortran_code"
                / "info"
                / "LITHO1.0.nc"
            },
        )
        == velmodel
    )

import pathlib
import shutil
import tempfile

from wasp.inversion_chen_new import automatic_usgs, set_directory_structure
from wasp.management import default_dirs

from .testutils import DATA_DIR, END_TO_END_DIR, get_tensor_info, get_velmodel_data

TENSOR = get_tensor_info()
VELMODEL = get_velmodel_data()


def test_automatic_usgs():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        set_directory_structure(TENSOR, directory=tempdir)
        shutil.copy2(END_TO_END_DIR / "data", tempdir)
        ddirs = default_dirs(config_path=DATA_DIR / "config.ini")
        print(ddirs)
        automatic_usgs(
            tensor_info=TENSOR,
            data_type=[
                "cgps",
                "gps",
                "insar",
                "strong_motion",
                "surf_tele",
                "tele_body",
            ],
            default_dirs=updated_default_dirs,
        )
    finally:
        shutil.rmtree(tempdir)

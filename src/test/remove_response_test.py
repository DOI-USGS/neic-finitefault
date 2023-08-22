import pathlib
import shutil
import tempfile

from obspy.io.sac import SACTrace

from wasp.remove_response import convert_response_acc, get_sacpz_file

from .testutils import RESULTS_DIR


def test_get_sacpz_file():
    header = SACTrace.read(RESULTS_DIR / "data" / "P" / "final_IU_TSUM_BHZ.sac")
    assert get_sacpz_file(header) == "SAC_PZs_IU_TSUM_BHZ_00"
    assert get_sacpz_file(header, data="IRIS") == "SACPZ.IU.TSUM.00.BHZ"


def test_convert_response_acc():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        # test remove 2
        with open(tempdir / "testresponse", "w") as f:
            f.write(
                """TESTLINE
         ZEROS   2
         0.0000E+00   0.0000E+00
         0.0000E+00   0.0000E+00"""
            )
        convert_response_acc(tempdir / "testresponse")
        with open(tempdir / "testresponse") as f:
            d = f.read()
        assert d == "TESTLINE\nZEROS\t 0\n"
        # test remove 3
        with open(tempdir / "testresponse", "w") as f:
            f.write(
                """TESTLINE
         ZEROS   3
         0.0000E+00   0.0000E+00
         0.0000E+00   0.0000E+00
         0.0000E+00   0.0000E+00"""
            )
        convert_response_acc(tempdir / "testresponse")
        with open(tempdir / "testresponse") as f:
            d = f.read()
        assert d == "TESTLINE\nZEROS\t 0\n"
    finally:
        shutil.rmtree(tempdir)

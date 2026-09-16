import json
import pathlib
import shutil
import tempfile

from ffm.tsunami_plot import setplot

from .testutils import RESULTS_DIR

# TODO include checks for amplitude_time plot


def test_setplot():
    tempdir = tempfile.mkdtemp()
    tempdir = pathlib.Path(tempdir)
    try:
        geoclaw_parameter_file = RESULTS_DIR / "NP1" / "geoclaw_parameter_file.json"
        f = open(geoclaw_parameter_file)
        parameter = json.load(f)
        plotdata = setplot(parameter_file=geoclaw_parameter_file)
        assert plotdata.format == "ascii"
        assert plotdata.getaxes("pcolor", figname="Surface").xlimits == [
            parameter["model_region"]["west"],
            parameter["model_region"]["east"],
        ]
        assert plotdata.getaxes("pcolor", figname="Surface").ylimits == [
            parameter["model_region"]["south"],
            parameter["model_region"]["north"],
        ]
    finally:
        shutil.rmtree(tempdir)

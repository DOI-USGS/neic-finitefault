import pathlib
import shutil
import tempfile

from obspy import UTCDateTime  # type:ignore
import numpy as np

from wasp.seismic_tensor import get_tensor, planes_from_tensor, write_tensor

from .testutils import END_TO_END_DIR, RESULTS_DIR, get_tensor_info


def test_get_tensor():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        target = {
            "mrr": 1.95e28,
            "mtt": -4.36e26,
            "mpp": -1.91e28,
            "mrt": 7.42e27,
            "mrp": -2.48e28,
            "mtp": 9.42e26,
            "date_origin": UTCDateTime(2015, 9, 16, 22, 54, 32),
            "lat": -31.57,
            "lon": -71.67,
            "depth": 22.4,
            "time_shift": 49.98,
            "half_duration": 33.4,
            "centroid_lat": -31.13,
            "centroid_lon": -72.09,
            "centroid_depth": 17.35,
            "datetime": "2015-09-16T22:54:32",
            "moment_mag": 3.2291872286240194e28,
        }
        # test cmt file
        cmt_file = END_TO_END_DIR / "Illapel_CMTSOLUTION"
        cmt_tensor = get_tensor(cmt_file=cmt_file)
        del cmt_tensor["timedelta"]
        assert cmt_tensor == target

        # test updating tensor
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tensor_info.json", tempdir / "tensor_info.json"
        )
        tensor = get_tensor(directory=tempdir)
        del tensor["timedelta"]
        assert tensor == target

        # test writing the tensor
        write_tensor(tensor, tempdir)

        # test xml file
        # TODO, xml reader doesn't appear to work. Fix so that this passes
        xml_file = END_TO_END_DIR / "info" / "quakeml.xml"
        xml_tensor = get_tensor(quake_file=xml_file)
        del xml_tensor["timedelta"]
        # TODO investigate why the xml and cmt are different
        # (collected at different points in time?)
        assert xml_tensor == {
            "mrr": 1.915e28,
            "mtt": 3.5e26,
            "mpp": -1.949e28,
            "mrt": -3.6e26,
            "mrp": -2.536e28,
            "mtp": 1.19e27,
            "date_origin": UTCDateTime(2015, 9, 16, 22, 54, 32, 860000),
            "lat": -31.5729,
            "lon": -71.6744,
            "depth": 22.44,
            "time_shift": 5.0,
            "half_duration": 38.05789881045027,
            "centroid_lat": -31.637,
            "centroid_lon": -71.741,
            "centroid_depth": 23.3,
            "datetime": "2015-09-16T22:54:32.860000",
            "moment_mag": 3.1905117287303885e28,
        }
    finally:
        shutil.rmtree(tempdir)


def test_planes_from_tensor():
    tensor = get_tensor_info()
    p1, p2 = planes_from_tensor(tensor)
    p1_target = {
        "strike": 6.613912311529926,
        "dip": 19.280827965117993,
        "rake": 109.27817171619564,
    }
    p2_target = {
        "strike": 166.28167341901923,
        "dip": 71.83929826714792,
        "rake": 83.41183885645165,
    }
    for key, target_values in p1["plane_info"].items():
        np.testing.assert_almost_equal(p1["plane_info"][key], p1_target[key], decimal=5)
    for key, target_values in p2["plane_info"].items():
        np.testing.assert_almost_equal(p2["plane_info"][key], p2_target[key], decimal=5)

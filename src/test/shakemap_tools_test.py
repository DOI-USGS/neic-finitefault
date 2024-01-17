import pathlib
import shutil
import tempfile

import numpy as np

from wasp.fault_plane import point_sources_param
from wasp.get_outputs import read_solution_static_format
from wasp.shakemap_tools import (
    equivalent_slip_length,
    locate_equivalent_slip,
    translate_xy_to_latlondep,
)

from .testutils import RESULTS_DIR, get_segments_data, get_tensor_info

TENSOR = get_tensor_info()
SEGMENTS = get_segments_data()


def get_strike():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solucion.txt", tempdir / "Solucion.txt")
        solution = read_solution_static_format(
            SEGMENTS["segments"],
            tempdir,
        )
    finally:
        shutil.rmtree(tempdir)
    slip = solution["slip"]
    slip_seg = slip[0]
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    slip_threshold = 0.1 * (max_slip / 100.0)
    if slip_threshold < 1.0:
        slip_threshold = 1.0
    if max_slip / 100.0 < 3.0:
        slip_threshold = 0.5
    if max_slip / 100.0 < 1.0:
        slip_threshold = 0.2
    cumslip_seg = slip_seg
    idxslip = np.where(cumslip_seg < 100.0 * slip_threshold)
    cumslip_seg[idxslip] = 0.0
    along_strike = np.sum(cumslip_seg, axis=0)
    along_dip = np.sum(cumslip_seg, axis=1)
    delta_strike = SEGMENTS["segments"][0]["delta_strike"]
    delta_dip = SEGMENTS["segments"][0]["delta_dip"]
    return along_strike, along_dip, delta_strike, delta_dip


def test_shakemap_tools():
    along_strike, along_dip, delta_strike, delta_dip = get_strike()
    eq_len_AS = equivalent_slip_length(along_strike / 100, delta_strike)
    eq_len_AD = equivalent_slip_length(along_dip / 100, delta_dip)
    assert eq_len_AS == 275.4
    assert eq_len_AD == 89.1
    left_edge_AS = locate_equivalent_slip(along_strike / 100, delta_strike, eq_len_AS)
    left_edge_AD = locate_equivalent_slip(along_dip / 100, delta_dip, eq_len_AD)
    assert left_edge_AS == 113.28652173913042
    assert left_edge_AD == -7.239116979408531
    corner_1, corner_2, corner_3, corner_4 = translate_xy_to_latlondep(
        SEGMENTS["segments"][0],
        TENSOR["lon"],
        TENSOR["lat"],
        TENSOR["depth"],
        eq_len_AS,
        eq_len_AD,
        left_edge_AS,
        left_edge_AD,
    )
    assert corner_1 == "-71.61 -30.55 20.01"
    assert corner_2 == "-71.29 -28.08 20.01"
    assert corner_3 == "-70.41 -28.16 49.43"
    assert corner_4 == "-70.72 -30.63 49.43"

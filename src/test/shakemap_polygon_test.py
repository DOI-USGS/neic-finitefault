#!/usr/bin/env python

import pathlib
import tempfile

from wasp.shakemap_polygon import ShakeRupture

SINGLE_SEGMENT_POLYGON = """# Source: USGS NEIC Rapid Finite Fault
# Eventid: us70004bn0
# Model created: 2025-04-08 07:39:31
-117.33 35.49 25.6
-117.31 35.50 1.1
-118.07 36.21 1.1
-118.09 36.20 25.6
-117.33 35.49 25.6
>
"""

MULTI_SEGMENT_POLYGON = """# Source: USGS NEIC Rapid Finite Fault
# Eventid: us7000pn9s
# Model created: 2025-04-08 08:13:48
96.02 21.11 1.1
96.05 21.11 18.9
95.99 22.64 18.9
95.96 22.64 1.1
96.02 21.11 1.1
>
96.07 20.49 1.1
96.10 20.49 18.9
96.05 20.94 18.9
96.03 20.93 1.1
96.07 20.49 1.1
>
96.26 19.24 1.1
96.28 19.24 18.9
96.11 20.40 18.9
96.08 20.40 1.1
96.26 19.24 1.1
>
96.47 18.09 1.1
96.49 18.09 18.9
96.29 19.15 18.9
96.27 19.15 1.1
96.47 18.09 1.1
>"""


def _test_rupture(fspfile, eventid, cmpstring):
    with tempfile.TemporaryDirectory() as tdir:
        tmpdir = pathlib.Path(tdir)
        rupture = ShakeRupture(eventid, fspfile)
        polygon_file = rupture.write_rupture(tmpdir)
        with open(polygon_file, "rt") as fobj:
            polylines = [
                line.rstrip() for line in fobj.readlines() if not line.startswith("#")
            ]
            seglines = [
                line for line in cmpstring.split("\n") if not line.startswith("#")
            ]
            for i, polyline in enumerate(polylines):
                segline = seglines[i]
                assert polyline == segline


def test_shake_rupture():
    single_fspfile = (
        pathlib.Path(__file__).parent / "data" / "single_segment_inversion.fsp"
    )
    single_eventid = "us70004bn0"
    multi_fspfile = (
        pathlib.Path(__file__).parent / "data" / "multi_segment_inversion.fsp"
    )
    multi_eventid = "us7000pn9s"
    _test_rupture(single_fspfile, single_eventid, SINGLE_SEGMENT_POLYGON)
    _test_rupture(multi_fspfile, multi_eventid, MULTI_SEGMENT_POLYGON)


if __name__ == "__main__":
    test_shake_rupture()

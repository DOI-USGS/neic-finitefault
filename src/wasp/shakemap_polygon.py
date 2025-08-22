#!/usr/bin/env python

import datetime
import pathlib
import sys
import typing

import numpy as np
import numpy.typing as npt

SM_PROFILE = pathlib.Path.home() / ".shakemap" / "profiles.conf"
# 2021-10-25 17:04:20
TIMESTR = "%Y-%m-%d %H:%M:%S"

FLOATPAT = r"[-+]?[0-9]*\.?[0-9]*"
COLUMNS = [
    "latitude",
    "longitude",
    "ew_coord",
    "ns_coord",
    "depth",
    "slip",
    "rake",
    "rupture_time",
    "rise_time",
    "moment",
]


def counterclockwise_sort(points: list) -> list:
    """
    Sorts a list of 3D points in counterclockwise order around their centroid.

    Args:
        points: A list of tuples, where each tuple represents a point (x, y, z).

    Returns:
        list:
            A new list of points sorted in counterclockwise order.
    """
    if not points:
        return []

    # Calculate the centroid of the points
    center_x = sum(x for x, y, z in points) / len(points)
    center_y = sum(y for x, y, z in points) / len(points)

    # Calculate angles relative to the centroid
    def angle_to_center(point):
        x, y, z = point
        return np.arctan2(y - center_y, x - center_x)

    # Sort points by angle
    return sorted(points, key=angle_to_center)


class ShakeRupture:
    def __init__(self, eventid: str, fsp_file: pathlib.Path) -> None:
        """
        Trim the edges of possibly many segments of a Finite Fault.

        Args:
            eventid (str):
                ComCat event id.
            fsp_file (pathlib.Path):
                Path to existing FSP file.
        Returns:
           None
        """
        self.eventid = eventid
        with open(fsp_file, "rt") as fobj:
            self.segments = self.read_fsp(fobj)

    def read_fsp(self, fobj: typing.TextIO) -> list:
        """
        Read an FSP format file and return segment data/metadata.

        Args:
            fobj (File-like object):
                Open fsp file.
        Returns:
           list:
                List of segment dictionaries, including metadata and data.
        """
        # read the header information for the segments
        segments = self.read_fsp_headers(fobj)
        # now populate the segments with arrays of data
        for segment in segments:
            fobj.seek(0)
            nrows = segment["nx"] * segment["nz"]
            data = np.loadtxt(fobj, skiprows=segment["nlines"], max_rows=nrows)
            data_array = self.organize_data(data, segment["nx"], segment["nz"])
            segment["data"] = data_array
        return segments

    def read_fsp_headers(self, fobj: typing.TextIO) -> list[dict[str, typing.Any]]:
        """
        Read FSP format headers and return segment metadata.

        Args:
            fobj (File-like object):
                Open fsp file.
        Returns:
           list:
                List of segment dictionaries with metadata.
        """
        segments = []
        segment: dict[str, typing.Any] = {}
        dx = None
        dz = None
        nlines = 0
        for line in fobj.readlines():
            nlines += 1
            if not line.startswith("%"):  # end of segment
                if len(segment):
                    segment["dx"] = dx
                    segment["dz"] = dz
                    segment["nlines"] = nlines - 1
                    segment["nx"] = int(segment["length"] / segment["dx"])
                    segment["nz"] = int(segment["width"] / segment["dz"])
                    segments.append(segment.copy())
                    segment = {}
                    continue
                else:
                    continue

            commline = line.strip("%").strip()
            if commline.startswith("Invs : Dx"):
                parts = commline.split(":")[1].split()
                dx = float(parts[2])
                dz = float(parts[6])
                continue
            if commline.startswith("Size"):
                parts = commline.split()
                segment["segment"] = 1
                segment["length"] = float(parts[4])
                segment["width"] = float(parts[8])
                continue
            if commline.startswith("SEGMENT"):
                parts = commline.split(":")
                segment["segment"] = int(parts[0].split()[-1])
                angleparts = parts[1].strip().split()
                segment["strike"] = float(angleparts[2])
                segment["dip"] = float(angleparts[6])
                continue
            if commline.startswith("LEN"):
                parts = commline.split()
                segment["length"] = float(parts[2])
                segment["width"] = float(parts[6])
                continue
            if commline.startswith("hypocenter"):
                parts = commline.split(":")[1].split(",")
                xstring = parts[0].split("=")[1].strip()
                ystring = parts[1].split("=")[1].strip()
                segment["hypo_x"] = float(xstring)
                segment["hypo_y"] = float(ystring)
                continue
        return segments

    def organize_data(
        self, data: npt.NDArray, nx: int, ny: int
    ) -> npt.NDArray[np.float32]:
        """
        Arrange 2D data array into labeled numpy record arrays.

        Args:
            data (np.array):
                2D array of files read from FSP files.
            nx (int):
                Number of cells in the strike direction.
            nx (int):
                Number of cells in the dip direction.
        Returns:
           np.recarray:
                numpy labeled record array of values (lat/lon/slip/etc).
        """
        dtypes = list(zip(COLUMNS, [np.float32] * len(COLUMNS)))
        data_array = np.empty((ny, nx), dtypes)
        for icol, column in enumerate(COLUMNS):
            data_array[column] = data[:, icol].reshape((ny, nx))
        return data_array

    def write_rupture(self, event_path: pathlib.Path) -> pathlib.Path:
        """
        Write ShakeMap compatible polygon text file with trimmed edges.

        Args:
            event_path (pathlib.Path):
                Path to directory where finite fault event products should be written.
        Returns:
           pathlib.Path:
                Path to polygon text file.
        """
        polygon_file = pathlib.Path(event_path) / "shakemap_polygon.txt"
        tnow = datetime.datetime.now()
        ordered_points = self.get_segment_corners()
        with open(polygon_file, "wt") as fobj:
            fobj.write("# Source: USGS NEIC Rapid Finite Fault\n")
            fobj.write(f"# Eventid: {self.eventid}\n")
            fobj.write(f"# Model created: {tnow.strftime(TIMESTR)}\n")
            for segment, points in ordered_points.items():
                for ptuple in points:
                    fobj.write(f"{ptuple[0]:.2f} {ptuple[1]:.2f} {ptuple[2]:.1f}\n")
                fobj.write(
                    f"{points[0][0]:.2f} {points[0][1]:.2f} {points[0][2]:.1f}\n"
                )
                fobj.write(">\n")

        return polygon_file

    def get_segment_corners(self) -> dict:
        """
        Write ShakeMap compatible polygon text file with trimmed edges.

        Returns:
           list:
                lon/lat/depth coordinates denoting trimmed corners of fault segments.
        """
        segdict = {}
        for segment in self.segments:
            slip = segment["data"]["slip"]
            latitude = segment["data"]["latitude"]
            longitude = segment["data"]["longitude"]
            depth = segment["data"]["depth"]
            total_slip = slip.sum()

            slip_row_sum = slip.sum(axis=0)
            trim_sum = slip_row_sum.sum()
            colidx = 0
            while trim_sum / total_slip > 0.9:
                colidx += 1
                leftedge = colidx
                rightedge = -1 - colidx
                trim_sum = slip_row_sum[leftedge:rightedge].sum()
            leftedge -= 1
            rightedge += 1

            slip_col_sum = slip.sum(axis=1)
            trim_sum = slip_col_sum.sum()
            rowidx = 0
            while trim_sum / total_slip > 0.9:
                rowidx += 1
                bottomedge = rowidx
                topedge = -1 - rowidx
                trim_sum = slip_col_sum[bottomedge:topedge].sum()
            bottomedge -= 1
            topedge += 1

            lat_corner1 = latitude[bottomedge, leftedge]
            lat_corner2 = latitude[bottomedge, rightedge]
            lat_corner3 = latitude[topedge, rightedge]
            lat_corner4 = latitude[topedge, leftedge]

            lon_corner1 = longitude[bottomedge, leftedge]
            lon_corner2 = longitude[bottomedge, rightedge]
            lon_corner3 = longitude[topedge, rightedge]
            lon_corner4 = longitude[topedge, leftedge]

            dep_corner1 = depth[bottomedge, leftedge]
            dep_corner2 = depth[bottomedge, rightedge]
            dep_corner3 = depth[topedge, rightedge]
            dep_corner4 = depth[topedge, leftedge]

            ordered_points = counterclockwise_sort(
                [
                    (lon_corner1, lat_corner1, dep_corner1),
                    (lon_corner2, lat_corner2, dep_corner2),
                    (lon_corner3, lat_corner3, dep_corner3),
                    (lon_corner4, lat_corner4, dep_corner4),
                ]
            )
            segdict[segment["segment"]] = ordered_points
        return segdict

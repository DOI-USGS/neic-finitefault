import io
import pathlib
from functools import cached_property
from typing import List, Tuple

import pandas as pd  # type: ignore [import-untyped]
from pydantic import BaseModel, Field

from ffm.model.result.segment import Segment, SegmentData
from ffm.model.setup.connection import Connection
from ffm.model.setup.risetime import RiseTime
from ffm.model.setup.segmentmetadata import SegmentMetadata


class Solution(BaseModel):
    """
    Class to handle model solution files, formats, and data
    """

    connections: List[Connection] = Field(
        [], description="Coordinates linking/connecting two segments together"
    )
    rise_time: RiseTime
    segments: List[Segment]

    @classmethod
    def from_static_solution_file(
        cls, segments_data: dict, solution_filepath: pathlib.Path
    ):
        """
        Instantiate Solution from segments_data and the static solution file

        :param cls: Solution
        :param segments_data: The segments information (from segments_data.json)
        :param solution_filepath: The path to the static solution file (Solution.txt)
        """
        if not solution_filepath.exists():
            raise FileNotFoundError(
                f"Provided file {solution_filepath} does not exist!"
            )
        segments_metadata = segments_data["segments"]
        segments = Solution._parse_static_solution_file(
            segments_metadata, solution_filepath
        )

        solution = cls(
            rise_time=segments_data["rise_time"],
            segments=segments,
            connections=segments_data.get("connections", []),
        )
        return solution

    @cached_property
    def num_segments(self) -> int:
        """
        Count the number of segments

        :param self: Solution instance
        :return: The number of segments
        """
        return len(self.segments)

    @staticmethod
    def __parse_static_data(
        lines: List[str], column_names: List[str]
    ) -> Tuple[pd.DataFrame, List[str]]:
        """
        Parse the fixed width data for a segment in static solution file (Solution.txt)

        :param lines: A subset of lines from the file
        :param column_names: List of the names of the data columns
        :return: (data, an updated subset of lines from the file with used lines removed)
        """
        data_io = io.StringIO()
        while True:
            data_io.write(lines.pop(0))
            if len(lines) == 0 or "#" in lines[0] or lines[0] in ["", "\n"]:
                break
        data_io.seek(0)
        data = pd.read_fwf(
            data_io,
            names=column_names,
        )
        return data, lines

    @staticmethod
    def __parse_static_header(
        lines: List[str],
    ) -> Tuple[float, float, int, int, List[tuple]]:
        """
        Parse the header for a segment in static solution file (Solution.txt)

        :param lines: A subset of lines from the file
        :return: (
                    delta in the x direction,
                    delta in the y direction,
                    number of subfaults in the x direction,
                    number of subfaults in the y direction,
                    The polygon defining the extents of the segment
                )
        """
        segment_metadata = lines[0]
        dx = float(
            segment_metadata.split("Dx =")[-1].split("ny")[0].replace("km", "").strip()
        )
        dy = float(segment_metadata.split("Dy =")[-1].split("km")[0].strip())
        nx = int(segment_metadata.split("nx(Along-strike)=")[-1].split("Dx")[0].strip())
        ny = int(segment_metadata.split("ny(downdip)=")[-1].split("Dy")[0].strip())
        # removed unused boundar of fault segment and polygon header self._lines
        # and add polygon
        polygon = []
        for i in range(8):
            header_line = lines[i]
            if "#" not in header_line and header_line:
                polygon += [
                    tuple(
                        [col for col in header_line.split(" ") if col not in ["", "\n"]]
                    )
                ]
        return dx, dy, nx, ny, polygon

    @staticmethod
    def _parse_static_solution_file(
        segments_metadata: List[dict], solution_filepath: pathlib.Path
    ) -> List[Segment]:
        """
        Parse the segments from a static solution file (Solution.txt)

        :param segments_metadata: The segments information (from segments_data.json)
        :param solution_filepath: The path to the static solution file (Solution.txt)
        :return: List of Segments
        """
        # read file
        with open(solution_filepath) as f:
            lines = f.readlines()
        num_subfaults = int(lines.pop(0).split("fault_segments=")[-1].strip())
        segments = []
        for segment_number in range(num_subfaults):
            cols = [
                "lat",
                "lon",
                "depth",
                "slip",
                "rake",
                "strike",
                "dip",
                "trup",
                "trise",
                "tfall",
                "moment",
            ]
            metadata = SegmentMetadata(**segments_metadata[segment_number])
            # first line has metadata
            dx, dy, nx, ny, polygon = Solution.__parse_static_header(
                [lines.pop(0) for _ in range(9)]
            )
            # now parse the data
            data, lines = Solution.__parse_static_data(lines, cols)

            # move data to corner
            flat_data = {}
            for col in cols:
                flat_data[col] = data[col].values.flatten()

            # setup segment
            segments += [
                Segment(
                    data=SegmentData(**flat_data),
                    dx=dx,
                    dy=dy,
                    metadata=metadata,
                    nx=nx,
                    ny=ny,
                    polygon=polygon,
                )
            ]
        return segments

    @cached_property
    def _old_format(self) -> dict:
        """
        Product the old format used in many other functions. Marked as internal as this will
        go away as other parts of the code are optimized to new the new Solution class

        :param self: Solution instance
        :return: Dictionary of segment data
        """
        return {
            "slip": [
                s.data.slip.reshape(s.metadata.dip_subfaults, s.metadata.stk_subfaults)
                for s in self.segments
            ],
            "rake": [
                s.data.rake.reshape(s.metadata.dip_subfaults, s.metadata.stk_subfaults)
                for s in self.segments
            ],
            "rupture_time": [
                s.data.trup.reshape(s.metadata.dip_subfaults, s.metadata.stk_subfaults)
                for s in self.segments
            ],
            "trise": [
                s.data.trise.reshape(s.metadata.dip_subfaults, s.metadata.stk_subfaults)
                for s in self.segments
            ],
            "tfall": [
                s.data.tfall.reshape(s.metadata.dip_subfaults, s.metadata.stk_subfaults)
                for s in self.segments
            ],
            "lat": [
                s.data.lat.reshape(s.metadata.dip_subfaults, s.metadata.stk_subfaults)
                for s in self.segments
            ],
            "lon": [
                s.data.lon.reshape(s.metadata.dip_subfaults, s.metadata.stk_subfaults)
                for s in self.segments
            ],
            "depth": [
                s.data.depth.reshape(s.metadata.dip_subfaults, s.metadata.stk_subfaults)
                for s in self.segments
            ],
            "moment": [
                s.data.moment.reshape(
                    s.metadata.dip_subfaults, s.metadata.stk_subfaults
                )
                for s in self.segments
            ],
        }

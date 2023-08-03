# -*- coding: utf-8 -*-
"""Management of fault planes files and solution model
"""


import os
from typing import Tuple


def __unpack_plane_data(
    plane_info: dict,
) -> Tuple[int, int, float, float, float, float]:
    """Extract information from a fault segment

    :param plane_info: The segment/plane information
    :type plane_info: dict
    :return: The number of strike subfaults, the number of dip subfaults,
            the increment of the strike, the increment of the dip, the strike of
            the hypocenter, the dip of the hypocenter
    :rtype: Tuple[int, int, float, float, float, float]
    """
    stk_subfaults = plane_info["stk_subfaults"]
    dip_subfaults = plane_info["dip_subfaults"]
    delta_strike = plane_info["delta_strike"]
    delta_dip = plane_info["delta_dip"]
    # TODO: Investigate why -1 is applied
    hyp_stk = plane_info["hyp_stk"] - 1
    hyp_dip = plane_info["hyp_dip"] - 1
    return stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip


if __name__ == "__main__":
    import os

    from wasp.get_outputs import read_solution_fsp_format

    #    segments, rise_time, point_sources = __read_planes_info()
    #    print(segments)
    os.chdir(
        "/home/pk/Inversion_Chen_Ji/Inversiones/20140401234647/paper_challa/pk.5/NP2/srcmod_solutions/"
    )
    tensor_info, solution, subfault = read_solution_fsp_format("hayes")
    print(solution)

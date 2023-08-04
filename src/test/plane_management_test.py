from wasp.plane_management import __unpack_plane_data


def test__unpack_plane_data():
    assert __unpack_plane_data(
        {
            "stk_subfaults": 1,
            "dip_subfaults": 2,
            "delta_strike": 0.2,
            "delta_dip": 0.1,
            "hyp_stk": 21,
            "hyp_dip": 11,
        }
    ) == (
        1,
        2,
        0.2,
        0.1,
        20,
        10,
    )

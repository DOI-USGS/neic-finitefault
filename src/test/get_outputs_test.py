import json
import os
import pathlib
import shutil
from tempfile import mkdtemp

import numpy as np
from obspy import read  # type: ignore

from wasp.get_outputs import (
    __get_channel,
    __is_number,
    get_data_dict,
    get_imagery,
    read_solution_fsp_format,
    read_solution_static_format,
    retrieve_gnss,
    synthetics_to_SAC,
)

from .testutils import (
    END_TO_END_DIR,
    RESULTS_DIR,
    get_cgnss_json,
    get_imagery_json,
    get_segments_data,
    get_static_json,
    get_tensor_info,
    update_manager_file_locations,
)

SEGMENTS = get_segments_data()
TENSOR = get_tensor_info()


def test_get_data_dict():
    tempdir = pathlib.Path(mkdtemp())
    try:
        cgnss_waves = get_cgnss_json()
        new_cgnss_waves = update_manager_file_locations(
            cgnss_waves, tempdir / "data", replace_dir=str(RESULTS_DIR / "data")
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "synthetics_cgnss.txt",
            tempdir / "synthetics_cgnss.txt",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "waveforms_cgnss.txt", tempdir / "waveforms_cgnss.txt"
        )
        os.mkdir(tempdir / "data")
        os.mkdir(tempdir / "data" / "cGNSS")
        for o, n in zip(cgnss_waves, new_cgnss_waves):
            shutil.copyfile(o["file"], n["file"])
        updated = get_data_dict(
            new_cgnss_waves,
            tempdir / "synthetics_cgnss.txt",
            True,
            tempdir / "waveforms_cgnss.txt",
        )
        for wform, target_min, target_max, target_mean in zip(
            updated,
            [-25.5165386, -15.0545864, -156.320694],
            [6.46761608, 7.3232789, 2.9956665],
            [-9.419298518062012, -4.989112241010887, -110.9082276890334],
        ):
            assert np.min(wform["synthetic"]) == target_min
            assert np.max(wform["synthetic"]) == target_max
            assert np.mean(wform["synthetic"]) == target_mean
        updated = get_data_dict(
            new_cgnss_waves, tempdir / "synthetics_cgnss.txt", False
        )
        for wform, target_min, target_max, target_mean in zip(
            updated,
            [-25.5165386, -15.0545864, -156.320694],
            [6.46761608, 7.3232789, 2.9956665],
            [-9.419298518062012, -4.989112241010887, -110.9082276890334],
        ):
            assert np.min(wform["synthetic"]) == target_min
            assert np.max(wform["synthetic"]) == target_max
            assert np.mean(wform["synthetic"]) == target_mean
    finally:
        shutil.rmtree(tempdir)


def test_get_channel():
    assert __get_channel("P") == ["BHZ"]
    assert __get_channel("BHE") == ["BHT"]
    assert __get_channel("HLZ") == ["HNZ", "HLZ", "BNZ"]
    assert __get_channel("HNE") == ["HNE", "HLE", "BNE"]
    assert __get_channel("HNN") == ["HNN", "HLN", "BNN"]
    assert __get_channel("LHZ") == ["LXZ", "LHZ", "LYZ"]
    assert __get_channel("LYE") == ["LXE", "LHE", "LYE"]
    assert __get_channel("LYN") == ["LXN", "LHN", "LYN"]
    assert __get_channel("dart") == ["dart"]


def test_get_imagery():
    ## just run without failing
    tempdir = pathlib.Path(mkdtemp())
    try:
        imagery_data = get_imagery_json()
        new_imagery = update_manager_file_locations(
            imagery_data,
            tempdir,
            replace_dir=str(RESULTS_DIR / "NP1"),
            file_key="name",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "imagery_synthetics.txt",
            tempdir / "imagery_synthetics.txt",
        )
        with open(tempdir / "imagery_data.json", "w") as i:
            json.dump(new_imagery, i)
        for key in ["insar_ascending", "insar_descending"]:
            for o, n in zip(imagery_data[key], new_imagery[key]):
                shutil.copyfile(o["name"], n["name"])
        get_imagery(tempdir)
    finally:
        shutil.rmtree(tempdir)


def test_is_number():
    assert __is_number("invalid") == False
    assert __is_number("20.3") == True


def test_read_solution_static_format():
    tempdir = mkdtemp()
    try:
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "Solution.txt",
            pathlib.Path(tempdir) / "Solution.txt",
        )
        solution = read_solution_static_format(
            SEGMENTS["segments"],
            tempdir,
        )
        for key, target in zip(
            [
                "slip",
                "rake",
                "rupture_time",
                "trise",
                "tfall",
                "lat",
                "lon",
                "depth",
                "moment",
            ],
            [
                573.014587,
                129.273148,
                118.995361,
                15.0,
                15.0,
                -29.326191,
                -70.841278,
                40.1408,
                6.50101e26,
            ],
        ):
            assert np.max([np.max(seg.flatten()) for seg in solution[key]]) == target
    finally:
        shutil.rmtree(tempdir)


def test_read_solution_fsp_format():
    tempdir = mkdtemp()
    try:
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "PublicationFiles" / "us20003k7a.fsp",
            pathlib.Path(tempdir) / "us20003k7a.fsp",
        )
        tensor, solution, subfault = read_solution_fsp_format(
            pathlib.Path(tempdir) / "us20003k7a.fsp",
        )
        assert tensor == {
            "lat": TENSOR["lat"],
            "lon": TENSOR["lon"],
            "depth": TENSOR["depth"],
        }
        seg_target = SEGMENTS["segments"][0]
        assert subfault == {
            "stk_subfaults": seg_target["stk_subfaults"],
            "dip_subfaults": seg_target["dip_subfaults"],
            "delta_strike": seg_target["delta_strike"],
            "delta_dip": seg_target["delta_dip"],
            "strike": seg_target["strike"],
        }
        for key, target in zip(
            [
                "slip",
                "rake",
                "trise",
                "lat",
                "lon",
                "depth",
            ],
            [
                6.5294,
                129.2683,
                30.0,
                -29.2679,
                -70.7734,
                42.112,
            ],
        ):
            assert np.max(solution[key]) == target
    finally:
        shutil.rmtree(tempdir)


def test_synthetics_to_SAC():
    tempdir = pathlib.Path(mkdtemp())
    try:
        waveform_type = ["body", "surf", "strong", "cgnss"]
        data_dir = ["BODY", "LONG", "STR", "cGNSS"]
        waveform_json = [
            "tele_waves.json",
            "surf_waves.json",
            "strong_motion_waves.json",
            "cgnss_waves.json",
        ]
        synthetic_txt = [
            "synthetics_body.txt",
            "synthetics_surf.txt",
            "synthetics_strong.txt",
            "synthetics_cgnss.txt",
        ]
        target_lon = [
            5.110109806060791,
            5.110109806060791,
            -30.172700881958008,
            -30.6037540435791,
        ]
        target_lat = [
            -52.64448165893555,
            -52.64448165893555,
            -70.79930114746094,
            -71.20388793945312,
        ]
        target_id = ["MPG", "MPG", "GO04", "ovll"]
        target_dt = [0.2, 4.0, 0.4, 0.4]
        target_min = [-329.32175, -1.349203, -2.30605, -25.51653]
        target_max = [147.96167, 1.2100694, 2.70080, 6.467616]
        for data_type in range(len(waveform_json)):
            with open(
                END_TO_END_DIR / "results" / "NP1" / waveform_json[data_type], "r"
            ) as t:
                waveform_dict = json.load(t)
            waveform_dict[0][
                "file"
            ] = f'{END_TO_END_DIR}/results/data/{waveform_dict[0]["file"]}'
            # copy only first station to test
            station_name = waveform_dict[0]["name"]
            station_comp = waveform_dict[0]["component"]
            with open(tempdir / waveform_json[data_type], "w") as f:
                json.dump(
                    [waveform_dict[0]],
                    f,
                    sort_keys=True,
                    indent=4,
                    separators=(",", ": "),
                    ensure_ascii=False,
                )
            out_file = tempdir / synthetic_txt[data_type]
            with open(out_file, "w") as outfile:
                with open(
                    END_TO_END_DIR / "results" / "NP1" / synthetic_txt[data_type], "r"
                ) as file:
                    for line in file.readlines()[0:1025]:  # copy only first station
                        outfile.write(line)
            outfile.close()
            # check targets
            synthetics_to_SAC(waveform_type[data_type], directory=tempdir)
            st = read(
                f"{tempdir}/forward_model/{data_dir[data_type]}/forward_{station_name}_{station_comp}.sac"
            )
            assert st[0].stats.sac.stla == target_lat[data_type]
            assert st[0].stats.sac.stlo == target_lon[data_type]
            assert st[0].stats.sac.kstnm == target_id[data_type]
            np.testing.assert_almost_equal(
                st[0].stats.sac.delta, target_dt[data_type], 7
            )
            np.testing.assert_almost_equal(np.min(st[0].data), target_min[data_type], 5)
            np.testing.assert_almost_equal(np.max(st[0].data), target_max[data_type], 5)
    finally:
        shutil.rmtree(tempdir)


def test_retrieve_gnss():
    tempdir = pathlib.Path(mkdtemp())
    try:
        static_data = get_static_json()
        shutil.copyfile(
            END_TO_END_DIR / "data" / "static_data.json", tempdir / "static_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "static_synthetics.txt",
            tempdir / "static_synthetics.txt",
        )

        names, lats, lons, observed, synthetic, error = retrieve_gnss(directory=tempdir)
        assert names == [
            "VALN",
            "ZAPA",
            "LSCH",
            "TOLO",
            "PEDR",
            "LVIL",
            "CERN",
            "CMBA",
            "SLMC",
            "PFRJ",
        ]
        assert lats == [
            -33.0279,
            -32.5528,
            -29.9082,
            -30.1699,
            -30.839,
            -31.9092,
            -32.5581,
            -31.1882,
            -31.777,
            -30.6748,
        ]
        assert lons == [
            -71.635,
            -71.4656,
            -71.246,
            -70.8061,
            -70.6891,
            -71.5138,
            -70.9289,
            -70.999,
            -70.9628,
            -71.6354,
        ]
        np.testing.assert_array_almost_equal(
            np.array(observed, dtype=float).flatten(),
            np.array(
                [
                    ["-1.26", "-0.42", "-0.966"],
                    ["-1.185", "1.1400000000000001", "-3.9800000000000004"],
                    ["-2.82", "-9.855", "-16.97"],
                    ["-0.9119999999999999", "-12.04", "-25.28"],
                    ["-3.64", "-10.02", "-53.31"],
                    ["-6.98", "6.614000000000001", "-35.3"],
                    ["-1.142", "3.8739999999999997", "-6.357"],
                    ["-12.41", "0.629", "-82.82000000000001"],
                    ["-5.12", "12.679000000000002", "-39.438"],
                    ["-24.69", "-24.07", "-143.53"],
                ],
                dtype=float,
            ).flatten(),
            decimal=5,
        )
        np.testing.assert_array_almost_equal(
            np.array(synthetic, dtype=float).flatten(),
            np.array(
                [
                    ["-0.391721070", "-2.05364060", "-3.88984251"],
                    ["-0.403379977", "-0.732790589", "-9.14662075"],
                    ["1.84858394", "-10.5170603", "-22.9566288"],
                    ["-1.66698122", "-10.3476543", "-23.8477077"],
                    ["-6.36605597", "-8.84063911", "-45.7251091"],
                    ["-3.55172944", "4.73673296", "-38.2507591"],
                    ["-3.08602810", "3.04987383", "-10.8689032"],
                    ["-13.3618250", "-2.44215941", "-80.3503571"],
                    ["-6.26744890", "8.48536301", "-39.5643044"],
                    ["-23.3768368", "-23.5915432", "-144.032288"],
                ],
                dtype=float,
            ).flatten(),
            decimal=5,
        )
        np.testing.assert_array_almost_equal(
            np.array(error, dtype=float).flatten(),
            np.array(
                [
                    [0.32, 0.13999999999999999, 0.116],
                    [0.42, 0.13, 0.13],
                    [0.31, 0.13999999999999999, 0.12],
                    [0.29, 0.13, 0.16],
                    [0.38999999999999996, 0.15, 0.18],
                    [0.4, 0.15, 0.15],
                    [0.52, 0.2, 0.16],
                    [0.35000000000000003, 0.128, 0.13],
                    [0.35000000000000003, 0.15, 0.15],
                    [0.45999999999999996, 0.18, 0.2],
                ],
                dtype=float,
            ).flatten(),
            decimal=5,
        )
    finally:
        shutil.rmtree(tempdir)

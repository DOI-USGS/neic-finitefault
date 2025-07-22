import json
import os
import pathlib
import shutil
import tempfile

import numpy as np

from wasp.input_files import (
    forward_model,
    from_synthetic_to_obs,
    input_chen_dart,
    input_chen_insar,
    input_chen_near_field,
    input_chen_static,
    input_chen_tele_body,
    input_chen_tele_surf,
    inputs_simmulated_annealing,
    model_space,
    plane_for_chen,
    write_green_file,
    write_velmodel,
)

from .testutils import (
    RESULTS_DIR,
    get_cgps_json,
    get_insar_json,
    get_sampling_filter,
    get_segments_data,
    get_static_json,
    get_strong_motion_json,
    get_surf_waves_json,
    get_tele_waves_json,
    get_tensor_info,
    get_velmodel_data,
    update_manager_file_locations,
)

CMT = get_tensor_info()
SAMPLE_FILTER = get_sampling_filter()
SEGMENTS = get_segments_data()


def test_forward_model():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        ## Don't have an example of this all the way through
        ## TODO: add baseline data
        segments = get_segments_data()
        velmodel = get_velmodel_data()
        rupt_vel = segments["segments"][0]["rupture_vel"]
        lambda_min = 0.4
        lambda_max = 1.25
        min_vel, max_vel = [lambda_min * rupt_vel, lambda_max * rupt_vel]
        forward_model(
            CMT,
            segments,
            {
                "slip": np.array([0] * 207).reshape(1, 9, 23),
                "rake": np.array([109.27817171619564] * 207).reshape(1, 9, 23),
                "trise": np.array([0] * 207).reshape(1, 9, 23),
                "tfall": np.array([1.5] * 207).reshape(1, 9, 23),
                "trup": np.array([1.5] * 207).reshape(1, 9, 23),
            },
            min_vel,
            max_vel,
            tempdir,
        )
        for f in ["fault&rise_time.txt"]:
            with open(tempdir / f, "r") as d:
                data = d.read()
            with open(RESULTS_DIR / "NP1" / f, "r") as t:
                target = t.read()
            data = "\n".join(data.split("\n")[:5])
            target = "\n".join(target.split("\n")[:5])
            assert data == target
    finally:
        shutil.rmtree(tempdir)


def test_input_chen_dart():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        ## dont have dart data so use surface waves as an example to just run
        ## TODO: Add target data
        new_surfwaves = update_manager_file_locations(
            get_surf_waves_json(all=True),
            tempdir,
            replace_dir=str(RESULTS_DIR / "data"),
        )
        with open(tempdir / "dart_waves.json", "w") as tw:
            json.dump(new_surfwaves, tw)
        os.mkdir(tempdir / "RAYLEIGH")
        os.mkdir(tempdir / "LOVE")
        for o, n in zip(get_surf_waves_json(all=True), new_surfwaves):
            shutil.copyfile(o["file"], n["file"])

        input_chen_dart(CMT, SAMPLE_FILTER, tempdir)
    finally:
        shutil.rmtree(tempdir)


def test_input_chen_insar():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_insar = update_manager_file_locations(
            get_insar_json(),
            tempdir,
            replace_dir=str(RESULTS_DIR / "NP1"),
            file_key="name",
        )
        with open(tempdir / "insar_data.json", "w") as tw:
            json.dump(new_insar, tw)
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_ascending.txt", tempdir / "insar_ascending.txt"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_descending.txt",
            tempdir / "insar_descending.txt",
        )

        input_chen_insar(tempdir)
        for f in ["insar_data.txt", "insar_weights.txt"]:
            with open(tempdir / f, "r") as nd:
                new_data = nd.read()
            with open(RESULTS_DIR / f, "r") as t:
                target = t.read()
            assert new_data == target
    finally:
        shutil.rmtree(tempdir)


def test_input_chen_near_field():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        # test for strong motion
        new_strong_motion = update_manager_file_locations(
            get_strong_motion_json(all=True),
            tempdir,
            replace_dir=str(RESULTS_DIR / "data"),
        )
        with open(tempdir / "strong_motion_waves.json", "w") as tw:
            json.dump(new_strong_motion, tw)
        os.mkdir(tempdir / "STR")
        for o, n in zip(get_strong_motion_json(all=True), new_strong_motion):
            shutil.copyfile(o["file"], n["file"])

        input_chen_near_field(CMT, SAMPLE_FILTER, "strong", tempdir)
        for f in [
            "filter_strong.txt",
            "channels_strong.txt",
            "wavelets_strong.txt",
            "waveforms_strong.txt",
        ]:
            with open(tempdir / f, "r") as nd:
                new_data = nd.read()
            with open(RESULTS_DIR / f, "r") as t:
                target = t.read()
            assert new_data == target

        # test for cgps
        new_cgps = update_manager_file_locations(
            get_cgps_json(all=True),
            tempdir,
            replace_dir=str(RESULTS_DIR / "data"),
        )
        with open(tempdir / "cgps_waves.json", "w") as tw:
            json.dump(new_cgps, tw)
        os.mkdir(tempdir / "cGPS")
        for o, n in zip(get_cgps_json(all=True), new_cgps):
            shutil.copyfile(o["file"], n["file"])

        input_chen_near_field(CMT, SAMPLE_FILTER, "cgps", tempdir)
        for f in [
            "filter_cgps.txt",
            "channels_cgps.txt",
            "wavelets_cgps.txt",
            "waveforms_cgps.txt",
        ]:
            with open(tempdir / f, "r") as nd:
                new_data = nd.read()
            with open(RESULTS_DIR / f, "r") as t:
                target = t.read()
            assert new_data == target

    finally:
        shutil.rmtree(tempdir)


def test_input_chen_static():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_static = get_static_json(all=True)
        with open(tempdir / "static_data.json", "w") as tw:
            json.dump(new_static, tw)
        shutil.copyfile(RESULTS_DIR / "NP1" / "gps_data", tempdir / "gps_data")

        input_chen_static(tempdir)
        for f in ["static_data.txt"]:
            with open(tempdir / f, "r") as nd:
                new_data = nd.read()
            with open(RESULTS_DIR / f, "r") as t:
                target = t.read()
            assert new_data == target
    finally:
        shutil.rmtree(tempdir)


def test_input_chen_surf_body():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_surfwaves = update_manager_file_locations(
            get_surf_waves_json(all=True),
            tempdir,
            replace_dir=str(RESULTS_DIR / "data"),
        )
        with open(tempdir / "surf_waves.json", "w") as tw:
            json.dump(new_surfwaves, tw)
        with open(
            pathlib.Path(__file__).parent / "data" / "config.ini",
            "r",
        ) as c:
            config = c.read()
        with open(tempdir / "config.ini", "w") as c:
            c.write(config.replace("/home/user/neic-finitefault", str(tempdir)))
        os.mkdir(tempdir / "RAYLEIGH")
        os.mkdir(tempdir / "LOVE")
        for o, n in zip(get_surf_waves_json(all=True), new_surfwaves):
            shutil.copyfile(o["file"], n["file"])

        input_chen_tele_surf(
            CMT, SAMPLE_FILTER, directory=tempdir, config_path=tempdir / "config.ini"
        )
        for f in [
            "surf_filter.txt",
            "channels_surf.txt",
            "wavelets_surf.txt",
            "Wavelets_surf_tele.txt",
        ]:
            with open(tempdir / f, "r") as nd:
                new_data = nd.read()
            with open(RESULTS_DIR / f, "r") as t:
                target = t.read()
            if f == "wavelets_surf.txt":
                new_data = "\n".join(new_data.split("\n")[2:])
                target = "\n".join(target.split("\n")[1:])
            assert new_data == target
    finally:
        shutil.rmtree(tempdir)


def test_input_chen_tele_body():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        target_dir = RESULTS_DIR / "NP1"
        new_tele_waves = update_manager_file_locations(
            get_tele_waves_json(all=True),
            tempdir,
            replace_dir=str(RESULTS_DIR / "data"),
        )
        with open(tempdir / "tele_waves.json", "w") as tw:
            json.dump(new_tele_waves, tw)
        os.mkdir(tempdir / "P")
        os.mkdir(tempdir / "SH")
        for o, n in zip(get_tele_waves_json(all=True), new_tele_waves):
            shutil.copyfile(o["file"], n["file"])

        input_chen_tele_body(CMT, SAMPLE_FILTER, tempdir)
        for f in [
            "Wavelets_tele_body.txt",
            "body_wave_weight.txt",
            "channels_body.txt",
            "filter_tele.txt",
            "instrumental_response.txt",
            "wavelets_body.txt",
        ]:
            with open(tempdir / f, "r") as nd:
                new_data = nd.read()
            with open(RESULTS_DIR / f, "r") as t:
                target = t.read()
            assert new_data == target
    finally:
        shutil.rmtree(tempdir)


def test_model_space():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        with open(RESULTS_DIR / "NP1" / "model_space.json", "r") as m:
            segments = json.load(m)
        model_space(segments, tempdir)
        for f in [
            "model_space.txt",
            "regularization_borders.txt",
            "special_model_space.txt",
            "special_regularization_borders.txt",
        ]:
            with open(tempdir / f, "r") as d:
                data = d.read()
            with open(RESULTS_DIR / "NP1" / f, "r") as t:
                target = t.read()
            assert data == target
    finally:
        shutil.rmtree(tempdir)


def test_plane_for_chen():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        segments = get_segments_data()
        velmodel = get_velmodel_data()
        rupt_vel = segments["segments"][0]["rupture_vel"]
        lambda_min = 0.4
        lambda_max = 1.25
        min_vel, max_vel = [lambda_min * rupt_vel, lambda_max * rupt_vel]
        plane_for_chen(CMT, segments, min_vel, max_vel, velmodel, tempdir)
        for f in ["fault&rise_time.txt", "point_sources.txt", "shear_model.txt"]:
            print(f)
            with open(tempdir / f, "r") as d:
                data = d.read()
            with open(RESULTS_DIR / "NP1" / f, "r") as t:
                target = t.read()
            assert data == target
    finally:
        shutil.rmtree(tempdir)


def test_from_synthetic_to_obs():
    ## Don't have an example of this all the way through
    ## TODO: add baseline test data
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        os.makedirs(tempdir / "cGPS")
        os.makedirs(tempdir / "P")
        os.makedirs(tempdir / "SH")
        os.makedirs(tempdir / "RAYLEIGH")
        os.makedirs(tempdir / "LOVE")
        os.makedirs(tempdir / "STR")
        for data_type, method in zip(
            ["cgps", "strong", "surf", "body"],
            [
                get_cgps_json,
                get_strong_motion_json,
                get_surf_waves_json,
                get_tele_waves_json,
            ],
        ):
            d = update_manager_file_locations(
                method(),
                tempdir,
                replace_dir=str(RESULTS_DIR / "data"),
            )
            for a, b in zip(method(), d):
                shutil.copyfile(a["file"], b["file"])
            files = d
            for f in [
                "synthetics_strong.txt",
                "synthetics_cgps.txt",
                "synthetics_body.txt",
                "synthetics_surf.txt",
            ]:
                shutil.copyfile(RESULTS_DIR / "NP1" / f, tempdir / f)
            from_synthetic_to_obs(
                files=files,
                data_type=data_type,
                data_prop=SAMPLE_FILTER,
                add_error=False,
                directory=tempdir,
            )
        for f in [
            "waveforms_strong.txt",
            "waveforms_cgps.txt",
            "waveforms_body.txt",
            "waveforms_surf.txt",
        ]:
            assert (tempdir / f).exists()

    finally:
        shutil.rmtree(tempdir)


def test_write_files_wavelet_observed():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        inputs_simmulated_annealing(
            {
                "seismic_moment": "test_moment",
                "moment_weight": "test_moment_weight",
                "slip_weight": "test_slip_weight",
                "time_weight": "test_time_weight",
                "max_source_dur": "test_source",
                "iterations": "test_iter",
                "cooling_rate": "test_cool",
                "initial_temperature": "test_temp",
            },
            tempdir,
        )
        with open(tempdir / "annealing.txt") as a:
            data = a.read()
        target = """
        test_iter -7 1 test_moment 90
        test_temp test_cool 4e-06 0.1 test_moment_weight test_slip_weight test_time_weight
        0 0.0001 0 test_source
        1"""
        data == target
    finally:
        shutil.rmtree(tempdir)


def test_write_green_file():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        with open(RESULTS_DIR / "strong_motion_gf.json") as s:
            sm_gf = json.load(s)
        with open(RESULTS_DIR / "cgps_gf.json") as c:
            cgps_gf = json.load(c)

        # test strong motion
        write_green_file(sm_gf, False, tempdir)
        with open(tempdir / "Green_strong.txt") as gs:
            green_strong = "\n".join(gs.read().split("\n")[0:-1])
        with open(RESULTS_DIR / "Green_strong.txt") as gst:
            green_strong_target = gst.read()
        assert green_strong == green_strong_target

        # test cgps
        write_green_file(cgps_gf, True, tempdir)
        with open(tempdir / "Green_cgps.txt") as gc:
            green_cgps = "\n".join(gc.read().split("\n")[0:-1])
        with open(RESULTS_DIR / "Green_cgps.txt") as gct:
            green_cgps_target = gct.read()
        assert green_cgps == green_cgps_target
    finally:
        shutil.rmtree(tempdir)


def test_write_velmodel():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        write_velmodel(get_velmodel_data(), directory=tempdir)
        with open(RESULTS_DIR / "NP1" / "vel_model.txt", "r") as v:
            target = v.read()
        with open(tempdir / "vel_model.txt", "r") as v:
            vmodel = v.read()
        assert vmodel == target
    finally:
        shutil.rmtree(tempdir)

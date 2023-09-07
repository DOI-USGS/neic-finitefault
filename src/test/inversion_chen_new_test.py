import glob
import json
import logging
import os
import pathlib
import shutil
import subprocess
import tempfile

import numpy as np
import pytest
from obspy import read

from wasp.data_management import filling_data_dicts
from wasp.green_functions import fk_green_fun1
from wasp.input_files import write_green_file, write_velmodel
from wasp.inversion_chen_new import (
    _automatic2,
    automatic_usgs,
    processing,
    set_directory_structure,
    writing_inputs0,
)
from wasp.management import default_dirs
from wasp.seismic_tensor import get_tensor, planes_from_tensor
from wasp.traces_properties import properties_json
from wasp.velocity_models import select_velmodel

from .testutils import (
    DATA_DIR,
    END_TO_END_DIR,
    HOME,
    RESULTS_DIR,
    get_tensor_info,
    get_velmodel_data,
)

TENSOR = get_tensor_info()
VELMODEL = get_velmodel_data()


def _handle_lowin():
    with open(HOME / "fortran_code" / "gfs_nm" / "long" / "low.in", "r") as f:
        file_data = f.read()
    with open(HOME / "fortran_code" / "gfs_nm" / "long" / "low.in", "a") as f:
        if "fd_bank" not in file_data:
            f.write(
                f"\n{str((HOME  / 'fortran_code' / 'gfs_nm'/ 'long'/'fd_bank').resolve())}"
            )


def _compare(f1, ftarget, tempdir):
    with open(f1) as td:
        target = json.load(td)
    with open(ftarget) as d:
        data = json.load(d)
    if isinstance(target, list):
        for idx in range(len(target)):
            match = False
            for idy in range(len(data)):
                if (
                    data[idy]["name"] == target[idx]["name"]
                    and data[idy]["component"] == target[idx]["component"]
                ):
                    data_dict = data[idy]
                    target_dict = target[idx]
                    del data_dict["file"]
                    del target_dict["file"]
                    for key, target_item in target_dict.items():
                        if data_dict[key] == None and target_item == []:
                            continue
                        assert data_dict[key] == target_item
                    match = True
                    break
            if not match:
                raise Exception(
                    f"No match for {target[idx]['name']} {target[idx]['component']}!"
                )
    else:
        for key, target_item in target.items():
            if key == "location":
                continue
            if isinstance(target_item, list):
                assert json.dumps(data[key]).replace(
                    f"{str(tempdir / '20150916225432' / 'ffm.0')}/", ""
                ) == json.dumps(target_item)
            else:
                assert data[key] == target_item


def _end_to_end(
    tensor_info,
    data_type,
    default_dirs,
    velmodel=None,
    dt_cgps=1.0,
    st_response=True,
    config_path=None,
    directory=pathlib.Path(),
):
    directory = pathlib.Path(directory)
    if "gps" in data_type:
        if os.path.isfile(os.path.join(directory, "data", "gps_data")):
            shutil.copy2(os.path.join(directory, "data", "gps_data"), directory)
    if "insar" in data_type:
        insar_files = glob.glob(os.path.join(directory, "data", "insar_a*txt"))
        insar_files = insar_files + glob.glob(
            os.path.join(directory, "data", "insar_d*txt")
        )
        for file in insar_files:
            if os.path.isfile(file):
                shutil.copy2(file, directory)
    data_dir = directory / "data"
    data_prop = properties_json(tensor_info, dt_cgps=dt_cgps, data_directory=directory)
    processing(
        tensor_info, data_type, data_prop, st_response=st_response, directory=data_dir
    )
    data_folder = os.path.join(directory, "data")
    insar_asc = glob.glob(str(directory) + "/insar_asc*txt")
    insar_asc = None if len(insar_asc) == 0 else insar_asc  # type:ignore
    insar_desc = glob.glob(str(directory) + "/insar_desc*txt")
    insar_desc = None if len(insar_desc) == 0 else insar_desc  # type:ignore
    filling_data_dicts(
        tensor_info,
        data_type,
        data_prop,
        data_folder,
        insar_asc=insar_asc,  # type:ignore
        insar_desc=insar_desc,  # type:ignore
        working_directory=directory,
    )
    writing_inputs0(
        tensor_info, data_type, config_path=config_path, directory=directory
    )
    if not velmodel:
        velmodel = select_velmodel(tensor_info, default_dirs)
    write_velmodel(velmodel, directory=directory)
    gf_bank_str = os.path.join(directory, "GF_strong")
    gf_bank_cgps = os.path.join(directory, "GF_cgps")
    get_gf_bank = default_dirs["strong_motion_gf_bank2"]
    if "cgps" in data_type:
        green_dict = fk_green_fun1(
            data_prop, tensor_info, gf_bank_cgps, cgps=True, directory=directory
        )
        write_green_file(green_dict, cgps=True, directory=directory)
        with open(os.path.join(directory, "logs", "GF_cgps_log"), "w") as out_gf_cgps:
            p1 = subprocess.Popen(
                [get_gf_bank, "cgps", f"{(directory)}/"], stdout=out_gf_cgps
            )
        p1.wait()
    if "strong_motion" in data_type:
        green_dict = fk_green_fun1(
            data_prop, tensor_info, gf_bank_str, directory=directory
        )
        write_green_file(green_dict, directory=directory)
        with open(
            os.path.join(directory, "logs", "GF_strong_log"), "w"
        ) as out_gf_strong:
            p2 = subprocess.Popen(
                [get_gf_bank, "strong", f"{(directory)}/"],
                stdout=out_gf_strong,
            )
        p2.wait()
    files = [
        directory / "Green_strong.txt",
        directory / "Green_cgps.txt",
        directory / "modelling_stats.json",
        directory / "gps_data",
        directory / "strong_motion_gf.json",
        directory / "cgps_gf.json",
        directory / "sampling_filter.json",
    ]
    files2 = glob.glob(str(directory) + "/channels_*txt")
    files3 = glob.glob(str(directory) + "/wavelets_*txt")
    files4 = glob.glob(str(directory) + "/waveforms_*txt")
    files5 = glob.glob(str(directory) + "/*waves.json")
    files6 = glob.glob(str(directory) + "/static*")
    files7 = glob.glob(str(directory) + "/filtro*") + glob.glob(
        str(directory) + "/surf_filter*"
    )
    files8 = [
        directory / "instrumental_response.txt",
        directory / "body_wave_weight.txt",
    ]
    files9 = glob.glob(str(directory) + "/insar*")
    files = (
        files
        + files2
        + files3
        + files4
        + files5
        + files6
        + files7
        + files8
        + files9  # type:ignore
    )
    folders = [directory / "NP1", directory / "NP2"]
    for folder in folders:
        for file in files:  # type:ignore
            if os.path.isfile(file):
                shutil.copy2(file, folder)
    info_np1, info_np2 = planes_from_tensor(tensor_info)
    plane1_folder = directory / "NP1"
    _automatic2(
        tensor_info=tensor_info,
        plane_data=info_np1,
        data_type=data_type,
        data_prop=data_prop,
        default_dirs=default_dirs,
        logger=logging.Logger("testlogger"),
        velmodel=velmodel,
        directory=plane1_folder,
    )


@pytest.mark.skipif(
    os.getenv("CI_REGISTRY") is not None or os.getenv("RUN_END_TO_END", False) == False,
    reason="Takes 2+ hours to run",
)
def test_automatic_usgs():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    _handle_lowin()
    try:
        shutil.copytree(END_TO_END_DIR / "data", tempdir / "data")
        set_directory_structure(TENSOR, directory=tempdir)
        for file in os.listdir(tempdir / "data"):
            if os.path.isfile(os.path.join(tempdir / "data", file)):
                shutil.copy2(
                    os.path.join(tempdir / "data", file),
                    tempdir / "20150916225432" / "ffm.0" / "data",
                )
        ddirs = json.dumps(default_dirs(config_path=DATA_DIR / "config.ini"))
        with open(DATA_DIR / "config.ini") as f:
            config = f.read().replace("/home/user/neic-finitefault", str(HOME))
        with open(tempdir / "config.ini", "w") as wf:
            wf.write(config)
        updated_default_dirs = json.loads(
            ddirs.replace("/home/user/neic-finitefault", str(HOME))
        )
        _end_to_end(
            tensor_info=TENSOR,
            data_type=[
                "cgps",
                "gps",
                "insar",
                "strong_motion",
                "surf_tele",
                "tele_body",
            ],
            dt_cgps=None,
            default_dirs=updated_default_dirs,
            config_path=tempdir / "config.ini",
            directory=tempdir / "20150916225432" / "ffm.0",
        )
        for f in [
            "tele_waves.json",
            "strong_motion_waves.json",
            "outlier_strong_motion_waves.json",
            "surf_waves.json",
            "static_data.json",
            "sampling_filter.json",
            "strong_motion_gf.json",
            "insar_data.json",
            "cgps_waves.json",
            "cgps_gf.json",
        ]:
            with open(RESULTS_DIR / "NP1" / f) as td:
                target = json.load(td)
            with open(tempdir / "20150916225432" / "ffm.0" / f) as d:
                data = json.load(d)
            if isinstance(target, list):
                for idx in range(len(target)):
                    match = False
                    for idy in range(len(data)):
                        if (
                            data[idy]["name"] == target[idx]["name"]
                            and data[idy]["component"] == target[idx]["component"]
                        ):
                            data_dict = data[idy]
                            target_dict = target[idx]
                            del data_dict["file"]
                            del target_dict["file"]
                            for key, target_item in target_dict.items():
                                if data_dict[key] == None and target_item == []:
                                    continue
                                assert data_dict[key] == target_item
                            match = True
                            break
                    if not match:
                        raise Exception(
                            f"No match for {target[idx]['name']} {target[idx]['component']}!"
                        )
            else:
                for key, target_item in target.items():
                    if key == "location":
                        continue
                    if isinstance(target_item, list):
                        assert json.dumps(data[key]).replace(
                            f"{str(tempdir / '20150916225432' / 'ffm.0')}/", ""
                        ) == json.dumps(target_item)
                    else:
                        assert data[key] == target_item
        with open(tempdir / "20150916225432" / "ffm.0" / "NP1" / "Solucion.txt") as f:
            solucion = f.read()
        with open(RESULTS_DIR / "NP1" / "Solucion.txt") as t:
            target = t.read()
        assert solucion == target
        # compare processed cGPS waveforms
        data_dir = RESULTS_DIR / "data"
        waveforms = glob.glob(str(data_dir / "cGPS") + "/*.sac")
        for target_file in waveforms:
            basename = os.path.basename(target_file)
            stream = read(
                str(tempdir / "20150916225432" / "ffm.0" / "data" / "cGPS" / basename)
            )
            target_stream = read(str(target_file))
            np.testing.assert_array_equal(stream[0].data, target_stream[0].data)
        # compare processed strong motion waveforms
        data_dir = RESULTS_DIR / "data"
        waveforms = glob.glob(str(data_dir / "STR") + "/*.sac")
        for target_file in waveforms:
            basename = os.path.basename(target_file)
            stream = read(
                str(tempdir / "20150916225432" / "ffm.0" / "data" / "STR" / basename)
            )
            target_stream = read(str(target_file))
            np.testing.assert_array_equal(stream[0].data, target_stream[0].data)
        # compare processed teleseismic waveforms
        data_dir = RESULTS_DIR / "data"
        for tp in ["P", "SH", "LONG"]:
            waveforms = glob.glob(str(data_dir / tp) + "/*.sac")
            for target_file in waveforms:
                basename = os.path.basename(target_file)
                stream = read(
                    str(tempdir / "20150916225432" / "ffm.0" / "data" / tp / basename)
                )
                target_stream = read(str(target_file))
                # TODO: investigate why KOWA surface waves are different
                if tp == "LONG" and "KOWA" in str(target_file):
                    continue
                np.testing.assert_array_equal(stream[0].data, target_stream[0].data)
    finally:
        shutil.rmtree(tempdir)


@pytest.mark.skipif(
    os.getenv("CI_REGISTRY") is not None or os.getenv("RUN_ALL", False) == False,
    reason="Takes 25+ minutes to run",
)
def test_automatic_cgps():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    _handle_lowin()
    try:
        shutil.copytree(END_TO_END_DIR / "data", tempdir / "data")
        set_directory_structure(TENSOR, directory=tempdir)
        for file in os.listdir(tempdir / "data"):
            if os.path.isfile(os.path.join(tempdir / "data", file)):
                shutil.copy2(
                    os.path.join(tempdir / "data", file),
                    tempdir / "20150916225432" / "ffm.0" / "data",
                )
        ddirs = json.dumps(default_dirs(config_path=DATA_DIR / "config.ini"))
        with open(DATA_DIR / "config.ini") as f:
            config = f.read().replace("/home/user/neic-finitefault", str(HOME))
        with open(tempdir / "config.ini", "w") as wf:
            wf.write(config)
        updated_default_dirs = json.loads(
            ddirs.replace("/home/user/neic-finitefault", str(HOME))
        )
        automatic_usgs(
            tensor_info=TENSOR,
            data_type=["cgps"],
            dt_cgps=None,
            default_dirs=updated_default_dirs,
            config_path=tempdir / "config.ini",
            directory=tempdir / "20150916225432" / "ffm.0",
        )
        # compare json files
        for f in [
            "cgps_waves.json",
            "cgps_gf.json",
        ]:
            _compare(
                RESULTS_DIR / "NP1" / f,
                tempdir / "20150916225432" / "ffm.0" / f,
                tempdir,
            )
        # compare solucion
        with open(tempdir / "20150916225432" / "ffm.0" / "NP1" / "Solucion.txt") as f:
            solucion = f.read()
        with open(RESULTS_DIR / "NP1" / "Solucion_cgps.txt", "r") as f:
            target_solucion = f.read()
        assert solucion == target_solucion
        # compare processed waveforms
        data_dir = RESULTS_DIR / "data"
        waveforms = glob.glob(str(data_dir / "cGPS") + "/*.sac")
        for target_file in waveforms:
            basename = os.path.basename(target_file)
            stream = read(
                str(tempdir / "20150916225432" / "ffm.0" / "data" / "cGPS" / basename)
            )
            target_stream = read(str(target_file))
            np.testing.assert_array_equal(stream[0].data, target_stream[0].data)
    finally:
        shutil.rmtree(tempdir)


def test_automatic_gps():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    _handle_lowin()
    try:
        shutil.copytree(END_TO_END_DIR / "data", tempdir / "data")
        set_directory_structure(TENSOR, directory=tempdir)
        for file in os.listdir(tempdir / "data"):
            if os.path.isfile(os.path.join(tempdir / "data", file)):
                shutil.copy2(
                    os.path.join(tempdir / "data", file),
                    tempdir / "20150916225432" / "ffm.0" / "data",
                )
        ddirs = json.dumps(default_dirs(config_path=DATA_DIR / "config.ini"))
        with open(DATA_DIR / "config.ini") as f:
            config = f.read().replace("/home/user/neic-finitefault", str(HOME))
        with open(tempdir / "config.ini", "w") as wf:
            wf.write(config)
        updated_default_dirs = json.loads(
            ddirs.replace("/home/user/neic-finitefault", str(HOME))
        )
        automatic_usgs(
            tensor_info=TENSOR,
            data_type=[
                "gps",
            ],
            dt_cgps=None,
            default_dirs=updated_default_dirs,
            config_path=tempdir / "config.ini",
            directory=tempdir / "20150916225432" / "ffm.0",
        )
        # compare json files
        for f in [
            "static_data.json",
        ]:
            _compare(
                RESULTS_DIR / "NP1" / f,
                tempdir / "20150916225432" / "ffm.0" / f,
                tempdir,
            )
        # compare solucion
        with open(tempdir / "20150916225432" / "ffm.0" / "NP1" / "Solucion.txt") as f:
            solucion = f.read()
        with open(RESULTS_DIR / "NP1" / "Solucion_gps.txt", "r") as f:
            target_solucion = f.read()
        assert solucion == target_solucion
    finally:
        shutil.rmtree(tempdir)


@pytest.mark.skipif(
    os.getenv("CI_REGISTRY") is not None or os.getenv("RUN_ALL", False) == False,
    reason="Takes 18+ minutes to run",
)
def test_automatic_insar():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    _handle_lowin()
    try:
        shutil.copytree(END_TO_END_DIR / "data", tempdir / "data")
        set_directory_structure(TENSOR, directory=tempdir)
        for file in os.listdir(tempdir / "data"):
            if os.path.isfile(os.path.join(tempdir / "data", file)):
                shutil.copy2(
                    os.path.join(tempdir / "data", file),
                    tempdir / "20150916225432" / "ffm.0" / "data",
                )
        ddirs = json.dumps(default_dirs(config_path=DATA_DIR / "config.ini"))
        with open(DATA_DIR / "config.ini") as f:
            config = f.read().replace("/home/user/neic-finitefault", str(HOME))
        with open(tempdir / "config.ini", "w") as wf:
            wf.write(config)
        updated_default_dirs = json.loads(
            ddirs.replace("/home/user/neic-finitefault", str(HOME))
        )
        automatic_usgs(
            tensor_info=TENSOR,
            data_type=["insar"],
            dt_cgps=None,
            default_dirs=updated_default_dirs,
            config_path=tempdir / "config.ini",
            directory=tempdir / "20150916225432" / "ffm.0",
        )
        # compare json files
        for f in [
            "insar_data.json",
        ]:
            _compare(
                RESULTS_DIR / "NP1" / f,
                tempdir / "20150916225432" / "ffm.0" / f,
                tempdir,
            )
        # compare solucion
        with open(tempdir / "20150916225432" / "ffm.0" / "NP1" / "Solucion.txt") as f:
            solucion = f.read()
        with open(RESULTS_DIR / "NP1" / "Solucion_insar.txt", "r") as f:
            target_solucion = f.read()
        assert solucion == target_solucion
    finally:
        shutil.rmtree(tempdir)


@pytest.mark.skipif(
    os.getenv("CI_REGISTRY") is not None or os.getenv("RUN_ALL", False) == False,
    reason="Takes 26+ minutes to run",
)
def test_automatic_strong_motion():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    _handle_lowin()
    try:
        shutil.copytree(END_TO_END_DIR / "data", tempdir / "data")
        set_directory_structure(TENSOR, directory=tempdir)
        for file in os.listdir(tempdir / "data"):
            if os.path.isfile(os.path.join(tempdir / "data", file)):
                shutil.copy2(
                    os.path.join(tempdir / "data", file),
                    tempdir / "20150916225432" / "ffm.0" / "data",
                )
        ddirs = json.dumps(default_dirs(config_path=DATA_DIR / "config.ini"))
        with open(DATA_DIR / "config.ini") as f:
            config = f.read().replace("/home/user/neic-finitefault", str(HOME))
        with open(tempdir / "config.ini", "w") as wf:
            wf.write(config)
        updated_default_dirs = json.loads(
            ddirs.replace("/home/user/neic-finitefault", str(HOME))
        )
        automatic_usgs(
            tensor_info=TENSOR,
            data_type=[
                "strong_motion",
            ],
            dt_cgps=None,
            default_dirs=updated_default_dirs,
            config_path=tempdir / "config.ini",
            directory=tempdir / "20150916225432" / "ffm.0",
        )
        # compare json files
        for f in [
            "strong_motion_waves.json",
            "outlier_strong_motion_waves.json",
        ]:
            _compare(
                RESULTS_DIR / "NP1" / f,
                tempdir / "20150916225432" / "ffm.0" / f,
                tempdir,
            )
        # compare solucion
        with open(tempdir / "20150916225432" / "ffm.0" / "NP1" / "Solucion.txt") as f:
            solucion = f.read()
        with open(RESULTS_DIR / "NP1" / "Solucion_strong_motion.txt", "r") as f:
            target_solucion = f.read()
        assert solucion == target_solucion
        # compare processed waveforms
        data_dir = RESULTS_DIR / "data"
        waveforms = glob.glob(str(data_dir / "STR") + "/*.sac")
        for target_file in waveforms:
            basename = os.path.basename(target_file)
            stream = read(
                str(tempdir / "20150916225432" / "ffm.0" / "data" / "STR" / basename)
            )
            target_stream = read(str(target_file))
            np.testing.assert_array_equal(stream[0].data, target_stream[0].data)
    finally:
        shutil.rmtree(tempdir)


def test_automatic_tele():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    _handle_lowin()
    try:
        shutil.copytree(END_TO_END_DIR / "data", tempdir / "data")
        set_directory_structure(TENSOR, directory=tempdir)
        for file in os.listdir(tempdir / "data"):
            if os.path.isfile(os.path.join(tempdir / "data", file)):
                shutil.copy2(
                    os.path.join(tempdir / "data", file),
                    tempdir / "20150916225432" / "ffm.0" / "data",
                )
        ddirs = json.dumps(default_dirs(config_path=DATA_DIR / "config.ini"))
        with open(DATA_DIR / "config.ini") as f:
            config = f.read().replace("/home/user/neic-finitefault", str(HOME))
        with open(tempdir / "config.ini", "w") as wf:
            wf.write(config)
        updated_default_dirs = json.loads(
            ddirs.replace("/home/user/neic-finitefault", str(HOME))
        )
        automatic_usgs(
            tensor_info=TENSOR,
            data_type=[
                "surf_tele",
                "tele_body",
            ],
            dt_cgps=None,
            default_dirs=updated_default_dirs,
            config_path=tempdir / "config.ini",
            directory=tempdir / "20150916225432" / "ffm.0",
        )
        # compare json files
        for f in [
            "surf_waves.json",
            "sampling_filter.json",
        ]:
            _compare(
                RESULTS_DIR / "NP1" / f,
                tempdir / "20150916225432" / "ffm.0" / f,
                tempdir,
            )
        # compare solucion
        with open(tempdir / "20150916225432" / "ffm.0" / "NP1" / "Solucion.txt") as f:
            solucion = f.read()
        with open(RESULTS_DIR / "NP1" / "Solucion_tele.txt", "r") as f:
            target_solucion = f.read()
        assert solucion == target_solucion
        # compare processed waveforms
        data_dir = RESULTS_DIR / "data"
        for tp in ["P", "SH", "LONG"]:
            waveforms = glob.glob(str(data_dir / tp) + "/*.sac")
            for target_file in waveforms:
                basename = os.path.basename(target_file)
                stream = read(
                    str(tempdir / "20150916225432" / "ffm.0" / "data" / tp / basename)
                )
                target_stream = read(str(target_file))
                # TODO: investigate why KOWA surface waves are different
                if tp == "LONG" and "KOWA" in str(target_file):
                    continue
                np.testing.assert_array_equal(stream[0].data, target_stream[0].data)
    finally:
        shutil.rmtree(tempdir)

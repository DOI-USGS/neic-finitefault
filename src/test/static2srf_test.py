import pathlib
import shutil
import tempfile

from wasp import get_outputs
from wasp.static2srf import static_to_srf

from .testutils import (
    RESULTS_DIR,
    get_segments_data,
    get_tensor_info,
    get_velmodel_data,
)


def test_static_to_srf():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solucion.txt", tempdir / "Solucion.txt")
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "strong_motion_waves.json",
            tempdir / "strong_motion_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "cgps_waves.json",
            tempdir / "cgps_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "static_data.json",
            tempdir / "static_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tele_waves.json",
            tempdir / "tele_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "surf_waves.json",
            tempdir / "surf_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.json",
            tempdir / "insar_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.txt",
            tempdir / "insar_data.txt",
        )
        segments = get_segments_data()
        solution = get_outputs.read_solution_static_format(
            segments["segments"], tempdir
        )
        static_to_srf(
            get_tensor_info(),
            segments,
            ["cgps", "gps", "insar", "strong_motion", "surf_tele", "tele_body"],
            get_velmodel_data(),
            solution,
            directory=tempdir,
        )

        with open(tempdir / "srf_sol_file.txt") as f:
            data = f.readlines()[:40]
        target = [
            "2.0\n",
            "#\n",
            "# Data :\tBODY\tSURF\tSTRONG\tcGPS\tGPS\tInSAR\tDART\tTRIL\tLEVEL\tOTHER\n",
            "# NoS  :\t19\t19\t3\t3\t10\t2166\t0\t0\t0\t0\n",
            "# PHImx :\t186.81\t186.81\t246.51\t312.71\t183.06\t2\t0.00\t0.0\t0.0\t0.0\n",
            "# Rmin :\t40.80\t40.80\t1.11\t0.89\t0.36\t--\t0.00\t0.0\t0.0\t0.0\n",
            "# ---------------------------------- FINITE-SOURCE RUPTURE MODEL ---------------------------------- \n",
            "#\n",
            "# Event : NEAR COAST OF CENTRAL CHILE 2015-09-16T22:54:32 NEIC\n",
            "# EventTAG: 2015-09-16T22:54:32\n",
            "#\n",
            "# Loc  : LAT = -31.57  LON = -71.67  DEP = 22.4\n",
            "#\n",
            "# ----------------------------------  ---------------------------------- \n",
            "#\n",
            "# VELOCITY-DENSITY STRUCTURE\n",
            "# No. of layers = 6\n",
            "#\n",
            "# DEPTH P_VEL S_VEL DENS QP QS\n",
            "# [km] [km/s] [km/s] [g/cm^3]\n",
            "# 0.00 3.35 1.44 2.04  1200.0  600.0\n",
            "# 0.07 6.23 3.61 2.71  1200.0  600.0\n",
            "# 12.15 6.75 3.87 2.83  1200.0  600.0\n",
            "# 25.09 7.65 4.36 2.97  1200.0  600.0\n",
            "# 40.90 8.08 4.47 3.38  1200.0  500.0\n",
            "# 236.90 8.59 4.66 3.45  360.0  140.0\n",
            "#\n",
            "# ----------------------------------  ---------------------------------- \n",
            "# 21/8/2023 created by degoldberg@usgs.gov\n",
            "#\n",
            "# SOURCE MODEL PARAMETERS\n",
            "PLANE 1\n",
            "-72.1958 \t -31.0308 \t 23 \t 9 \t 412.33 \t 134.32 \n",
            "6.61 \t 19.28 \t 0.22 \t -53.67 \t 67.16 \n",
            "POINTS 207\n",
            "-72.4357 -32.7937 2.69 6.61 19.28 267.56 123.94 1.0 -1 -1\n",
            "121.47 67.7587 31 0 0 0 0 \n",
            "  0.000000e+00  4.935634e-02  1.952683e-01  4.313587e-01  7.473094e-01  1.129312e+00\n",
            "  1.560670e+00  2.022533e+00  2.494714e+00  2.956577e+00  3.387935e+00  3.769938e+00\n",
            "  4.085888e+00  4.321979e+00  4.467891e+00  4.517247e+00  4.467891e+00  4.321979e+00\n",
        ]
        for idx, target_line in enumerate(target):
            if "@usgs.gov" in target_line:
                continue
            assert data[idx] == target_line
    finally:
        shutil.rmtree(tempdir)

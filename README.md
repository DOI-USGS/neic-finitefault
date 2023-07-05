Wavelet and simulated Annealing SliP inversion (WASP)

This code uses a wavelet based method to ivnert for slip on a fault using both regional and tele-seismic data. The code is absed on the original tele-seismic approach of Ji et al. (2002) and modified to include regional data (GNSS and strong motion) and for real-time operations by the Centro Sismologico Nacional of Chile. Details of the implementation can be found in Koch et al. (2019) and Goldberg et al. (2022).

Users of this code should consider citing the following relevant publications:

- Ji, C., D. J. Wald, and D. V. Helmberger (2002). Source description of the 1999 Hector Mine, California, earthquake, Part I: Wavelet domain inversion theory and resolution analysis, Bulletin of the Seismological Society of America, 92, no. 4, 1192–1207, https://doi.org/10.1785/0120000916.
- Koch, P., F. Bravo, S. Riquelme, and J. G. F. Crempien (2019). Near-real-time finite-fault inversions for large earthquakes in Chile using strong-motion data, Seismological Research Letters, 90, no. 5, 1971–1986, https://doi.org/10.1785/0220180294.
  -Goldberg, D. E., P. Koch, D. Melgar, S. Riquelme, and W. L. Yeck (2022). Beyond the Teleseism: Introducing Regional Seismic and Geodetic Data into Routine USGS FiniteFault Modeling, Seismological Research Letters, 93, 3308–3323, https://doi.org/10.1785/0220220047.
  -Zhu, L., & Rivera, L. A. (2002). A note on the dynamic and static displacements from a point source in multilayered media: A note on the dynamic and static displacements from a point source. Geophysical Journal International, 148(3), 619–627. https://doi.org/10.1046/j.1365-246X.2002.01610.x.

The surface wave GF bank (fd_bank) is now available in Zenodo, where the user needs to search for "Surface wave Green's functions for WASP slip inversion code" or go directly to: https://zenodo.org/record/7236739#.Y9q4BOzMKDV
Put the fd_bank file in fortran_code/gfs_nm/long/.

Example json files for the modelling of multi-segment faults are available now in multi_segment_example. With these files, the user can re-create the geometry for the solution of the 2011 Tohoku event published by the USGS.

## Installation

## Prior Dependencies

The following dependencies are not handled by any install scripts and/or package managers and must be available on your system prior to installing the dependencies and code:

1. [Python](https://www.python.org/downloads/). Currently supporting versions with a [security status](https://devguide.python.org/versions/).
2. [Poetry](https://python-poetry.org/)

## Install Scripts and

Automated installation of the dependencies and fortran code has been provided in the form of an [install script](./install.sh). Currently this install script only supports installation on linux systems (specifically Ubuntu for the system packages). Installation of Python dependencies and code is managed with the provided Poetry environment setupy by pyproject.toml and package-lock.json. To install the dependencies and code run the two commands:

1. `sudo ./install.sh` (with optional configurations)
2. `sudo ./environment.d/<operating_system>.sh`
   1. If you want the configurations in to be loaded automatically consider adding them to your .bashrc or .bash_profile: `echo "source /home/user/neic-finitefault/environment.d/ubuntu.sh" >> ~/.bashrc`
3. `poetry install`

> Note 1: the installation of system packages, GEOS, GMT, and PROJ, requires that the install script be run as root. A full list of configurations can be found by running `sudo ./install.sh --help`
> Note 2: the install scripts can also be run individually if some of the dependencies (e.g. proj, gmt, etc) are already satisfied on your system. The usage for each individual script can be accessed using the `--help` flag (e.g. `./install.d/proj.sh --help`).
> Note 3: if packages cannot initially be found when running the packages script (ubuntu_packages.sh), then your system package manager may need to be updated (`apt update -y`) and/or upgraded (`apt upgrade -y`).

The following documents provide more information about the installation process:

- [Data Dependencies](./docs/data-dependencies.md): Provides a list of data required to run the code
- [Code Dependencies](./docs/code-dependecies.md): Provides a list of dependencies required to run the code
- [Manual Installation](./docs/code-dependecies.md): Provides a list of steps to manually install dependencies and code without reference to a specific operating system.

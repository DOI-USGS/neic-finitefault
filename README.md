- [Wavelet and simulated Annealing SliP inversion (WASP)](#wavelet-and-simulated-annealing-slip-inversion-wasp)
  - [Authors](#authors)
  - [References](#references)
- [Installation](#installation)
  - [Prior Dependencies](#prior-dependencies)
  - [Install Scripts](#install-scripts)
  - [Using the Poetry Environment](#using-the-poetry-environment)
- [Local Testing](#local-testing)

# Wavelet and simulated Annealing SliP inversion (WASP)

---

This code uses a nonlinear simulated annealing inversion method to model slip amplitude, rake, rupture time, and rise time on a discretized fault plane, finding the solution that best fits the observations in the wavelet domain.

WASP currently accomodates observations from (1) teleseismic broadband stations, (2) regional strong-motion accelerometer stations, (3) static and high-rate Global Navigation Satellite Systems stations, and (4) Interferometric Synthetic Aperture Radar.

The code is based on the approach of Ji et al. (2002). Regional Green's functions are calculated using the method of Zhu & Rivera (2002). Details of the implementation can be found in Koch et al. (2019) and Goldberg et al. (2022).

## Authors

- **[Pablo Koch](https://www.cmm.uchile.cl/?cmm_people=pablo-koch)** - [National Seismological Center, University of Chile](https://www.sismologia.cl/)
- **[Dara Goldberg](https://www.usgs.gov/staff-profiles/dara-e-goldberg)** - [National Earthquake Information Center (NEIC) of the Geologic Hazards Science center](https://www.usgs.gov/centers/geohazards)
- **Heather Hunsinger** - [Geologic Hazards Science Center](https://www.usgs.gov/centers/geohazards)
- **[Diego Melgar](https://earthsciences.uoregon.edu/profile/dmelgarm/)** - [Department of Earth Sciences, University of Oregon](https://earthsciences.uoregon.edu)
- **[Sebastian Riquelme](http://www.dgf.uchile.cl/academicas-y-academicos/profesores-expertos)** - [National Seismological Center, University of Chile](https://www.sismologia.cl/)
- **[William Yeck](https://www.usgs.gov/staff-profiles/william-l-yeck)** - [National Earthquake Information Center (NEIC) of the Geologic Hazards Science Center](https://www.usgs.gov/centers/geohazards)
- **[Kirstie Haynie](https://www.usgs.gov/staff-profiles/kirstie-l-haynie)** - [Geologic Hazards Science Center](https://www.usgs.gov/centers/geohazards)

## References

Users of this code should consider citing the following relevant publications:

- Ji, C., D. J. Wald, and D. V. Helmberger (2002). Source description of the 1999 Hector Mine, California, earthquake, Part I: Wavelet domain inversion theory and resolution analysis, Bulletin of the Seismological Society of America, 92, no. 4, 1192–1207, https://doi.org/10.1785/0120000916.
- Koch, P., F. Bravo, S. Riquelme, and J. G. F. Crempien (2019). Near-real-time finite-fault inversions for large earthquakes in Chile using strong-motion data, Seismological Research Letters, 90, no. 5, 1971–1986, https://doi.org/10.1785/0220180294.
- Goldberg, D. E., P. Koch, D. Melgar, S. Riquelme, and W. L. Yeck (2022). Beyond the Teleseism: Introducing Regional Seismic and Geodetic Data into Routine USGS FiniteFault Modeling, Seismological Research Letters, 93, 3308–3323, https://doi.org/10.1785/0220220047.
- Zhu, L., & Rivera, L. A. (2002). A note on the dynamic and static displacements from a point source in multilayered media: A note on the dynamic and static displacements from a point source. Geophysical Journal International, 148(3), 619–627. https://doi.org/10.1046/j.1365-246X.2002.01610.x.

# Installation

## Prior Dependencies

The following dependencies are not handled by any install scripts and/or package managers and must be available on your system prior to installing the dependencies and code:

1. [Python](https://www.python.org/downloads/). Currently supporting versions with a [security status](https://devguide.python.org/versions/).
2. [Poetry](https://python-poetry.org/)

## Install Scripts

Automated installation of the dependencies and fortran code has been provided in the form of an [install script](./install.sh). Currently this install script only supports installation on linux systems (specifically Ubuntu for the system packages). Installation of Python dependencies and code is managed with the provided Poetry environment setup by pyproject.toml and package-lock.json. To install the dependencies and code run the two commands:

1. `sudo bash install.sh <path to the local neic-finitefault repository> <the os (currently only "ubuntu")> <Python major version (e.g. 3.9 or 3.10)>` (with other optional configurations available, run `sudo bash install.sh -h` for the help information)
   1. > NOTE: The scripts in [./install.d](./install.d/) may be run individually to suit the individuals needs. For example, if GMT is already installed, there may be no need to run [gmt.sh](./install.d/gmt.sh).
2. `sudo ./environment.d/<operating_system>.sh`
   1. If you want the configurations in to be loaded automatically consider adding them to your .bashrc or .bash_profile: `echo "source /home/user/neic-finitefault/environment.d/ubuntu.sh" >> ~/.bashrc`
3. `poetry install`
4. Manually install okada_wrapper after other Python dependencies: `pip install okada_wrapper`

> Note 1: the installation of system packages, GEOS, GMT, and PROJ, requires that the install script be run as root. A full list of configurations can be found by running `sudo ./install.sh --help`
> Note 2: the install scripts can also be run individually if some of the dependencies (e.g. proj, gmt, etc) are already satisfied on your system. The usage for each individual script can be accessed using the `--help` flag (e.g. `./install.d/proj.sh --help`).
> Note 3: if packages cannot initially be found when running the packages script (ubuntu_packages.sh), then your system package manager may need to be updated (`apt update -y`) and/or upgraded (`apt upgrade -y`).

The following documents provide more information about the installation process:

- [Data Dependencies](./docs/data-dependencies.md): Provides a list of data required to run the code
- [Code Dependencies](./docs/code-dependecies.md): Provides a list of dependencies required to run the code
- [Manual Installation](./docs/code-dependecies.md): Provides a list of steps to manually install dependencies and code without reference to a specific operating system.

## Using the Poetry Environment

After running `poetry install`, you will need to activate your environment. This can be done a number of ways:

- `poetry shell` (from within the project at the same location of the pyproject.toml)
- `source <path to poetry virtual environments>/<environment name>/bin/activate`

The following commands may also be useful:

- `poetry config --list`: shows your poetry configuration including the path to poetry virtual environments
- `poetry env list`: Shows the virtual environments associated with the project you are in
- `poetry env info`: Shows the information about the currently activated virtual environment

You can also skip activating the environment by prefixing commands run with `poetry run` (Example: `poetry run wasp --help`). These commands must be run from within the project at the same location of the pyproject.toml.

See official [Poetry documentation](https://python-poetry.org/docs/managing-environments/) for a full description of managing environments.

# Local Testing

Tests and linting can both be run locally:

1. To run all python unit tests: `poetry run poe test`
2. To run python linting: `poetry run poe lint`

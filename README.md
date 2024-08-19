- [Wavelet and simulated Annealing SliP inversion (WASP)](#wavelet-and-simulated-annealing-slip-inversion-wasp)
  - [Authors](#authors)
  - [References](#references)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Wasp Installation Scripts](#wasp-installation-scripts)
- [Local Testing](#local-testing)

# Wavelet and simulated Annealing SliP inversion (WASP)

---

This code uses a nonlinear simulated annealing inversion method to model slip amplitude, rake, rupture time, and rise time on a discretized fault plane, finding the solution that best fits the observations in the wavelet domain.

WASP currently accomodates observations from (1) teleseismic broadband stations, (2) regional strong-motion accelerometer stations, (3) static and high-rate Global Navigation Satellite Systems stations, and (4) Interferometric Synthetic Aperture Radar.

The code is based on the approach of Ji et al. (2002). Regional Green's functions are calculated using the method of Zhu & Rivera (2002). Details of the implementation can be found in Koch et al. (2019) and Goldberg et al. (2022).

## Authors

- **[Pablo Koch](https://www.cmm.uchile.cl/?cmm_people=pablo-koch)** - [National Seismological Center, University of Chile](https://www.sismologia.cl/)
- **[Dara Goldberg](https://www.usgs.gov/staff-profiles/dara-e-goldberg)** - [National Earthquake Information Center (NEIC) of the Geologic Hazards Science Center](https://www.usgs.gov/centers/geohazards)
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

## Prerequisites

In order to compile and/or install the source code there are a number of prerequisite requirements:

1. gfortran: To compile the code in [fortran_code](./fortran_code/)
2. cmake: To compile the code in [fortran_code](./fortran_code/)
3. gcc: To provide support to miniconda for compiling c code
4. miniconda/anaconda: To install python dependencies. Conda can be installed using the provided script: [conda_install.sh](./conda_install.sh)

## Wasp Installation Scripts

Automated installation of the dependencies and fortran code has been provided in the form of the install script [install.sh](./install.sh). Currently this install script only supports installation on linux systems as the fortran code cannot be compiled on MacOS. To instal the code please ensure that all of the [prerequisites](#prerequisites) are available and miniconda/anaconda has been initialized

1. `source install.sh <path to the local neic-finitefault repository>` (with other optional configurations available, run `sudo bash user_install.sh -h` for the help information)
   1. > NOTE: The scripts in [./install.d](./install.d/) may be run individually to suit the individuals needs. For example, to only rerun compilation of the fortran you can singularly run [wap.sh](./install.d/wasp.sh).
2. `conda activate ff-env`

The following documents provide more information about the installation process:

- [Data Dependencies](./docs/data-dependencies.md): Provides a list of data required to run the code
- [Code Dependencies](./docs/code-dependencies.md): Provides a list of dependencies required to run the code
- [Manual Installation](./docs/manual-installation.md): Provides a list of steps to manually install dependencies and code without reference to a specific operating system.

# Local Testing

Tests and linting can both be run locally:

1. To run all python unit tests: `poe test`
   1. The full end to end inversion tests take a consideral amount of time to run. As a result, they are skipped by default and can be enabled by setting the following environment variables to "True"
       - RUN_ALL
       - RUN_END_TO_END
2. To run python linting: `poe lint`

# Manual Installation

This document provides manual installation instructions for users that do not wish to use the provided install scripts and/or want to use other package management tools.

- [Manual Installation](#manual-installation)
- [Install Dependencies](#install-dependencies)
  - [Installing Packages on Your System](#installing-packages-on-your-system)
  - [Installing Dependencies with C code](#installing-dependencies-with-c-code)
  - [Installing Python Dependencies](#installing-python-dependencies)
    - [Okada Wrapper](#okada-wrapper)
- [Get the Data Dependencies](#get-the-data-dependencies)
- [Add Configurations](#add-configurations)
- [Compile the Fortran Code](#compile-the-fortran-code)

# Install Dependencies
## Installing Packages on Your System
In order to add the required dependencies on your system you need to be able to compile Fortran code and (depending on the package manager that you use) compile C code. Depending on how some of the dependencies that require C code are installed you may also need their packages (ghostscript, netCDF, gdal, python developer tools, sqlite, etc). For example, on Linux (Ubuntu), the packages can be installed with apt:

```
sudo apt install -y \
  build-essential \
  cmake \
  gcc \
  gfortran \
  ghostscript \
  libnetcdf-dev \
  libgdal-dev \
  "python3.10-dev" \
  sqlite3;
```

## Installing Dependencies with C code
GMT, PROJ, and GEOS all include C code and are required for some of the Python dependencies below (Cartopy, PyGMT, etc). If you choose to use a package manager like [Anaconda](https://www.anaconda.com/), these may be install and compiled for you. However, Python dependencies installed with pip and/or Poetry may require installing these dependencies manually. Links to instructions for installing these dependencies are linked in [./docs/code-dependencies](./code-dependecies.md#other-dependencies). The [install scripts](../install.d/) install these from source.

## Installing Python Dependencies
The full list of Python dependencies can be found in the [pyproject.toml](../pyproject.toml) file and in [./docs/code-dependencies](./code-dependecies.md#python-dependencies). The provided [Poetry](https://python-poetry.org/) environment (in pyproject.toml) can be used to create a virtual environment with these dependencies: `poetry install`. If you prefer to use [Anaconda](https://www.anaconda.com/) to manage your Python dependencies by creating a finite fault environment (`conda create -n ff-env python=3.10`) and then adding the packages (`conda activate ff-env && conda install <package(s)>`).

### Okada Wrapper
The Okada wrapper is sometimes difficult to install along with all other dependencies since the install (setup.py) requires that NumPy is already available. This can create a causality dilemma. As a result, it is recommended that Okada wrapper be installed separately. With the provided [poetry](https://python-poetry.org/) environment, Okada can be installed with a [task](../pyproject.toml#L49): `poetry run poe okada`. Otherwise, Okada can be installed with Pip: `pip install okada-wrapper`.

# Get the Data Dependencies
Download the [required data dependencies](./data-dependecies.md) and put them in the following locations within the source code.
1. fd_bank: `./fortran_code/gfs_nm/long/fd_bank`
2. LITHO1.0.nc: `./fortran_code/info/LITHO1.0.nc`
3. tectonicplates: `./tectonicplates`

# Add Configurations
1. Setup your configuration file at `./config.ini`. An example [config.ini](./examples/config.ini) file has been provided.
2. Add a line at the bottom of [./fortran_code/gfs_nm/long/low.in](../fortran_code/gfs_nm/long/low.in) that points to the location of the fd_bank file. Example: `/home/user/neic-finitefault/fortran_code/gfs_nm/long/fd_bank`

# Compile the Fortran Code
The code in each Fortran directory must be compiled. In Linux (Ubuntu), for example:

```
cd "/home/user/neic-finitefault/fortran_code" \
    && cd bin_inversion_gfortran_f95 \
    && make clean \
    && make \
    && cd .. \
    && cd bin_str_f95 \
    && make clean \
    && make \
    && cd .. \
    && cd src_dc_f95 \
    && make clean \
    && make;
```

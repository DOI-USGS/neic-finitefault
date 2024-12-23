# Manual Installation

This document provides manual installation instructions for users that do not wish to use the provided install scripts and/or want to use other package management tools.

- [Manual Installation](#manual-installation)
- [Install Dependencies](#install-dependencies)
  - [Installing Packages on Your System](#installing-packages-on-your-system)
  - [Installing Python Dependencies](#installing-python-dependencies)
- [Get the Data Dependencies](#get-the-data-dependencies)
- [Add Configurations](#add-configurations)
- [Compile the Fortran Code](#compile-the-fortran-code)

# Install Dependencies
## Installing Packages on Your System
In order to add the required dependencies on your system you need to be able to compile Fortran code and (depending on the package manager that you use) compile C code. Depending on how some of the dependencies that require C code are installed you may also need their packages (ghostscript, netCDF, gdal, python developer tools, sqlite, etc). For example, on Linux (Ubuntu), the packages can be installed with apt:

```
sudo apt install -y \
  cmake \
  curl \
  gcc \
  gfortran \
  && apt clean;
```

## Installing Python Dependencies
The full list of Python dependencies can be found in the provided miniforge/conda environment file [environment.yml](../install.d/environment.yml). These can be installed individually with your preferred package manager or directly with miniforge/conda using the following commands: 

1. `conda env create -f ./install.d/environment.yml`
2. `conda activate ff-env`
3. `pip install -e <path to neic-finitefault directory>`

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

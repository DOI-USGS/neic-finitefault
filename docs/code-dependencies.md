# Code Dependencies

This document outlines the dependencies required to run the code. This may be helpful for those that do not want to use the install script or Poetry for dependency management.

## Python Dependencies

1. Python: Recommend a [stable/security](https://devguide.python.org/versions/) version.
2. [Cartopy](https://scitools.org.uk/cartopy/docs/latest/installing.html): For managing/manipulating spatial shapes and for plotting
3. [matplotlib](https://matplotlib.org/): For plotting
4. [netCDF4](https://unidata.github.io/netcdf4-python/): Provides the `Dataset` object
5. [NumPy](https://numpy.org/): For array computation
6. [ObsPy](https://docs.obspy.org/): For waveform computation
7. [pyGMT](https://www.pygmt.org/dev/index.html): For plotting maps
8. [SciPy](https://scipy.org/): For fundamental algorithms
9. [shapely](https://shapely.readthedocs.io/en/stable/manual.html): For spatial analysis (in ObsPy)

> A list of dependencies and versions are available in conda environment file [environment.yml](../install.d/environment.yml)

## Fortran Dependencies

1. [GFortran](https://fortran-lang.org/en/learn/os_setup/install_gfortran/): For compiling the FORTRAN code
2. make

## C Dependencies

1. gcc


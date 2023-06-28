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

> Some of these dependencies may have their own secondary dependencies. A package/environment manager such as [Poetry](https://python-poetry.org/) or [Anaconda](https://www.anaconda.com/) will help resolve installing all dependencies.

## Fortran Dependencies

1. [GFortran](https://fortran-lang.org/en/learn/os_setup/install_gfortran/): For compiling the FORTRAN code

## C Dependencies

1. build-essential
2. gcc

## Other dependencies

1. [GMT](https://github.com/GenericMappingTools/gmt/blob/master/INSTALL.md): Must be installed separately pyGMT is installed with pip
2. [GEOS](https://geos.readthedocs.io/en/latest/users.html): Must be installed separately if Cartopy is installed with pip
3. [PROJ](https://proj.org/en/9.2/install.html): Must be installed separately if Cartopy is installed with pip

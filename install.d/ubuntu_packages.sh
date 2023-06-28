#! /bin/bash -ex

PYTHON_VERSION=${1:-"3.10"}

# install packages
apt install -y \
  build-essential \
  cmake \
  curl \
  gcc \
  gfortran \
  ghostscript \
  libnetcdf-dev \
  libgdal-dev \
  python${PYTHON_VERSION}-dev \
  sqlite3 \
  && apt clean;



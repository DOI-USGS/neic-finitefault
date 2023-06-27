#! /bin/bash -ex

PYTHON_VERSION=${PYTHON_VERSION:-"3.10"}

# must run as root or use sudo
if [ "$(whoami)" != "root" ]; then
  echo "Run this script as root:"
  echo "sudo ${0}"
  exit 1;
fi

# install packages
apt update -y \
  && apt upgrade -y \
  && apt install -y \
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



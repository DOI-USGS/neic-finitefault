#! /bin/bash -ex
# ==============================================================================
# Define usage
#     Sets up usage function and parses arguments
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Install package with Ubuntu package manager 'apt'. "
    echo "apt update/upgrade may be required if packages can't be found."
    echo "----------------------------------------------------------------"
}
## parse arguments
REQUIRED_ARGS=()
while [[ $# -gt 0 ]]; do
    # shellcheck disable=SC2221,SC2222
    case $1 in
        -h|--help)
            usage
            exit 0
            ;;
        -*|--*)
            echo "Invalid option $1"
            usage
            exit 1
            ;;
        *)
            REQUIRED_ARGS+=("$1")
            shift
            ;;
    esac
done

# ==============================================================================
# Download and install packages with apt
# ==============================================================================
# install packages
echo "Installing packages with apt"
apt install -y \
  build-essential \
  cmake \
  curl \
  gcc \
  gfortran \
  ghostscript \
  git \
  libssl-dev \
  libffi-dev \
  libnetcdf-dev \
  libgdal-dev \
  libgeos-dev \
  python3-dev \
  sqlite3 \
  && apt clean;

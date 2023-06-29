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
    echo ""
    echo "  PYTHON_VERSION string   major version of python available"
    echo "                          default=3.10"
    echo "                          (example: 3.9)"
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
### set positional argument's defaults
set -- "${REQUIRED_ARGS[@]}" # restore positional parameters
PYTHON_VERSION=${1:-"3.10"}

# ==============================================================================
# Download and install packages with apt
# ==============================================================================
# install packages
apt install -y \
  build-essential \
  cmake \
  curl \
  gcc \
  gfortran \
  ghostscript \
  git \
  libnetcdf-dev \
  libgdal-dev \
  "python${PYTHON_VERSION}-dev" \
  sqlite3 \
  && apt clean;



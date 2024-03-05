#!/bin/bash

# ==============================================================================
# Define usage
#     Install conda (if needed) and create the conda environment
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Runs install scripts in install.d and sets the "
    echo "environment with an environment.d file"
    echo "----------------------------------------------------------------"
    echo "  -p,--python string      The python version to install"
    echo "                          default=3.10.13"
    echo "                          (example: -p 3.11)"
}

## parse arguments
REQUIRED_ARGS=()
while [[ $# -gt 0 ]]; do
    # shellcheck disable=SC2221,SC2222
    case $1 in
        -p|--python)
            PYTHON_VERSION="$2"
            shift
            shift
            ;;
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
### set defaults
PYTHON_VERSION=${PYTHON_VERSION:-"3.10.13"}

# Get system information and conda install url
system="$(uname)"
if [ "$system" == 'Linux' ]; then
    profile=~/.bashrc
    mini_conda_url=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    # shellcheck disable=SC1090
    source $profile;
elif [ "$system" == 'FreeBSD' ] || [ "$system" == 'Darwin' ]; then
    profile=~/.bash_profile
    mini_conda_url=https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    # shellcheck disable=SC1090
    source $profile;
else
    echo "Currently only Linux and Unix systems are supported. Exiting."
    exit
fi

# Download and install conda if not available already
conda --version
# shellcheck disable=SC2181
if [ $? -ne 0 ]; then
    echo "Conda cannot be found. Installing..."
    sha256sum ./miniconda.sh;
    curl -L $mini_conda_url -o ./miniconda.sh;

    bash ./miniconda.sh -fbp "${HOME}"/miniconda;
    # shellcheck disable=SC1090
    . "${HOME}"/miniconda/etc/profile.d/conda.sh
    # remove the shell script
    rm ./miniconda.sh;
    # send source to profile
    echo ". ${HOME}/miniconda/etc/profile.d/conda.sh" >> $profile;
    # shellcheck disable=SC1090
    source $profile;
else
    echo "Conda already installed. No need to download/install."
fi

# Create finite fault env
echo "Creating a finite fault environment called 'ff-env'"
conda create -n ff-env python="${PYTHON_VERSION}" geos poetry pygmt;
echo "Environment created. Use 'conda activate ff-env' to start the environment"
conda init;
# shellcheck disable=SC1090
source $profile;

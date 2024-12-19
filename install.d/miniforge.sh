#!/bin/bash

# ==============================================================================
# Define usage
#     Install miniforge (if needed)
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Install Miniforge"
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

# Get system information and conda install url
system="$(uname)"
mini_forge_url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
if [ "$system" == 'Linux' ]; then
    profile=~/.bashrc
    # shellcheck disable=SC1090
    source $profile;
elif [ "$system" == 'FreeBSD' ] || [ "$system" == 'Darwin' ]; then
    profile=~/.bash_profile
    # shellcheck disable=SC1090
    source $profile;
else
    echo "Currently only Linux and Unix systems are supported. Exiting."
    exit
fi

# Download and install miniforge if not available already
conda --version
# shellcheck disable=SC2181
if [ $? -ne 0 ]; then
    echo "Miniforge cannot be found. Installing..."
    sha256sum ./miniforge.sh;
    curl -L "${mini_forge_url}" -o ./miniforge.sh;

    bash ./miniforge.sh -fbp "${HOME}"/miniforge;
    # shellcheck disable=SC1090
    . "${HOME}"/miniforge/etc/profile.d/conda.sh
    # remove the shell script
    rm ./miniforge.sh;
    # send source to profile
    echo ". ${HOME}/miniforge/etc/profile.d/conda.sh" >> $profile;
    # shellcheck disable=SC1090
    source $profile;
else
    echo "Miniforge already installed. No need to download/install."
fi

# shellcheck disable=SC1090
source $profile;
conda init;

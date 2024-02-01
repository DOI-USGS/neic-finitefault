#!/bin/bash

# ==============================================================================
# Define usage
#     Sets up usage function, parses and validates arguments
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Runs install scripts in install.d and sets the "
    echo "environment with an environment.d file"
    echo "----------------------------------------------------------------"
    echo ""
    echo "  FINITEFAULT_DIR string  full path to the local repository"
    echo "                          (example: /home/user/neic-finitefault)"
    echo "  -p,--python number      The python version to install"
    echo "                          default=3.10.13"
    echo "                          (example: -p 3.11)"
    echo "  --fd-bank string        the location of the fd_bank file"
    echo "                          if not specified it will be downloaded"
    echo "                          default=download"
    echo "                          (example: --fd-bank /home/user/fd_bank)"
    echo "  --lith string        the location of the fd_bank file"
    echo "                          if not specified it will be downloaded"
    echo "                          default=download"
    echo "                          (example: --lith /home/user/fd_bank)"
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
        --fd-bank)
            FD_BANK="$2"
            shift
            shift
            ;;
        --lith)
            LITHO1="$2"
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
### set positional arguments
set -- "${REQUIRED_ARGS[@]}" # restore positional parameters
FINITEFAULT_DIR=$1
### set defaults
PYTHON_VERSION=${PYTHON_VERSION:-"3.10.13"}
FD_BANK=${FD_BANK:-"download"}
LITHO1=${LITHO1:-"download"}

# Running install scripts in install.d/
INSTALL_DIR="${FINITEFAULT_DIR}/install.d"
# shellcheck source=./install.d/conda.sh
source "${INSTALL_DIR}/conda.sh" -p "${PYTHON_VERSION}";
# shellcheck source=./install.d/poetry.sh
source "${INSTALL_DIR}/poetry.sh";
# shellcheck source=./install.d/data_dependencies.sh
source "${INSTALL_DIR}/data_dependencies.sh" "${FINITEFAULT_DIR}" --fd-bank "${FD_BANK}" --lith "${LITHO1}";
# shellcheck source=./install.d/wasp.sh
source "${INSTALL_DIR}/wasp.sh" "${FINITEFAULT_DIR}";
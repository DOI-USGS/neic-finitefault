#!/bin/bash

# must run as root or use sudo
if [ "${USER}" != "root" ]; then
  echo "Run this script as root:"
  echo "sudo ${0}"
  exit 1;
fi

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
    echo "  OPERATING_SYSTEM string the operating system to define install"
    echo "                          (example: ubuntu)"
    echo "                          (note: currently only ubuntu supported)"
    echo "  -x,--skip-packages bool skip system specific install"
    echo "                          default=false"
    echo "                          (example: -x, skips ubuntu_packages.sh)"
    echo "  -s,--skip-env bool      skips sourcing environment files"
    echo "                          default=false"
    echo "                          (example: -s, skips environment.d/*)"
    echo "  --dcw string            version of DCW data to support GMT"
    echo "                          default=2.1.1"
    echo "                          (example: --dcw 2.1.1)"
    echo "  --gmt string            version of GMT to install"
    echo "                          default=6.4.0"
    echo "                          (example: --gmt 6.4.0)"
    echo "  --gshhg string          version of GSHHG data to support GMT"
    echo "                          default=2.3.7"
    echo "                          (example: 2.3.7)"
    echo "  --geos string           version of GEOS to install"
    echo "                          default=3.11.2"
    echo "                          (example: --geos 3.11.2)"
    echo "  --proj string           version of PROJ to install"
    echo "                          default=9.2.0"
    echo "                          (example: --proj 9.2.0)"
    echo "  --fd-bank string        the location of the fd_bank file"
    echo "                          if not specified it will be downloaded"
    echo "                          default=download"
    echo "                          (example: --fd-bank /home/user/fd_bank)"
    echo "  --lith string        the location of the lith file"
    echo "                          if not specified it will be downloaded"
    echo "                          default=download"
    echo "                          (example: --lith /home/user/LITHO1.0.nc)"
    echo "  -c,--cleanup bool       remove gmt, proj, and geos source code"
    echo "                             from /opt after they're installed"
    echo "                          default=false"
    echo "                          (example: -c, deletes files from /opt)"
}

## parse arguments
REQUIRED_ARGS=()
while [[ $# -gt 0 ]]; do
    # shellcheck disable=SC2221,SC2222
    case $1 in
        -x|--skip-packages)
            SKIP_SYSTEM_INSTALL=true
            shift
            ;;
        -s|--skip-env)
            SKIP_ENV=true
            shift
            ;;
        --dcw)
            DCW_VERSION="$2"
            shift
            shift
            ;;
        --gmt)
            GMT_VERSION="$2"
            shift
            shift
            ;;
        --gshhg)
            GSHHG_VERSION="$2"
            shift
            shift
            ;;
        --geos)
            GEOS_VERSION="$2"
            shift
            ;;
        --proj)
            PROJ_VERSION="$2"
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
        -c|--cleanup)
            CLEANUP=true
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
OPERATING_SYSTEM=$2
### set defaults
SKIP_SYSTEM_INSTALL=${SKIP_SYSTEM_INSTALL:false}
SKIP_ENV=${SKIP_ENV:false}
DCW_VERSION=${DCW_VERSION:-"2.1.1"}
GMT_VERSION=${GMT_VERSION:-"6.4.0"}
GSHHG_VERSION=${GSHHG_VERSION:-"2.3.7"}
GEOS_VERSION=${GEOS_VERSION:-"3.11.2"}
PROJ_VERSION=${PROJ_VERSION:-"9.2.0"}
FD_BANK=${FD_BANK:-"download"}
LITHO1=${LITHO1:-"download"}
### validate operating system choice
case $OPERATING_SYSTEM in
  ubuntu)
    ;;

  *)
    echo -n "Invalid option for OPERATING_SYSTEM: ${OPERATING_SYSTEM}" && exit
    ;;
esac

# ==============================================================================
# Execute scripts
#     Run install scripts and source environment files
# ==============================================================================


# Running install scripts in install.d/
INSTALL_DIR="${FINITEFAULT_DIR}/install.d"
if [ "$SKIP_SYSTEM_INSTALL" = true ] ; then
    echo "Skipping running ${INSTALL_DIR}/${OPERATING_SYSTEM}_packages.sh";
else
    # shellcheck source=./install.d/ubuntu_packages.sh
    source "${INSTALL_DIR}/${OPERATING_SYSTEM}_packages.sh";
fi
# shellcheck source=./install.d/gmt.sh
source "${INSTALL_DIR}/gmt.sh" "${CLEANUP}" "${DCW_VERSION}" "${GMT_VERSION}" "${GSHHG_VERSION}";
# shellcheck source=./install.d/libgeos.sh
source "${INSTALL_DIR}/libgeos.sh" "${CLEANUP}"  "${GEOS_VERSION}";
# shellcheck source=./install.d/proj.sh
source "${INSTALL_DIR}/proj.sh" "${CLEANUP}"  "${PROJ_VERSION}";
# shellcheck source=./install.d/data_dependencies.sh
source "${INSTALL_DIR}/data_dependencies.sh" "${FINITEFAULT_DIR}" --fd-bank "${FD_BANK}" --lith "${LITHO1}";
# shellcheck source=./install.d/wasp.sh
source "${INSTALL_DIR}/wasp.sh" "${FINITEFAULT_DIR}"

# Source environment.d file for environment variables
ENV_DIR="${FINITEFAULT_DIR}/environment.d"
if [ "$SKIP_ENV" = true ] ; then
    echo "Skipping sourcing env ${ENV_DIR}/${OPERATING_SYSTEM}.sh";
else
    # shellcheck source=./environment.d/ubuntu.sh
    source "${ENV_DIR}/${OPERATING_SYSTEM}.sh";    
fi

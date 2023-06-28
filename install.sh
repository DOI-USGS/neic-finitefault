#! /bin/bash -ex

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
    echo "  PYTHON_VERSION string   the major version of python available"
    echo "                          (example: 3.10)"
    echo "  -x,--skip-packages bool skip system specific install"
    echo "                          default=false"
    echo "                          (example: -x, skips ubuntu.sh)"
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
    echo "  -c,--cleanup bool       remove gmt, proj, and geos source code"
    echo "                             from /opt after they're installed"
    echo "                          default=false"
    echo "                          (example: -c, deletes files from /opt)"
}

## parse arguments
REQUIRED_ARGS=()
while [[ $# -gt 0 ]]; do
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
PYTHON_VERSION=$3
### set defaults
SKIP_SYSTEM_INSTALL=${SKIP_SYSTEM_INSTALL:false}
SKIP_ENV=${SKIP_ENV:false}
DCW_VERSION=${DCW_VERSION:-"2.1.1"}
GMT_VERSION=${GMT_VERSION:-"6.4.0"}
GSHHG_VERSION=${GSHHG_VERSION:-"2.3.7"}
GEOS_VERSION=${GEOS_VERSION:-"3.11.2"}
PROJ_VERSION=${PROJ_VERSION:-"9.2.0"}
CLEANUP=${CLEANUP:-false}
### validate python version choice (in security status)
#### See https://devguide.python.org/versions/ for statuses of versions
case $PYTHON_VERSION in
  3.8)
    ;;
  3.9)
    ;;
  3.10)
    ;;
  *)
    echo -n "Invalid option for PYTHON_VERSION: ${PYTHON_VERSION}" && exit
    ;;
esac
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
    source "${INSTALL_DIR}/${OPERATING_SYSTEM}_packages.sh" "${PYTHON_VERSION}";    
fi
source "${INSTALL_DIR}/gmt.sh.sh" $CLEANUP "${DCW_VERSION}" "${GMT_VERSION}" "${GSHHG_VERSION}";
source "${INSTALL_DIR}/libgeos.sh" $CLEANUP  "${GEOS_VERSION}";
source "${INSTALL_DIR}/proj.sh" $CLEANUP  "${PROJ_VERSION}";

# Source environment.d file for environment variables
ENV_DIR="${FINITEFAULT_DIR}/environment.d"
if [ "$SKIP_ENV" = true ] ; then
    echo "Skipping sourcing env ${ENV_DIR}/${OPERATING_SYSTEM}.sh";
else
    source "${ENV_DIR}/${OPERATING_SYSTEM}.sh";    
fi

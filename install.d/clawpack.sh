#! /bin/bash

set -ex

# Install Clawpack/GeoClaw if requested to facilitate DART observations
# ==============================================================================
# Define usage
#     Sets up usage function, parses and validates arguments
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Install Clawpack/GeoClaw for DART usage"
    echo "----------------------------------------------------------------"
    echo ""
    echo "  CLAW_DIR string  the directory where clawpack will be installed"
    echo "                   (example: /home/user)"
    echo "  --yes install without prompt"
   }


bathymetry_download () {
    TOPO_FILE="${CLAW_DIR}/clawpack/topo30.grd";
    if [ -e "$TOPO_FILE" ]; then
        echo -e "\nBathymetry already exists at ${TOPO_FILE}. Skipping download.";
    else    
        echo -e "\nCopying TOPO30 (Becker et al., 2009)"; 
        curl -o "${TOPO_FILE}" -L "https://topex.ucsd.edu/pub/srtm30_plus/topo30/topo30.grd";
    fi
}

clawpack_install () {
    echo -e "\nInstalling Clawpack";
    git clone https://github.com/clawpack/clawpack.git "${CLAW_DIR}/clawpack";
    cd "${CLAW_DIR}/clawpack";
    git checkout v5.14.0;
    git submodule init;
    git submodule update;

    cd "${OLDPWD}";
    clawpack_profile;

    # topo/bathy data
    bathymetry_download;
}


clawpack_profile () {
    system="$(uname)";
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
    # Check CLAW setup in profile and add if necessary
    if [[ "$CLAW" == "${CLAW_DIR}/clawpack" ]]; then
        echo "CLAW already specified in profile"
    else
        echo "Adding CLAW info to profile"
        echo "{ ; # CLAWPACK; export CLAW=${CLAW_DIR}/clawpack; export OMP_NUM_THREADS=nproc/4 }" >> $profile;
    fi

    if [[ "$PYTHONPATH" == *"${CLAW_DIR}/clawpack"* ]]; then
        echo "PythonPath already contains Clawpack"
    else
        echo "Adding CLAW to PythonPath"
        echo "export PYTHONPATH=${CLAW_DIR}/clawpack:\$PYTHONPATH" >> $profile;
    fi
    if [[ "$FC" == "gfortran" ]]; then
        echo "FC already set to gfortran"
    else
        echo "export FC=gfortran" >> $profile;
    fi
    source $profile;
}


clawpack_feedback () {
  while IFS= read -n1 -r -p "Do you want to install Clawpack? Required for DART observations. [y/n] " && [[ $REPLY != d ]]; do
  case $REPLY in
    y) clawpack_install && break;;
    n) echo -e "\nSkipping Clawpack installation.\n" && break;;
    *) echo -e "\nInvalid input (valid choices: y, n)\n";;
  esac
  done
}


## parse arguments
SILENT=1
REQUIRED_ARGS=()
while [[ $# -gt 0 ]]; do
    # shellcheck disable=SC2221,SC2222
    case $1 in
        -y|--yes)
            SILENT=0
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
CLAW_DIR=${1%%/}
if [ -z ${1+x} ]
then
    echo "Argument CLAW_DIR must be set";
    exit 1;
fi

echo -e "\nChecking for existing Clawpack installation";
if [ -d "${CLAW_DIR}/clawpack" ]; then
  echo "Clawpack found";
  bathymetry_download;
  clawpack_profile;
elif (( SILENT == 0 )); then
  clawpack_install;
else
  echo "Clawpack installation not found";
  clawpack_feedback;
fi


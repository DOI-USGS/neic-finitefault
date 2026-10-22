#! /bin/bash

set -e

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
    echo "  --t, --topo           downlaod the full global topography file"
    echo "                          (local topography can be queried for dynamically using 'ffm get-data topo')"
   }


bathymetry_download () {
    TOPO_FILE="${CLAW_DIR}/clawpack/topo30.grd";
    if [ -e "$TOPO_FILE" ]; then
        echo -e "\nBathymetry already exists at ${TOPO_FILE}. Skipping download.";
    else    
        echo -e "\nCopying TOPO30 (Becker et al., 2009)"; 
        curl -o "${TOPO_FILE}" -L "https://topex.ucsd.edu/pub/srtm30_plus/topo30/topo30.grd";
    fi

    echo GLOBAL_TOPOGRAPHY=${TOPO_FILE} >> ${CLAW_DIR}/.env
}

clawpack_install () {
    # activate ff-env
    
    system="$(uname)"
    if [ "$system" == 'Linux' ]; then
        source ~/.bashrc
    elif [ "$system" == 'FreeBSD' ] || [ "$system" == 'Darwin' ]; then
        source ~/.bash_profile
    fi

    status=0
    conda activate ff-env || status=$?
    if [ "$status" -eq "1" ]; then
        echo "Activating the ff-env conda environment failed."
        echo "Ensure that it has been created and conda init has run."
        exit 1;
    fi
    
    conda install -y meson-python;
    echo "Installing Clawpack";
    pip install --src="${CLAW_DIR}" --user --no-build-isolation -e \
        git+https://github.com/clawpack/clawpack.git@v5.14.0#egg=clawpack

    clawpack_profile;
}


clawpack_profile () {
    system="$(uname)";
    if [ "$system" == 'Linux' ]; then 
        total_cpus=$(nproc)
        # shellcheck disable=SC1090
    elif [ "$system" == 'FreeBSD' ] || [ "$system" == 'Darwin' ]; then 
        total_cpus=$(sysctl -n hw.ncpu)
        # shellcheck disable=SC1090
    else    
        echo "Currently only Linux and Unix systems are supported. Exiting."
        exit    
    fi
    # Check CLAW setup in profile and add if necessary
    OMP_NUM_THREADS=$(( total_cpus / 4 ))
    if [[ "${OMP_NUM_THREADS}" -eq 0 ]]; then
        OMP_NUM_THREADS=1
    fi

    cat << EOF > ${CLAW_DIR}/.env
CLAW=${CLAW_DIR}/clawpack
OMP_NUM_THREADS=${OMP_NUM_THREADS}
PYTHONPATH=${CLAW_DIR}/clawpack:\${PYTHONPATH}
FC=gfortran
EOF
}


clawpack_feedback () {
  while IFS= read -n1 -r -p "Do you want to install Clawpack and/or install global topography? Required for DART observations. [c: clawpack, t: global topography, b: both clawpack and global topography, n: no install] " && [[ $REPLY != d ]]; do
  case $REPLY in
    c) clawpack_install && break;;
    t) bathymetry_download && break;;
    b) clawpack_install && bathymetry_download && break;;
    n) echo -e "\nSkipping Clawpack installation.\n" && break;;
    *) echo -e "\nInvalid input (valid choices: c, t, b, n)\n";;
  esac
  done
}


## parse arguments
TOPO=1
SILENT=1
REQUIRED_ARGS=()
while [[ $# -gt 0 ]]; do
    # shellcheck disable=SC2221,SC2222
    case $1 in
        -y|--yes)
            SILENT=0
            shift
            ;;
        -t|--topo)
            TOPO=0
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
if (( SILENT == 0 )); then
  clawpack_install;
  if (( TOPO == 0 )); then
    bathymetry_download
  fi
else
  echo "Clawpack installation not found";
  clawpack_feedback;
fi


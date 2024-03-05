#!/bin/bash

# ==============================================================================
# Define usage
#     Sets up usage function and parses arguments
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Get data dependencies."
    echo "----------------------------------------------------------------"
    echo ""
    echo "  FINITEFAULT_DIR string  the directory where the finite fault"
    echo "                          repository exists"
    echo "                          (example: /home/user/neic-finitefault)"
    echo "  --fd-bank string        the location of the fd_bank file"
    echo "                          if not specified it will be downloaded"
    echo "                          default=download"
    echo "                          (example: --fd-bank /home/user/fd_bank)"
    echo "  --lith string        the location of the lith file"
    echo "                          if not specified it will be downloaded"
    echo "                          default=download"
    echo "                          (example: --lith /home/user/LITHO1.0.nc)"
   }
## parse arguments
REQUIRED_ARGS=()
while [[ $# -gt 0 ]]; do
    # shellcheck disable=SC2221,SC2222
    case $1 in
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
FD_BANK=${FD_BANK:-"download"}
LITHO1=${LITHO1:-"download"}

# ==============================================================================
# Download required data and compile FORTRAN code
# ==============================================================================
# assume running from top level of the project directory
if [ -z ${FINITEFAULT_DIR+x} ]
then 
    echo "Argument FINITEFAULT_DIR must be set";
    exit 1; 
fi
FORTRAN_DIR="${FINITEFAULT_DIR}/fortran_code"
# Get the required data files 
## fd_bank
FD_FILE="${FORTRAN_DIR}/gfs_nm/long/fd_bank"
if [[ "$FD_BANK" == *"download" ]]; then
    if [ -f "${FD_FILE}" ]; then
        echo "fd_bank already exists at ${FD_FILE}. Skipping download.";
    else
        echo "Downloading the fd_bank file"
        curl -o "${FD_FILE}" -L \
            "https://zenodo.org/record/7236739/files/fd_bank?download=1"
    fi
else
    echo "Copying the provided fd_bank file to the required location"
    cp "${FD_BANK}" "${FD_FILE}"
fi
fdminsize=875050000
fdsize=$(wc -c <"$FD_FILE")
if [ "$fdsize" -ge $fdminsize ]; then
    echo "fd_bank size verified.";
else
    echo "fd_bank size does not match the expected value";
    exit 1;
fi
## LITHO1.0.nc
LITHO1_FILE="${FORTRAN_DIR}/info/LITHO1.0.nc"
if [[ "$LITHO1" == *"download" ]]; then
    if [ -f "${LITHO1_FILE}" ]; then
        echo "LITHO1.0.nc already exists at ${LITHO1_FILE}. Skipping download.";
    else
        echo "Downloading the LITHO1.0.nc file"
        curl -o "${LITHO1_FILE}" -L \
            "https://ds.iris.edu/files/products/emc/emc-files/LITHO1.0.nc"; 
    fi
else
    echo "Copying the provided LITHO1.0.nc file to the required location"
    cp "${LITHO1}" "${LITHO1_FILE}"
fi
lithminsize=44350184
lithsize=$(wc -c <"$LITHO1_FILE")
if [ "$lithsize" -ge $lithminsize ]; then
    echo "LITHO1.0.nc size verified.";
else
    echo "LITHO1.0.nc size does not match the expected value";
    exit 1;
fi
# tectonicsplates master branch
TECTONICS_DIR="${FINITEFAULT_DIR}/tectonicplates"
TECTONICS_URL="https://github.com/fraxen/tectonicplates.git"
echo "Getting a copy of the tectonicplates data from master branch"
git clone --single-branch --branch master "${TECTONICS_URL}" "${TECTONICS_DIR}"


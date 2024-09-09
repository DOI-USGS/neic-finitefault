#!/bin/bash

# ==============================================================================
# Define usage
#     Sets up usage function and parses arguments
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Compile the wasp fortran code."
    echo "----------------------------------------------------------------"
    echo ""
    echo "  FINITEFAULT_DIR string  the directory where the finite fault"
    echo "                          repository exists"
    echo "                          (example: /home/user/neic-finitefault)"
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
### set positional arguments
set -- "${REQUIRED_ARGS[@]}" # restore positional parameters
FINITEFAULT_DIR=${1%%/}
if [ -z ${1+x} ]
then 
    echo "Argument FINITEFAULT_DIR must be set";
    exit 1;
fi

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
# build fortran
echo "Compiling the FORTRAN code with make"
cd "${FORTRAN_DIR}" \
    && cd bin_inversion_gfortran_f95 \
    && make clean \
    && make \
    && cd .. \
    && cd bin_str_f95 \
    && make clean \
    && make \
    && cd .. \
    && cd src_dc_f95 \
    && make clean \
    && make \
    && cd ../..;

# update configuration for fd_bank location
LOWIN_FILE="${FORTRAN_DIR}/gfs_nm/long/low.in"
FD_FILE="${FORTRAN_DIR}/gfs_nm/long/fd_bank"
echo "Updating the location of the fd_bank file in ${LOWIN_FILE}"
if grep -Fq "fd_bank" "${LOWIN_FILE}"
then
    echo "There appears to be a path to fd_bank in ${LOWIN_FILE} already."
    echo "Please check that it is correct. Skipping adding the file line."
else
    echo -e "\n$FD_FILE" >> "$LOWIN_FILE"
fi

# write the config file
echo "Writing the config file"
CONFIG_FILE="${FINITEFAULT_DIR}/config.ini"

cat <<EO_CONFIG > "${CONFIG_FILE}"

[PATHS]
code_path = ${FINITEFAULT_DIR}
surf_gf_bank = %(code_path)s/fortran_code/gfs_nm/long/low.in
modelling = %(code_path)s/fortran_code/bin_inversion_gfortran_f95
get_near_gf = %(code_path)s/fortran_code/bin_str_f95
compute_near_gf = %(code_path)s/fortran_code/src_dc_f95
info = %(code_path)s/fortran_code/info
cartopy_files = %(code_path)s/fortran_code/tectonicplates

EO_CONFIG
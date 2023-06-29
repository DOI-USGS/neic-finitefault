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
FINITEFAULT_DIR=$1

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
echo "HERE ${FINITEFAULT_DIR}"
echo "THERE ${FORTRAN_DIR}"
# Get the required data files 
## fd_bank
FD_FILE="${FORTRAN_DIR}/gfs_nm/long/fd_bank"
curl -o "${FD_FILE}" -L \
    "https://zenodo.org/record/7236739/files/fd_bank?download=1"
## LITHO1.0.nc
LITHO1_FILE="${FORTRAN_DIR}/info/LITHO1.0.nc"
curl -o "${LITHO1_FILE}" -L \
  "https://ds.iris.edu/files/products/emc/emc-files/LITHO1.0.nc"; 
# tectonicsplates master branch
TECTONICS_DIR="${FINITEFAULT_DIR}/tectonicplates"
TECTONICS_URL="https://github.com/fraxen/tectonicplates.git"
git clone --single-branch --branch master "${TECTONICS_URL}" "${TECTONICS_DIR}"

# build fortran
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
grep -qF -- "$FD_FILE" "$LOWIN_FILE" || echo "$FD_FILE" >> "$LOWIN_FILE"

# write the config file
CONFIG_FILE="${FINITEFAULT_DIR}/config.ini"

cat <<EO_CONFIG > "${CONFIG_FILE}"

[PATHS]
code_path = ${FINITEFAULT_DIR}
surf_gf_bank = %(code_path)s/fortran_code/gfs_nm/long/low.in
modelling = %(code_path)s/fortran_code/bin_inversion_gfortran_f95
get_near_gf = %(code_path)s/fortran_code/bin_str_f95
compute_near_gf = %(code_path)s/fortran_code/src_dc_f95
info = %(code_path)s/fortran_code/info
cartopy_files = %(code_path)s/fortran_code/tectonicplates-master

EO_CONFIG
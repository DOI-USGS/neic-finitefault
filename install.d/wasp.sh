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
    echo "Downloading the fd_bank file"
    curl -o "${FD_FILE}" -L \
      "https://zenodo.org/record/7236739/files/fd_bank?download=1"
else
    echo "Copying the provided fd_bank file to the required location"
    cp "${FD_BANK}" "${FD_FILE}"
fi
## LITHO1.0.nc
LITHO1_FILE="${FORTRAN_DIR}/info/LITHO1.0.nc"
if [[ "$LITHO1" == *"download" ]]; then
    echo "Downloading the LITHO1.0.nc file"
    curl -o "${LITHO1_FILE}" -L \
        "https://ds.iris.edu/files/products/emc/emc-files/LITHO1.0.nc"; 
else
    echo "Copying the provided LITHO1.0.nc file to the required location"
    cp "${LITHO1}" "${LITHO1_FILE}"
fi
# tectonicsplates master branch
TECTONICS_DIR="${FINITEFAULT_DIR}/tectonicplates"
TECTONICS_URL="https://github.com/fraxen/tectonicplates.git"
echo "Getting a copy of the tectonicplates data from master branch"
git clone --single-branch --branch master "${TECTONICS_URL}" "${TECTONICS_DIR}"

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
echo "Updating the location of the fd_bank file in ${LOWIN_FILE}"
grep -qF -- "$FD_FILE" "$LOWIN_FILE" || echo "$FD_FILE" >> "$LOWIN_FILE"

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
cartopy_files = %(code_path)s/fortran_code/tectonicplates-master

EO_CONFIG
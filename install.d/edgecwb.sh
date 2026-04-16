#!/bin/bash

# ==============================================================================
# Define usage
#     Sets up usage function and parses arguments
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Get the EdgeCwb CWBQuery.jar Java client."
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
CWB_QUERY_ZIP="${FINITEFAULT_DIR}/CWBQuery.zip";
CWB_QUERY="${FINITEFAULT_DIR}/CWBQuery";
curl -o "${CWB_QUERY_ZIP}" -L \
            "https://code.usgs.gov/ghsc/neic/edgecwb/edgecwbgroup/edgecwbfiles/-/package_files/18579/download";
unzip "${CWB_QUERY_ZIP}" -d "${CWB_QUERY}";
rm "${CWB_QUERY_ZIP}";
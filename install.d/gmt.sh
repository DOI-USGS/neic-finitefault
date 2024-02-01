#!/bin/bash

# ==============================================================================
# Define usage
#     Sets up usage function and parses arguments
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Runs scripts to download GMT source code as well as supporting "
    echo "data. The source code is then compiled with cmake."
    echo "----------------------------------------------------------------"
    echo ""
    echo "  CLEANUP bool            Whether to delete build directories after global install"
    echo "                          default=false"
    echo "                          (example: true)"
    echo "  DCW_VERSION string      version of DCW data to support GMT"
    echo "                          default=2.1.1"
    echo "                          (example: 2.1.1)"
    echo "  GMT_VERSION string      version of GMT to install"
    echo "                          default=6.4.0"
    echo "                          (example: 6.4.0)"
    echo "  GSHHG_VERSION string    version of GSHHG data to support GMT"
    echo "                          default=2.3.7"
    echo "                          (example: 2.3.7)"
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
### set positional argument's defaults
set -- "${REQUIRED_ARGS[@]}" # restore positional parameters
CLEANUP=${1:-false}
DCW_VERSION=${2:-"2.1.1"}
GMT_VERSION=${3:-"6.4.0"}
GSHHG_VERSION=${4:-"2.3.7"}

# ==============================================================================
# Download and install GMT
# ==============================================================================
# Download GMT source code and supporting data
echo "Downloading GMT supporting data DCW ${DCW_VERSION}"
GMT_URL="https://github.com/GenericMappingTools"
curl -o "/opt/dcw-gmt-${DCW_VERSION}.tar.gz" -L \
  "${GMT_URL}/dcw-gmt/releases/download/${DCW_VERSION}/dcw-gmt-${DCW_VERSION}.tar.gz";
cd /opt && tar -xvf "/opt/dcw-gmt-${DCW_VERSION}.tar.gz";
echo "Downloading GMT ${GMT_VERSION}"
curl -o "/opt/gmt-${GMT_VERSION}-src.tar.gz" -L \
  "${GMT_URL}/gmt/releases/download/${GMT_VERSION}/gmt-${GMT_VERSION}-src.tar.gz";
cd /opt && tar -xvf "/opt/gmt-${GMT_VERSION}-src.tar.gz";
echo "Downloading GMT supporting data GSHHG ${GSHHG_VERSION}"
curl -o "/opt/gshhg-gmt-${GSHHG_VERSION}.tar.gz" -L \
  "${GMT_URL}/gshhg-gmt/releases/download/${GSHHG_VERSION}/gshhg-gmt-${GSHHG_VERSION}.tar.gz";
cd /opt && tar -xvf "/opt/gshhg-gmt-${GSHHG_VERSION}.tar.gz";
mv "/opt/gshhg-gmt-${GSHHG_VERSION}" "/opt/gmt-${GMT_VERSION}/share/gshhg-gmt";
mv "/opt/dcw-gmt-${DCW_VERSION}" "/opt/gmt-${GMT_VERSION}/share/dcw-gmt";

# Build the GMT source code
echo "Installing GMT with cmake"
cd "/opt/gmt-${GMT_VERSION}" \
  && mkdir -p build \
  && cd build \
  && cmake .. \
  && cmake --build . \
  && cmake --build . --target install;

# cleanup
if [ "$CLEANUP" = true ] ; then
    echo "Cleaning up source code for GMT.";
    rm -r /opt/dcw-gmt-"${DCW_VERSION}"*;
    rm -r /opt/gshhg-gmt-"${GSHHG_VERSION}"*;
    rm -r /opt/gmt-"${GMT_VERSION}"-src*;
fi
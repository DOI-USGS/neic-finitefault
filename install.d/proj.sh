#!/bin/bash

# ==============================================================================
# Define usage
#     Sets up usage function and parses arguments
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Runs scripts to download PROJ source code. "
    echo "The source code is then compiled with cmake."
    echo "----------------------------------------------------------------"
    echo ""
    echo "  CLEANUP bool            Whether to delete build directories after global install"
    echo "                          default=false"
    echo "                          (example: true)"
    echo "  PROJ_VERSION string     version of PROJ to install"
    echo "                          default=3.11.2"
    echo "                          (example: 3.11.2)"
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
PROJ_VERSION=${2:-"9.2.0"}

# ==============================================================================
# Download and install PROJ
# ==============================================================================
# Download PROJ source
echo "Downloading PROJ ${PROJ_VERSION} source code"
curl -o "/opt/proj-${PROJ_VERSION}.tar.gz" -L \
  "https://download.osgeo.org/proj/proj-${PROJ_VERSION}.tar.gz"; 
cd /opt && tar -xvf "/opt/proj-${PROJ_VERSION}.tar.gz";

# Build the PROJ source code
echo "Installing PROJ with cmake"
cd "/opt/proj-${PROJ_VERSION}"\
    && mkdir -p build \
    && cd build \
    && cmake -DBUILD_TESTING=OFF .. \
    && cmake --build . \
    && cmake --build . --target install;

# cleanup
if [ "$CLEANUP" = true ] ; then
    echo 'Cleaning up source code for PROJ.';
    rm -r /opt/proj-"${PROJ_VERSION}"*;
fi
#! /bin/bash -ex

# ==============================================================================
# Define usage
#     Sets up usage function and parses arguments
# ==============================================================================
## define usage
function usage {
    echo "----------------------------------------------------------------"
    echo "Runs scripts to download GEOS source code. "
    echo "The source code is then compiled with cmake."
    echo "----------------------------------------------------------------"
    echo ""
    echo "  CLEANUP bool            version of DCW data to support GMT"
    echo "                          default=false"
    echo "                          (example: true)"
    echo "  GEOS_VERSION string     version of GEOS to install"
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
GEOS_VERSION=${2:-"3.11.2"}

# ==============================================================================
# Download and install GEOS
# ==============================================================================
# Download GEOS source
curl -o "/opt/geos-${GEOS_VERSION}.tar.bz2" -L \
  "https://download.osgeo.org/geos/geos-${GEOS_VERSION}.tar.bz2"; 
cd /opt && tar xvfj "/opt/geos-${GEOS_VERSION}.tar.bz2";

# Build the gmt source code
cd "/opt/geos-${GEOS_VERSION}" \
    && mkdir -p build \
    && cd build \
    && cmake .. \
    && make \
    && make install;

# cleanup
if [ "$CLEANUP" = true ] ; then
    echo "Cleaning up source code for GEOS.";
    rm -r "/opt/geos-${GEOS_VERSION}*";
fi
#! /bin/bash -ex

CLEANUP=${1:-false}
GEOS_VERSION=${2:-"3.11.2"}

# Download GEOS source
curl -o "/opt/geos-${GEOS_VERSION}.tar.bz2" -L \
  "https://download.osgeo.org/geos/geos-${GEOS_VERSION}.tar.bz2"; 
cd /opt && tar xvfj "/opt/geos-${GEOS_VERSION}.tar.bz2";

# Build the gmt source code
cd /opt/geos-${GEOS_VERSION} \
    && mkdir -p build \
    && cd build \
    && cmake .. \
    && make \
    && make install;

# cleanup
if [ "$CLEANUP" = true ] ; then
    echo 'Cleaning up source code for geos.';
    rm -r /opt/geos-${GEOS_VERSION}*;
fi
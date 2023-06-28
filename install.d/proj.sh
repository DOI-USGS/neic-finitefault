#! /bin/bash -ex

CLEANUP=${1:-false}
PROJ_VERSION=${2:-"9.2.0"}

# Download proj source
curl -o "/opt/proj-${PROJ_VERSION}.tar.gz" -L \
  "https://download.osgeo.org/proj/proj-${PROJ_VERSION}.tar.gz"; 
cd /opt && tar -xvf "/opt/proj-${PROJ_VERSION}.tar.gz";

# Build the proj source code
cd /opt/proj-${PROJ_VERSION}\
    && mkdir -p build \
    && cd build \
    && cmake .. \
    && cmake --build . \
    && cmake --build . --target install;

# cleanup
if [ "$CLEANUP" = true ] ; then
    echo 'Cleaning up source code for proj.';
    rm -r /opt/proj-${PROJ_VERSION}*;
fi
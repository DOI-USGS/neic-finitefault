#! /bin/bash -ex

CLEANUP=${1:-false}
DCW_VERSION=${2:-"2.1.1"}
GMT_VERSION=${3:-"6.4.0"}
GSHHG_VERSION=${4:-"2.3.7"}

# Download GMT source and supporting data
GMT_URL="https://github.com/GenericMappingTools"
curl -o "/opt/dcw-gmt-${DCW_VERSION}.tar.gz" -L \
  "${GMT_URL}/dcw-gmt/releases/download/${DCW_VERSION}/dcw-gmt-${DCW_VERSION}.tar.gz";
cd /opt && tar -xvf "/opt/dcw-gmt-${DCW_VERSION}.tar.gz";
curl -o "/opt/gmt-${GMT_VERSION}-src.tar.gz" -L \
  "${GMT_URL}/gmt/releases/download/${GMT_VERSION}/gmt-${GMT_VERSION}-src.tar.gz";
cd /opt && tar -xvf "/opt/gmt-${GMT_VERSION}-src.tar.gz";
curl -o /opt/gshhg-gmt-${GSHHG_VERSION}.tar.gz -L \
  "${GMT_URL}/gshhg-gmt/releases/download/${GSHHG_VERSION}/gshhg-gmt-${GSHHG_VERSION}.tar.gz";
cd /opt && tar -xvf "/opt/gshhg-gmt-${GSHHG_VERSION}.tar.gz";
mv "/opt/gshhg-gmt-${GSHHG_VERSION}" "/opt/gmt-${GMT_VERSION}/share/gshhg-gmt";
mv "/opt/dcw-gmt-${DCW_VERSION}" "/opt/gmt-${GMT_VERSION}/share/dcw-gmt";

# Build the gmt source code
cd /opt/gmt-${GMT_VERSION} \
  && mkdir -p build \
  && cd build \
  && cmake .. \
  && cmake --build . \
  && cmake --build . --target install;

# cleanup
if [ "$CLEANUP" = true ] ; then
    echo 'Cleaning up source code for gmt.';
    rm -r /opt/dcw-gmt-${DCW_VERSION}*;
    rm -r /opt/gshhg-gmt-${GSHHG_VERSION}*;
    rm -r /opt/gmt-${GMT_VERSION}-src*;
fi
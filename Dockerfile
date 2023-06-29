ARG FROM_IMAGE=code.usgs.gov:5001/devops/images/usgs/python:3.10

# ubuntu packages
FROM ${FROM_IMAGE} as packages

USER root

ARG PYTHON="3.10"
COPY install.d/ubuntu_packages.sh ./install.d/
## packages
RUN apt update -y && apt upgrade -y
RUN bash install.d/ubuntu_packages.sh $PYTHON

USER usgs-user

# base packages
FROM packages as dependencies
ARG DCW_VERSION
ARG GMT_VERSION
ARG GSHHG_VERSION
ARG GEOS_VERSION
ARG PROJ_VERSION

USER root

## gmt
COPY install.d/gmt.sh ./install.d/
RUN bash install.d/gmt.sh true "${DCW_VERSION}" "${GMT_VERSION}" "${GSHHG_VERSION}"
ENV PATH "/usr/local/bin/gmt:${PATH}"
## geos
COPY install.d/libgeos.sh ./install.d/
RUN bash install.d/libgeos.sh true "${GEOS_VERSION}"
ENV PATH "/usr/local/bin/geosop:${PATH}"
ENV GEOS_INCLUDE_PATH "/usr/local/include"
ENV GEOS_LIBRARY_PATH "/usr/local/lib"
## proj
COPY install.d/proj.sh ./install.d/
RUN bash install.d/proj.sh true "${PROJ_VERSION}"
## environment variables
ENV PATH "/usr/local/bin/gmt:${PATH}"
ENV GEOS_INCLUDE_PATH "/usr/local/include"
ENV GEOS_LIBRARY_PATH "/usr/local/lib"
ENV PATH "/usr/local/bin/proj:${PATH}"
## cleanup scripts
RUN rm -rf install.d

USER usgs-user



# Add and install/compile code
FROM dependencies as finitefault
ARG FD_BANK
ARG LITHO1

USER root

## copy code
RUN mkdir /home/usgs-user/finitefault
COPY . /home/usgs-user/finitefault/
ENV FINITEFAULT_DIR /home/usgs-user/finitefault

## compile fortran code
RUN cd $FINITEFAULT_DIR \
    && bash install.d/wasp.sh \
    "${FINITEFAULT_DIR}" \
    --fd-bank "${FINITEFAULT_DIR}/${FD_BANK:-download}" \
    --lith "${FINITEFAULT_DIR}/${LITHO1:-download}"
## install python dependencies
RUN cd $FINITEFAULT_DIR && poetry build 
RUN pip install $FINITEFAULT_DIR/dist/neic_finitefault*.whl
RUN pip install okada-wrapper==18.12.07.3
## update permissions
RUN chown -R usgs-user:usgs-user "${FINITEFAULT_DIR}"
RUN chmod -R 777 "${FINITEFAULT_DIR}"
## cleanup
RUN rm -rf "${FINITEFAULT_DIR}/install.d"
RUN if [ "${FD_BANK}" != "fortran_code/gfs_nm/long/fd_bank" ]; then rm -f "${FINITEFAULT_DIR}/${FD_BANK}"; fi
RUN if [ "${LITHO1}" != "fortran_code/info/LITHO1.0.nc" ]; then rm -f "${FINITEFAULT_DIR}/${LITHO1}"; fi
RUN apt remove -y \
    cmake \
    curl \
    gcc \
    gfortran \
    git;
USER usgs-user

# test code
FROM finitefault as test

USER usgs-user

# (placeholder for unit tests, currently just try importing)
RUN python -c "from wasp import *"
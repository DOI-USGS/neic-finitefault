ARG PYTHON_VERSION
ARG PYTHON_VERSION=${PYTHON_VERSION:-3.10}
ARG FROM_IMAGE=code.usgs.gov:5001/devops/images/usgs/python:${PYTHON_VERSION}

# ubuntu packages
FROM ${FROM_IMAGE} as packages

USER root

COPY install.d/ubuntu_packages.sh ./install.d/
RUN apt update -y && apt upgrade -y
RUN bash install.d/ubuntu_packages.sh ${PYTHON_VERSION}

USER usgs-user

# base packages
FROM packages as dependencies
ARG DCW_VERSION
ARG GMT_VERSION
ARG GSHHG_VERSION
ARG GEOS_VERSION
ARG PROJ_VERSION
ARG FD_BANK
ARG LITHO1
ARG FD_BANK="${FD_BANK:-download}"
ARG LITHO1="${LITHO1:-download}"

USER root

COPY . /home/usgs-user/finitefault/
## GMT
COPY install.d/* ./install.d/
RUN bash install.d/gmt.sh true "${DCW_VERSION}" "${GMT_VERSION}" "${GSHHG_VERSION}"
ENV PATH "/usr/local/bin/gmt:${PATH}"
## GEOS
RUN bash install.d/libgeos.sh true "${GEOS_VERSION}"
ENV PATH "/usr/local/bin/geosop:${PATH}"
ENV GEOS_INCLUDE_PATH "/usr/local/include"
ENV GEOS_LIBRARY_PATH "/usr/local/lib"
## PROJ
RUN bash install.d/proj.sh true "${PROJ_VERSION}"
## environment variables
ENV PATH "/usr/local/bin/gmt:${PATH}"
ENV GEOS_INCLUDE_PATH "/usr/local/include"
ENV GEOS_LIBRARY_PATH "/usr/local/lib"
ENV PATH "/usr/local/bin/proj:${PATH}"
# get data dependencies
ENV FINITEFAULT_DIR /home/usgs-user/finitefault
RUN bash ./install.d/data_dependencies.sh \
    "${FINITEFAULT_DIR}" \
    --fd-bank "${FINITEFAULT_DIR}/${FD_BANK}" \
    --lith "${FINITEFAULT_DIR}/${LITHO1}"
## cleanup scripts
RUN rm -rf install.d fortran_code src poetry.lock pyproject.toml
RUN if [ "${FD_BANK}" != "fortran_code/gfs_nm/long/fd_bank" ]; then rm -f "${FINITEFAULT_DIR}/${FD_BANK}"; fi
RUN if [ "${LITHO1}" != "fortran_code/info/LITHO1.0.nc" ]; then rm -f "${FINITEFAULT_DIR}/${LITHO1}"; fi

USER usgs-user

# Add and install/compile code
FROM dependencies as finitefault

USER root

## copy code
RUN mkdir /home/usgs-user/finitefault
COPY . /home/usgs-user/finitefault/
ENV FINITEFAULT_DIR /home/usgs-user/finitefault

## compile fortran code
RUN cd $FINITEFAULT_DIR \
    && bash install.d/wasp.sh \
    "${FINITEFAULT_DIR}"
## install python dependencies
RUN cd $FINITEFAULT_DIR && poetry build 
RUN pip install $FINITEFAULT_DIR/dist/neic_finitefault*.whl
RUN pip install okada-wrapper==18.12.07.3
## update permissions
RUN chown -R usgs-user:usgs-user "${FINITEFAULT_DIR}"
RUN chmod -R 777 "${FINITEFAULT_DIR}"
## cleanup
RUN rm -rf "${FINITEFAULT_DIR}/install.d"
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

ARG PYTHON_VERSION
ARG PYTHON_VERSION=${PYTHON_VERSION:-3.10}
ARG FROM_IMAGE=code.usgs.gov:5001/devops/images/usgs/python:"${PYTHON_VERSION}"-pygmt

# ubuntu packages
FROM ${FROM_IMAGE} as packages

USER root

COPY install.d/ubuntu_packages.sh ./install.d/
RUN apt update -y && apt upgrade -y
RUN bash install.d/ubuntu_packages.sh

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
WORKDIR /home/usgs-user/finitefault/
## PROJ
RUN bash install.d/proj.sh true "${PROJ_VERSION}"
ENV PATH "/usr/local/bin/proj:${PATH}"
# get data dependencies
ENV FINITEFAULT_DIR /home/usgs-user/finitefault
RUN bash ./install.d/data_dependencies.sh \
    "${FINITEFAULT_DIR}" \
    --fd-bank "${FINITEFAULT_DIR}/${FD_BANK}" \
    --lith "${FINITEFAULT_DIR}/${LITHO1}"
## cleanup files
RUN rm -rf install.d src poetry.lock pyproject.toml pb2002_boundaries.gmt
RUN if [ "${FD_BANK}" != "fortran_code/gfs_nm/long/fd_bank" ]; then rm -f "${FINITEFAULT_DIR}/${FD_BANK}"; fi
RUN if [ "${LITHO1}" != "fortran_code/info/LITHO1.0.nc" ]; then rm -f "${FINITEFAULT_DIR}/${LITHO1}"; fi
## set permissions
RUN chown usgs-user:usgs-user "fortran_code/gfs_nm/long/fd_bank"
RUN chmod 777 "fortran_code/gfs_nm/long/fd_bank"
RUN chown usgs-user:usgs-user "fortran_code/info/LITHO1.0.nc"
RUN chmod 777 "fortran_code/info/LITHO1.0.nc"
RUN chown -R usgs-user:usgs-user "tectonicplates"
RUN chmod -R 777 "tectonicplates"
USER usgs-user

# Add and install and compile code
FROM dependencies as finitefault

USER root

## copy code
COPY . /home/usgs-user/finitefault/
ENV FINITEFAULT_DIR /home/usgs-user/finitefault
WORKDIR /home/usgs-user/finitefault/

## compile fortran code
RUN cd $FINITEFAULT_DIR \
    && bash install.d/wasp.sh \
    "${FINITEFAULT_DIR}"
## install python packages
RUN cd $FINITEFAULT_DIR && poetry build 
RUN pip install $FINITEFAULT_DIR/dist/neic_finitefault*.whl
RUN pip install okada-wrapper==18.12.07.3
## update permissions
RUN chown -R usgs-user:usgs-user "${FINITEFAULT_DIR}"
RUN chmod -R 777 "${FINITEFAULT_DIR}"
## cleanup files and packages
RUN rm -rf "${FINITEFAULT_DIR}/install.d" src pyproject.toml poetry.lock 
RUN apt remove -y \
    cmake \
    curl \
    git
USER usgs-user

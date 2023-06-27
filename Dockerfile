ARG FROM_IMAGE=code.usgs.gov:5001/devops/images/usgs/python:3.10

# ubuntu packages
FROM ${FROM_IMAGE} as packages

USER root

ENV PYTHON "3.10"
COPY install.d/1_ubuntu_packages.sh ./install.d/
## packages
RUN bash install.d/1_ubuntu_packages.sh

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
COPY install.d/2_gmt.sh ./install.d/
RUN bash install.d/2_gmt.sh true
ENV PATH "/usr/local/bin/gmt:${PATH}"
## geos
COPY install.d/3_libgeos.sh ./install.d/
RUN bash install.d/3_libgeos.sh true
ENV PATH "/usr/local/bin/geosop:${PATH}"
ENV GEOS_INCLUDE_PATH "/usr/local/include"
ENV GEOS_LIBRARY_PATH "/usr/local/lib"
## proj
COPY install.d/4_proj.sh ./install.d/
RUN bash install.d/4_proj.sh true
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

USER root

## copy code
RUN mkdir /home/usgs-user/finitefault
COPY . /home/usgs-user/finitefault
ENV FINITEFAULT_DIR /home/usgs-user/finitefault

## compile fortran code
RUN apt install -y git
RUN cd $FINITEFAULT_DIR && bash install.d/5_wasp.sh
RUN apt remove -y git
## install python dependencies
RUN cd $FINITEFAULT_DIR && poetry build
RUN pip install $FINITEFAULT_DIR/dist/neic_finitefault*.whl
## update permissions
RUN chown -R usgs-user:usgs-user /home/usgs-user/finitefault
RUN chmod -R 777 /home/usgs-user/finitefault
## cleanup
RUN rm -rf /home/usgs-user/finitefault/install.d

USER usgs-user

# test code
FROM finitefault as test

USER usgs-user

# (placeholder for unit tests, currently just try importing)
RUN python -c "from wasp import *"
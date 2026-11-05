ARG FROM_IMAGE=code.usgs.gov:5001/devops/images/usgs/ubuntu

### ubuntu packages
FROM ${FROM_IMAGE} AS packages
USER root

## Add deprecation warning
RUN echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
RUN echo "The Dockerfile was a temporary stop gap until Mac compatible Fortran compiling was achieved."
RUN echo "As of !214, the Fortran code is now Mac compatible."
RUN echo "The Dockerfile is deprecated and no longer supported as of version 2.0.0 and will be removed after version 2.1.0."
RUN echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
RUN sleep 5;

RUN mkdir /neic-finitefault
WORKDIR /neic-finitefault
COPY install.d /neic-finitefault/install.d
RUN apt install -y \
    cmake \
    curl \
    gcc \
    git \
    gfortran \
    && apt clean;

### ffm
FROM packages AS ffm-python

COPY . /neic-finitefault/
RUN chmod -R 777 /neic-finitefault
RUN chown usgs-user /neic-finitefault

USER usgs-user
WORKDIR /neic-finitefault
SHELL ["/bin/bash", "-i", "-c"]
RUN bash install.sh /neic-finitefault
RUN echo ". ${HOME}/miniforge/etc/profile.d/conda.sh" >> "${HOME}/.bashrc"
RUN echo "conda activate ff-env" >> "${HOME}/.bashrc"
RUN echo "source ${HOME}/.bashrc" >> "${HOME}/.bash_profile"
ENV PATH="${HOME}/miniforge/bin:$PATH"
ENV PATH="${HOME}/miniforge/envs/ff-env/bin:$PATH"
RUN bash install.d/clawpack.sh /neic-finitefault --yes
RUN conda activate ff-env && pip install -e .



FROM ffm-python AS ffm-test
ARG CI_REGISTRY
ARG RUN_END_TO_END=False
ENV CI_REGISTRY=$CI_REGISTRY
ENV RUN_END_TO_END=$RUN_END_TO_END

RUN conda activate ff-env && ffm get-data topo -74.5 -70 -33.5 -20 -d /neic-finitefault/clawpack/ -f netcdf
ENV GLOBAL_TOPOGRAPHY=/neic-finitefault/clawpack/topo.nc
CMD ["/bin/bash", "-i", "-c", "poe test"]

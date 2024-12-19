ARG FROM_IMAGE=code.usgs.gov:5001/devops/images/usgs/ubuntu

### ubuntu packages
FROM ${FROM_IMAGE} as packages
USER root

RUN mkdir /neic-finitefault
WORKDIR /neic-finitefault
COPY install.d /neic-finitefault/install.d
RUN apt install -y \
    cmake \
    curl \
    gcc \
    gfortran \
    && apt clean;

### miniforge
FROM packages as wasp-python

WORKDIR /root
SHELL ["/bin/bash", "-i", "-c"]
RUN bash /neic-finitefault/install.d/miniforge.sh
ENV PATH=/root/miniforge/bin:$PATH
RUN bash /neic-finitefault/install.d/ff-env.sh /neic-finitefault
RUN echo "conda activate ff-env" >> /root/.bashrc
COPY pyproject.toml /neic-finitefault/pyproject.toml

### wasp-fortran
FROM wasp-python as wasp-fortran

WORKDIR /neic-finitefault
COPY fortran_code /neic-finitefault/fortran_code
RUN install.d/wasp.sh /neic-finitefault
RUN chmod 777 -R /neic-finitefault

### data dependencies
FROM wasp-fortran as wasp-dependencies

RUN apt install -y git
RUN bash /neic-finitefault/install.d/data_dependencies.sh /neic-finitefault
COPY pb2002_boundaries.gmt /neic-finitefault/pb2002_boundaries.gmt
# cleanup
RUN apt remove -y git
RUN rm -rf /neic-finitefault/install.d

### add python source code and install
FROM wasp-dependencies as wasp
ARG RUN_ALL
ARG CI_REGISTRY
ENV RUN_ALL "${RUN_ALL:-True}"
ENV CI_REGISTRY "${CI_REGISTRY}"


COPY src /neic-finitefault/src
RUN pip install -e .
ENV PATH /root/miniforge/envs/ff-env/bin:$PATH

## by default the command is to run the test
CMD ["conda", "run", "--no-capture-output", "-n", "ff-env", "poe", "end_to_end_test"]
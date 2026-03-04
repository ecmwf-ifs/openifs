# (C) Copyright 2011- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Tested with 11.2 and 12.2, both work with OpenIFS 48r1
# 11.2 is closer to ATOS default gnu (at time of writing)
FROM docker.io/library/gcc:13.2.0-bookworm

## Force container shell to be bash
SHELL ["/bin/bash", "-l", "-c"]

ARG OPENIFS_DIR=
ARG OPENIFS_EXPT_DATA=openifs-expt
ARG OPENIFS_EXPT_URL=
ARG SCM_URL=
ARG OPENIFS_REPO_URL=
ARG OPENIFS_BRANCH=
ARG IFSDATA_URL=
ARG CLIMATE_URL=
ARG CLIMATE_VERSION=
ARG RTABLES_URL=

## Install required compilers and libraries, build OpenMPI, and cleanup in one layer
RUN apt update && \
    apt install -y git cmake python3 python3-ruamel.yaml python3-yaml python3-venv \
        libomp-dev libboost-dev libboost-date-time-dev libboost-filesystem-dev \
        libboost-serialization-dev libboost-program-options-dev netcdf-bin \
        libnetcdf-dev libnetcdff-dev liblapack-dev libeigen3-dev vim emacs \
        wget bc ca-certificates && \
    update-ca-certificates && \
    cd /tmp && \
    wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.10.tar.gz && \
    tar -xvf openmpi-5.0.10.tar.gz && \
    cd openmpi-5.0.10 && \
    ./configure --disable-libxml2 --prefix=/usr/local && \
    make -j4 all && \
    make install && \
    ldconfig && \
    cd / && \
    rm -rf /tmp/openmpi-5.0.10* && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create user and directory structure
RUN groupadd --gid 1000 openifs && \
    useradd --uid 1000 --gid openifs --shell /bin/bash --create-home openifs && \
    mkdir -p /home/openifs/${OPENIFS_DIR} && \
    mkdir -p /home/openifs/${OPENIFS_EXPT_DATA}/ && \
    chown -R openifs:openifs /home/openifs

USER openifs

# Download and extract experiment data in one layer
WORKDIR /home/openifs/${OPENIFS_EXPT_DATA}
RUN wget ${SCM_URL} && \
    tar -xvf $(basename ${SCM_URL}) && \
    rm $(basename ${SCM_URL}) && \
    wget ${OPENIFS_EXPT_URL} && \
    tar -xvf $(basename ${OPENIFS_EXPT_URL}) && \
    rm $(basename ${OPENIFS_EXPT_URL})

# Copy OpenIFS directory
WORKDIR /home/openifs/${OPENIFS_DIR}
COPY --chown=openifs:openifs ${OPENIFS_DIR} /home/openifs/${OPENIFS_DIR}

# Update OIFS_HOME in config file, copy source command to bashrc
# source the config and copy 
RUN sed -i 's|export OIFS_HOME="${HOME}/.*"|export OIFS_HOME="${HOME}/'${OPENIFS_DIR}'"|g' oifs-config.edit_me.sh && \
    echo "source /home/openifs/${OPENIFS_DIR}/oifs-config.edit_me.sh">> /home/openifs/.bashrc && \
    source /home/openifs/${OPENIFS_DIR}/oifs-config.edit_me.sh && \
    cd $OIFS_EXPT/ab7z/2016092500 && \
    cp $OIFS_RUN_SCRIPT/oifs-run .

# Download and extract all OpenIFS data files in one layer
WORKDIR /home/openifs/${OPENIFS_DIR}/openifs-data
RUN mkdir -p ifsdata rtables ${CLIMATE_VERSION} && \
    cd rtables && \
    wget ${RTABLES_URL} && \
    tar -xvzf $(basename ${RTABLES_URL}) && \
    rm $(basename ${RTABLES_URL}) && \
    cd ../ifsdata && \
    wget ${IFSDATA_URL} && \
    tar -xvzf $(basename ${IFSDATA_URL}) && \
    rm $(basename ${IFSDATA_URL}) && \
    cd ../${CLIMATE_VERSION} && \
    wget ${CLIMATE_URL} && \
    tar -xvzf $(basename ${CLIMATE_URL}) && \
    rm $(basename ${CLIMATE_URL})

WORKDIR /home/openifs/${OPENIFS_DIR}
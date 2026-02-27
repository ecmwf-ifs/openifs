# Tested with 11.2 and 12.2, both work with OpenIFS 48r1
# 11.2 is closer to ATOS default gnu (at time of writing)
FROM docker.io/library/gcc:13.2.0-bookworm

## Force container shell to be bash
SHELL ["/bin/bash", "-l", "-c"]

ARG OPENIFS_DIR=
ARG OPENIFS_DATA=openifs-expt
ARG SCM_URL=
ARG OPENIFS_REPO_URL=
ARG OPENIFS_BRANCH=

## Install required compilers and libraries
RUN apt update && \
    apt install -y git && \
    apt install -y cmake && \
    apt install -y python3 python3-ruamel.yaml python3-yaml python3-venv && \
    apt install -y libomp-dev && \
    apt install -y libboost-dev libboost-date-time-dev libboost-filesystem-dev libboost-serialization-dev libboost-program-options-dev&& \
    apt install -y netcdf-bin libnetcdf-dev libnetcdff-dev && \
    apt install -y liblapack-dev && \
    apt install -y libeigen3-dev && \
    # apt install -y --no-install-recommends libopenmpi-dev openmpi-bin && \
    apt install -y vim emacs && \
    apt install -y wget bc && \
    apt install -y ca-certificates && \
    update-ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR "/tmp"
RUN wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.10.tar.gz && \
    tar -xvf openmpi-5.0.10.tar.gz
WORKDIR "/tmp/openmpi-5.0.10"
RUN ./configure --disable-libxml2 --prefix=/usr/local && make -j4 all && \
    make install && \
    ldconfig

# Create a separate user to run as nonroot
RUN groupadd --gid 1000 openifs && \
     useradd --uid 1000 --gid openifs --shell /bin/bash --create-home openifs

USER openifs

# Create working directories
RUN mkdir ~/${OPENIFS_DIR} && \
    mkdir -p ~/${OPENIFS_DATA}/

WORKDIR /home/openifs/${OPENIFS_DATA}
RUN wget ${SCM_URL} && \
    tar -xvf $(basename ${SCM_URL}) && \
    rm $(basename ${SCM_URL})

# Clone OpenIFS repository directly in the container
# WORKDIR /home/openifs
# RUN git clone --depth 1 --branch ${OPENIFS_BRANCH} --single-branch ${OPENIFS_REPO_URL} ${OPENIFS_DIR}

WORKDIR /home/openifs/${OPENIFS_DIR}

COPY --chown=openifs:openifs ${OPENIFS_DIR} /home/openifs/${OPENIFS_DIR} 

# Update OIFS_HOME in config file
RUN sed -i 's|export OIFS_HOME="${HOME}/.*"|export OIFS_HOME="${HOME}/'${OPENIFS_DIR}'"|g' oifs-config.edit_me.sh

RUN echo "source /home/openifs/${OPENIFS_DIR}/oifs-config.edit_me.sh">> /home/openifs/.bashrc
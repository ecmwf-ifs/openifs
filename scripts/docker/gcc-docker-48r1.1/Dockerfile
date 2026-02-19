# Tested with 11.2 and 12.2, both work with OpenIFS 48r1
# 11.2 is closer to ATOS default gnu (at time of writing)
FROM docker.io/library/gcc:11.2.0-bullseye
#FROM docker.io/library/gcc:12.3.0-bullseye

## Force container shell to be bash
SHELL ["/bin/bash", "-l", "-c"]

ARG OPENIFS_DIR=openifs-48r1
ARG OPENIFS_DATA=openifs-expt
ARG SCM=scm_openifs 

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
    apt install -y vim emacs && \
    apt install -y wget bc && \
## Install python3-pip and the dependencies to run 
## plotscm.py. Note, even if the script runs the plots 
## need to be moved back to the host machine to view 
## (could just move all the output instead)
    apt install -y python3-pip && \
    pip3 install numpy && \
    pip3 install netcdf4 && \
    pip3 install matplotlib && \
    pip3 install pandas && \
    pip3 install xarray

WORKDIR "/tmp"
RUN wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.5.tar.gz && \
    tar -xvf openmpi-4.1.5.tar.gz
WORKDIR "/tmp/openmpi-4.1.5"
RUN ./configure --prefix=/usr/local && make -j4 all && \
    make install && \
    ldconfig

# Create a separate user to run as nonroot
RUN groupadd --gid 1000 openifs && \
    useradd --uid 1000 --gid openifs --shell /bin/bash --create-home openifs

USER openifs

# Note that $HOME is /home/openifs user..
# create working directory of $HOME/$OPENIFS_DIR so 
# that the OIFS env (in oifs-config.edit_me.sh) points to 
# the correct directory without having to change the 
# config  
RUN mkdir ~/${OPENIFS_DIR} && \
    mkdir -p ~/${OPENIFS_DATA}/${SCM}

WORKDIR /home/openifs/${OPENIFS_DIR}

COPY --chown=openifs:openifs ${OPENIFS_DIR} /home/openifs/${OPENIFS_DIR} 
COPY --chown=openifs:openifs ${SCM} /home/openifs/${OPENIFS_DATA}/${SCM}

RUN echo "source /home/openifs/${OPENIFS_DIR}/oifs-config.edit_me.sh">> /home/openifs/.bashrc




FROM ubuntu:20.04

# set an encoding to make things work smoothly.
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
# set timezone
ENV TZ=Europe/London

ARG DEBIAN_FRONTEND=noninteractive
RUN apt --yes -qq update \
 && apt --yes -qq upgrade \
 && apt --yes -qq install \
                  bzip2 \
                  cmake \
                  cpio \
                  curl \
                  g++ \
                  gcc \
                  gfortran \
                  git \
                  libblas-dev \
                  libhdf5-dev \
                  libidn11-dev \
                  liblapack-dev \
                  liblua5.3-dev \
                  libopenmpi-dev \
                  libpython3-all-dev \
                  libscalapack-openmpi-dev \
                  lua5.3 \
                  openmpi-bin \
                  pkg-config \
                  python3-dev \
                  python3-pip \
                  python3-tk \
                  uuid-dev \
 && rm -rf /var/lib/apt/lists/*

# run common commands
COPY run_commands /opt/docker/bin/run_commands
RUN /opt/docker/bin/run_commands

ENV PATH $PATH:/root/.local/bin

ENV LD_LIBRARY_PATH /root/.local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

SHELL ["/bin/bash"]

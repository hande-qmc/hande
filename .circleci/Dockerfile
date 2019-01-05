FROM ubuntu:18.04

# Set an encoding to make things work smoothly.
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# Add a timestamp for the build. Also, bust the cache.
ADD https://now.httpbin.org/when/now /opt/docker/etc/timestamp

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get --yes -qq update \
 && apt-get --yes -qq upgrade \
 && apt-get --yes -qq install \
                      bzip2 \
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
                      liblua5.3 \
                      libopenmpi-dev \
                      libpython3-all-dev \
                      libscalapack-openmpi-dev \
                      lua5.3 \
                      openmpi-bin \
                      pkg-config \
                      python-dev \
                      python-pip \
                      python-tk \
                      python3-dev \
                      python3-pip \
                      python3-tk \
                      uuid-dev \
 && rm -rf /var/lib/apt/lists/*

# Run common commands
COPY run_commands /opt/docker/bin/run_commands
RUN /opt/docker/bin/run_commands

ENV PATH $PATH:/root/.local/bin:/opt/cmake/bin

ENV LD_LIBRARY_PATH /root/.local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

SHELL ["/bin/bash"]

FROM continuumio/miniconda3:25.3.1-1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    autoconf \
    automake \
    bison \
    flex \
    git \
    libfreetype6-dev \
    libxml2-dev \
    libtool \
    m4 \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    zlib1g-dev \
    curl \
 && rm -rf /var/cache/apt/*

## Python packages
COPY docker/squidpy/environment.yml /opt/environment.yml
RUN conda env create -f /opt/environment.yml

ENV CONDA_ENV_NAME=squidpy
ENV CONDA_PREFIX=/opt/conda/envs/${CONDA_ENV_NAME}
ENV PATH=${CONDA_PREFIX}/bin:$PATH

COPY bin/common /opt
COPY bin/analysis /opt

WORKDIR /opt
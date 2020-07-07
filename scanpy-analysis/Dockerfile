FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    autoconf \
    automake \
    bison \
    fastqc \
    flex \
    git \
    libfreetype6-dev \
    libigraph0-dev \
    libxml2-dev \
    libtool \
    m4 \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    zlib1g-dev \
 && rm -rf /var/cache/apt/*

## Python packages
COPY requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

COPY data /opt/data
COPY bin /opt

WORKDIR /opt

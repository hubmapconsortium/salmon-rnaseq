
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    openjdk-17-jre-headless \
    python3 \
    python3-pip \
    unzip \
 && rm -rf /var/cache/apt/*

ADD https://downloads.openmicroscopy.org/bio-formats/6.10.0/artifacts/bftools.zip /opt/bftools.zip

COPY docker/analysis/requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

WORKDIR /opt
RUN unzip bftools.zip \
 && rm bftools.zip

COPY bin/ome-tiff/convert_ome_tiff.py /opt
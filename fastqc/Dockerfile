#FROM biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1

FROM ubuntu:18.04

RUN apt-get update

RUN apt-get install -y python3

RUN apt-get install -y python3-pip

RUN apt-get install -y fastqc

WORKDIR /out

COPY fastqc_wrapper.py /opt

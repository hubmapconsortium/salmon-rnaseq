FROM pypy:3.8

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    gcc \
    git \
    make \
    python3 \
    python3-pip \
    wget \
    zlib1g-dev \
 && rm -rf /var/cache/apt/*

## Python packages
COPY docker/trim_reads/requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

WORKDIR /opt

ENV COMMIT=7c04ce7898ad5909bd309c6ba3cd9c3bd0651f0e

RUN wget -O seqtk.tar.gz https://github.com/lh3/seqtk/archive/$COMMIT.tar.gz \
 && tar -xf seqtk.tar.gz \
 && cd seqtk-$COMMIT \
 && make \
 && cp seqtk /opt \
 && cd .. \
 && rm -rf seqtk.tar.gz seqtk-$COMMIT

COPY bin/common /opt
COPY bin/trim_reads /opt

FROM mruffalo/salmon-custom:1.3.0-hbm1-5-g6e3ad33829a9

WORKDIR /opt

RUN curl https://giygas.compbio.cs.cmu.edu/gencode.v39.intron-exon.sidx.tar | tar -xf - \
 && curl -O https://giygas.compbio.cs.cmu.edu/gencode.v39.annotation.expanded.tx2gene.tsv

RUN apt-get update \
 && apt-get -y install \
    python3-pip \
 && rm -rf /var/cache/apt/* /var/lib/apt/lists/*

COPY docker/salmon/requirements.txt /opt/requirements.txt
RUN python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

RUN mkdir /out

COPY bin/common /opt
COPY bin/salmon /opt

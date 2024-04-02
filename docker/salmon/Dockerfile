FROM mruffalo/salmon-custom:1.3.0-hbm1-5-g6e3ad33829a9

WORKDIR /opt

RUN curl https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/gencode-v35-salmon-index.tar | tar -xf -
RUN curl https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/v1_index.tar.gz | tar -xzf -
RUN curl https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/v2_index.tar.gz | tar -xzf -
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visiumv1.tx2gene.tsv
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visiumv2.tx2gene.tsv
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v1.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v2.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v3.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v4.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v5.txt


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
COPY data/multiome10x /opt

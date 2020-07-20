FROM mruffalo/salmon-sciseq:1.2.1

WORKDIR /opt

RUN curl https://giygas.compbio.cs.cmu.edu/salmon/grch38_index.tar | tar -xf -
RUN curl -O https://giygas.compbio.cs.cmu.edu/salmon/Homo_sapiens.GRCh38.cdna.all.fa.gz.map

RUN mkdir /out

COPY salmon_wrapper.py /opt
COPY bulk_salmon_wrapper.py /opt

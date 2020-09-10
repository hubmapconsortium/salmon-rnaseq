FROM mruffalo/salmon-snare:1.2.1

WORKDIR /opt

RUN curl http://quark/gencode-v35-salmon-index.tar | tar -xf -

RUN mkdir /out

COPY salmon_wrapper.py /opt
COPY bulk_salmon_wrapper.py /opt

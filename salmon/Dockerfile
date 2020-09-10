FROM mruffalo/salmon-snare:1.2.1

WORKDIR /opt

RUN curl https://giygas.compbio.cs.cmu.edu/salmon/gencode-v35-salmon-index.tar | tar -xf -

RUN apt-get update \
 && apt-get -y install \
    python3-pip \
 && rm -rf /var/cache/apt/* /var/lib/apt/lists/*

COPY requirements.txt /opt/requirements.txt
RUN python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

RUN mkdir /out

COPY salmon_wrapper.py /opt
COPY bulk_salmon_wrapper.py /opt

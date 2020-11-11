FROM pypy:3

## Python packages
COPY docker/barcode_adj/requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

COPY data /opt/data
COPY bin/common /opt
COPY bin/analysis /opt

WORKDIR /opt

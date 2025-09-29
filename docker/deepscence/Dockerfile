FROM continuumio/miniconda3:25.3.1-1

ENV DEBIAN_FRONTEND=noninteractive

## Python packages
COPY docker/deepscence/environment.yml /opt/environment.yml
RUN conda env create -f /opt/environment.yml

ENV CONDA_ENV_NAME=deepscence
ENV CONDA_PREFIX=/opt/conda/envs/${CONDA_ENV_NAME}
ENV PATH=${CONDA_PREFIX}/bin:$PATH

COPY data /opt/data
COPY bin/common /opt
COPY bin/analysis /opt

WORKDIR /opt
RUN chmod -R +x /opt
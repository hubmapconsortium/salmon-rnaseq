<h1 align="center">
  <br>
</h1>

<div class="flex-container" align="center">
    <a href="https://img.shields.io/badge/Python-3.%7C3.8%7C3.9%7C3.10%7C3.11-blue">
    <img src="https://img.shields.io/badge/Python-3.7%7C3.8%7C3.9%7C3.10%7C3.11-blue"
        alt="Python Version">
    <a href="https://github.com/psf/black">
    <img src="https://img.shields.io/badge/code%20style-black-000000.svg"
        alt="Format Version">
    <a href="https://github.com/pre-commit/pre-commit">
    <img src="https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white"
        alt="pre commit">
    </br>
    <a href="https://www.docker.com">
    <img src="https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white"
        alt="Docker">
    <a href="https://www.commonwl.org">
    <img src="https://img.shields.io/badge/cwltool-red?style=for-the-badge&logo=cwltool&logoColor=white"
        alt="cwltool">
</div>

<p align="center" style="color:green">
  <a href="#about">About</a> •
  <a href="#requirements">Requirements</a> •
  <a href="#installation">Installation</a> •
  <a href="#usage">Usage</a> •
  <a href="#license">License</a>
</p>

# About

The HuBMAP scRNA-seq pipeline is built on Salmon, Scanpy, and scVelo, and is
implemented as a CWL workflow wrapping command-line tools encapsulated in
Docker containers.

# Requirements

We require [Docker](https://www.docker.com/) to run the pipeline.

Once Docker is installed run the Docker Daemon. On Linux, this is typically
done by running ``sudo dockerd``. On Mac, this is done by clicking the Docker GUI app. If you are using a M1 or newer silicon Mac, add ``export DOCKER_DEFAULT_PLATFORM=linux/amd64`` to your ``.zshrc`` file to avoid docker warnings and possible errors.

# Installation

Clone the repository and install the requirements. The ``master`` branch and ``latest`` published Docker images may not always
be in sync; checking out a version like ``v2.0.6`` is *highly* recommended
before running the pipeline, unless building Docker images locally.

```bash
git clone https://github.com/tmsincomb/salmon-rnaseq.git
cd salmon-rnaseq
pip install -e .
```

# Usage
```bash
salmon-rnaseq --help
salmon-rnaseq --assay ASSAY --fastq_dir FASTQ_DIR --threads THREADS -o OUTPUT_DIR
```

Supported assays:

* ``10x_v2`` (single-cell)
* ``10x_v2_sn`` (single-nucleus)
* ``10x_v3`` (single-cell)
* ``10x_v3_sn`` (single-nucleus)
* ``snareseq``
* ``sciseq``
* ``slideseq``
  
# License

[![License](https://img.shields.io/github/license/hubmapconsortium/salmon-rnaseq)](https://github.com/hubmapconsortium/salmon-rnaseq/blob/main/LICENSE)

Copyright © HuBMAP

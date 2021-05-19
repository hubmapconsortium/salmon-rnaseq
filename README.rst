.. image:: https://travis-ci.com/hubmapconsortium/salmon-rnaseq.svg?branch=master
    :target: https://travis-ci.com/hubmapconsortium/salmon-rnaseq
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

HuBMAP scRNA-seq pipeline: Salmon, Scanpy, scVelo
=================================================

Overview
--------

The HuBMAP scRNA-seq pipeline is built on Salmon, Scanpy, and scVelo, and is
implemented as a CWL workflow wrapping command-line tools encapsulated in
Docker containers.

Requirements
------------

Running the pipeline requires a CWL workflow execution engine and container
runtime; we recommend Docker and the ``cwltool`` reference implementation.
``cwltool`` is written in Python and can be installed into a sufficiently
recent Python environment with ``pip install cwltool``. Afterward, clone this
repository, check out a tag, and invoke the pipeline as::

  cwltool pipeline.cwl --assay ASSAY --fastq_dir FASTQ_DIR --threads THREADS

(The ``master`` branch and ``latest`` Docker images may not always be in sync;
checking out a version like ``v2.0`` is *highly* recommended.)

Supported assays:

* ``10x_v2`` (single-cell)
* ``10x_v2_sn`` (single-nucleus)
* ``10x_v3`` (single-cell)
* ``10x_v3_sn`` (single-nucleus)
* ``snareseq``
* ``sciseq``
* ``slideseq``

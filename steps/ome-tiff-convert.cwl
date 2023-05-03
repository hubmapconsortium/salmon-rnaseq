#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: v1.0
label: Conversion of tiff file to OME-TIFF
requirements:
  DockerRequirement:
    dockerPull: hubmap/visium-ome-tiff:latest
baseCommand: /opt/convert_ome_tiff.py

inputs:
  img_dir:
    label: "Directory containing tiff file"
    type: Directory
    inputBinding:
      position: 1
outputs:
  ome_tiff_file:
    type: File
    outputBinding:
      glob: "*.ome.tiff"
    label: "Visium image in ome-tiff format"


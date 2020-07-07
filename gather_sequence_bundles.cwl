#!/usr/bin/env cwl-runner
class: CommandLineTool
id: gather_sequence_bundles
label: gather sequence bundles
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: mruffalo/scanpy-snare:latest
  InlineJavascriptRequirement: {}

  InitialWorkDirRequirement:
    listing:
      - $(inputs.sequence_directory)

inputs:
  sequence_directory:
    type: Directory
    inputBinding:
      position: 1
    doc: The directory with sample fastq or fastq.gz files.

outputs:
  fastq_r1:
    type: File[]
    outputBinding:
      glob: "bundles.json"
      loadContents: true
      outputEval: |
                ${
                    var bundle_array_str = self[0].contents;
                    var bundle_array = JSON.parse(bundle_array_str);
                    var file_array = [];
                    for (var i = 0; i < bundle_array.length; i++) {
                        var file = bundle_array[i].fastq_r1;
                        file_array.push(file);
                    }
                    return file_array;
                 }

  fastq_r2:
    type: File[]
    outputBinding:
      glob: "bundles.json"
      loadContents: true
      outputEval: |
                ${
                    var bundle_array_str = self[0].contents;
                    var bundle_array = JSON.parse(bundle_array_str);
                    var file_array = [];
                    for (var i = 0; i < bundle_array.length; i++) {
                        var file = bundle_array[i].fastq_r2;
                        file_array.push(file);
                    }
                    return file_array;
                 }

  dataset_ids:
    type: string[]
    outputBinding:
      glob: "bundles.json"
      loadContents: true
      outputEval: |
                ${
                    var bundle_array_str = self[0].contents;
                    var bundle_array = JSON.parse(bundle_array_str);
                    var label_array = [];
                    for (var i = 0; i < bundle_array.length; i++) {i
                        var file = bundle_array[i].dataset_id;
                        label_array.push(file);
                    }
                    return label_array;
                 }

baseCommand: [/opt/gather_sequence_files.py]

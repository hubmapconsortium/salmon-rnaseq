cwlVersion: v1.0
class: CommandLineTool
label: Correct SNARE-seq barcodes
hints:
  DockerRequirement:
    dockerPull: mruffalo/scanpy-snare:latest
baseCommand: /opt/correct_snareseq_barcodes.py

inputs:
  fastq_r1:
    type: File
    inputBinding:
      position: 1
  fastq_r2:
    type: File
    inputBinding:
      position: 2
outputs:
  barcode_umi_fastq:
    type: File
    outputBinding:
      glob: 'barcode_umi.fastq'
  transcript_fastq:
    type: File
    outputBinding:
      glob: 'transcript.fastq'

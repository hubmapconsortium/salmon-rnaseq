cwlVersion: v1.0
class: CommandLineTool
label: Correct SNARE-seq barcodes
hints:
  DockerRequirement:
    dockerPull: hubmap/scanpy-snare:1.3-snare
baseCommand: /opt/correct_snareseq_barcodes.py

inputs:
  fastq_dir:
    type: Directory[]
    inputBinding:
      position: 1
outputs:
  barcode_umi_fastq:
    type: File
    outputBinding:
      glob: 'barcode_umi.fastq'
  transcript_fastq:
    type: File
    outputBinding:
      glob: 'transcript.fastq'

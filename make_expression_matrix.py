#Takes directory containing a quant-sf file for each sample
#These quant-sf files are effectively TSV files with the following fields
#Name, Length, Effective Length, TPM, and NumReads
#Creates a hdf file containing matrices for TPM and NumReads

import pandas as pd
import sys, os

def get_sample_id(quant_file_name):
    #Get the sample id from the file name
    return quant_file_name[:-8]

def initialize_matrix_dfs(source_df):

    #Get transcript names
    transcript_names = source_df["Name"].tolist()

    #Turn them into column headings for each data frame
    tpm_columns = [transcript_name + "_tpm" for transcript_name in transcript_names]
    num_reads_columns = [transcript_name + "_num_reads" for transcript_name in transcript_names]

    #Add sample id to each set of columns
    tpm_columns.insert(0, "sample_id")
    num_reads_columns.insert(0, "sample_id")

    #Create data frames
    tpm_df = pd.DataFrame(columns=tpm_columns)
    num_reads_df = pd.DataFrame(columns=num_reads_columns)

    #Index by sample_id
    indexed_tpm_df = tpm_df.set_index("sample_id", drop=False)
    indexed_num_reads_df = num_reads_df.set_index("sample_id", drop=False)

    return indexed_tpm_df, indexed_num_reads_df

def find_quant_sf_files(directory: Path) -> Iterable[Path]:
    pattern = '**/*quant.sf'
    for extension in FASTQ_EXTENSIONS:
        yield from directory.glob(pattern.format(extension=extension))

def main(directory: Path):

    quant_sfs = find_quant_sf_file(directory)

    tpm_df = None
    num_reads_df = None

    #For each sample
    for quant_sf in quant_sfs:
        sample_id = get_sample_id(quant_sf)#Get the sample_id
        source_df = pd.read_csv(quant_sf_dir + "/" + quant_sf, delimiter='\t')#And open the source file

        if tpm_df == None:#If we haven't initialized our matrices yet
            tpm_df, num_reads_df = initialize_matrix_dfs(source_df)#Do it now

        #Append the appropriate column as a row to each of the matrices
        tpm_df.append(source_df["TPM"].tolist().insert(0, sample_id))
        num_reads_df.append(source_df["NumReads"].tolist().insert(0, sample_id))

    #Write out to hdf5 file
    tpm_df.to_hdf('expression_matrices.h5', key='tpm', mode='w')
    num_reads_df.to_hdf('expression_matrices.h5', key='num_reads')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.directory)

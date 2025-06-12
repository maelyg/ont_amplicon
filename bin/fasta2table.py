#!/usr/bin/env python
import argparse
import pandas as pd
import os.path
from Bio import SeqIO


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")
    parser.add_argument("--fasta", type=str, required=True, help='provide fasta file')
    parser.add_argument("--sample", type=str, required=True, help='provide sample name')
    parser.add_argument("--tophits", type=str, required=True, help='provide blast top hits')
    args = parser.parse_args()
    fasta_file = args.fasta
    sample_name = args.sample
    blast = args.tophits
    if os.path.getsize(blast) > 0:
        fasta_df = fasta_to_dataframe(fasta_file)
    else:
        fasta_df = pd.DataFrame(columns=["qseqid", "consensus_seq"])

    if os.path.getsize(blast) > 0:
        blastn_results = pd.read_csv(blast, sep="\t", header=0)
        blastn_results.drop(['sample_name'], axis=1, inplace=True)
        merged_df = pd.merge(fasta_df, blastn_results, on = ['qseqid'], how = 'outer')
        merged_df.insert(0, "sample_name", sample_name)
        merged_df['n_read_cont_cluster'] = merged_df['qseqid'].str.split('_').str[2]
        merged_df['n_read_cont_cluster'] = merged_df['n_read_cont_cluster'].str.replace("RC","").astype(int)
        merged_df = merged_df.sort_values(["n_read_cont_cluster"], ascending=[False])

        # merged_df.to_csv(str(sample_name) + "_blastn_top_hits.txt", index=None, sep="\t")
        merged_df.to_csv(os.path.basename(blast).replace("_top_hits_tmp.txt", "_top_hits.txt"), index=None, sep="\t")

    else:
        print("DataFrame is empty!")
        for col in ['sgi', 'sgi', 'sacc', 'length', 'nident', 'pident', 'mismatch', 
                    'gaps', 'gapopen', 'qstart', 'qend', 'qlen', 'sstart', 
                    'send', 'slen', 'sstrand', 'evalue', 'bitscore', 
                    'qcovhsp', 'stitle', 'staxids', 'qseq', 'sseq', 
                    'sseqid', 'qcovs', 'qframe', 'sframe', 'species', 
                    'broad_taxonomic_category', 'FullLineage', 'target_organism_match', 'n_read_cont_cluster']:
            if col not in fasta_df.columns:
                fasta_df[col] = None
                fasta_df.to_csv(os.path.basename(blast).replace("_top_hits_tmp.txt", "_top_hits.txt"), index=None, sep="\t")
# Function to convert FASTA file to DataFrame
def fasta_to_dataframe(fasta_file):
    records = SeqIO.parse(fasta_file, "fasta")
    
    # List to hold the sequence data
    data = []
    
    for record in records:
        # Append ID and sequence to the list
        data.append([record.id, str(record.seq)])
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=["qseqid", "consensus_seq"])
    return df


if __name__ == "__main__":
    main()
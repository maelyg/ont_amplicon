#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from functools import reduce
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections


def parse_args():
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")
    parser.add_argument("--sample", type=str, required=True, help='Provide sample name')
    parser.add_argument("--blastn_results", type=str, required=True)
    parser.add_argument("--nanostat", type=str, required=True)
    parser.add_argument("--bed", type=str, required=True)
    parser.add_argument("--coverage", type=str, required=True)
    parser.add_argument("--target_size", type=str, required=True)
    parser.add_argument("--contig_seqids", type=str, help="Path to mapping file")
    parser.add_argument("--reads_fasta", type=str, help="Reads fasta file")
    parser.add_argument("--consensus", type=str, help="Reads fasta file")
    parser.add_argument("--mapping_quality", type=str, required=True)
    return parser.parse_args()


def read_filtered_read_count(nanostat_path):
    with open(nanostat_path) as f:
        for line in f:
            if "number_of_reads" in line:
                return int(float(line.strip().split("\t")[1]))
    raise ValueError("number_of_reads not found in NanoStat file")


def load_and_prepare_data(blast_path, coverage_path, bed_path, mapping_quality, filtered_read_counts):
    blast_df = pd.read_csv(blast_path, sep="\t", header=0)
    blast_df.rename(columns={"length": "alignment_length"}, inplace=True)

    samtools_cov = pd.read_csv(coverage_path, sep="\t", usecols=["#rname", "endpos", "numreads", "meandepth"], header=0)
    samtools_cov.rename(columns={
        "#rname": "qseqid",
        "endpos": "query_match_length",
        "numreads": "qseq_mapping_read_count",
        "meandepth": "qseq_mean_depth"
    }, inplace=True)
    samtools_cov['qseq_pc_mapping_read'] = samtools_cov['qseq_mapping_read_count'] / filtered_read_counts * 100

    mosdepth = pd.read_csv(bed_path, sep="\t", header=0)
    mosdepth.columns = ["qseqid", "start", "end", "region", "base_counts_at_depth_30X"]
    mosdepth['qseq_pc_depth_30X'] = np.where(
        mosdepth['base_counts_at_depth_30X'] > mosdepth['end'],
        100,
        mosdepth['base_counts_at_depth_30X'] / mosdepth['end'] * 100
    )
    mosdepth['qseq_pc_depth_30X'] = mosdepth['qseq_pc_depth_30X'].round(1)

    mq = pd.read_csv(mapping_quality, sep="\t", header=None)
    mq.columns = ["qseqid", "mean_MQ"]
    return blast_df, samtools_cov, mosdepth[["qseqid", "qseq_pc_depth_30X"]],mq


def merge_dataframes(blast_df, samtools_cov, mosdepth_df, mq_df, df_passes_90_5):
    return reduce(
        lambda left, right: pd.merge(left, right, on="qseqid", how='outer').fillna(0),
        [blast_df, samtools_cov, mosdepth_df, mq_df, df_passes_90_5]
    )


def apply_qc_flags(df, target_size):
    #######30X_DEPTH_FLAG#######
    #Conditions:
    #GREEN: If sgi != 0 and the qseq_pc_depth_30X is >=90
    #ORANGE: If sgi != 0 and the qseq_pc_depth_30X is between 75 and 90.
    #RED: If sgi != 0 and the qseq_pc_depth_30X is <75.
    #GREY: If sgi == 0.

    df['30X_DEPTH_FLAG'] = np.select(
        [
            (df['sacc'] != 0) & (df['qseq_pc_depth_30X'] >= 90),
            (df['sacc'] != 0) & (df['qseq_pc_depth_30X'] >= 75) & (df['qseq_pc_depth_30X'] < 90),
            (df['sacc'] != 0) & (df['qseq_pc_depth_30X'] < 75),
            (df['sacc'].isin([0, None, '', '0', '-']))
        ],
        ['GREEN', 'ORANGE', 'RED', 'GREY'],
        default=""
    )

    #######MAPPED_READ_COUNT_FLAG#######
    #Conditions:
    #GREEN: If sgi != 0 and the qseq_mapping_read_count is >=1000
    #ORANGE: If sgi != 0 and the qseq_mapping_read_count is between 200 and 1000.
    #RED: If sgi != 0 and the qseq_mapping_read_count is <200.
    #GREY: If sgi == 0.
    df['MAPPED_READ_COUNT_FLAG'] = np.select(
        [
            (df['sacc'] != 0) & (df['qseq_mapping_read_count'] >= 1000),
            (df['sacc'] != 0) & (df['qseq_mapping_read_count'] >= 200) & (df['qseq_mapping_read_count'] < 1000),
            (df['sacc'] != 0) & (df['qseq_mapping_read_count'] < 200),
            (df['sacc'].isin([0, None, '', '0', '-']))
        ],
        ['GREEN', 'ORANGE', 'RED', 'GREY'],
        default=""
    )

    # Mean coverage flag
    df['MEAN_COVERAGE_FLAG'] = np.select(
        [
            (df['sacc'] != 0) & (df['qseq_mean_depth'] >= 500),
            (df['sacc'] != 0) & (df['qseq_mean_depth'] >= 100) & (df['qseq_mean_depth'] < 500),
            (df['sacc'] != 0) & (df['qseq_mean_depth'] < 100),
            (df['sacc'].isin([0, None, '', '0', '-']))
        ],
        ['GREEN', 'ORANGE', 'RED', 'GREY'],
        default=""
    )

    #######TARGET_ORGANISM_FLAG#######
    #Conditions:
    #GREEN: If sgi != 0 and the target_organism_match is Y
    #RED: If sgi != 0 and the target_organism_match is N
    #GREY: If sgi == 0.
    df['TARGET_ORGANISM_FLAG'] = np.select(
        [
            (df['sacc'] != 0) & (df['target_organism_match'] == 'Y') & (df['pident'] >= 90),
            (df['sacc'] != 0) & (df['target_organism_match'] == 'Y') & (df['pident'] < 90),
            (df['sacc'] != 0) & (df['target_organism_match'] == 'N'),
            (df['sacc'].isin([0, None, '', '0', '-']))
        ],
        ['GREEN', 'ORANGE', 'RED', 'GREY'],
        default=""
    )

    #######TARGET_SIZE_FLAG#######
    #Conditions:
    #GREEN: If sgi != 0 and the query_match_length is within ±20% of the target_size.
    #ORANGE: If sgi != 0 and the query_match_length is between:
    #    Target size + 20% to 40%, OR
    #    Target size - 20% to -40%.
    #RED: If sgi != 0 and the query_match_length is outside the range of ±40% of the target_size.
    #GREY: If sgi == 0.
    target_size = float(target_size)
    df['TARGET_SIZE_FLAG'] = np.select(
        [
            (df['sacc'] != 0) &
            (df['query_match_length'].between(target_size * 0.8, target_size * 1.2)),
            (df['sacc'] != 0) &
            (df['query_match_length'].between(target_size * 0.6, target_size * 0.8)) |
            (df['query_match_length'].between(target_size * 1.2, target_size * 1.4)),
            (df['sacc'] != 0) &
            ((df['query_match_length'] < target_size * 0.6) |
             (df['query_match_length'] > target_size * 1.4)),
            (df['sacc'].isin([0, None, '', '0', '-']))
        ],
        ['GREEN', 'ORANGE', 'RED', 'GREY'],
        default=""
    )
#     df['PC_READ_LENGTH_FLAG_70_15'] = np.where(
#         (df['sacc'] != 0) & 
# #        (df['TARGET_SIZE_FLAG'] == 'GREEN') &
#         (df['qseq_mapping_read_count'] >= 200) &
#         (df['pc_read_length_passes_70_15'] == True),
#         "GREEN",
#         np.where((df['sacc'] != 0) & 
#                 (df['qseq_mapping_read_count'] >= 200) & 
# #                (df['TARGET_SIZE_FLAG'] == 'GREEN') &
#                 (df['pc_read_length_passes_70_15'] == False),
#                 "RED",
#             np.where((df['sacc'] == 0) |
#                 (df['qseq_mapping_read_count'] < 200), 
#                 "GREY",
#                 ""
#             )
#         )
#     )

#     df['PC_READ_LENGTH_FLAG_80_5'] = np.where(
#         (df['sacc'] != 0) & 
# #        (df['TARGET_SIZE_FLAG'] == 'GREEN') &
#         (df['qseq_mapping_read_count'] >= 200) &
#         (df['pc_read_length_passes_80_5'] == True),
#         "GREEN",
#         np.where((df['sacc'] != 0) & 
#                 (df['qseq_mapping_read_count'] >= 200) & 
# #                (df['TARGET_SIZE_FLAG'] == 'GREEN') &
#                 (df['pc_read_length_passes_80_5'] == False),
#                 "RED",
#             np.where((df['sacc'] == 0) |
#                 (df['qseq_mapping_read_count'] < 200), 
#                 "GREY",
#                 ""
#             )
#         )
#     )
    df['READ_LENGTH_FLAG'] = np.where(
        (df['sacc'] != 0) & 
        (df['num_passing_90'] >= 200),
        "GREEN",
        np.where((df['sacc'] != 0) & 
                (df['num_passing_90'] < 200) & 
                (df['num_passing_90'] >= 50),
                "ORANGE",
            np.where((df['sacc'] != 0) &
                (df['num_passing_90'] < 50), 
                "RED",
                np.where(df['sacc'].isin([0, None, '', '0', '-']),
                    "GREY",
                    ""
                )
            )
        )
    )
    # df['READ_LENGTH_FLAG'] = np.where(
    #     (df['sacc'] != 0) & 
    #     (df['num_passing_90'] >= 250),
    #     "GREEN",
    #     np.where((df['sacc'] != 0) & 
    #             (df['num_passing_90'] < 250) & 
    #             (df['num_passing_90'] >= 100),
    #             "ORANGE",
    #         np.where((df['sacc'] != 0) &
    #             (df['num_passing_90'] < 100), 
    #             "RED",
    #             np.where(df['sacc'].isin([0, None, '', '0', '-']),
    #                 "GREY",
    #                 ""
    #             )
    #         )
    #     )
    # )

    df['MEAN_MQ_FLAG'] = np.where(
        (df['sacc'] != 0) & 
        (df['mean_MQ'] >= 30),
        "GREEN",
        np.where((df['sacc'] != 0) & 
                (df['mean_MQ'] < 30) & 
                (df['mean_MQ'] >= 10),
                "ORANGE",
            np.where((df['sacc'] != 0) &
                (df['mean_MQ'] < 10),
                "RED",
                np.where(df['sacc'].isin([0, None, '', '0', '-']),
                    "GREY",
                    ""
                )
            )
        )
    )

    flag_score_map = {
        'GREEN': 2,
        'ORANGE': 1,
        'RED': 0,
        'GREY': 0
    }
    
    flag_columns = [
        '30X_DEPTH_FLAG',
        'MAPPED_READ_COUNT_FLAG',
        'MEAN_COVERAGE_FLAG',
        'TARGET_ORGANISM_FLAG',
        'TARGET_SIZE_FLAG',
        'READ_LENGTH_FLAG',
        'MEAN_MQ_FLAG'
    ]

    # Convert flag values to scores
    for col in flag_columns:
        df[col + '_SCORE'] = df[col].map(flag_score_map)

    # Total score for each cluster
    df['TOTAL_CONF_SCORE'] = df[[col + '_SCORE' for col in flag_columns]].sum(axis=1)

    # Optionally normalize: score out of 14 (7 flags × max score of 2)
    df['NORMALISED_CONF_SCORE'] = df['TOTAL_CONF_SCORE'] / (2 * len(flag_columns))  # Result: 0 to 1 scale

     #df['qseq_pc_mapping_read'] = df['qseq_pc_mapping_read'].round(1)
    df['qseq_pc_mapping_read'] = df['qseq_pc_mapping_read'].apply(lambda x: float("{:.1f}".format(x)))
    df['qseq_mean_depth'] = df['qseq_mean_depth'].apply(lambda x: float("{:.1f}".format(x)))
    df['NORMALISED_CONF_SCORE'] = df['NORMALISED_CONF_SCORE'].apply(lambda x: float("{:.3f}".format(x)))
    return df


def save_summary(df, sample_name):
    df = df.sort_values(["qseq_pc_mapping_read", "target_organism_match"], ascending=[False, False])
    df.drop("pc_read_length_passes_90_5" , axis=1, inplace=True)
    df.drop("30X_DEPTH_FLAG_SCORE" , axis=1, inplace=True)
    df.drop("MAPPED_READ_COUNT_FLAG_SCORE" , axis=1, inplace=True)
    df.drop("MEAN_COVERAGE_FLAG_SCORE" , axis=1, inplace=True)
    df.drop("TARGET_ORGANISM_FLAG_SCORE" , axis=1, inplace=True)
    df.drop("TARGET_SIZE_FLAG_SCORE" , axis=1, inplace=True)
    df.drop("READ_LENGTH_FLAG_SCORE" , axis=1, inplace=True)
    df.drop("MEAN_MQ_FLAG_SCORE" , axis=1, inplace=True)

    output_file = f"{sample_name}_top_blast_with_cov_stats.txt"
    df.to_csv(output_file, index=False, sep="\t")
    print(f"Saved final summary to {output_file}")
    print(df.head())

def parse_mapping_file(mapping_path):
    """Parse tab-delimited file of read_id and reference."""
    mapping = collections.defaultdict(list)
    with open(mapping_path, 'r') as f:
        for line in f:
            read_id, reference = line.strip().split()
            if reference not in mapping[read_id]:  # Manual deduplication
                mapping[read_id].append(reference)
    return mapping

def get_read_lengths(fasta_path):
    """Return dictionary of read_id: length from fasta file."""
    read_lengths = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        read_lengths[record.id] = len(record.seq)
    return read_lengths

def get_reference_lengths(fasta_path):
    """Return dictionary of reference_id: length from consensus fasta."""
    ref_lengths = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        ref_lengths[record.id] = len(record.seq)
    return ref_lengths

def group_lengths_by_reference(mapping, read_lengths):
    """Return reference: list of lengths."""
    reference_lengths = collections.defaultdict(list)
    for read_id, refs in mapping.items():
        length = read_lengths.get(read_id)
        if length:
            # Only use the first reference assigned (or apply a priority rule)
            ref = refs[0] if isinstance(refs, list) else refs
            reference_lengths[ref].append(length)
    #print(reference_lengths)
    return reference_lengths

def plot_coverage_bar(results_dict, output_path="ref_coverage_bar.png"):
    """
    Bar plot showing % of reads ≥80% reference length for each reference.
    """
    refs = list(results_dict.keys())
    fractions = [v["fraction"] for v in results_dict.values()]

    colors = ['green' if v["passes"] else 'red' for v in results_dict.values()]

    plt.figure(figsize=(12, 6))
    bars = plt.bar(refs, fractions, color=colors)
    plt.axhline(y=0.05, color='blue', linestyle='--', label="5% threshold")

    plt.xticks(rotation=90)
    plt.ylim(0, 1)
    plt.ylabel("Fraction of reads ≥ 80% of reference")
    plt.title("Read coverage across references (5/80 rule)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Bar chart saved as {output_path}")

def analyze_read_lengths_against_reference(reference_lengths, grouped_read_lengths, crl, rpc):
    """
    Analyze read lengths to determine if >=rpc% of reads are >=crl% of the reference length.
    For each reference, check if 5% of reads are ≥80% of reference length
    Parameters:
        read_lengths (dict): read_id -> length of read
        reference_lengths (dict): reference name -> reference length
        crl (int): cutoff percent of reference length
        rpc (int): required percent of reads passing

    Returns:
        df (DataFrame): DataFrame with reference and pass/fail status
    """
    
    results = {}
    for ref, lengths in grouped_read_lengths.items():
        if ref not in reference_lengths:
            print(f"Warning: {ref} not found in consensus fasta.")
            continue

        num_reads = len(lengths)
        ref_len = reference_lengths[ref]
        threshold = (crl / 100) * ref_len

        num_passing = sum(1 for l in lengths if l >= threshold)
        fraction = num_passing / num_reads if num_reads > 0 else 0
        passes = fraction >= (rpc / 100)


        results[ref] = {
            "ref_len": ref_len,
            "num_reads": num_reads,
            "num_passing": num_passing,
            "fraction": fraction,
            "passes": passes
        }

        print(f"{ref}: RefLen={ref_len}, Reads={num_reads}, >={crl}%Ref={num_passing}")
        print(f"Passes {rpc}/{crl}? {'YES' if passes else 'NO'}\n")
    # Save results to file
    #plot_coverage_bar(results)
    df = pd.DataFrame([
    { "qseqid": ref, f"num_passing_{crl}": res["num_passing"], f"pc_read_length_passes_{crl}_{rpc}": res["passes"] }
    for ref, res in results.items()
    ])
    print(df)
    return df

def main():
    args = parse_args()
    filtered_read_counts = read_filtered_read_count(args.nanostat)
    blast_df, samtools_cov, mosdepth_df, mq_df = load_and_prepare_data(
        args.blastn_results,
        args.coverage,
        args.bed,
        args.mapping_quality,
        filtered_read_counts
    )
    mapping = parse_mapping_file(args.contig_seqids)
    read_lengths = get_read_lengths(args.reads_fasta)
    reference_lengths = get_reference_lengths(args.consensus)
    grouped_read_lengths = group_lengths_by_reference(mapping, read_lengths)
    
    #df_passes_70_15 = analyze_read_lengths_against_reference(reference_lengths, grouped_read_lengths, 70, 15)
    df_passes_90_5 = analyze_read_lengths_against_reference(reference_lengths, grouped_read_lengths, 90, 5)

    #df_passes_70_15.to_csv("rpc_read_length_passes_70_15.csv", index=False)
    #df_passes_90_5.to_csv("rpc_read_length_passes_80_5.csv", index=False)
    merged_df = merge_dataframes(blast_df, samtools_cov, mosdepth_df, mq_df, df_passes_90_5)
    flagged_df = apply_qc_flags(merged_df, args.target_size)
    save_summary(flagged_df, args.sample)

if __name__ == "__main__":
    main()

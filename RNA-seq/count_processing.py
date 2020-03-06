#! /bin/python
"""This script will compile a series of RNA-seq HTseq-counts counts from single samples into a multi-sample expression
matrix suitable for use in DeSeq2. The script will label each column in the matrix with the appropriate sample name by
referring to a file containing sample metadata.
Additionally, a single csv containing the count stats for all samples will be produced."""
import os
import pandas as pd
import re


def merge_counts(folder, suffix, regex, meta_location, save_location):
    """Merges all of the htseq-count output files in folder into a single expression matrix."""
    meta = pd.read_csv(meta_location)
    for file in os.listdir(folder):
        if file.endswith(suffix):
            sample_number = re.search(regex, file)[1]
            sample = meta[meta["sample_number"] == sample_number]["sample"].values[0]
            if 'dataset' not in locals():
                dataset = pd.read_csv(os.path.join(folder, file), sep="\t", header=None, names=['gene', sample])
            # print(os.path.join(folder, file))
            else:
                counts = pd.read_csv(os.path.join(folder, file), sep="\t", header=None, names=['gene', sample])
                dataset = pd.merge(dataset, counts, on='gene')
    all_counts = dataset.iloc[:-5]
    stats = dataset.iloc[-5:].rename({'gene': 'count'}, axis='columns')
    all_counts.to_csv(save_location, sep="\t")
    stats.to_csv("data/chr_count_stats.tsv", sep="\t")


# Usage:
merge_counts(folder="rna-seq_bams",
             suffix=".counts",
             regex="HI\.4184\.00.\.Index_\d{1,2}\.(\S{3,4}).counts",
             meta_location="data/rna-seq_meta.csv",
             save_location="data/counts.tsv")

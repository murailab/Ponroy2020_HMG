#!/bin/bash

mkdir trimGalore

find fastq -name '*.fastq.gz*' | cut -d _ -f -2 | uniq | \
parallel -j 12 "trim_galore --output_dir trimGalore --fastqc_args '--outdir trimGalore' \
--gzip --paired {}_R1.fastq.gz {}_R2.fastq.gz"
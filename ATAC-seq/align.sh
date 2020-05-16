#!/bin/bash

# usage: bash align.sh | parallel -j 12

samples=`ls fastq/*.fastq.gz | cut -d / -f 2 | cut -d . -f 1-5 | cut -d _ -f 1,2 | uniq`
genome="reference_data/Ensembl/Homo_sapiens/GRCh37/Sequence/BWAIndex/genome.fa"
ncores="8"

for s in $samples; do

	if [ -f bam/${s}.bam ]; then
		continue
	fi

	echo "bwa mem -t 2 ${genome} trimGalore/${s}_R1_val_1.fq.gz trimGalore/${s}_R2_val_2.fq.gz | samtools sort -@ 2 -O bam -T bam/${s}.tmp -o bam/${s}.bam -"

done
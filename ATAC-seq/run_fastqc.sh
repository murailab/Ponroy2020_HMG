#!/bin/bash

mkdir fastqc

fastqc fastq/*.fastq.gz --threads 12 --outdir fastqc 
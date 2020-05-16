#!/bin/bash

mkdir bigwig

ls rmdup/*.bam | \
parallel -j 12  "bamCoverage --bam {} --outFileName "bigwig/{/.}.bw" --outFileFormat bigwig --normalizeUsing RPKM --binSize 1"
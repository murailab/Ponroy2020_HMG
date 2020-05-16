#!/bin/bash

# usage: bash 05-call_peaks.sh | parallel -j 4

# Alex's command:
# macs2 callpeak -t sample.bam -c input.bam --name sample_vs_input -f BAMPE \
# 	--bw 250 -g mm --mfold 10 30 -p 1e-5

# We'll use the defaults for all other settings for now

#for bam in splitBams/*_openRegions.bam; do

#sample=`echo $bam | cut -d '/' -f 2 | cut -d . -f 1-5`

#macs2 callpeak -t ${bam} -n ${sample} -f BAMPE --outdir 05-peaks

#done


#Todd's version using GNU parallel

ls splitBams/*.bam | parallel -j 6 macs2 callpeak -t {} -n {/.} -f BAMPE --outdir 05-peaks



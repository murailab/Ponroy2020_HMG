#!/bin/bash

computeMatrix reference-point \
  --referencePoint TSS \
  -R output/RNAseq_de_genes.down_lfc_threshold.bed \
  -S \
  bigwig/HI.4713.003.N701---N502.1_GM168_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N701---N502.1_GM168_openRegions.rmdup.bw \
  bigwig/HI.4713.003.N702---N502.2_GM191_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N702---N502.2_GM191_openRegions.rmdup.bw \
  bigwig/HI.4713.003.N703---N502.3_GM260_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N703---N502.3_GM260_openRegions.rmdup.bw \
  bigwig/HI.4713.003.N704---N502.4_ATCC_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N704---N502.4_ATCC_openRegions.rmdup.bw \
  bigwig/HI.4713.003.N705---N502.5_GM04_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N705---N502.5_GM04_openRegions.rmdup.bw \
  bigwig/HI.4713.003.N706---N502.6_AG38_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N706---N502.6_AG38_openRegions.rmdup.bw \
  --upstream 1000 \
  --downstream 0 \
  --binSize 1000 \
  --numberOfProcessors 20 \
  --outFileName output/down_de_genes_matrix.gz \
  --outFileNameMatrix output/down_de_genes_matrix.tab \
  --outFileSortedRegions output/down_de_genes_matrix.bed
  
computeMatrix reference-point \
  --referencePoint TSS \
  -R output/RNAseq_de_genes.up_lfc_threshold.bed \
  -S \
  bigwig/HI.4713.003.N701---N502.1_GM168_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N701---N502.1_GM168_openRegions.rmdup.bw \
  bigwig/HI.4713.003.N702---N502.2_GM191_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N702---N502.2_GM191_openRegions.rmdup.bw \
  bigwig/HI.4713.003.N703---N502.3_GM260_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N703---N502.3_GM260_openRegions.rmdup.bw \
  bigwig/HI.4713.003.N704---N502.4_ATCC_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N704---N502.4_ATCC_openRegions.rmdup.bw \
  bigwig/HI.4713.003.N705---N502.5_GM04_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N705---N502.5_GM04_openRegions.rmdup.bw \
  bigwig/HI.4713.003.N706---N502.6_AG38_openRegions.rmdup.bw \
  bigwig/HI.4713.004.N706---N502.6_AG38_openRegions.rmdup.bw \
  --upstream 1000 \
  --downstream 0 \
  --binSize 1000 \
  --numberOfProcessors 20 \
  --outFileName output/up_de_genes_matrix.gz \
  --outFileNameMatrix output/up_de_genes_matrix.tab \
  --outFileSortedRegions output/up_de_genes_matrix.bed

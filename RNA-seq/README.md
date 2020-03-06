# Pipeline for processing RNA-seq data for Ponroy et. al, 2020

## Sequencng
Sequencing was performed by Genome Quebec ...

## Alignment  
FASTQ files were aligned against Ensembl GRCh37 using STAR_2.5.1b with the following parameters:  
     --runMode alignReads  
     --genomeDir REF_GENOME(Ensemble GRCh37)   
     --outSAMunmapped Within  
     --outSAMtype BAM SortedByCoordinate  
The complete command used to generate each BAM file can be found at the @PG or the @CO flag in the bam header.  
  
## Feature Counting

Read counts at the gene level were generated using htseq with the following command with the Ensembl GRCh37 gene annotation gtf as the reference.
```bash
#! /bin/bash
find rna-seq_bams  -name "*bam" | \ 
parallel -j 24 "htseq-count -f bam -r pos -s reverse -i gene_name {} genes.gtf > {.}.counts" 
```
  
## Gene Expression Matrix for DeSeq2
The resulting counts for each sample were compiled into a single csv
for downstream analysis with DeSeq2 using the python script [count_processing.py](count_processing.py) contained in this repository. 

## Differential Gene Expression Analysis  
Differential gene expression analysis was performed using DeSeq2. The code that produced the figures
of the publication can be found in [RNA-seq_analysis.md](RNA-seq_analysis.md)

# ATAC-seq Pipeline
  
The following series of scripts were used to process the raw reads into the results presented in Ponroy, et. al, 2020. Preparation of DNA for ATAC-seq was performed in accordance with (ref). The raw sequence reads are available at the SRA database. 

## QC and Preprocessing

### Quality Control with FastQC  
Execute run_fastqc.sh. 

```bash
#!/bin/bash

mkdir fastqc

fastqc fastq/*.fastq.gz --threads 12 --outdir fastqc  
```
### Trimming Reads with [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
Execute run_trimGalore.sh.  
```bash
#!/bin/bash

mkdir trimGalore

find fastq -name '*.fastq.gz*' | cut -d _ -f -2 | uniq | \
parallel -j 12 "trim_galore --output_dir trimGalore --fastqc_args '--outdir trimGalore' \
--gzip --paired {}_R1.fastq.gz {}_R2.fastq.gz"
```
## Alignment with bwa  
Execute align.sh. 
```bash
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
```

## Generate Indexes for Bam Files using [Samtools](http://www.htslib.org/)
```bash
find bam -name '*.bam*' | parallel -j 12  samtools index {}
```  
  
## Isolate Open Chromatin Regions  
To generate bam files containing only reads that correspond to regions of free of nucleosomes (open chromatin), the aligned reads were filtered for paired reads that are less than 100bp apart. This was done using R with the code found in [nucFree.md](nucFree.md).


## Peak Calling with [MACS](https://taoliu.github.io/MACS/)
```bash
#!/bin/bash

mkdir peaks

ls nucFreeBams/*.bam | parallel -j 6 macs2 callpeak -t {} -n {/.} -f BAMPE --outdir peaks
```

## Mark Duplicated Reads with Picard
Execute mark_duplicates.sh.  
```bash
#!/bin/bash

mkdir rmdup

ls nucFreeBams/*.bam | \
parallel -j 6 "java -jar ${picard} MarkDuplicates I={} O=rmdup/{/.}.rmdup.bam M=rmdup/{/.}.dups.txt REMOVE_DUPLICATES=true"
find rmdup -name '*.bam*' | parallel -j 12  samtools index {}
```  

## Make Bigwig files with [deepTools](http://gensoft.pasteur.fr/docs/deepTools/2.3.1/index.html)
Execute make_bw.sh.  
```bash
#!/bin/bash

mkdir bigwig

ls rmdup/*.bam | \
parallel -j 12  "bamCoverage --bam {} --outFileName "bigwig/{/.}.bw" --outFileFormat bigwig --normalizeUsing RPKM --binSize 1"
```

## Differential Accessibility Analysis  

[ATAC-seq_analysis.md](ATAC-seq_analysis.md]) contains the analysis of the ATAC-seq data.

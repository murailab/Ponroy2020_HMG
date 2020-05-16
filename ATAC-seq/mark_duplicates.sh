#!/bin/bash

picard=~/picard.jar

mkdir rmdup

ls nucFreeBams/*.bam | \
parallel -j 6 "java -jar ${picard} MarkDuplicates I={} O=rmdup/{/.}.rmdup.bam M=rmdup/{/.}.dups.txt REMOVE_DUPLICATES=true"
find rmdup -name '*.bam*' | parallel -j 12  samtools index {}
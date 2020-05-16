---
title: "Isolate Nucleosome-Free Regions from ATAC-seq Bam Files"
author: "Todd Farmer"
date: "11/25/2018"
output: 
    html_document:
      code_folding: show
      df_print: paged
      keep_md: yes
      number_sections: yes
      theme: flatly
      toc: yes
      toc_float: true
          
---


The primary purpose of the code contained in this document is to extract nucleosome-free regions from ATAC-seq alignments using R. The code is adopted from [Analysis of ATAC-seq data in R and Bioconductor](https://rockefelleruniversity.github.io/RU_ATAC_Workshop.html) by Thomas Carrol.  
  
## Libraries 

```r
library(Rsubread)
library(Rsamtools)
library(ggplot2)
library(magrittr)
library(GenomicAlignments)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(rtracklayer)
library(knitr)
```


## Bam files to be filtered

```r
path <- "bam"
files <- list.files(path, pattern = "bam$")
files <- file.path(path, files)
kable(files, format = "markdown", caption="ATAC-seq bam files")
```



|x                                       |
|:---------------------------------------|
|bam/HI.4713.003.N701---N502.1_GM168.bam |
|bam/HI.4713.003.N702---N502.2_GM191.bam |
|bam/HI.4713.003.N703---N502.3_GM260.bam |
|bam/HI.4713.003.N704---N502.4_ATCC.bam  |
|bam/HI.4713.003.N705---N502.5_GM04.bam  |
|bam/HI.4713.003.N706---N502.6_AG38.bam  |
|bam/HI.4713.004.N701---N502.1_GM168.bam |
|bam/HI.4713.004.N702---N502.2_GM191.bam |
|bam/HI.4713.004.N703---N502.3_GM260.bam |
|bam/HI.4713.004.N704---N502.4_ATCC.bam  |
|bam/HI.4713.004.N705---N502.5_GM04.bam  |
|bam/HI.4713.004.N706---N502.6_AG38.bam  |

## Function to scan a Bam file and save the fragments from nucleosome-free regions.  
  
The nucFreeFilter function produces a bam file containing only fragments with a length of < 100bp that is  saved in the specified folder. If the plot argument is set to "TRUE" then the function produces  
1) a plot of the chromosomal location of the fragments in the bam file,  
2) a plot of the distribution of fragment sizes on a linear y-scale and  
3) a plot of the distribution of fragment sizes on a log y-scale. 
  
Note: this function is memory intensive.

```r
nucFreeFilter <- function(bam, save_dir, plot=TRUE){
  
  ##Read bam files 
  param <- ScanBamParam(mapqFilter = 1,
                        what=c("qname", "flag", "mapq", "isize","seq", "qual", "mrnm"),
                        flag = scanBamFlag(isSecondaryAlignment = FALSE,
                                           isUnmappedQuery=FALSE,
                                           isNotPassingQualityControls = FALSE,
                                           isPaired = TRUE, isProperPair = TRUE)
                                           )
  
  atacReads <- readGAlignmentPairs(bam, 
                             param, 
                             index =paste(bam,".bai",sep=""),
                             use.names=FALSE)

  #Save BAM files of nucleosome-free regions
  atacReads_read1 <- GenomicAlignments::first(atacReads)
  insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
  atacReads_Open <- atacReads[insertSizes < 100 , ]
  openRegionBam <- gsub("\\.bam", "_openRegions\\.bam", basename(bam))
  openRegionBamLoc <- file.path(save_dir, openRegionBam)
  export(atacReads_Open, openRegionBamLoc, format = "bam")
  
  if(plot==TRUE){
    
    #plot of chromosomal read location
    print(idxstatsBam(bam) %>%
          ggplot(aes(seqnames, mapped, fill = seqnames)) + 
          geom_bar(stat = "identity") +
          coord_flip() +
          ggtitle(bam))
  
    
  #plots of distribution of fragment sizes

  fragLenPlot <- table(insertSizes) %>%
                 data.frame %>%
                 rename(InsertSize = insertSizes, Count = Freq) %>%
                 mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count))) %>%
                 ggplot(aes(x = InsertSize, y = Count)) + 
                 geom_line() +
                 ggtitle(bam)
  
  print(fragLenPlot +
          geom_vline(xintercept = c(180, 247), colour = "red") +
          geom_vline(xintercept = c(315, 437), colour = "darkblue") +
          geom_vline(xintercept = c(100), colour = "darkgreen") + 
          theme_bw())

  print(fragLenPlot + scale_y_continuous(trans = "log2") +
          geom_vline(xintercept = c(180, 247), colour = "red") +
          geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
          geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw())
  }
}
```

## Iterating Over All Files Without Plots

```r
for(bam in files){
  nucFreeFilter(bam, "nucFreeBams", plot = FALSE)
} 
```

## Example of Output  with "plot = TRUE"

```r
nucFreeFilter(files[[1]], "nucFreeBams2", plot = TRUE)
```

![](nucleoFreeBAM_figures/example-1.png)<!-- -->![](nucleoFreeBAM_figures/example-2.png)<!-- -->![](nucleoFreeBAM_figures/example-3.png)<!-- -->
  


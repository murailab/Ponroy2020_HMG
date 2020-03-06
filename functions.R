# Helper functions for the RNA-seq analysis

#' Retrieve the most variable genes in a DESeq analysis
#' 
#' @param transformation Object of class DESeqTransform as returned by 
#' varianceStabilizingTransformation or rlog
#' @param n Numeric, specifying the number genes to retrieve. Default: 1000

getHVG <- function(transformation, n = 1000) {
  
  gene_variance <- rowVars(assay(transformation))
  rownames(assay(transformation))[head(order(gene_variance, decreasing = TRUE), n)]
  
}

#' A minimalist ggplot theme
theme_pub <- function(base_size = 11, base_family = "") {
  
  theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black", size = rel(1.2)),
      strip.text.y = element_text(colour = "black", size = rel(1.2)),
      title = element_text(size = rel(0.9)),
      axis.text = element_text(colour = "black", size = rel(0.8)),
      axis.title = element_text(colour = "black", size = rel(1.2)),
      legend.title = element_text(colour = "black", size = rel(0.9)),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", size = rel(0.7)),
      axis.ticks = element_line(colour = "black", size = rel(0.7))
    )
}

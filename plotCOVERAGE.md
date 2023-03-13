## plotCOVERAGE.R

I have used this script to plot multiple mapping files of reads to binned viral contigs viromes.

To run this script requires three `.tsv` files

#### file 1 : postbin_contigcheck_FINAL1.tsv

| sample      | totalNumReads |
| ----------- | ----------- |
| metagenome2 | 12345       |
| metagenome2 | 67891       |

#### file 2 : postbin_contigcheck_FINAL2.tsv

| ctg      | ctgLength |
| ----------- | ----------- |
| contig000001 | 12345       |
| contig000002 | 67891       |

#### file 3 : postbin_contigcheck_FINAL3.tsv
This infomation is acquired by mapping the reads to the contigs of interest (in my case the binned contigs). Thenrunning `samtools depth` on the `BAM` files. 

| sample      | ctg           | locus | depth |
| ----------- | ------------- | ----- | ----- |
| metagenome1 | contig000001  | 1     | 34    |
| metagenome2 | contig000002  | 1     | 456   |

in these example I show the headers, but in the script the headers are added in.

```r
#!/usr/bin/env Rscript

# Script use: 
## Visualising the BAM files mapped to contigs using ggplot
## BAM files must be converted to a tsv numeric file using `samtools depth`
## check the wrangle datasheets header parts for infromation for what should be contained within
## each tsv file in order for this script to work.

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(svglite)
library(reshape)
library(forcats)

theme_set(theme_classic())
set.seed(1234)

################################################
############## IMPORT DATASHEETS ###############
################################################
mapping.df <- read.csv("postbin_contigcheck_FINAL1.tsv", header=FALSE, sep="\t")
samples.df <- read.csv("postbin_contigcheck_FINAL2.tsv", header=FALSE, sep="\t")
contigs.df <- read.csv("postbin_contigcheck_FINAL3.tsv", header=FALSE, sep="\t")

################################################
############## WRANGLE DATASHEETS ##############
################################################
# rename headers
mapping.df <- rename(mapping.df, c(V1="sample", V2="ctg", V3="locus", V4="numReads"))
samples.df <- rename(samples.df, c(V1="sample", V2="totalNumReads", V3="bam_path"))
contigs.df <- rename(contigs.df, c(V1="ctg", V2="ctgLength"))

# extend the mapping.df to contain the ctg legnth
mapping.df$totalNumReads <- samples.df$totalNumReads[match(mapping.df$sample, samples.df$sample)]
mapping.df$ctgLength <- contigs.df$ctgLength[match(mapping.df$ctg, contigs.df$ctg)]

# RPKM caluclation to data
mapping.df <- mapping.df %>% mutate(RPKM=paste0(numReads/((ctgLength/1e3)*(totalNumReads/1e6))))
mapping.df <- mapping.df %>% mutate(minRPKM=paste0(10/((ctgLength/1e3)*(totalNumReads/1e6))))            
mapping.df$RPKM <- as.numeric(mapping.df$RPKM)
mapping.df$minRPKM <- as.numeric(mapping.df$minRPKM)
columns <- c("sample")

################################################
################ RPKM Plotting ################# 
################################################

for (i in 1:nrow(samples.df)) {
  
  # # establish all the variable names for use in the plotting and naming:
  samp <- as.character(samples.df[i, columns[1]])
  
  # # Subset the data for the data of interest and pull out mapping to the two DWV genomes:
  samp.mapping.data = mapping.df %>% subset(sample == samp)

  # # Plot the line-graph of the coverage across the genome:
  p1 <- ggplot() +
    geom_line(data=samp.mapping.data, aes(x=locus, y=RPKM), color="#041562", linewidth = 0.8) +
    geom_line(data=samp.mapping.data, aes(x=locus, y=minRPKM), color="red", linewidth = 0.8, alpha = 0.2) +
    scale_y_log10(limits = c(NA,10000)) +
    scale_x_continuous(limits=c(1,NA), expand = c(0, 0), breaks=c(1,2500,5000,7500,10000)) +
    labs(x="",y="RPKM") + scale_color_manual(values = colors) +
    ggtitle(paste0("sampleID: ",samp))
  
  # If the reads mapped to more than one genome/contig performed a facet_wrap:
  if( length(unique(samp.mapping.data$ctg)) > 1 ){ p2 <- p1 + facet_wrap(vars(ctg)) } else { p2 <- p1 }
  p2
  # generate each figure as its own plot
  assign(paste("plot_", samp, sep=""), p2)
  ggsave(plot = last_plot(),paste0("figures/plot-RPKMs_sample-", samp, format(Sys.time(), "_%Y-%m-%d"), ".png"), dpi=600)
  #ggsave(plot = last_plot(),paste0("figures/plot-RPKM_sample-", samp, format(Sys.time(), "_%Y-%m-%d"), ".svg"), dpi=600)
  
}

################################################
############## numReads Plotting ############### 
################################################

#for (i in 1:nrow(samples.df)){
  
  # # establish all the variable names for use in the plotting and naming:
  samp <- as.character(samples.df[i, columns[1]])
  
  # # Subset the data for the data of interest and pull out mapping to the two DWV genomes:
  samp.mapping.data = mapping.df %>% subset(sample == samp)
  
  # # Plot the line-graph of the coverage across the genome:
  p1 <- ggplot() +
    geom_hline(yintercept = 10, linewidth = 1, colour = "red", alpha = .2,) +
    #annotate("rect", xmin = 1, xmax = 10000, ymin = 0, ymax = 1, alpha = .2, fill = "red") +
    geom_line(data=samp.mapping.data, aes(x=locus, y=numReads), color="#041562", linewidth = 0.8) +
    scale_y_log10(limits = c(NA,10000)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x="",y="RPKM") + scale_color_manual(values = colors) +
    ggtitle(paste0("sampleID: ",samp))
  
  # If the reads mapped to more than one genome/contig performed a facet_wrap:
  if( length(unique(samp.mapping.data$ctg)) > 1 ){ p2 <- p1 + facet_wrap(vars(ctg)) } else { p2 <- p1 }
  p2
  
  # generate each figure as its own plot
  assign(paste("plot_", samp, sep=""), p2)
  ggsave(plot = last_plot(),paste0("figures/plot-numReads_sample-", samp, format(Sys.time(), "_%Y-%m-%d"), ".png"), dpi=600)
  #ggsave(plot = last_plot(),paste0("figures/plot-numReads_sample-", samp, format(Sys.time(), "_%Y-%m-%d"), ".svg"), dpi=600)
  
}

```

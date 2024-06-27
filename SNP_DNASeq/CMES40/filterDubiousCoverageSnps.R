#!/usr/bin/env Rscript

.libPaths(c("/shared/projects/xci/r_packages/","/shared/software/conda/envs/r-4.1.1/lib/R/library/"))
# .libPaths(c("/shared/projects/xci/r_packages/","/shared/software/conda/envs/r-4.0.3/lib/R/library/"))

#Identify the SNPs called with a dubious coverage from DNA-Seq reads. That are SNPs with 0.3<AR>0.7 since AR should be ~0.5

library(limma)
library(tidyverse)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)

#Load SNPs in a tab format generated with VariantsToTable
snp.path <- args[1]
bed.out <- args[2]

snp.tab <- read.table(file = snp.path, header = T, sep = "\t",dec = ".")
colnames(snp.tab) <- c("chr","pos","id","ref","alt","qual","filter","allele.count","allele.freq","depth","MQ","FS","SOR", "MQRankSum","ReadPosRankSum","GT","AD","DP","GQ","PL")
snp.tab$source <- ifelse(snp.tab$id==".","rhESC","dbsnp")
snp.tab$id <- ifelse(snp.tab$id==".",paste(snp.tab$chr, snp.tab$pos, sep = "_"), snp.tab$id)



#Remove scaffolds
chrom <- unique(snp.tab$chr)
chrom <- chrom[grepl(pattern = "chr([1-9]|[1-9][0-9]|[XY])$", x = chrom, perl = T)]
snp.tab <- subset(snp.tab, chr %in% chrom)



#Compute allelic ratio
ar.df <- data.frame(snp.tab$chr, snp.tab$pos, snp.tab$id, read.alt=as.numeric(strsplit2(snp.tab$AD,",")[,2]), dp= as.numeric(snp.tab$DP))
ar.df$ratio <- ar.df$read.alt/ar.df$dp


#Remove the SNPs called with a dubious coverage
ar.filter <- subset(ar.df, !((ratio < 0.3) | (ratio > 0.7)))
#Prepapre the table for the output
bed.filter <- subset(snp.tab, id %in% ar.filter$snp.tab.id, select = c("chr","pos","pos"))
#Bed is a 0-based format
bed.filter$pos <- bed.filter$pos-1

#Convert to GRanges to sort the table
bed.gr <- GRanges(bed.filter$chr, IRanges(bed.filter$pos, bed.filter$pos.1))
bed.gr.sorted <- sort(bed.gr)

#Write the coordinates of the SNPs to retain to bed file
write.table(x = data.frame(seqnames(bed.gr.sorted),
                           as.integer(IRanges::start(bed.gr.sorted)-1), 
                           as.integer(IRanges::end(bed.gr.sorted))),file = bed.out, col.names = F, row.names = F, quote = F, sep = "\t")

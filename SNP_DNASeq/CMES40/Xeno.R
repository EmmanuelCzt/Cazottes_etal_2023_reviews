#!/usr/bin/env Rscript

.libPaths(c("/shared/projects/xci/r_packages/","/shared/software/conda/envs/r-4.1.1/lib/R/library/"))
# .libPaths(c("/shared/projects/xci/r_packages/","/shared/software/conda/envs/r-4.0.3/lib/R/library/"))

args = commandArgs(trailingOnly=TRUE)

rhesus <- args[1]
mouse <- args[2]
sample.name <- args[3]
MM.threshold <- args[4]

parts <- unlist(strsplit(rhesus, .Platform$file.sep))
outfolder.dest <- paste(do.call(file.path, as.list(parts[1:length(parts) - 1])), "XenoFilteR",sample.name, sep = "/")

if (file.exists(paste(outfolder.dest,"Filtered_bams", sep = "/"))){
  unlink(paste(outfolder.dest,"Filtered_bams", sep = "/"), recursive = TRUE)
}else{
  dir.create(paste(outfolder.dest, sep = "/"), recursive = T)
}
 
library(XenofilteR)
library(BiocParallel)

sample.list <- data.frame(rhesus, mouse)

bp.param <- MulticoreParam() #use MulticoreParam so that all workers know about global functions. DO NOT use SnowParam, otherwise biocparallel ERROR: cannot find load functions


XenofilteR(sample.list, destination.folder = outfolder.dest, bp.param = bp.param, output.names = sample.name,MM_threshold = MM.threshold)

#!/usr/bin/env Rscript

options(echo=FALSE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

cat("Loading bambu library.\n")
suppressMessages(library("bambu"))

cat("Loading bam files and annotation files.\n")

bam_files <- args[1]
print(bam_files)
  
fa_file <- args[2]
print(fa_file)

gtf_file <- args[3]
print(gtf_file)

prefix_name <- args[4]
print(prefix_name)

gff_name = paste(prefix_name, "gff", sep=".")
print(gff_name)


bambuAnnotations <- prepareAnnotations(gtf_file)

se.discoveryOnly <- bambu(reads = bam_files, annotations = bambuAnnotations,genome = fa_file, ncore = 1, quant = FALSE)

writeToGTF(se.discoveryOnly, gff_name)

se <- bambu(reads = bam_files, annotations = bambuAnnotations,genome = fa_file, ncore = 1)

writeBambuOutput(se, path = "./bambu/")
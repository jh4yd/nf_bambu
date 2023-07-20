#!/usr/bin/env Rscript

cat("Loading bambu library.\n")
suppressMessages(library(optparse))
suppressMessages(library("bambu"))
suppressMessages(library("Rsamtools"))

# options(echo=FALSE) # if you want see commands in output file
# args <- commandArgs(trailingOnly = TRUE)

my_callback <- function(s4, lflag, val, ps4) {
    return(strsplit(val, ', '))
}

option_list <- list (
    make_option ("--bamfiles", type = 'character', action = 'store', default=NA),
    make_option ("--ref_genome", type = 'character', action = 'store', default=NA),
    make_option ("--ref_annotation", type = 'character', action = 'store'),
    make_option ("--ncore", type = 'integer', action = 'store', default=1)

)

parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser)

# print(arguments$bamfiles)
# print(arguments$ref_genome)
# print(arguments$ref_annotation)
# print(arguments$ncore)

cat("Loading bam files and annotation files.\n")


bam_files <- arguments$bamfiles
print(bam_files)
  
fa_file <- arguments$ref_genome
print(fa_file)

gtf_file <- arguments$ref_annotation
print(gtf_file)

# prefix_name <- args[4]
# print(prefix_name)

ncore <- arguments$ncore
print(ncore)

# gff_name = paste(prefix_name, "gff", sep=".")
# print(gff_name)

cat("processing annotation file")
bambuAnnotations <- prepareAnnotations(gtf_file)
cat("processing bam files")
bam_files <- Rsamtools::BamFileList(strsplit(bam_files, ", ")[[1]])
print(bam_files)
se.multiSample <- bambu(reads = bam_files, annotations = bambuAnnotations,genome = fa_file, ncore = ncore)
writeBambuOutput(se.multiSample, path = "./bambu/")
writeToGTF(rowRanges(se.multiSample), file = tempfile(fileext = ".gff"))


# writeToGTF(se.discoveryOnly, gff_name)

# bam_files <- Rsamtools::BamFileList(strsplit(bam_files, ", ")[[1]])
# print(bam_files)
# se <- bambu(reads = bamFiles_list, annotations = bambuAnnotations,genome = fa_file, ncore = 1)


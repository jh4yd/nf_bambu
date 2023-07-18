#!/usr/bin/env Rscript

suppressMessages(library("DRIMSeq"))
suppressMessages(library("GenomicFeatures"))

cat("Loading counts, conditions and parameters.\n")
cts <- as.matrix(read.csv("merged/all_counts.tsv", sep="\t", row.names="Reference", stringsAsFactors=FALSE))

# Set up sample data frame:
#changed this to sample_id
coldata <- read.csv("de_analysis/coldata.tsv", row.names="sample_id", sep=",", stringsAsFactors=TRUE)

coldata$sample_id <- rownames(coldata)
coldata$condition <- factor(coldata$condition, levels=rev(levels(coldata$condition)))

de_params <- read.csv("de_analysis/de_params.tsv", sep="\t", stringsAsFactors=FALSE)

cat("Loading annotation database.\n")

annotationtype <- de_params$annotation_type[[1]]
txdb <- makeTxDbFromGFF("annotation.gtf",  format = annotationtype)
txdf <- select(txdb, keys(txdb,"GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx<- tab[match(txdf$GENEID, names(tab))]

strip_version<-function(x) {
    tmp<-data.frame(strsplit(x,".", fixed=TRUE), stringsAsFactors=FALSE)
    tmp<-as.vector(tmp[1,])
    colnames(tmp) <- c()
    rownames(tmp) <- c()
    return(tmp)
}

#rownames(cts) <- strip_version(rownames(cts))

cts <- cts[rownames(cts) %in% txdf$TXNAME, ] # FIXME: filter for transcripts which are in the annotation. Why they are not all there? 

# Reorder transcript/gene database to match input counts:
txdf <- txdf[match(rownames(cts), txdf$TXNAME), ]
rownames(txdf) <- NULL

# Create counts data frame:
counts<-data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts)

cat("Filtering counts using DRIMSeq.\n")

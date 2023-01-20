#!/usr/bin/env Rscript
library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
library(NMF)
library(gridExtra)
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
outputdf <- args[2]
sample_name <- args[3]

vcfs <- read_vcfs_as_granges(input, sample_name, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

### cosmic
cosmic <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
sp_url <- paste(cosmic, sep = "")
ref_signatures <- read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order <- match(row.names(mut_mat), ref_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe> 
ref_signatures <- ref_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names>
row.names(ref_signatures) = ref_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
ref_signatures <- as.matrix(ref_signatures[,4:33])
##

ff <- fit_to_signatures(mut_mat, ref_signatures)

cos_sim_ori_rec <- cos_sim_matrix(mut_mat, ff$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
colnames(cos_sim_ori_rec) = "cos_sim"

data <- as.matrix(ff$contribution)
data <- t(data)
data <- data/rowSums(data)

res <- merge(data, cos_sim_ori_rec, by="row.names")
rownames(res) <- res$Row.names
res$Row.names <- NULL

write.table(t(res), file=outputdf, sep="\t", quote=FALSE)
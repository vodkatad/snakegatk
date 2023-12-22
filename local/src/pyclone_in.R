library(vcfR)

vcf_f  <- snakemake@input[['vcf']]
out_f  <- snakemake@output[['outtsv']]
sample <- snakemake@wildcards[['sample']]

vcf <- read.vcfR(vcf_f)

A <- extract.info(vcf, element = "A", as.numeric = TRUE)
B <- extract.info(vcf, element = "B", as.numeric = TRUE)

AD <- extract.gt(vcf, element = "AD", as.numeric = FALSE)
ADs <- strsplit(AD[,sample], split=',', fixed=TRUE)
ref <- sapply(ADs, function(x) {as.numeric(x[[1]][1])})
alt <- sapply(ADs, function(x) {as.numeric(x[[2]][1])})

pyclone_in <- data.frame(mutation_id=rownames(AD), sample_id=sample, ref_counts=ref, alt_counts=alt, major_cn=A,
                         minor_cn=B, normal_cn=2, tumour_content=1, error_rate=0.001)

pyclone_in <- pyclone_in[!is.na(pyclone_in$major_cn) & !is.na(pyclone_in$minor_cn),]
pyclone_in <- pyclone_in[pyclone_in$minor_cn!=0,] # remove LOH
#pyclone_in <- pyclone_in[!grepl('chrY', pyclone_in$mutation_id),]
#pyclone_in <- pyclone_in[!grepl('chrX', pyclone_in$mutation_id),]
pyclone_in[grepl('chrY', pyclone_in$mutation_id), 'normal_cn'] <- 1
#pyclone_in$tot <- pyclone_in$ref_counts+ pyclone_in$alt_counts
#pyclone_in <- pyclone_in[pyclone_in$tot > 50,]
#pyclone_in$tot <- NULL
#try with a thr on VAF
#pyclone_in <- pyclone_in[(pyclone_in$alt / (pyclone_in$ref+pyclone_in$alt)) > 0.05,]
pyclone_in <- pyclone_in[(pyclone_in$minor_cn+pyclone_in$major_cn) < 5,]
write.table(pyclone_in, file=out_f, quote=FALSE, sep="\t")

save.image('beh.Rdata')
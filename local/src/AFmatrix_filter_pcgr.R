#!/opt/R/R-3.6.3/bin/Rscript

# steps for targeted annotation:
#- load merged.table_nomultiallele
#- for each column open mutect/CRC0322LMX0A02201TUMD13000.pcgr_acmg.grch38.snvs_indels.tiers.tsv, keep only TIER 3, put 0 for other muts
#- remove rows with all 0

# produces parallel matrix with genes level info

# TODO get 'mutect' as a parameter from rule to work also for paired mutect
library(getopt)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'merged_af', 'm', 1, 'character',
  'output_gene', 'g', 1, 'character',
  'output_af', 'o', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (!is.null(opt$help) | is.null(opt$merged_af) | is.null(opt$output_af) | is.null(opt$output_gene)) {
  cat(getopt(opts, usage=TRUE))
  stop('-m/-o/-g are mandatory arguments!')
}

af_table_f <- opt$merged_af
af_table_out_f <- opt$output_af
genes_table_out_f <- opt$output_gene

af <- read.table(af_table_f, sep="\t", header=TRUE, row.names = 1)
#genes <- as.data.frame(matrix('', nrow=nrow(af), ncol=1), stringsAsFactors=FALSE)
genes <- data.frame()
samples <- colnames(af)
samples <- gsub('.','-', samples, fixed=TRUE) # so sad, for CRC0282-
#colnames(genes) <- 'gene'
#rownames(genes) <- rownames(af)


load_filter_pcgr <- function(sample) {
  d <- read.table(paste0('mutect/', sample, '.pcgr_acmg.grch38.snvs_indels.tiers.tsv'), sep="\t", header=T, stringsAsFactors = FALSE, quote='')
  d$id <- paste0(d$CHROM, ":", d$POS, ":", d$REF, ":", d$ALT)
  W_TIERS <- c('TIER 1', 'TIER 2', 'TIER 3')
  d <- d[d$TIER %in% W_TIERS,]
  list(ids=paste0('chr', d$id), genes=d$SYMBOL)
}
save.image('pippo.Rdata')
for (i in seq(1, length(samples))) {
  keep <- load_filter_pcgr(samples[i])
  keepid <- intersect(keep[['ids']], rownames(af)) # our af matrix has out of target variants filtered, pcgr do not
  # the order of the first one is kept:
  #> intersect(c(1,2,3), c(3,2))
  #[1] 2 3
  #> intersect(c(3,2), c(1,2,3))
  #[1] 3 2
  igenes <- keep[['genes']]
  igenes <- igenes[keep[['ids']] %in% keepid]
  kg <- data.frame(row.names=keepid, genes=igenes) # is order mantained this way? TODO check for an example where we remove smt
  af[!rownames(af) %in% keepid, i] <- 0
  #genes[rownames(genes) %in% keepid,'gene'] <- igenes
  #tgenes <- merge(genes, kg, by='ids', all.x=TRUE)
  genes <- rbind(genes, kg)
}
genes <- unique(genes)

r <- rowSums(af) == 0
print(table(r))
af <- af[!r,]
genes <- genes[!r,, drop=FALSE]

write.table(af, file=af_table_out_f, sep="\t", quote=FALSE)
write.table(genes, file=genes_table_out_f, sep="\t", quote=FALSE)
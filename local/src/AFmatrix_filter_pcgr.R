#!/opt/R/R-3.6.3/bin/Rscript

# steps for targeted annotation:
#- load merged.table_nomultiallele
#- for each column open mutect/CRC0322LMX0A02201TUMD13000.pcgr_acmg.grch38.snvs_indels.tiers.tsv, keep only TIER 3, put 0 for other muts
#- remove rows with all 0


# TODO get 'mutect' as a parameter from rule to work also for paired mutect
library(getopt)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'merged_af', 'm', 1, 'character',
  'output_af', 'o', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (!is.null(opt$help) | is.null(opt$merged_af) | is.null(opt$output_af)) {
  cat(getopt(opts, usage=TRUE))
  stop('-m/-o are mandatory arguments!')
}

af_table_f <- opt$merged_af
af_table_out_f <- opt$output_af

af <- read.table(af_table_f, sep="\t", header=TRUE, row.names = 1)

samples <- colnames(af)
samples <- gsub('.','-', samples, fixed=TRUE) # so sad, for CRC0282-

load_filter_pcgr <- function(sample) {
  d <- read.table(paste0('mutect/', sample, '.pcgr_acmg.grch38.snvs_indels.tiers.tsv'), sep="\t", header=T, stringsAsFactors = FALSE, quote='')
  d$id <- paste0(d$CHROM, ":", d$POS, ":", d$REF, ":", d$ALT)
  W_TIERS <- c('TIER 1', 'TIER 2', 'TIER 3')
  d <- d[d$TIER %in% W_TIERS,]
  c(paste0('chr', d$id))
}

save.image('pippo.Rdata')
for (i in seq(1, length(samples))) {
  keep <- load_filter_pcgr(samples[i])
  af[rownames(af) %in% keep, i] <- 0
}

r <- rowSums(af) == 0
print(table(r))
af <- af[!r,]

write.table(af, file=af_table_out_f, sep="\t", quote=FALSE)
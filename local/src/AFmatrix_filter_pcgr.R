#!/opt/R/R-3.6.3/bin/Rscript
library(reshape)
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
  'output_long', 'l', 1, 'character',
  'wanted_tiers', 'w', 1, 'character',
  'mutect_dir', 'u', 1, 'character',
  'output_af', 'o', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (!is.null(opt$help) | is.null(opt$merged_af) | is.null(opt$output_af) | is.null(opt$output_gene) | is.null(opt$output_long) |  is.null(opt$wanted_tiers)) {
  cat(getopt(opts, usage=TRUE))
  stop('-m/-g/-l/-w/-o are mandatory arguments!')
}

mutect_dir <- 'mutect' # we keep mutect as the default to avoid having to add variables to already run projects, we'll add a parameter only to the paired new ones
if (!is.null(opt$mutect_dir)) {
  mutect_dir <- opt$mutect_dir
}

W_TIERS <- unlist(strsplit(opt$wanted_tiers, ',', fixed=TRUE))

af_table_f <- opt$merged_af
af_table_out_f <- opt$output_af
genes_table_out_f <- opt$output_gene
long_af_genes <- opt$output_long

af <- read.table(af_table_f, sep="\t", header=TRUE, row.names = 1)
#genes <- as.data.frame(matrix('', nrow=nrow(af), ncol=1), stringsAsFactors=FALSE)
genes <- data.frame(stringsAsFactors=FALSE)
samples <- colnames(af)
samples <- gsub('.','-', samples, fixed=TRUE) # so sad, for CRC0282-
#colnames(genes) <- 'gene'
#rownames(genes) <- rownames(af)


load_filter_pcgr <- function(sample) {
  d <- read.table(paste0(mutect_dir, '/', sample, '.pcgr_acmg.grch38.snvs_indels.tiers.tsv'), sep="\t", header=T, stringsAsFactors = FALSE, quote='')
  d$id <- paste0(d$CHROM, ":", d$POS, ":", d$REF, ":", d$ALT)
  #W_TIERS <- c('TIER 1', 'TIER 2', 'TIER 3')
  d <- d[d$TIER %in% W_TIERS,]
  list(ids=paste0('chr', d$id), genes=d$SYMBOL, cds=d$CDS_CHANGE)
}

keepids <- c()
for (i in seq(1, length(samples))) {
  keep <- load_filter_pcgr(samples[i])
  keepid <- intersect(keep[['ids']], rownames(af)) 
  # our af matrix has out of target variants filtered, pcgr do not ??
  keepids <- c(keepids, keepid)
  # the order of the first one is kept:
  #> intersect(c(1,2,3), c(3,2))
  #[1] 2 3
  #> intersect(c(3,2), c(1,2,3))
  #[1] 3 2
  igenes <- keep[['genes']]
  igenes <- igenes[keep[['ids']] %in% keepid]
  icds <- keep[['cds']]
  icds <- icds[keep[['ids']] %in% keepid]
  kg <- data.frame(symbol=keepid, genes=igenes, cds=icds) # is order mantained this way? TODO check for an example where we remove smt
  
  #af[!rownames(af) %in% keepid, i] <- 0 # no! # YES! we need to remove the variants not in the right tiers! TODO
  
  #genes[rownames(genes) %in% keepid,'gene'] <- igenes
  #tgenes <- merge(genes, kg, by='ids', all.x=TRUE)
  if (i != 1) {
    genes <- rbind(genes, kg, stringsAsFactors=FALSE)
  } else {
    genes <- kg
  }
}
genes <- unique(genes) # what's the difference
#genes <- genes[genes$symbol %in% unique(genes$symbol),, drop=FALSE]
print(length(unique(keepids)))
# We never remove anything cause it was not in af to start with...
af[!rownames(af) %in% keepids,] <- rep(0, ncol(af))
r <- rowSums(af) == 0
print(table(r))
af <- af[!r,]
#genes <- genes[!r,, drop=FALSE]

write.table(af, file=af_table_out_f, sep="\t", quote=FALSE)
write.table(genes, file=genes_table_out_f, sep="\t", quote=FALSE, row.names = FALSE)

print('Checking if filtering of af based on TIERs and muts info based on AF worked')
stopifnot(nrow(genes) == nrow(af))
m <- merge(af, genes, by.x='row.names', by.y='symbol')
stopifnot(nrow(genes) == nrow(m))

cols <- colnames(m)
cols <- cols[grepl('CRC', cols, fixed=TRUE)]
long <- melt(m, id.vars=c('Row.names', 'genes', 'cds'), measure.vars=cols)
colnames(long) <- c('mut_id','gene','cds','lgenealogy','af')
write.table(long, file=long_af_genes, sep="\t", quote=FALSE)

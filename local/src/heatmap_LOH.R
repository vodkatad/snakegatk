### pheatmap for loh

setwd("/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_subs/loh/final/")

df1_f <- read.table("CRC0065PRX0A01201TUMD03000V2_CRC0065LMX0B02205TUMD02000V2_loh_seg_dim_10000000.tsv",
                    quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

df2_f <- read.table("CRC1169PRX0A01001TUMD08000V2_CRC1169LMX0A02001TUMD08000V2_loh_seg_dim_10000000.tsv",
                    quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

df3_f <- read.table("CRC1473PRX0B01001TUMD03000V2_CRC1473LMX0B02002TUMD03000V2_loh_seg_dim_10000000.tsv",
                    quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

df4_f <- read.table("CRC1599PRX0A02002TUMD03000V2_CRC1599LMX0A02001TUMD03000V2_loh_seg_dim_10000000.tsv",
                    quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

df1_f$model <- "CRC0065"
df2_f$model <- "CRC1169"
df3_f$model <- "CRC1473"
df4_f$model <- "CRC1599"

df <- rbind(df1_f, df2_f, df3_f, df4_f)
df$start <- NULL
df$end <- NULL
df$event <- as.factor(df$event)

dfnochr <- df
dfnochr$chr <- NULL

result <- as.data.frame(dfnochr %>%
  pivot_wider(names_from = segment_id, values_from = event))
rownames(result) <- result$model
result$model <- NULL
result[] <- lapply(result, function(x) {
  x <- gsub("LOH_early", 1, x)
  x <- gsub("NO", 0, x)
  x <- gsub("LOH_diff", 2, x)
  as.numeric(x)
})

an_col <- as.data.frame(matrix(ncol = 1, nrow = length(colnames(result))))
rownames(an_col) <- colnames(result)
an_col$V1 <- df1_f$chr
names(an_col)[names(an_col)=="V1"] <- "Chr"
#names(an_col)[names(an_col)=="V1"] <- "Enzyme"
an_col$Chr <- factor(an_col$Chr, 
                               levels = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 
                                          'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 
                                          'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 
                                          'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'))

gappini <- df1_f
gappini2 <- gappini %>% group_by(chr)
gappini2 <- gappini2 %>% summarise(
  segment_id = max(segment_id)
)

righini <- c(1,2,3,4)


pheatmap(result, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = an_col, 
         show_colnames = FALSE, legend_breaks = c(0,1,2), color = colorRampPalette(c("grey", "black", "blue"))(3),
         legend_labels = c("NO", "LOH_early", "LOH_diff"), gaps_col = gappini2$segment_id, gaps_row = righini)


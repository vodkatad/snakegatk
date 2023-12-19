library(tidyverse)

mixed_f <- snakemake@input[["tsv"]]
res <- snakemake@output[["loh"]]

#mixed_f <- "/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_subs/loh/mixed_bed/CRC1599PRX0A02002TUMD03000V2_CRC1599LMX0A02001TUMD03000V2.bed"
mixed <- read.table(mixed_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(mixed) <- c("chr_e", "start_e", "end_e", "A_e","B_e", "chr_l", "start_l", "end_l", "A_l","B_l", "intersection")

mixed$loh <- ifelse((mixed$B_e != 0 & mixed$B_l == 0), "LOH", "NONE")

loh <- mixed %>% filter(loh == "LOH")

loh$start <- pmin(loh$start_e, loh$start_l)
loh$end <- pmax(loh$end_e, loh$end_l)

loh <- loh[,c(1, 11, 13, 14)]
loh <- loh[,c(1, 3, 4, 2)]
names(loh)[names(loh)== "chr_e"] <- "chr"
loh <- loh %>% filter(intersection > 100000)

write.table(loh, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

#se overlappa con due sceglie intersezione con segmento maggiore 
#se overlappo con solo una delle due 
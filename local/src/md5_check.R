links <- "/mnt/cold1/snaketree/prj/snakegatk_real/local/share/data/Pri_Met_pairs/aws_links.tsv"
links <- read.table(links, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colonne <- colnames(links)
colnames(links)[colnames(links) == "R1.MD5"] <- "R1.MD5_old"
colnames(links)[colnames(links) == "R2.MD5"] <- "R2.MD5_old"
colnames(links)[colnames(links) == "Sample"] <- "Sample_R1"
links <- cbind(links, links$Sample_R1)
colnames(links)[colnames(links) == "links$Sample_R1"] <- "Sample_R2"
links$Sample_R1_complete <- NA
links$Sample_R2_complete <- NA

for (i in seq(length(links$Sample_R1))) {
  links[i, "Sample_R1_complete"] <- paste0(links[i, "Sample_R1"], "_SA_L001_R1_001.fastq.gz")
}
for (i in seq(length(links$Sample_R2))) {
  links[i, "Sample_R2_complete"] <- paste0(links[i, "Sample_R2"], "_SA_L001_R2_001.fastq.gz")
}

links_fin_R1 <- as.data.frame(cbind(links$Sample_R1_complete, links$R1.MD5_old))
colnames(links_fin_R1) <- c("genealogy", "md5")
links_fin_R2 <- as.data.frame(cbind(links$Sample_R2_complete, links$R2.MD5_old))
colnames(links_fin_R2) <- c("genealogy", "md5")
links_fin <- as.data.frame(rbind(links_fin_R1, links_fin_R2))

links3 <- "/mnt/cold1/snaketree/prj/snakegatk_real/local/share/data/Pri_Met_pairs/aws_links_3.tsv"
links3 <- read.table(links3, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(links3) <- colonne
colnames(links3)[colnames(links3) == "R1.MD5"] <- "R1.MD5_old"
colnames(links3)[colnames(links3) == "R2.MD5"] <- "R2.MD5_old"
colnames(links3)[colnames(links3) == "Sample"] <- "Sample_R1"
links3 <- cbind(links3, links3$Sample_R1)
colnames(links3)[colnames(links3) == "links3$Sample_R1"] <- "Sample_R2"

links3$Sample_R1_complete <- NA
links3$Sample_R2_complete <- NA

for (i in seq(length(links3$Sample_R1))) {
  links3[i, "Sample_R1_complete"] <- paste0(links3[i, "Sample_R1"], "_SA_L001_R1_001.fastq.gz")
}
for (i in seq(length(links3$Sample_R2))) {
  links3[i, "Sample_R2_complete"] <- paste0(links3[i, "Sample_R2"], "_SA_L001_R2_001.fastq.gz")
}

links3_fin_R1 <- as.data.frame(cbind(links3$Sample_R1_complete, links3$R1.MD5_old))
colnames(links3_fin_R1) <- c("genealogy", "md5")
links3_fin_R2 <- as.data.frame(cbind(links3$Sample_R2_complete, links3$R2.MD5_old))
colnames(links3_fin_R2) <- c("genealogy", "md5")
links3_fin <- as.data.frame(rbind(links3_fin_R1, links3_fin_R2))

finale <- as.data.frame(rbind(links_fin, links3_fin))
names(finale)[names(finale) == "md5"] <- "md5_old"

md5 <- "/mnt/cold1/snaketree/prj/snakegatk_real/local/share/data/Pri_Met_pairs/new_md5_our.txt"
md5 <- read.table(md5, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
md5[c('md5', 'genealogy')] <- str_split_fixed(md5$V1, ' ', 2)
md5$genealogy <- trimws(md5$genealogy, which = c("left"))

md5$V1 <- NULL

merged <- merge(finale, md5, by = "genealogy")

# > setdiff(md5$genealogy, finale$genealogy)
# character(0)
# > setdiff(finale$genealogy, md5$genealogy)
# character(0)

merged$same <- NA

for (i in seq(length(merged$genealogy))) {
  if (merged[i, "md5_old"] == merged[i, "md5"]) {
    merged[i, "same"] <- TRUE
  } else {
    merged[i, "same"] <- FALSE
  }
}

write.table(merged, "md5_check.tsv", quote = FALSE, sep = "\t", col.names = TRUE)

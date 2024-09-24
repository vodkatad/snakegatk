library(stringr)
library(readxl)

md5new <- "/mnt/cold1/snaketree/prj/snakegatk_real/local/share/data/WES_2023/us_defall_md5.txt"
md5new <- read.table(md5new, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

md5new <- as.data.frame(str_split_fixed(md5new$V1, "  ", 2))
md5new <- md5new[,c("V2","V1")]
colnames(md5new) <- c("Sample", "md5_new")
#md5new$ssample <- substr(md5new$Sample, 1, 15)
#md5new$ssample <- gsub("_l", "", md5new$ssample)
#md5new$ssample <- sub("_$", "", md5new$ssample)
md5new$R <- str_sub(md5new$Sample, -14, -14)
#md5new$number <- substr(md5new$ssample, 13, 15)
#md5new$number <- as.numeric(md5new$number)
#md5new <- md5new[order(md5new$number),]
md5R1 <- md5new %>% filter(R == "1")
md5R1$ssample <- substr(md5R1$Sample, 1, 26)
names(md5R1)[colnames(md5R1) == "md5_new"] <- "md5R1_new"
#md5R1 <- md5R1[,c(3,2)]
md5R2 <- md5new %>% filter(R == "2")
md5R2$ssample <- substr(md5R2$Sample, 1, 26)
names(md5R2)[colnames(md5R2) == "md5_new"] <- "md5R2_new"
#md5R2 <- md5R2[,c(3,2)]

md5original <- "/mnt/cold1/snaketree/prj/snakegatk_real/local/share/data/WES_2023/Sequencing_Report_WES_finale.xlsx"
md5original <- read_excel(md5original)
md5original <- as.data.frame(cbind(md5original$`Sample ID`, md5original$R1_MD5_original, md5original$R2_MD5_original))
colnames(md5original) <- c("ssample", "md5R1_ori", "md5R2_ori")

mergedR1 <- merge(md5original, md5R1, by = "ssample")
merged <- merge(mergedR1, md5R2, by = "ssample")
merged$Sample.x <- NULL
merged$R.x <- NULL
merged$Sample.y <- NULL
merged$R.y <- NULL

for (i in seq(merged$ssample)) {
  if (merged[i, "md5R1_ori"] == merged[i, "md5R1_new"]) {
    merged[i, "test_R1"] <- TRUE
  } else {
    merged[i, "test_R1"] <- FALSE
  }
}

for (i in seq(merged$ssample)) {
  if (merged[i, "md5R2_ori"] == merged[i, "md5R2_new"]) {
    merged[i, "test_R2"] <- TRUE
  } else {
    merged[i, "test_R2"] <- FALSE
  }
}

write.table(merged, file = "/mnt/cold1/snaketree/prj/snakegatk_real/local/share/data/WES_2023/md5_final_test.tsv", 
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

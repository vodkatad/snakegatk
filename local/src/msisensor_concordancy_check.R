#!/usr/bin/env Rscript

library(getopt)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'pdo', 'o', 1, 'character',
  'pdx', 'x', 1, 'character',
  'output', 'r', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

pdo <- read.table(opt$pdo, sep='\t', header=FALSE)
names(pdo) <- "sample"
pdo$passage <- as.numeric(substr(pdo$sample,16,17))
pdo$time <- ifelse(pdo$passage==3, "early", "late")
early <- pdo[pdo$time=="early",]
late <- pdo[pdo$time=="late",]

pdx <- read.table(opt$pdx, sep='\t', header=FALSE)
names(pdx) <- "sample"

files_pdo <- list.files(path="msipro", pattern="*.txt$",full.names = TRUE)
files_pdx <- list.files(path="../biobanca_earlylate_xeno/msipro", pattern="*.txt$",full.names = TRUE)

model <- pdo$sample
new <- data.frame(row.names=unique(substr(model,1,7)), stringsAsFactors = FALSE)

for( m in 1:length(model) ) {
  tmp <- grep(model[m],files_pdo,value=TRUE)
  if( model[m] %in% early$sample ) {
    new[substr(model[m],1,7),"pdo_early"] <- model[m]
    res <- read.table(tmp, sep='\t', header=TRUE)
    new[substr(model[m],1,7),"msipro_early"] <- res[,3]
  } else {
    new[substr(model[m],1,7),"pdo_late"] <- model[m]
    res <- read.table(tmp, sep='\t', header=TRUE)
    new[substr(model[m],1,7),"msipro_late"] <- res[,3]
  }
}

model <- pdx$sample

for( m in 1:length(model) ) {
  tmp <- grep(model[m],files_pdx,value=TRUE)
  new[substr(model[m],1,7),"pdx"] <- model[m]
  res <- read.table(tmp, sep='\t', header=TRUE)
  new[substr(model[m],1,7),"msipro_pdx"] <- res[,3]
}

new$msi_pdo_early <- ifelse(new$msipro_early>5,"MSI","MSS")
new$msi_pdo_late <- ifelse(new$msipro_late>5,"MSI","MSS")
new$earlylate_concordancy <- ifelse(new$msi_pdo_early==new$msi_pdo_late, "YES","NO")
new$msi_pdx <- ifelse(new$msipro_pdx>5,"MSI","MSS")
new$pdxo_concordancy <- ifelse(new$msi_pdo_early==new$msi_pdx, "YES","NO")


write.table(new, opt$output, sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)

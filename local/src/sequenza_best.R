cptable_f <- snakemake@input[["cptable"]]
sample <- snakemake@wildcards[['sample']]
log_f <- snakemake@log[['vlog']]

sink(log_f)
sessionInfo()
sink()

load(cptable_f)

#which(CRC0740PRX0A02001TUMD08000V2_sequenza_cp_table$lpp == max(CRC0740PRX0A02001TUMD08000V2_sequenza_cp_table$lpp), arr.ind=T)
#    row col
#3.6  27  91
#> rownames(CRC0740PRX0A02001TUMD08000V2_sequenza_cp_table$lpp)[27]
#[1] "3.6"
#> colnames(CRC0740PRX0A02001TUMD08000V2_sequenza_cp_table$lpp)[91]
#[1] "1"

cp_table_obj <- eval(as.symbol(paste0(sample, '_sequenza_cp_table')))
# da riga 60 qui: https://bitbucket.org/sequenzatools/sequenza/src/master/R/graphics.R
# vedo che la tabella non è ancora scaled, lo fa per il plot

max_lpp <- max(cp_table_obj$lpp)

ind <- which(cp_table_obj$lpp == max_lpp, arr.ind=TRUE)
res <- data.frame(matrix(NA, nrow=nrow(ind), ncol=3))
colnames(res) <- c('purity', 'ploidy', 'LPP')
for (i in seq(1, nrow(ind))) {
    ind_rows <- ind[i,1]
    ind_cols <- ind[i,2]
    res[i,] <- c(colnames(cp_table_obj$lpp)[ind_cols], rownames(cp_table_obj$lpp)[ind_rows], max_lpp)
} 

write.table(res, file=snakemake@output[[1]], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

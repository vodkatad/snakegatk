library(ggplot2, quietly=TRUE)

idata <- snakemake@input[["bed"]]
debug <- snakemake@params[["debug"]]
output <- snakemake@output[["pdf"]]

if (debug == "yes") {
  save.image(file=paste0(output,'.debug','.RData'))
}


data <- read.table(gzfile(idata), sep="\t")
colnames(data) <- c("chr", "begin", "end", "depth")
wanted_chr <- paste0("chr", seq(1,22))
wanted_chr <- c(wanted_chr, "chrY", "chrX")
wdata <- data[data$chr %in% wanted_chr,]
# TODO parametrize
wdata$depth <- factor(wdata$depth, levels= c("0:1","1:5","5:10","10:50","50:100","100:150","150:inf"))
m <- data.frame(depth=c("0:1","1:5","5:10","10:50","50:100","100:150","150:inf"), callability=c("NOCOV","LOW_COVERAGE","CALLABLE","CALLABLE","CALLABLE","CALLABLE","HIGHCOV"))
mdata <- merge(wdata, m, by="depth")
ggplot(mdata, aes(x=depth))+ geom_bar(aes(fill=callability))+theme_bw()
ggsave(output)

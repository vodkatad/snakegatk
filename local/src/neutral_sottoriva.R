loweraf <- as.numeric(snakemake@wildcards[["loweraf"]])
higheraf <- as.numeric(snakemake@wildcards[["higheraf"]])
data <- snakemake@input[["afmatrix"]]
afcolumn <- snakemake@params[["afcolumn"]]
fit <- snakemake@output[["fit"]]
histo <- snakemake@output[["hist"]]
debug <- snakemake@params[["debug"]]
r2 <- snakemake@output[["r2"]]

if (debug == "yes") {
  save.image(file=paste0(fit,'.debug','.RData'))
}

data <- read.table(gzfile(data), sep="\t", header=TRUE)
af <- data[,afcolumn, drop=FALSE]

#https://github.com/andreasottoriva/neutral-tumor-evolution/blob/master/Identification%20of%20neutral%20tumor%20evolution%20across%20cancer%20types%20-%20%20Simulation%20Results.ipynb
exsubcl <- af[af<higheraf & af>loweraf]
exsubcl_nohigh <- af[af >= 0] # < 0.6 for may 2019 used plots (always better and easy reproducibility...)
if (length(exsubcl) != 0) {
    excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]<=exsubcl[1:length(exsubcl)]))
    invf <- 1/exsubcl - 1/higheraf # subtract to use fit without intercept
    maxi <- length(invf)
    labels <-  c(1, floor(maxi/5), floor(maxi/2), floor(maxi/0.5), maxi)
    model <- lm(excum~invf+0)
    sfit <- summary(model)
    print(sfit)
    coeffs <- coefficients(model)
    beta <- unname(coeffs[2])
    int <- unname(coeffs[1])
    dr2 <- data.frame(r=sfit$r.squared, intercept=int, slope=beta, subcl=length(exsubcl), all=length(exsubcl_nohigh))
    write.table(dr2, file=r2, sep="\t", quote=F)
    pdf(fit)
    plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", main=paste0("R2=", round(sfit$r.squared, digits=3)))
    oi <- invf[order(invf)]
    oex <- exsubcl[order(-exsubcl)]
    axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
    abline(model, col="red")
    graphics.off()
} else {
    dr2 <- data.frame(r=NA, intercept=NA, slope=NA, subcl=length(exsubcl), all=length(exsubcl_nohigh))
    write.table(dr2, file=r2, sep="\t", quote=F)
    pdf(fit)
    graphics.off()
}
pdf(histo)
hist(exsubcl_nohigh, breaks=50, cex=1.5, xlab="Allelic frequency (f)", ylab="Number of muts", border="white", col="black", main="")
graphics.off()

if (debug == "yes") {
  save.image(file=paste0(fit,'.debug','.RData'))
}


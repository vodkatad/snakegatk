loweraf <- as.numeric(snakemake@wildcards[["loweraf"]])
higheraf <- as.numeric(snakemake@wildcards[["higheraf"]])
data <- snakemake@input[["afmatrix"]]
afcolumn <- snakemake@params[["afcolumn"]]
fit <- snakemake@output[["fit"]]
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

v <- loweraf
#u=max_range
u <- higheraf

#steps=u:-0.001:v
steps <- seq(u, v, -0.001)

#cumsum=Array(Int64,0)
#v=Array(Float64,0)

#for i in steps
#push!(cumsum,length(filter((x)->x>=i,VAF)))
#push!(v,i)
#end
#cumsum=cumsum-cumsum[1]
cumsum <- c()
v <- c()
for (i in steps) {
      cumsum <- c(cumsum, sum(exsubcl >= i))
  v <- c(v, i)
}
cumsum <- cumsum-cumsum[1] #? # qui non so che cazzo faccia ??? # itś Julia like therefore 1 based counts
# not needed cause we filter muts!
# toglie a tutte le comulative le più grandi del max freq u, ok

#DFcumsum = DataFrame(cumsum=cumsum,v=v)
DFcumsum <- data.frame(cumsum=cumsum,v=v)
DFcumsum$invvaf_sub <- 1/DFcumsum$v  - 1/u # rimappa a zero x -1/fmax

##x,y=hist(VAF,0.0:0.01:1)

#fit constrained fit using GLM fit function
#lmfit=fit(LinearModel, cumsum ~ invVAF + 0 , DFcumsum)
model <- lm(cumsum~invvaf_sub+0, data=DFcumsum)

maxi <- nrow(DFcumsum)
labels <-  c(1, floor(maxi/5), floor(maxi/2), floor(maxi/0.5), maxi)

sfit <- summary(model)
print(sfit)
coeffs <- coefficients(model)
beta <- unname(coeffs[2])
int <- unname(coeffs[1])
dr2 <- data.frame(r=sfit$r.squared, intercept=int, slope=beta, subcl=length(exsubcl), all=length(exsubcl_nohigh))
write.table(dr2, file=r2, sep="\t", quote=F)
pdf(fit)
invf <- DFcumsum$invvaf_sub
excum <- DFcumsum$cumsum
plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", main=paste0("R2=", round(sfit$r.squared, digits=3)))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
abline(model, col="red")
graphics.off()

if (debug == "yes") {
  save.image(file=paste0(fit,'.debug','.RData'))
}


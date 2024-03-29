library(vcfR)

setwd('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot')

vcf <- read.vcfR('CRC1599LMX0A02001TUMD03000V2_nolohcn3.vcf.gz') #CRC1599PRX0A02002TUMD03000V2_nolohcn3.vcf.gz
vaf <- extract.gt(vcf, element = "AF", as.numeric = TRUE)

loweraf <- 0.1
higheraf <- 0.2

af <- vaf[,1]
exsubcl <- af[af<higheraf & af>loweraf]
exsubcl_nohigh <- af[af >= 0]
excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]<=exsubcl[1:length(exsubcl)]))
invf <- 1/exsubcl - 1/higheraf # subtract to use fit without intercept
maxi <- length(invf)
labels <-  c(1, floor(maxi/5), floor(maxi/2), floor(maxi/0.5), maxi)
model <- lm(excum~invf+0)
sfit <- summary(model)
coeffs <- coefficients(model)
plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", main=paste0("R2=", round(sfit$r.squared, digits=3)))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
abline(model, col="red")


length(vaf)
length(exsubcl)

vcf <- read.vcfR('CRC1599PRX0A02002TUMD03000V2_nolohcn3.vcf.gz')
vaf <- extract.gt(vcf, element = "AF", as.numeric = TRUE)

af <- vaf[,2]
exsubcl <- af[af<higheraf & af>loweraf]
exsubcl_nohigh <- af[af >= 0]
excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]<=exsubcl[1:length(exsubcl)]))
invf <- 1/exsubcl - 1/higheraf # subtract to use fit without intercept
maxi <- length(invf)
labels <-  c(1, floor(maxi/5), floor(maxi/2), floor(maxi/0.5), maxi)
model2 <- lm(excum~invf+0)
sfit2 <- summary(model2)
coeffs2 <- coefficients(model2)
plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", main=paste0("R2=", round(sfit2$r.squared, digits=3)))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
abline(model2, col="red")


length(vaf)
length(exsubcl)

#egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot$ cat CRC1599PRX0A02002TUMD03000V2.wlength3.txt 
#39620343
#egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot$ cat CRC1599LMX0A02001TUMD03000V2.wlength3.txt 
#44902560


coeffs/44902560

coeffs2/39620343

0.9*(coeffs/44902560)
### cn2-3
vcf <- read.vcfR('CRC1599LMX0A02001TUMD03000V2_nolohcnint2_4.vcf.gz') #CRC1599PRX0A02002TUMD03000V2_nolohcn3.vcf.gz
vaf <- extract.gt(vcf, element = "AF", as.numeric = TRUE)

loweraf <- 0.12
higheraf <- 0.24

af <- vaf[,1]
exsubcl <- af[af<higheraf & af>loweraf]
exsubcl_nohigh <- af[af >= 0]
excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]<=exsubcl[1:length(exsubcl)]))
invf <- 1/exsubcl - 1/higheraf # subtract to use fit without intercept
maxi <- length(invf)
labels <-  c(1, floor(maxi/5), floor(maxi/2), floor(maxi/0.5), maxi)
model <- lm(excum~invf+0)
sfit <- summary(model)
coeffs <- coefficients(model)
plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", main=paste0("R2=", round(sfit$r.squared, digits=3)))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
abline(model, col="red")


vcf <- read.vcfR('CRC1599PRX0A02002TUMD03000V2_nolohcnint2_4.vcf.gz')
vaf <- extract.gt(vcf, element = "AF", as.numeric = TRUE)

af <- vaf[,2]
exsubcl <- af[af<higheraf & af>loweraf]
exsubcl_nohigh <- af[af >= 0]
excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]<=exsubcl[1:length(exsubcl)]))
invf <- 1/exsubcl - 1/higheraf # subtract to use fit without intercept
maxi <- length(invf)
labels <-  c(1, floor(maxi/5), floor(maxi/2), floor(maxi/0.5), maxi)
model2 <- lm(excum~invf+0)
sfit2 <- summary(model2)
coeffs2 <- coefficients(model2)
plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", main=paste0("R2=", round(sfit2$r.squared, digits=3)))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
abline(model2, col="red")


length(vaf)
length(exsubcl)

#==> CRC1599LMX0A02001TUMD03000V2.wlength234.txt <==
#  101868154

#==> CRC1599PRX0A02002TUMD03000V2.wlength234.txt <==
#  94830349

m <- coeffs/101868154

p <- coeffs2/94830349

m/p

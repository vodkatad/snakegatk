targ <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdx/mutect/mut_burden.tsv', sep="\t", header=T, stringsAsFactors = F)
wes <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/mutect_paired/mut_burden.tsv', sep="\t", header=T, stringsAsFactors = F)

wes$lmodel <- substr(wes$X, 0,9)
targ$lmodel <- substr(targ$X, 0,9)

m <- merge(wes, targ, by="lmodel")
plot(m$burden.x, m$burden.y, xlab="paired WES", ylab="unpaired targeted")

## phyper X-M

d <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/phyper_0.05_0.2.tsv', sep="\t", header=T)
library(melt)
library(ggplot2)
d$id <- rownames(d)
dl <- melt(d)
dl[dl$value == 0, 'value'] <- 1e-230
dl$log10 <- -log(dl$value)/log(10)
ggplot(data=dl, aes(x=value, fill=variable))+geom_histogram(position="dodge")+theme_bw()+theme(text=element_text(size = 18))
ggplot(data=dl, aes(x=value, fill=variable))+geom_histogram(position="dodge")+theme_bw()+xlim(-0.01, 0.1)+theme(text=element_text(size = 18))

o <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/o', sep="\t", header=F, stringsAsFactors = F)

o$loweraf <- as.numeric(sapply(strsplit(o$V1, c("_"), fixed=TRUE), `[[`, 2))
o$higheraf <-  as.numeric(gsub('\\.tsv:.+', '', sapply(strsplit(o$V1, c("_"), fixed=TRUE), `[[`, 3)))
#o[is.na(o$V3), 'V3'] <- 1
o <- o[o$higheraf > o$loweraf,]
ggplot(data=o[o$loweraf==0.12,], aes(x=higheraf, y=-log10(V3)))+geom_point()+theme_bw()+ylab('-log10(Psubcl)')+geom_hline(yintercept=2)+theme(text=element_text(size = 18))
ggplot(data=o, aes(x=higheraf, y=-log10(V3), color=as.factor(loweraf)))+geom_point(alpha=0.8, size=5)+theme_bw()+ylab('-log10(Psubcl)')+geom_hline(yintercept=2)+theme(text=element_text(size = 18))+scale_colour_brewer(palette='Spectral')

ggplot(data=o, aes(x=higheraf, y=-log10(V2)))+geom_point(alpha=0.8, size=5)+theme_bw()+ylab('-log10(Pcl)')+geom_hline(yintercept=2)+theme(text=element_text(size = 18))


ggplot(data=o, aes(x=loweraf, y=-log10(V3), color=as.factor(higheraf)))+geom_point(alpha=0.8, size=5)+theme_bw()+ylab('-log10(Psubcl)')+geom_hline(yintercept=2)+theme(text=element_text(size = 18))+scale_colour_brewer(palette='Spectral')

ggplot(data=o, aes(x=loweraf, y=-log10(V2), color=as.factor(higheraf)))+geom_point(alpha=0.8, size=5)+theme_bw()+ylab('-log10(Pcl)')+geom_hline(yintercept=2)+theme(text=element_text(size = 18))+scale_colour_brewer(palette='Spectral')



ggplot(data=o, aes(x=higheraf, y=-log10(V2), color=loweraf))+geom_point()+theme_bw()+ylab('-log10(Pcl)')+geom_hline(yintercept=2)+theme(text=element_text(size = 18))


o <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/oo', sep="\t", header=F, stringsAsFactors = F)
o$smodel <- sapply(strsplit(o$V1, c(":"), fixed=TRUE), `[[`, 2)
o$loweraf <- as.numeric(sapply(strsplit(o$V1, c("_"), fixed=TRUE), `[[`, 2))
o$higheraf <-  as.numeric(gsub('\\.tsv:.+', '', sapply(strsplit(o$V1, c("_"), fixed=TRUE), `[[`, 3)))
#o[is.na(o$V3), 'V3'] <- 1
o <- o[o$higheraf > o$loweraf,]

ggplot(data=o[o$loweraf==0.05,], aes(x=higheraf, y=-log10(V3), color=as.factor(loweraf)))+
  geom_point(alpha=0.7)+theme_bw()+ylab('-log10(Psubcl)')+geom_hline(yintercept=2)+
  geom_line(aes(group=smodel))+
  theme(text=element_text(size = 18))+scale_color_manual(values="brown")


ggplot(data=o[o$loweraf==0.15,], aes(x=higheraf, y=-log10(V3), color=as.factor(loweraf)))+
  geom_point(alpha=0.7)+theme_bw()+ylab('-log10(Psubcl)')+geom_hline(yintercept=2)+
  geom_line(aes(group=smodel))+
  theme(text=element_text(size = 18))+scale_color_manual(values="dodgerblue2")


ggplot(data=o[o$loweraf==0.05,], aes(x=higheraf, y=-log10(V2), color=as.factor(loweraf)))+
  geom_point(alpha=0.7)+theme_bw()+ylab('-log10(Psubcl)')+geom_hline(yintercept=2)+
  geom_line(aes(group=smodel))+
  theme(text=element_text(size = 18))+scale_color_manual(values="brown")


ggplot(data=o[o$loweraf==0.15,], aes(x=higheraf, y=-log10(V2), color=as.factor(loweraf)))+
  geom_point(alpha=0.7)+theme_bw()+ylab('-log10(Psubcl)')+geom_hline(yintercept=2)+
  geom_line(aes(group=smodel))+
  theme(text=element_text(size = 18))+scale_color_manual(values="dodgerblue2")


ggplot(data=o[o$loweraf==0.05,], aes(x=higheraf, y=-log10(V2), color=as.factor(loweraf)))+
  geom_point(alpha=0.7)+theme_bw()+ylab('-log10(Psubcl)')+geom_hline(yintercept=2)+
  geom_line(aes(group=smodel))+
  theme(text=element_text(size = 18))+scale_color_manual(values="brown")


ggplot(data=o[o$loweraf==0.15,], aes(x=higheraf, y=-log10(V2), color=as.factor(loweraf)))+
  geom_point(alpha=0.7)+theme_bw()+ylab('-log10(Psubcl)')+geom_hline(yintercept=2)+
  geom_line(aes(group=smodel))+
  theme(text=element_text(size = 18))+scale_color_manual(values="dodgerblue2")



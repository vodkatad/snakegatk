library(ggplot2)
library(reshape)


d <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/cnvkit/cmp.tsv', sep= "\t", header=F)
colnames(d) <-  c('file', 'chr', 'b', 'e', 'cnvkit', 'seq')

alls <- unique(d$file)

count_class <- function(s, mydata, thr_delta) {
  myd <- mydata[mydata$file == s, ]
  myd <- myd[myd$cnvkit!=-1 | myd$seq!=-1,]
  miss_cnvkit <- nrow(myd[myd$cnvkit==-1,])
  miss_seq <- nrow(myd[myd$seq==-1,])
  myd <- myd[myd$cnvkit!=-1 & myd$seq!=-1,]
  same <- nrow(myd[myd$cnvkit == myd$seq,])
  myd <- myd[myd$cnvkit != myd$seq,]
  delta <- abs(myd$cnvkit - myd$seq)
  different <- length(delta)
  larger <- sum(delta > thr_delta)
  return(c(miss_cnvkit, miss_seq, same, different-larger, larger))
}

tt <- t(as.data.frame(sapply(alls, count_class, d, 3)))
colnames(tt) <- c('miss_cnvkit', 'miss_seq', 'same', 'different', 'larger3')


pd <- melt(tt)
colnames(pd) <- c('sample', 'class', 'N')
pd$sample <- gsub('cnvkit/', '', pd$sample, fixed=T)
pd$sample <- gsub('_unionbedg.bed', '', pd$sample, fixed=T)

ggplot(data=pd, aes(x=sample, y=N, fill=class))+geom_col(position='stack')+theme_bw(base_size = 20)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+scale_fill_manual(values=rev(c('red', 'orange', 'darkgreen', 'grey', 'violet')))


count_class_len <- function(s, mydata, thr_delta) {
  myd <- mydata[mydata$file == s, ]
  myd <- myd[myd$cnvkit!=-1 | myd$seq!=-1,]
  miss_cnvkit <- sum(myd[myd$cnvkit==-1,'len'])
  miss_seq <- sum(myd[myd$seq==-1,'len'])
  myd <- myd[myd$cnvkit!=-1 & myd$seq!=-1,]
  same <- sum(myd[myd$delta == 0,'len'])
  different <- sum(myd[myd$delta <= thr_delta & myd$delta !=0,'len'])
  larger <- sum(myd[myd$delta > thr_delta,'len'])
  return(c(miss_cnvkit, miss_seq, same, different-larger, larger))
}

d$len <- d$e - d$b
d$delta <- abs(d$cnvkit - d$seq)
tt <- t(as.data.frame(sapply(alls, count_class_len, d, 3)))
colnames(tt) <- c('miss_cnvkit', 'miss_seq', 'same', 'different', 'larger3')


pd <- melt(tt)
colnames(pd) <- c('sample', 'class', 'len')
pd$sample <- gsub('cnvkit/', '', pd$sample, fixed=T)
pd$sample <- gsub('_unionbedg.bed', '', pd$sample, fixed=T)

ggplot(data=pd, aes(x=sample, y=len, fill=class))+geom_col(position='stack')+theme_bw(base_size = 20)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+scale_fill_manual(values=rev(c('red', 'orange', 'darkgreen', 'grey', 'violet')))



#/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot$ cat cnvkit/*comp.txt > oo

dd <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/oo')
dd$tool <- ifelse(grepl('sequenza', dd$V2), 'sequenza', 'cnvkit')
colnames(dd) <- c('N_segments', 'V2', 'tool')


ggplot(data=dd, aes(x=tool, y=N_segments))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0, size=1)+theme_bw(base_size = 20)

# egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot$ cat cnvkit/*_len.txt  > o

dd2 <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/o')
dd2$tool <- ifelse(grepl('sequenza', dd2$V1), 'sequenza', 'cnvkit')
colnames(dd2) <- c('V1', 'tot_len', 'tool')


dd3 <- merge(dd, dd2, by.x='V2', by.y='V1')
dd3$ave_len <- dd3$tot_len / dd3$N_segments
ggplot(data=dd3, aes(x=tool.x, y=ave_len))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0, size=1)+theme_bw(base_size = 20)

#egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot$ for f in sequenza/*/*_confints_CP.txt; do echo -en "$f\t"; head -n2 $f | tail -n1 | cut -f 3; done |tr "/" "\t"| cut -f 2,4

#https://groups.google.com/g/sequenza-user-group/c/y9mQogdBuh8?pli=1

pl <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/cnvkit/ploidy.txt', sep="\t")
colnames(pl) <- c('V1', 'ploidy')
pl$ploidy <- round(pl$ploidy, 2)
pp <- merge(pd, pl, by.x="sample", by.y="V1")

pp1 <- pp[abs(3-pp$ploidy) < 0.5,]
pp2 <- pp[abs(3-pp$ploidy) >= 0.5,]

ppp1 <- unique(pp1[, c('sample', 'ploidy')])
ppp2 <- unique(pp2[, c('sample', 'ploidy')])

ggplot(data=pp1, aes(x=sample, y=len, fill=class))+geom_col(position='stack')+theme_bw(base_size = 20)+
scale_fill_manual(values=rev(c('red', 'orange', 'darkgreen', 'grey', 'violet')))+
  scale_x_discrete(label = ppp1$ploidy)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(data=pp2, aes(x=sample, y=len, fill=class))+geom_col(position='stack')+theme_bw(base_size = 20)+
  scale_fill_manual(values=rev(c('red', 'orange', 'darkgreen', 'grey', 'violet')))+
  scale_x_discrete(label = ppp2$ploidy)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

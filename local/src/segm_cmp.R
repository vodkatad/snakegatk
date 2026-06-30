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

### early late xeno ##############################################

d <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_earlylate_xeno/cnvkit/cmp.tsv', sep= "\t", header=F)
colnames(d) <-  c('file', 'chr', 'b', 'e', 'cnvkit', 'seq')

alls <- unique(d$file)

d$len <- d$e - d$b
d$delta <- abs(d$cnvkit - d$seq)
tt <- t(as.data.frame(sapply(alls, count_class_len, d, 3)))
colnames(tt) <- c('miss_cnvkit', 'miss_seq', 'same', 'different', 'larger3')


pd <- melt(tt)
colnames(pd) <- c('sample', 'class', 'len')
pd$sample <- gsub('cnvkit/', '', pd$sample, fixed=T)
pd$sample <- gsub('_unionbedg.bed', '', pd$sample, fixed=T)

pl <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_earlylate_xeno/cnvkit/ploidy.txt', sep="\t")
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

### segm eval
#library(ggplot2)
library(dplyr)
d <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/cnvkit/eval_ploidies_all.tsv.gz'), sep="\t", header=F)
colnames(d) <- c('sample', 'ploidy', 'chr', 'b', 'e', 'n_ok', 'n_ov')
dd <- d[d$n_ov != 0,] # TODO evaluate n.
dd$frac <- dd$n_ok / dd$n_ov

ddave <- dd |> 
  dplyr::group_by(sample, ploidy) |>
  dplyr::summarise(
    avefrac = mean(frac), minfrac=min(frac), maxfrac=max(frac), n = n()
  )


ggplot(data=ddave, aes(x=sample, y=avefrac, fill=ploidy))+geom_col(position='dodge')+
  theme_bw(base_size=18)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())



ggplot(data=ddave, aes(x=ploidy, y=avefrac))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(height=0, aes(color=sample))+
  theme_bw(base_size=18)+
  theme(legend.position='none')



ggplot(data=ddave, aes(x=ploidy, y=minfrac))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(height=0, aes(color=sample))+
  theme_bw(base_size=18)+
  theme(legend.position='none')


pl <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/cnvkit/ploidy.txt', sep="\t")

colnames(pl) <- c('sample', 'pl')
pl$realploidy <- cut(pl$pl, breaks=c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5))
# pl$s <- substr(pl$sample, 0,10)
# ddave$s <- substr(ddave$sample, 0,10)
m <- merge(ddave, pl, by='sample')

ggplot(data=m, aes(x=ploidy, y=avefrac))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(height=0, aes(color=realploidy))+
  theme_bw(base_size=18)

  mean(dd[dd$sample=='CRC0053LMX0A02204TUMD07000V2', 'frac'])
  
  m$p <- as.numeric(gsub('p', '', m$ploidy))
  m$frac_pen <- m$avefrac/(m$p/6)

  ggplot(data=m, aes(x=ploidy, y=frac_pen))+geom_boxplot(outlier.shape=NA)+
    geom_jitter(height=0, aes(color=realploidy))+
    theme_bw(base_size=18)
  
  ##
  ggplot(data=dd, aes(x=frac))+geom_histogram()+facet_wrap(~ploidy)+theme_bw(base_size=18)
  
  dd$len <- dd$e-dd$b
  
  ddlen_fragmentsok <- dd |> 
    dplyr::group_by(sample, ploidy) |>
    dplyr::filter(frac > 0.9) |>
    dplyr::summarise(
      totlen = sum(len), n = n()
    )
  
  m2 <- merge(ddlen_fragmentsok, pl, by='sample')
  
  ggplot(data=m2, aes(x=ploidy, y=totlen))+geom_boxplot(outlier.shape=NA)+
    geom_jitter(height=0, aes(color=realploidy))+
    theme_bw(base_size=18)
  
  
  ggplot(data=m2, aes(x=ploidy, y=n))+geom_boxplot(outlier.shape=NA)+
    geom_jitter(height=0, aes(color=realploidy))+
    theme_bw(base_size=18)
  
  
  
  ggplot(data=dd, aes(x=n_ov))+geom_histogram()+facet_wrap(~ploidy)+theme_bw(base_size=18)+scale_x_log10()
  ggplot(data=dd, aes(x=n_ok))+geom_histogram()+facet_wrap(~ploidy)+theme_bw(base_size=18)+scale_x_log10()
  
  
  ggplot(data=d, aes(x=n_ov))+geom_histogram()+facet_wrap(~ploidy)+theme_bw(base_size=18)+scale_x_log10()
  ggplot(data=d, aes(x=n_ok))+geom_histogram()+facet_wrap(~ploidy)+theme_bw(base_size=18)+scale_x_log10()
  
#
d <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/cnvkit/eval_seqploidies_all.tsv.gz'), sep="\t", header=F)
colnames(d) <- c('sample', 'chr', 'b', 'e', 'cn', 'n_ov', 'ok2', 'ok3', 'ok4', 'ok5', 'ok6')
dim(d)
dd <- d[d$n_ov != 0,] # TODO evaluate n.
dim(dd)

# 5 is cn, 6 n_ov, numerator is 
# 7 8 9 10 11
# 2 3 4 5  6
# x[5]+5 is > 11 then 11
map_cn <- function(i) {
  mapped <- i+5
  if (mapped > 11) {
    return(11)
  } else {
    return(mapped)
  }
}

#fracok <- apply(dd, 1, function(x) { x[map_cn(as.integer(x[5]))] /  as.integer(x[6])})
dd$frac <- NA 
for (i in 1:nrow(dd)) {
  dd[i, 'frac'] <- dd[i, map_cn(dd[i,'cn'])] / dd[i, 'n_ov']
}
hist(dd$frac)
dd$sel <- NA 
for (i in 1:nrow(dd)) {
  sel <- which.max(dd[i, c(7,8,9,10,11)])
  dd[i, 'sel'] <- sel+1
}
                                       
table(dd$sel == dd$cn)

dd6 <- dd[dd$cn <= 6,]
table(dd6$sel == dd6$cn)

library(reshape)
dd$id <- paste0(dd$sample, dd$chr, dd$b)
long <- melt(dd, id.vars='id', measure.vars=c("ok2", "ok3", "ok4", "ok5","ok6"))
             
mm <- merge(long, dd[, c('id', 'cn', 'n_ov')], by="id")  

mm$frac <- mm$value / mm$n_ov

library(dplyr)
mm$sample <- sapply(strsplit(mm$id, "X"), function(x) {x[[1]]})
ddave <- mm |> 
  dplyr::group_by(sample, variable) |>
  dplyr::summarise(
    avefrac = mean(frac), minfrac=min(frac), maxfrac=max(frac), n = n()
  )


ggplot(data=ddave, aes(x=variable, y=avefrac))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(height=0)+
  theme_bw(base_size=18)

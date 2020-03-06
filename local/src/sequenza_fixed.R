#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = T)
sample <- args[1]
seqz <- args[2]
outdir <- args[3]
cores <- as.numeric(args[4])
library("sequenza")

#extractData <- sequenza.extract(seqz,window = 10e5,min.reads=50,min.reads.normal=20) # def 10e6 40 10, this is sottoriva's our depth is 30x..
extractData <- sequenza.extract(seqz, min.reads=20, min.reads.normal=10)# parallel = cores) # def 10e6 40 10, this is sottoriva's our depth is 30x..
#extractData.CP <- sequenza.fit(extractData, segment.filter = 5e6)
#extractData.CP <- sequenza.fit(extractData, mc.cores = cores)
#sequenza.results(extractData, extractData.CP, out.dir = outdir ,sample.id = sample, cellularity=as.numeric(args[[5]]), ploidy=as.numeric(args[[6]])
sequenza.results(extractData, out.dir = outdir ,sample.id = sample, cellularity=as.numeric(args[[5]]), ploidy=as.numeric(args[6]) )
#save.image(paste0(sample,"_sequenza.RData"))

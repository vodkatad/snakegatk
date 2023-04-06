library(ggplot2)

fit_r2_f <- snakemake@input[["fit_r2_slope"]]
r2_plot_f <- snakemake@output[["r2_plot"]]
delta_plot_f <- snakemake@output[["delta_plot"]]
log_f <- snakemake@log[["log"]]
r2_thr <- as.numeric(snakemake@params[['r2_thr']])
subcl_thr <- as.numeric(snakemake@params[['subcl_thr']])

#r2_thr <- 0.9
#subcl_thr <- 10
#fit_r2_f <- '/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/bestbet_0.05_0.2.tsv'

printf <- function(...) cat(sprintf(...))

fit_r2 <- read.table(fit_r2_f, sep="\t", header=TRUE)

fit_r2$lmodel <- substr(rownames(fit_r2), 0, 10)
fit_r2$smodel <- substr(rownames(fit_r2), 0, 7)
fit_r2$mp <- substr(rownames(fit_r2), 8, 10)
fit_r2$mp <- factor(fit_r2$mp, levels=c('PRX', 'LMX'))

ggplot(data=fit_r2, aes(x=r, fill=mp))+
  geom_histogram(position="dodge", bins=25)+
  theme_bw()+theme(text=element_text(size = 18))+ 
  scale_fill_manual(values=c('#adacac', '#595959'))
ggsave(r2_plot_f)

fit_keep <- fit_r2[fit_r2$r > 0.90 & fit_r2$subcl > 10,]

pairs <- as.data.frame(table(fit_keep$smodel))
with_pair <- fit_keep[fit_keep$smodel %in% pairs[pairs$Freq == 2,'Var1'],]

sink(log_f)
printf('Total\t%d\nFiltered\t%d\nPairs\t%d\n', nrow(fit_r2), nrow(fit_keep), length(unique(with_pair$smodel)))
sink()

# order by delta
lar <- function(x, data) {
  sub <- data[data$smodel == x,]
  sub[sub$mp=="LMX",'intercept'] - sub[sub$mp=="PRX",'intercept']
}


# order by absolute slope of met
#meslo <- function(x, data) {
#  sub <- data[data$smodel == x,]
#  sub[sub$mp=="LMX",'intercept']
#}


with_pair_order <- as.data.frame(sapply(unique(with_pair$smodel), lar, with_pair))
with_pair_order <- with_pair_order[order(with_pair_order[,1]),, drop=FALSE]

d3 <- merge(with_pair_order, with_pair, by.x="row.names", by.y="smodel")
colnames(d3)[2] <- 'delta'
colnames(d3)[1] <- 'smodel'
ggplot(data=d3, aes(x=reorder(smodel,delta), y=intercept, fill=mp))+
  geom_col(position="dodge")+theme_bw()+
  theme(text=element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Patient')+ylab('Slope μ/β')+
  scale_fill_manual(values=c('#adacac', '#595959'))+
  guides(fill=guide_legend(title=""))

ggsave(delta_plot_f)

met <- d3[d3$mp == "LMX",]
pri <- d3[d3$mp == "PRX",]
if (!all(met$smodel==pri$smodel)) {
  stop('llama! Qualquadra non cosa in pri-met pairs')
}
ti <- t.test(met$intercept, pri$intercept, alternative="greater", paired=TRUE)
wi <- wilcox.test(met$intercept, pri$intercept, alternative="greater", paired=TRUE)

sink(log_f, append=TRUE)
printf('ttest\t%.4f\nwilcox\t%.4f\nPri_avg\t%.4f\nMet_avg\t%.4f\nMet_larger\t%.1f\n', 
       ti$p.value, wi$p.value, mean(pri$intercept), mean(met$intercept), sum(d3$delta > 0)/2)
summary(fit_r2[fit_r2$mp == "PRX", 'r'])
summary(fit_r2[fit_r2$mp == "LMX", 'r'])
sink()



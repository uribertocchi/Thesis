library(tidyverse)
library(ggpubr)
library(jtools)
library(stringr)
library(scales)
library(gridExtra)
library(psych)

EG_name <- read_csv("/specific/elkon/uribertocchi/Meth_RNA-seq/Roadmap data/Roadmap_WGBS/EG_name.csv")
EG_name <- EG_name[,-1]
EG <- as.data.frame(colnames(WGBS_data))
EG <- left_join(EG, EG_name,by = c("colnames(WGBS_data)" = 'number'))
EG <- EG$'name'

#correlate WGBS and RNAseq
WGBS_data <- read_csv("WGBS data with methylation avg per 1kb TSS.csv")
WGBS_data <- as.data.frame(WGBS_data)
rownames(WGBS_data) <- RNAseq_data$X1
WGBS_data <- WGBS_data[,-1]
colnames(WGBS_data) <- EG

RNAseq_data <- read_csv("/specific/elkon/uribertocchi/Meth_RNA-seq/Roadmap data/Roadmap_WGBS/RNAseq data.csv")
RNAseq_data <- as.data.frame(RNAseq_data)
rownames(RNAseq_data) <- RNAseq_data[,1]
RNAseq_data <- RNAseq_data[,-1]
colnames(RNAseq_data) <- EG

cor <- cor(as.data.frame(WGBS_data),as.data.frame(RNAseq_data),method = "spearman", use = 'complete.obs')

#scatterplot selected sample
data <- as.data.frame(cbind(WGBS_data$E071, RNAseq_data$E071))
rownames(data) <- rownames(RNAseq_data)
ggscatter(data = data, x = 'V1', y = 'V2',
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "DNA Methylation (%)", ylab = "Gene Expression (RPKM)") + yscale('log10')

data <- as.data.frame(cbind(WGBS_data$E071, RNAseq_data$E071))
rownames(data) <- rownames(RNAseq_data)

data$quantile <- case_when(data$V1 >= 0.75 ~ '75%-100%',
                           data$V1 <= 0.25 ~ '25%', 
                           TRUE ~ '26%-74%'
)  

data$quantile <- factor(data$quantile , levels=c("25%", "26%-74%", "75%-100%"))

give.n <- function(y) {
  return( 
    data.frame(
      y = 0.5+1.1*max(y),  #may need to modify this depending on your data
      label = paste('n =', length(y), '\n'), '\n')
  )}



#boxplot loop (working with new version of R  )
g <- lapply(1:33, function(i){
  df <- data.frame(x=WGBS_data[,i],y=RNAseq_data[,i])
  df$quantile <- case_when(df[,1] >= 0.75 ~ '75%-100%',
                           df[,1] <= 0.25 ~ '25%', 
                           TRUE ~ '26%-74%') 
  df$quantile <- factor(df$quantile , levels=c("25%", "26%-74%", "75%-100%"))
  ggplot(data = df, aes(x = as.factor(df$quantile), y = as.numeric(df$y))) + geom_boxplot() + labs(x = "DNA Methylation (%)", y = "Gene Exp. RPKM (log2)", title = EG[i]) + scale_y_continuous(trans = log2_trans(),
                                                                                                                                                                                               breaks = trans_breaks("log2", function(x) 2^x),
                                                                                                                                                                                               labels = trans_format("log2", math_format(2^.x)))  + stat_summary(fun.data = give.n, geom = "text", vjust = .8) + theme_apa()
})
grid.arrange(grobs = g, ncol = 6)

#Spearman correlation for WGBS and RNAseq in for each gene
#transpose data
t_WGBS <- t(WGBS_data)
rownames(t_WGBS) <- EG
colnames(t_WGBS) <- rownames(RNAseq_data)

t_RNA <- t(RNAseq_data)
rownames(t_RNA) <- EG
colnames(t_RNA) <- rownames(RNAseq_data)

t_cor <- cor(t_WGBS, t_RNA, method = "spearman", use = "pairwise.complete.obs")
t_cor <- corAndPvalue(t_WGBS, t_RNA, method = "spearman")

#Spearman correlation
M <- data.frame(matrix(ncol=3,nrow=ncol(t_RNA)))
M[,1] <- as.character()
G <- colnames(t_RNA)
for(i in 1:ncol(t_RNA)){
  M[i,1] <- G[i]     
  cor <- corr.test(as.data.frame(t_RNA[,i]),as.data.frame(t_WGBS[,i]),
                   method="spearman",adjust="none")
  M[i,2] <- cor$r
  M[i,3] <- cor$p
}
rownames(M) <- M$X1
M <- M[,-1]
M <- M[order(M$X2, decreasing = TRUE), ]  # Top N highest values by group
Spearman_correlation_WGBSXRNAseq <- M
write.csv(Spearman_correlation_WGBSXRNAseq, "Spearman_correlation_WGBSXRNAseq.csv")

#plot top 10
top10genes <- M[1:10,1:2]
top10genes$'names' <- rownames(top10genes)

s <- lapply(1:10, function(i){
  gene <- top10genes$'names'[i]
  df <- data.frame(x=t_WGBS[,gene],y= t_RNA[,gene])
  sp <- ggplot(data = df, aes(x = as.factor(df$x), y = as.numeric(df$y))) + stat_cor(method = "spearman", label.x = 3, label.y = 30) + geom_point() + labs(x = "DNA Methylation (%)", y = "Gene Exp. RPKM (log2)", title = gene) + theme_bw() + theme(axis.text.x=element_blank(),
                                                                                                                                                                                           axis.ticks.x=element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
  
                                                                                                                                                                                                     
})


grid.arrange(grobs = s, ncol = 2)

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggplotify)

WGBS_data <- read.csv("Data/WGBS_data.csv", row.names=1)
RNAseq_data <- read.csv("Data/RNAseq_data.csv", row.names=1)

#hierarchical clustering:
##Give columns tissue names
EG_name <- read.csv("~/GitHub/Thesis/data/EG_name.csv", row.names=1)
EG <- as.data.frame(colnames(WGBS_data))
EG <- left_join(EG, EG_name,by = c("colnames(WGBS_data)" = 'number'))
EG <- EG$'name'
EG_name <- EG_name[,-1]
colnames(RNAseq_data) <- EG
colnames(WGBS_data) <- EG

## make hierarchical clustering heatmap
###RNAseq
counts_rna <- as.matrix(RNAseq_data)
se_rna <- SummarizedExperiment(assays=list(counts=counts_rna))
rld_mat_rna <- assay(se_rna)
rld_cor_rna <- cor(rld_mat_rna,method='spearman')
head(rld_cor_rna) 
pheatmap(rld_cor_rna)
###WGBS
counts_WGBS <- as.matrix(WGBS_data)
se_WGBS <- SummarizedExperiment(assays=list(counts=counts_WGBS))
rld_mat_WGBS <- assay(se_WGBS)
rld_cor_WGBS <- cor(rld_mat_WGBS, use = 'complete.obs', method='pearson')
head(rld_cor_WGBS) 
pheatmap(rld_cor_WGBS)

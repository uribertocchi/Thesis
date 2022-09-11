.libPaths(new = "/specific/elkon/uribertocchi/RLibs/")

library(tidyverse)
library(rtracklayer)
library(GenomeInfoDb)
library(plyranges)
library(stringr)

#import all samples and merge them into one
setwd("/specific/elkon/uribertocchi/Meth_RNA-seq/Roadmap data/Roadmap_WGBS")
bw_merge = function(mypath){ #create an object as path to folder
  filenames=list.files(path=mypath, pattern = "*.bigwig", recursive = TRUE)
  datalist = lapply(filenames, function(x){import.bw(con=x)})
  Reduce(function(x,y) {join_overlap_left(x,y)}, datalist)
}

WGBS_path <- "/specific/elkon/uribertocchi/Meth_RNA-seq/Roadmap data/Roadmap_WGBS"
WGBS_joined <- bw_merge(WGBS_path) #join all samples
head(WGBS_joined)
save.image(file = "Roadmap_WGBS_280822.RData")

WGBS_samplenames <- EG_mnemonics_name$X1
head(WGBS_samplenames)

#change metadata columns names
granges <- as.data.frame(mcols(WGBS_joined))
colnames(granges) <- WGBS_samplenames
mcols(WGBS_joined) <- granges
head(WGBS_joined)
write_csv(as.data.frame(WGBS_joined), "/specific/elkon/uribertocchi/Meth_RNA-seq/Roadmap data/WGBS_joined.csv")

#Cut TSS Location and expand tss to 1kb
tss <- read_csv("/specific/elkon/uribertocchi/Meth_RNA-seq/Roadmap data/Roadmap_WGBS/TSS.csv")
numbers = as_tibble((0:1000))
tss$fake_col <- 1
numbers$fake_col <- 1
cj_tss <- full_join(tss, numbers, by = "fake_col") %>%
  mutate(TSS_start = TSS_start + value)%>%
  select(-c(fake_col, value))
tss_df_regions = data.frame(gene_id = cj_tss$Ensembl,chromosome = cj_tss$Chromosome,
                        start=cj_tss$TSS_start,
                        strand = cj_tss$Strand)
tss_granges <- as(tss_df_regions,'GRanges') # note that names have to match with GRanges slots

#join fraction methylation with tss
joined_tss_WGBS <- tss_df_regions %>% full_join(as.data.frame(WGBS_joined), 
                           by = c("start", 'chromosome' = 'seqnames'))
head(joined_tss_WGBS)
write.csv(joined_tss_WGBS, "joined_tss_WGBS.csv")

tss_WGBS_mean <- aggregate( joined_tss_WGBS[, 8:44], by = list(joined_tss_WGBS$gene_id), FUN = mean, na.rm = TRUE)
tss_WGBS_mean <- as.data.frame(tss_WGBS_mean)
tss_WGBS_mean2 <- tss_WGBS_mean[,-1]
rownames(tss_WGBS_mean2) <- tss_WGBS_mean[,1]
tss_WGBS_mean <- tss_WGBS_mean2
write.csv(tss_WGBS_mean, "tss_WGBS_mean.csv")

#Keep only samples we have data in RNA-seq and WGBS, make gene_id the row name
WGBS_names <- colnames(tss_WGBS_mean)
WGBS_names <- WGBS_names[-1]
RNAseq_names <- colnames(X57epigenomes_RPKM_pc)
RNAseq_names <- intersect(RNAseq_names, WGBS_names)
RNAseq_names <- c('gene_id', RNAseq_names)
RNAseq <- X57epigenomes_RPKM_pc[, RNAseq_names]
RNAseq <- as.data.frame(RNAseq)
RNAseq2 <- RNAseq[,-1]
rownames(RNAseq2) <- RNAseq[,1]
RNAseq <- RNAseq2


WGBS_names <- intersect(RNAseq_names, WGBS_names)
RNAseq_rownames <- rownames(RNAseq)
WGBS_tss_mean <- tss_WGBS_mean[RNAseq_rownames,WGBS_names]
WGBS_data <- WGBS_tss_mean

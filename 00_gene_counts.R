# Get gene counts from aligned RNAseq data
# Compare counts across treatments

library(Rsubread)
library(GenomicRanges)
library(DESeq2)
library(tidyr)
library(biomaRt)
library(systemPipeR)
library(mygene)

# these are the annotations used
hg19 <- getInBuiltAnnotation(annotation="hg19")
sam_files <- list.files(pattern="*.sam")
overlaps <- featureCounts(sam_files, annot.inbuilt="hg19")
countsGene <- overlaps$counts
# save counts
write.csv(countsGene, "data/counts_gene_featurecounts.csv")

# make annotation df
targets <- overlaps[["targets"]]
RNA <- factor(c("3P", "3P", "SS", "SS", "3P","3P", "SS", "SS"))
protein <- factor(c("KO", "KO","KO","KO","WT","WT","WT","WT"))
treatment <- paste0(protein, RNA)
coldata <- data.frame(targets, treatment, protein, RNA)

# WT vs KO
DESeq_WT_KO <- DESeqDataSetFromMatrix(countData = countsGene,
                                      colData = coldata,
                                      design = ~ protein)
DESeq_WT_KO <- DESeq(DESeq_WT_KO)
DESeq_WT_KO <- DESeq2::results(DESeq_WT_KO)
WT_KO <- DESeq_WT_KO@listData %>% as.data.frame()
WT_KO$GeneName <- DESeq_WT_KO@rownames

# 3P vs SS
DESeq_3P_SS <- DESeqDataSetFromMatrix(countData = countsGene,
                                      colData = coldata,
                                      design = ~ RNA)
DESeq_3P_SS <- DESeq(DESeq_3P_SS)
DESeq_3P_SS <- DESeq2::results(DESeq_3P_SS)
SS_3P <- DESeq_3P_SS@listData %>% as.data.frame()
SS_3P$GeneName <- DESeq_3P_SS@rownames

# WT 3P SS
DESeq_WT3P_WTSS <- DESeqDataSetFromMatrix(countData = countsGene[,5:8],
                                          colData = coldata[5:8,],
                                          design = ~ treatment)
DESeq_WT3P_WTSS <- DESeq(DESeq_WT3P_WTSS)
DESeq_WT3P_WTSS <- DESeq2::results(DESeq_WT3P_WTSS)
WT3P_WTSS <- DESeq_WT3P_WTSS@listData %>% as.data.frame()
WT3P_WTSS$GeneName <- DESeq_WT3P_WTSS@rownames

# WT SS KO SS
DESeq_KOSS_WTSS <- DESeqDataSetFromMatrix(countData = countsGene[,c(3,4,7,8)],
                                          colData = coldata[c(3,4,7,8),],
                                          design = ~ treatment)
DESeq_KOSS_WTSS <- DESeq(DESeq_KOSS_WTSS)
DESeq_KOSS_WTSS <- DESeq2::results(DESeq_KOSS_WTSS)
KOSS_WTSS <- DESeq_KOSS_WTSS@listData %>% as.data.frame()
KOSS_WTSS$GeneName <- DESeq_KOSS_WTSS@rownames

# KO SS 3P
DESeq_KO3P_KOSS <- DESeqDataSetFromMatrix(countData = countsGene[,1:4],
                                          colData = coldata[1:4,],
                                          design = ~ treatment)
DESeq_KO3P_KOSS <- DESeq(DESeq_KO3P_KOSS)
DESeq_KO3P_KOSS <- DESeq2::results(DESeq_KO3P_KOSS)
KO3P_KOSS <- DESeq_KO3P_KOSS@listData %>% as.data.frame()
KO3P_KOSS$GeneName <- DESeq_KO3P_KOSS@rownames

# WT 3P KO 3P
DESeq_KO3P_WT3P <- DESeqDataSetFromMatrix(countData = countsGene[,c(1,2,5,6)],
                                          colData = coldata[c(1,2,5,6),],
                                          design = ~ treatment)
DESeq_KO3P_WT3P <- DESeq(DESeq_KO3P_WT3P)
DESeq_KO3P_WT3P <- DESeq2::results(DESeq_KO3P_WT3P)
KO3P_WT3P <- DESeq_KO3P_WT3P@listData %>% as.data.frame()
KO3P_WT3P$GeneName <- DESeq_KO3P_WT3P@rownames

# add gene names to everything
gene_names<-queryMany(KO3P_KOSS$GeneName, 
                      scopes="entrezgene", 
                      fields="symbol",
                      return.as = "DataFrame", 
                      returnall=F) %>% 
    as.data.frame()
gene_names <- gene_names[,c(1,4)]
all_DESeq <- list(KO3P_KOSS=KO3P_KOSS, 
                  KO3P_WT3P=KO3P_WT3P, 
                  WT3P_WTSS=WT3P_WTSS, 
                  KOSS_WTSS=KOSS_WTSS)
all_DESeq <- lapply(all_DESeq, 
                  function(df) merge(df, gene_names, 
                                     by.x="GeneName", 
                                     by.y="query"))
# save work
saveRDS(all_DESeq, file="results/all_DESeq.rds")
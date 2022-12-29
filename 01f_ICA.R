library(tidyverse)
library(pheatmap)
library("RColorBrewer")
library(fastICA)

expression<-read.csv("PARP13 Table of Values.csv")[,c(1, 10:17)]
row.names(expression)<-expression[,1]
expression<-expression[,-1]
expression<-expression[rowMeans(expression[])>1,] # remove low-expressing genes
expression<-log(expression + 1) %>% as.matrix()
#boxplot(expression) # all samples have same distribution after transformation

ICA<-fastICA(expression, 3)

### create Figure 1F
samples_ICA<-t(scale(t(ICA[["A"]]), scale = F))
colnames(samples_ICA)<-c("KO3p1", "KO3p2", "KOss1", "KOss2", 
                         "WT3p1", "WT3p2", "WTss1", "WTss2")
rownames(samples_ICA) <-c("IC1", "IC2", "IC3")
breaksList <- seq(min(samples_ICA), max(samples_ICA), by = 0.001)
jpeg("ICA.jpg", width = 3, height = 1.5, res = 300, units = "in")
pheatmap(samples_ICA, fontsize = 6, angle_col = 90,
         fontsize_col = 8, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(length(breaksList)),
         show_rownames = T,
         breaks = breaksList,
         legend_breaks = c(-.1, 0, .1))
dev.off()

### pull genes that contribute to each pattern
genes<-ICA[["S"]] %>% as.data.frame() %>% mutate(index = row.names(.))
genes_ICA<-merge(table, genes)
write.csv(genes_ICA, 
          "results/ICA_RNAseq_genes.csv", 
          quote = F, row.names = F)

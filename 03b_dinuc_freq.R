library(tidyverse)
library(ggrepel)

# dinucleotide frequencies derived from compseq of the EMBOSS suit
# using fasta sequences of either CLIP peaks (Figure 3B) or 
# RNAseq reads (not shown in manuscript)

compseq_SS<-read.csv("results/SS_peaks_2.out", sep="\t")
compseq_SS[,2]<-compseq_SS[,1]
compseq_SS[,1]<-row.names(compseq_SS)
compseq_SS[,5]<-lapply(compseq_SS[,5], function(x){
    x<-as.character(x)
    return(substr(x, 1, nchar(x)-1) %>% as.numeric())
}) %>% unlist()
compseq_3P<-read.csv("results/3P_peaks_2.out", sep="\t")
compseq_3P[,2]<-compseq_3P[,1]
compseq_3P[,1]<-row.names(compseq_3P)
compseq_3P[,5]<-lapply(compseq_3P[,5], function(x){
    x<-as.character(x)
    return(substr(x, 1, nchar(x)-1) %>% as.numeric())
}) %>% unlist()

data<-merge(compseq_3P[,c(1,5)], compseq_SS[,c(1,5)], by="Word")
colnames(data)<-c("kmer", "WT3P", "WTSS")

ggplot(data, aes(x=WTSS, y=WT3P))+
    geom_abline(slope=1, color="red", alpha=.5)+
    geom_point(size = 3, alpha = .7)+
    theme_classic() +
    theme(text = element_text(size = 12),
          axis.text = element_text(size = 12)) +
    xlab("ssRNA enrichment")+
    ylab("3p-RNA enrichment") +
    ylim(.6, 2.7)+
    xlim(.6, 2.7) +
    geom_text_repel(aes(label = kmer),
                    max.overlaps = 10)
ggsave("2mer_CLIP_peak_comparison.jpeg", units = "in", height = 3, width = 3)



##### for comparison, how does the RNA-seq dinucleotide frequency look?
#
# compseq_SS_RNA<-read.csv("WTSS_seqs_2.out", sep="\t")
# compseq_SS_RNA[,2]<-compseq_SS_RNA[,1]
# compseq_SS_RNA[,1]<-row.names(compseq_SS_RNA)
# compseq_SS_RNA[,5]<-lapply(compseq_SS_RNA[,5], function(x){
#     x<-as.character(x)
#     return(substr(x, 1, nchar(x)-1) %>% as.numeric())
# }) %>% unlist()
# compseq_3P_RNA<-read.csv("WT3P_seqs_2.out", sep="\t")
# compseq_3P_RNA[,2]<-compseq_3P_RNA[,1]
# compseq_3P_RNA[,1]<-row.names(compseq_3P_RNA)
# compseq_3P_RNA[,5]<-lapply(compseq_3P_RNA[,5], function(x){
#     x<-as.character(x)
#     return(substr(x, 1, nchar(x)-1) %>% as.numeric())
# }) %>% unlist()
# 
# data_RNA<-merge(compseq_3P_RNA[,c(1,5)], compseq_SS_RNA[,c(1,5)], by="Word")
# colnames(data_RNA)<-c("kmer", "WT3P", "WTSS")
# 
# ggplot(data_RNA, aes(x=WTSS, y=WT3P))+
#     geom_abline(slope=1, color="red", alpha=.5)+
#     geom_point(size = 3, alpha = .7)+
#     theme_classic() +
#     theme(text = element_text(size = 12),
#           axis.text = element_text(size = 12)) +
#     xlab("ssRNA enrichment")+
#     ylab("3p-RNA enrichment") +
#     ylim(.3, 1.4)+
#     xlim(.3, 1.4) +
#     geom_text_repel(aes(label = kmer), max.overlaps = 17, 
#                     box.padding = 1)
# ggsave("2mer_RNA_comparison.jpeg")

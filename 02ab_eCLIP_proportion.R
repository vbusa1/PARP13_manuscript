setwd("/Users/veronica/Documents/Leung/PARP13/manuscript/Review")

library(tidyverse)
library(ggsignif)
library(ggforce)
library(Biostrings)

WTSS<-read.delim("../../CLIP\ analysis/WTSS_peaks.tsv", 
                 header = F, stringsAsFactors = F) %>%
    select(-V7, -V10)
WT3P<-read.delim("../../CLIP\ analysis/WT3P_peaks.tsv", 
                 header = F, stringsAsFactors = F) %>%
    select(-V7, -V10)

WTSS_region<-WTSS[,c("V8", "V9")] %>% unique()
WT3P_region<-WT3P[,c("V8", "V9")] %>% unique()

WTSS_region$sample<-"ssRNA"
WT3P_region$sample<-"3p-RNA"

regions<-rbind(WTSS_region, WT3P_region)
regions[regions$V9 %in% c("distintron500",
                          "distnoncoding_intron500",
                          "proxintron500",
                          "proxnoncoding_intron500"), "V9"]<-"intron"
regions<-regions %>% filter(V9 != "stop_codon", V9 != "miRNA")

ggplot(regions, aes(x = sample, 
                    fill = factor(V9, levels = c("intron", "noncoding_exon", "5utr", "CDS", "3utr")))) +
    geom_bar(position = "fill") +
    theme_classic() +
    labs(x = NULL, y = "Proportion", fill = "Gene\nRegion") +
    scale_fill_discrete(labels = c("Intron", "Noncoding exon", "5' UTR", "CDS", "3' UTR")) +
    theme(text = element_text(size = 18),
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 65, h = 1))
ggsave("PARP13_CLIP_gene_region_binding.png", units = "in", width = 4, height = 3)

chi<-regions %>% group_by(sample, V9) %>% 
    summarise(count = n()) %>% 
    spread(V9, count) %>%
    as.data.frame()
row.names(chi)<-chi[,1]
chi<-chi[,2:ncol(chi)]
chisq.test(chi)

## does region bound reflect upregulation in RNA-seq?
RNAseq<-read.csv("../PARP13 update/PARP13 Table of Values.csv")[,c("genename",
                                                                   "RNASEQ.SET.2.KO.vs.WT.SS.log2FoldChange", 
                                                                   "RNASEQ.SET.2.KO.vs.WT.SS.p.adjusted",
                                                                   "RNASEQ.SET.2.P3.vs.SS.WT.log2FoldChange",
                                                                   "RNASEQ.SET.2.P3.vs.SS.WT.p.adjusted",
                                                                   "RNASEQ.SET.2.KO.vs.WT.3P.log2FoldChange", 
                                                                   "RNASEQ.SET.2.KO.vs.WT.3P.p.adjusted",
                                                                   "RNASEQ.SET.2.P3.vs.SS.KO.log2FoldChange",
                                                                   "RNASEQ.SET.2.P3.vs.SS.KO.p.adjusted")]
RNAseq <- RNAseq[which(rowSums(RNAseq[,2:9], na.rm = T) != 0),]

RNAseq_bound <- merge(RNAseq, regions, 
                      by.x = "genename", by.y = "V8", 
                      all.x = T)
RNAseq_bound[is.na(RNAseq_bound$V9), "V9"] <- "none"

# ggplot(filter(RNAseq_bound, !is.na(V9)),
#        aes(x = V9,
#            y = RNASEQ.SET.2.KO.vs.WT.SS.log2FoldChange)) +
#     geom_signif(comparisons = list(c("3utr", "5utr"),
#                                    c("3utr", "CDS"),
#                                    c("3utr", "intron"),
#                                    c("3utr", "noncoding_exon"),
#                                    c("CDS", "noncoding_exon")), 
#                 map_signif_level=TRUE,
#                 y_position = c(1.7, 1.9, 2.1, 2.3, 1.7)) +
#     geom_sina(aes(color = V9)) +
#     geom_boxplot(color = "black", alpha = 0) +
#     theme_classic()


############################
WTSS<-read.delim("../../CLIP\ analysis/WTSS_peaks.tsv", 
                 header = F, stringsAsFactors = F) %>%
    select(-V7, -V10)

WTSS[WTSS$V9 %in% c("distintron500",
                    "distnoncoding_intron500",
                    "proxintron500",
                    "proxnoncoding_intron500"), "V9"] <- "intron"
WTSS<-WTSS %>% filter(V9 != "stop_codon", V9 != "miRNA")


## does region bound reflect upregulation in RNA-seq?
RNAseq<-read.csv("../../PARP13 update/PARP13 Table of Values.csv")[,c("genename",
                                                                   "RNASEQ.SET.2.KO.vs.WT.SS.log2FoldChange", 
                                                                   "RNASEQ.SET.2.KO.vs.WT.SS.p.adjusted")]
RNAseq <- RNAseq[which(rowSums(RNAseq[,2:3], na.rm = T) > 0),] %>% unique()

binding <- data.frame(genes = RNAseq$genename,
                      UTR3 = 0,
                      UTR5 = 0,
                      CDS = 0,
                      intron = 0,
                      noncoding_exon = 0)
for(gene in unique(binding$genes[binding$genes %in% WTSS$V8])){
    hold <- WTSS %>% filter(V8 == gene)
    binding[which(binding$genes == gene), "UTR3"] <- nrow(hold[which(hold$V9 == "3utr"),])
    binding[which(binding$genes == gene), "UTR5"] <- nrow(hold[which(hold$V9 == "5utr"),])
    binding[which(binding$genes == gene), "CDS"] <- nrow(hold[which(hold$V9 == "CDS"),])
    binding[which(binding$genes == gene), "intron"] <- nrow(hold[which(hold$V9 == "intron"),])
    binding[which(binding$genes == gene), "noncoding_exon"] <- nrow(hold[which(hold$V9 == "noncoding_exon"),])
}

RNAseq_bound <- merge(RNAseq, binding,
                      by.x = "genename", by.y = "genes")

plot <- RNAseq_bound %>%
    filter(((UTR3 > 0) + (UTR5 > 0) + (CDS > 0) + (intron > 0) + (noncoding_exon > 0)) == 1) %>%
    gather("site", "events", 
           -genename, 
           -RNASEQ.SET.2.KO.vs.WT.SS.log2FoldChange, 
           -RNASEQ.SET.2.KO.vs.WT.SS.p.adjusted) %>%
    filter(events > 0)

ggplot(plot,
       aes(x = factor(site, levels = c("UTR3", "CDS", "UTR5", "noncoding_exon", "intron")),
           y = RNASEQ.SET.2.KO.vs.WT.SS.log2FoldChange)) +
    geom_sina(aes(color = site),
              scale = "width") +
    geom_boxplot(alpha = 0, width = 0.5) +
    geom_signif(comparisons = list(c("CDS", "UTR3")), 
                map_signif_level=TRUE,
                y_position = c(1.65)) +
    theme_classic() +
    labs(x = "Gene Region", y = bquote(log[2]*FC~KO/WT~+ssRNA)) +
    scale_color_manual(values = c("#00B0F6", "#F8766D", "#A3A500", "#E76BF3", "#00BF7D")) + 
    theme(legend.position = "none", text = element_text(size = 15), axis.text.x=element_blank())
ggsave("CLIP_vs_FC_ss.jpeg", width = 3.5, height = 3.5, units = "in")



#### get unique peaks (e.g. not abutting)

WTSS_merge <- matrix(NA,
                     nrow = 1, 
                     ncol = 8) %>% as.data.frame()
colnames(WTSS_merge) <- colnames(WTSS)

for(g in unique(WTSS$V8)){
    gene <- WTSS %>% filter(V8 == g)
    interval <- c(gene$V2, gene$V3)
    dup <- interval[duplicated(interval)]
    if(length(dup) > 0){
        interval <- interval[!(interval %in% dup)] %>% sort()
        new_dat <- matrix(NA,nrow = length(interval)/2, 
                          ncol = 8) %>% as.data.frame()
        colnames(new_dat) <- colnames(WTSS)
        new_dat$V8 <- g
        new_dat$V1 <- gene[1, "V1"]
        new_dat$V2 <- interval[seq(1, length(interval), 2)]
        new_dat$V3 <- interval[seq(2, length(interval), 2)]
        for(i in 1:(length(interval)/2)){
            new_dat[i, "V9"] <- gene[which(gene$V2 == new_dat[i, "V2"]),"V9"]
            new_dat[i, "V6"] <- gene[which(gene$V2 == new_dat[i, "V2"]),"V6"]
        }
    }else{
        new_dat <- gene
    }
    WTSS_merge <- rbind(WTSS_merge, new_dat)
}

WT3P_merge <- matrix(NA,
                     nrow = 1, 
                     ncol = 8) %>% as.data.frame()
colnames(WT3P_merge) <- colnames(WT3P)

for(g in unique(WT3P$V8)){
    gene <- WT3P %>% filter(V8 == g)
    interval <- c(gene$V2, gene$V3)
    dup <- interval[duplicated(interval)]
    if(length(dup) > 0){
        interval <- interval[!(interval %in% dup)] %>% sort()
        new_dat <- matrix(NA,nrow = length(interval)/2, 
                          ncol = 8) %>% as.data.frame()
        colnames(new_dat) <- colnames(WT3P)
        new_dat$V8 <- g
        new_dat$V1 <- gene[1, "V1"]
        new_dat$V2 <- interval[seq(1, length(interval), 2)]
        new_dat$V3 <- interval[seq(2, length(interval), 2)]
        for(i in 1:(length(interval)/2)){
            new_dat[i, "V9"] <- gene[which(gene$V2 == new_dat[i, "V2"]),"V9"]
            new_dat[i, "V6"] <- gene[which(gene$V2 == new_dat[i, "V2"]),"V6"]
        }
    }else{
        new_dat <- gene
    }
    WT3P_merge <- rbind(WT3P_merge, new_dat)
}


WTSS_write <- WTSS_merge[complete.cases(WTSS_merge),] %>% select(-V8, -V9)
write_delim(WTSS_write,
            "WTSS_peaks_merge.bed",
            col_names = F,
            delim = "\t")

WTSS_fasta <- readDNAStringSet("WTSS_merge.fa") %>%
    as.data.frame()
WTSS_fasta$info <- row.names(WTSS_fasta)
WTSS_fasta <- WTSS_fasta %>% 
    separate(info, c("chr", "start", "end"))

WTSS_all <- merge(WTSS_fasta, WTSS_merge,
                  by.y = c("V1", "V2", "V3"),
                  by.x = c("chr", "start", "end"))
WTSS_all$CG <- grepl("CG", WTSS_all$x, 
                     ignore.case = T)
WTSS_all$TTCGAC <- grepl("TTCGAC", WTSS_all$x, 
                         ignore.case = T)
toMatch <- c("UUCG", "UACG", "AUCG", "AACG")
WTSS_all$motif <- grepl(paste(toMatch,collapse="|"), WTSS_all$x, 
                        ignore.case = T)

dup <- WTSS_all$V8[duplicated(WTSS_all$V8)]
WTSS_mult_binding <- WTSS_all[which(WTSS_all$V8 %in% dup),]

sum(WTSS_mult_binding$CG) # 590 of 689 = 86%
sum(WTSS_mult_binding$TTCGAC) # 6 1%
sum(WTSS_mult_binding$motif) # 85 12%

sum(WTSS_all$CG) # 1053 of 1190 = 88%
sum(WTSS_all$TTCGAC) # 11 1%
sum(WTSS_all$motif) # 148 12%




WT3P_write <- WT3P_merge[complete.cases(WT3P_merge),] %>% select(-V8, -V9)
write_delim(WT3P_write,
            "WT3P_peaks_merge.bed",
            col_names = F,
            delim = "\t")

WT3P_fasta <- readDNAStringSet("WT3P_merge.fa") %>%
    as.data.frame()
WT3P_fasta$info <- row.names(WT3P_fasta)
WT3P_fasta <- WT3P_fasta %>% 
    separate(info, c("chr", "start", "end"))

WT3P_all <- merge(WT3P_fasta, WT3P_merge,
                  by.y = c("V1", "V2", "V3"),
                  by.x = c("chr", "start", "end"))
WT3P_all$CG <- grepl("CG", WT3P_all$x, 
                     ignore.case = T)
WT3P_all$TTCGAC <- grepl("TTCGAC", WT3P_all$x, 
                         ignore.case = T)
toMatch <- c("UUCG", "UACG", "AUCG", "AACG")
WT3P_all$motif <- grepl(paste(toMatch,collapse="|"), WT3P_all$x, 
                        ignore.case = T)

dup <- WT3P_all$V8[duplicated(WT3P_all$V8)]
WT3P_mult_binding <- WT3P_all[which(WT3P_all$V8 %in% dup),]

sum(WT3P_mult_binding$CG) # 820 of 960 = 85%
sum(WT3P_mult_binding$TTCGAC) # 8 1%
sum(WT3P_mult_binding$motif) # 119 12%

sum(WT3P_all$CG) # 1486 of 1667 = 89%
sum(WT3P_all$TTCGAC) # 15 1%
sum(WT3P_all$motif) # 209 13%


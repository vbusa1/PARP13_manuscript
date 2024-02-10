library(tidyverse)
library(gridExtra)
library(ggrepel)

WT3P_WTSS<-readRDS("results/all_DESeq.rds")["WT3P_WTSS"] %>% 
    as.data.frame()
WT3P_WTSS<-WT3P_WTSS[complete.cases(WT3P_WTSS),c(8,3,7)]
colnames(WT3P_WTSS)<-c("gene", "log2fc", "p")
KO3P_WT3P<-readRDS("results/all_DESeq.rds")["KO3P_WT3P"] %>% 
    as.data.frame()
KO3P_WT3P<-KO3P_WT3P[complete.cases(KO3P_WT3P),c(8,3,7)]
colnames(KO3P_WT3P)<-c("gene", "log2fc", "p")
KOSS_WTSS<-readRDS("results/all_DESeq.rds")["KOSS_WTSS"] %>% 
    as.data.frame()
KOSS_WTSS<-KOSS_WTSS[complete.cases(KOSS_WTSS),c(8,3,7)]
colnames(KOSS_WTSS)<-c("gene", "log2fc", "p")
KO3P_KOSS<-readRDS("results/all_DESeq.rds")["KO3P_KOSS"] %>% 
    as.data.frame()
KO3P_KOSS<-KO3P_KOSS[complete.cases(KO3P_KOSS),c(8,3,7)]
colnames(KO3P_KOSS)<-c("gene", "log2fc", "p")

# # make volcano plots
# plot1<- ggplot(WT3P_WTSS, aes(x = -log2fc,
#                               y = -log10(p)))+
#     geom_point(data = subset(WT3P_WTSS, abs(log2fc) > 1 & p < .01), 
#                color = "red") +
#     geom_point(data = subset(WT3P_WTSS, abs(log2fc) <= 1 | p > .01), 
#                color = "gold") +
#     geom_point(data = subset(WT3P_WTSS, abs(log2fc) <= 1 & p > .01), 
#                color = "black") +
#     theme_classic() +
#     labs(x = bquote(log[2]*FC~3*p/ssRNA~WT),
#          y = bquote(-log[10]*p)) +
#     xlim(-10.1, 12.2) + ylim(0, 200) +
#     geom_text_repel(data = filter(WT3P_WTSS, 
#                                   gene %in% c("IFIT1", "IFIT2", "IFIT3", "IFNB1",
#                                               "OAS2", "OAS3", "DDX58", "CCL5",
#                                               "RSAD2", "ISG15", "STAT1", "IFI6")), 
#                     aes(label = gene),
#                     max.overlaps = 15,
#                     size = 2.5,
#                     fontface = 2,
#                     nudge_y = 7)
# plot2<- ggplot(KO3P_WT3P, aes(x = log2fc,
#                               y = -log10(p)))+
#     geom_point(data = subset(KO3P_WT3P, abs(log2fc) > 1 & p < .01), 
#                color = "red") +
#     geom_point(data = subset(KO3P_WT3P, abs(log2fc) <= 1 | p > .01), 
#                color = "gold") +
#     geom_point(data = subset(KO3P_WT3P, abs(log2fc) <= 1 & p > .01), 
#                color = "black") +
#     theme_classic() +
#     theme(axis.ticks.y = element_blank(),
#           axis.text.y = element_blank()) +
#     labs(x = bquote(log[2]*FC~3*pRNA~WT/KO), 
#          y = NULL) +
#     xlim(-10.1, 12.2) + ylim(0, 200) +
#     geom_text_repel(data = filter(KO3P_WT3P, 
#                                   -log10(p) > 80 | 
#                                       log2fc > 10 |
#                                       log2fc < -8),
#                     aes(label = gene),
#                     fontface = 2,
#                     size = 2.5,
#                     nudge_y = 7)
# plot3<- ggplot(KOSS_WTSS, aes(x = log2fc,
#                               y = -log10(p)))+
#     geom_point(data = subset(KOSS_WTSS, abs(log2fc) > 1 & p < .01), 
#                color = "red") +
#     geom_point(data = subset(KOSS_WTSS, abs(log2fc) <= 1 | p > .01), 
#                color = "gold") +
#     geom_point(data = subset(KOSS_WTSS, abs(log2fc) <= 1 & p > .01), 
#                color = "black") +
#     theme_classic() +
#     labs(x = bquote(log[2]*FC~ssRNA~WT/KO),
#          y = bquote(-log[10]*p)) +
#     xlim(-10.1, 12.2) + ylim(0, 200) +
#     geom_text_repel(data = filter(KOSS_WTSS, 
#                                       log2fc > 10 |
#                                       log2fc < -8 |
#                                   gene == "TNFRSF10D"),
#                     aes(label = gene),
#                     fontface = 2,
#                     size = 2.5,
#                     nudge_y = 15)
# plot4<- ggplot(KO3P_KOSS, aes(x = -log2fc,
#                               y = -log10(p)))+
#     geom_point(data = subset(KO3P_KOSS, abs(log2fc) > 1 & p < .01), 
#                color = "red") +
#     geom_point(data = subset(KO3P_KOSS, abs(log2fc) <= 1 | p > .01), 
#                color = "gold") +
#     geom_point(data = subset(KO3P_KOSS, abs(log2fc) <= 1 & p > .01), 
#                color = "black") +
#     theme_classic() +
#     theme(axis.ticks.y = element_blank(),
#           axis.text.y = element_blank()) +
#     labs(x = bquote(log[2]*FC~3*p/ssRNA~KO), 
#          y = NULL) +
#     xlim(-10.1, 12.2) + ylim(0, 200) +
#     geom_text_repel(data = filter(KO3P_KOSS, 
#                                   gene == "IFIT2"), 
#                     aes(label = gene),
#                     max.overlaps = 15,
#                     size = 2.5,
#                     fontface = 2)
# 
# jpeg("RNAseq_volcano.jpeg", width = 6, height = 4, quality = 100, res = 300, units = "in")
# grid.arrange(plot1, plot4, plot3, plot2, ncol = 2, nrow = 2)
# dev.off()
# 
# jpeg("RNAseq_volcano_3pss_WT.jpeg", width = 3, height = 2, quality = 100, res = 300, units = "in")
# plot1
# dev.off()
# 
# jpeg("RNAseq_volcano_3p_WTKO.jpeg", width = 3, height = 2, quality = 100, res = 300, units = "in")
# ggplot(KO3P_WT3P, aes(x = log2fc,
#                       y = -log10(p)))+
#     geom_point(data = subset(KO3P_WT3P, abs(log2fc) > 2 & p < .01), 
#                color = "red") +
#     geom_point(data = subset(KO3P_WT3P, abs(log2fc) <= 2 | p > .01), 
#                color = "gold") +
#     geom_point(data = subset(KO3P_WT3P, abs(log2fc) <= 2 & p > .01), 
#                color = "black") +
#     theme_classic() +
#     labs(x = bquote(log[2]*FC~3*pRNA~WT/KO)) +
#     xlim(-10.1, 12.2) + ylim(0, 200)
# dev.off()
# 
# jpeg("RNAseq_volcano_ss_WTKO.jpeg", width = 3, height = 2, quality = 100, res = 300, units = "in")
# plot3
# dev.off()
# 
# jpeg("RNAseq_volcano_3pss_KO.jpeg", width = 3, height = 2, quality = 100, res = 300, units = "in")
# ggplot(KO3P_KOSS, aes(x = -log2fc,
#                       y = -log10(p)))+
#     geom_point(data = subset(KO3P_KOSS, abs(log2fc) > 2 & p < .01), 
#                color = "red") +
#     geom_point(data = subset(KO3P_KOSS, abs(log2fc) <= 2 | p > .01), 
#                color = "gold") +
#     geom_point(data = subset(KO3P_KOSS, abs(log2fc) <= 2 & p > .01), 
#                color = "black") +
#     theme_classic() +
#     labs(x = bquote(log[2]*FC~3*p/ssRNA~KO)) +
#     xlim(-10.1, 12.2) + ylim(0, 200)
# dev.off()

### make FoldChange-FoldChange graphs
FC_viral <- merge(KO3P_KOSS, WT3P_WTSS,
                  by = "gene")
FC_viral$log2fc.x <- -FC_viral$log2fc.x
FC_viral$log2fc.y <- -FC_viral$log2fc.y


pdf("RNAseq_l2FC_compare_KO_WT.pdf", 
    width = 4.5, 
    height = 4.5)
ggplot(FC_viral,
       aes(x = log2fc.x,
           y = log2fc.y)) +
    geom_abline(slope = 1,
                color = "grey",
                linetype = 2) +
    geom_point(shape = 1) +
    labs(x = bquote(log[2]*FC~3*p/ssRNA~KO),
         y = bquote(log[2]*FC~3*p/ssRNA~WT)) +
    theme_classic() +
    lims(x = c(-3.8, 12.1),
         y = c(-3.8, 12.1)) +
    geom_point(data = filter(FC_viral,
                             gene %in% c("IFIT1", "IFIT2", "IFIT3", "IFNB1",
                                         "OAS1", "OAS2", "OASL", "DDX58", "CCL5",
                                         "RSAD2", "ISG15", "STAT1", "MX2")), 
               color = "red") +
    geom_text_repel(data = filter(FC_viral,
                                  gene %in% c("IFIT1", "IFIT2", "IFIT3", "IFNB1",
                                              "OAS1", "OAS2", "OASL", "DDX58", "CCL5",
                                              "RSAD2", "ISG15", "STAT1", "MX2")),
                    aes(label = gene),
                    max.overlaps = 10,
                    color = "red",
                    nudge_x = 0.5) +
    theme(text = element_text(size = 16))
dev.off()

FC_genotype <- merge(KOSS_WTSS, KO3P_WT3P,
                     by = "gene")

pdf("RNAseq_l2FC_compare_SS_3P.pdf", 
    width = 4.5, 
    height = 4.5)
ggplot(FC_genotype,
       aes(x = log2fc.x,
           y = log2fc.y)) +
    geom_abline(slope = 1,
                color = "grey",
                linetype = 2) +
    geom_point(shape = 1) +
    labs(x = bquote(log[2]*FC~ssRNA~WT/KO),
         y = bquote(log[2]*FC~3*pRNA~WT/KO)) +
    theme_classic() +
    lims(x = c(-9.2, 10.6),
         y = c(-9.2, 10.6)) +
    geom_point(data = filter(FC_genotype,
                             gene %in% c("IFIT1", "IFIT2", "IFIT3", "IFNB1",
                                         "TNFRSF10D", "UBE2L6", "USP41", "DDX58", "USP18",
                                         "RSAD2", "ISG15", "STAT1", "HERC5")),
               color = "red") +
    geom_text_repel(data = filter(FC_genotype,
                                  gene %in% c("IFIT1", "IFIT2", "IFIT3", "IFNB1",
                                              "TNFRSF10D", "UBE2L6", "USP41", "DDX58", "USP18",
                                              "RSAD2", "ISG15", "STAT1", "HERC5")),
                    aes(label = gene),
                    max.overlaps = 25,
                    color = "red") +
    theme(text = element_text(size = 16))
dev.off()


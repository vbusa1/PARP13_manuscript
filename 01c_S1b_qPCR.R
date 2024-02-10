library(tidyverse)
library(pheatmap)
library(ggpubr)
library(ggrepel)
library("RColorBrewer")

setwd("/Users/Veronica/Documents/Leung/PARP13/manuscript/Review/")

data<-read.csv("qPCR_collected.csv") %>% select(-Junk)

# fill in 40 for all undetermined-- possibly lower expression than that,
# but should be sufficient filler
for(c in 1:length(data$Ct)){
  if(data$Ct[c] == "Undetermined"){
    data$Ct[c]<-40
  }
}

# KO rep7 ssRNA-treated has substantially less cDNA than the other samples;
# this makes sense because < 1 ug RNA was loaded into rev txn reaction (7.5 uL max)

# since KO rep7 and KO rep8 were grown in parallel (it's just based on rows of the plate
# to assign 7 vs 8) I am setting KO rep8 ssRNA-treated as the control for KO rep7
# 3p-ssRNA-treated data

# names of controls
ctrl_cell <- "WT"
ctrl_transcript <- c("GAPDH", "TBP")
ctrl_treatment <- "ssRNA"

summary_data<-data %>% group_by(Cell.Line, Target, Treatment, Replicate) %>%
  dplyr::summarize(CT_mean=mean(as.numeric(Ct)), CT_sd=sd(as.numeric(Ct))) %>%
  mutate(dCT_mean=0, ddCT=0) %>% 
  as.data.frame()

for(i in 1:nrow(summary_data)){
  summary_data[i,"dCT_mean"]<-summary_data[i,"CT_mean"]-
    summary_data[which(summary_data$Cell.Line == summary_data[i,"Cell.Line"] &
                      summary_data$Target %in% ctrl_transcript &
                      summary_data$Replicate == summary_data[i,"Replicate"] &
                      summary_data$Treatment == summary_data[i,"Treatment"]),"CT_mean"]
}
for(i in 1:nrow(summary_data)){
  try(summary_data[i,"ddCT"] <- summary_data[i,"dCT_mean"]-
    (summary_data[which(summary_data$Cell.Line == ctrl_cell &
     summary_data$Target == summary_data[i,"Target"] &
     summary_data$Replicate == summary_data[i,"Replicate"] &
     summary_data$Treatment == ctrl_treatment),"dCT_mean"]),
    silent = T)
}

summary_data_plot<-summary_data %>% mutate(Expression=2^-ddCT,
                                      upper_error=((2^-(ddCT+CT_sd))-Expression),
                                      lower_error=(Expression-(2^-(ddCT+CT_sd))))
#summary_data_plot<-summary_data_plot[!(summary_data_plot$Target=="GAPDH"),]

ggplot(summary_data_plot, aes(x = Treatment, 
                         ymin = Expression + lower_error, 
                         ymax = Expression + upper_error, 
                         fill = Cell.Line)) +
  facet_grid(Target~Replicate, 
             scales="free") +
  geom_bar(position = position_dodge(), 
           aes(y = Expression), 
           stat = "identity") +
  geom_errorbar(position = position_dodge(width=0.9), 
                colour = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic() +
  xlab(NULL) +
  ylab("Relative Expression\n(Normalized to GAPDH)") +
  scale_y_log10()


heatmap <- summary_data_plot %>% 
  group_by(Target, Cell.Line, Treatment) %>%
  summarise(Expression = mean(Expression))
heatmap <- heatmap %>%
  ungroup() %>%
  mutate(Treatment = paste0(Cell.Line, substr(Treatment, 1, 2))) %>%
  select(-Cell.Line) %>%
  spread(Treatment, Expression) %>%
  filter(Target != "GAPDH",
         Target != "TBP",
         Target != "CERK") %>%
  as.data.frame()
row.names(heatmap) <- heatmap$Target
heatmap <- select(heatmap, -Target) %>%
  t() %>%
  scale(center = F) %>%
  t()
heatmap <- heatmap[c("TRAILR4", "WRB", "ZIC2",
                     "AGO2",
                     "SERPING1", "PPP1R15A", "ERAP2", "IFITM1", "IFIT1", "IFIT2",
                     "ISG15", "RSAD2","RIG-I", "OAS2", "OASL", "CCL4"),] %>% t()

breaksList <- seq(min(heatmap), max(heatmap), by = 0.01)
pdf("qPCR_confirm.pdf", width = 3, height = 1.5)
pheatmap(heatmap, 
         fontsize = 50, 
         angle_col = 90,
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(length(breaksList)),
         show_rownames = T,
         breaks = breaksList,
         legend = F,
         gaps_col = c(3, 4))
dev.off()

summary_data <- data %>% group_by(Cell.Line, Target, Treatment, Replicate) %>%
  mutate(CT_mean=mean(as.numeric(Ct)), 
         CT_sd=sd(as.numeric(Ct)),
         dCT_mean=0, ddCT=0) %>% 
  as.data.frame()

for(i in 1:nrow(summary_data)){
  summary_data[i,"dCT_mean"]<-as.numeric(summary_data[i,"Ct"])-
    summary_data[which(summary_data$Cell.Line == summary_data[i,"Cell.Line"] &
                         summary_data$Target %in% ctrl_transcript &
                         summary_data$Replicate == summary_data[i,"Replicate"] &
                         summary_data$Treatment == summary_data[i,"Treatment"]),"CT_mean"][1]
}
for(i in 1:nrow(summary_data)){
  try(summary_data[i,"ddCT"] <- summary_data[i,"dCT_mean"]-
        mean(summary_data[which(summary_data$Cell.Line == ctrl_cell &
                              summary_data$Target == summary_data[i,"Target"] &
                              summary_data$Treatment == ctrl_treatment),"dCT_mean"]),
      silent = T)
}

summary_data_plot2 <- summary_data %>% 
  select(-Replicate) %>%
  group_by(Cell.Line, Target, Treatment) %>%
  mutate(Expression=2^-mean(ddCT),
         upper_error=((2^-mean(ddCT+CT_sd))-Expression),
         lower_error=(Expression-(2^-mean(ddCT+CT_sd))))%>% 
  filter(Target %in% c("TRAILR4", "WRB", "ZIC2",
                       "AGO2",
                       "SERPING1", "PPP1R15A", "ERAP2", "IFITM1", "IFIT1", "IFIT2",
                       "ISG15", "RSAD2","RIG-I", "OAS2", "OASL", "CCL4")) %>% 
  mutate(Treatment = paste0(Cell.Line, substr(Treatment, 1, 2))) %>%
  select(-Cell.Line) %>%
  ungroup() %>%
  group_by(Target) %>%
  mutate(ylim_min = min(.33, min(2^-(ddCT))),
         ylim_max = max(3, max(2^-(ddCT))))

# ggplot(summary_data_plot2,
#        aes(x = Treatment,
#            ymin = Expression + lower_error, 
#            ymax = Expression + upper_error)) +
#   facet_wrap(~Target,
#              ncol = 2) +
#   geom_bar(position = position_dodge(), 
#            aes(y = Expression), 
#            stat = "identity",
#            fill = "white", color = "black") +
#   geom_errorbar(position = position_dodge(), 
#                 width = 0.3,
#                 colour = "black") +
#   geom_point(aes(y = 2^-(ddCT)),
#              shape = 1, size = 2.5) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme_classic() +
#   xlab(NULL) +
#   ylab("Relative Expression\n(Normalized to GAPDH)") +
#   scale_y_log10()

plots <- list()
for(t in c("TRAILR4", "WRB", "ZIC2",
           "AGO2",
           "SERPING1", "PPP1R15A", "ERAP2", "IFITM1", "IFIT1", "IFIT2",
           "ISG15", "RSAD2","RIG-I", "OAS2", "OASL", "CCL4")){
  plot_target <- filter(summary_data_plot2,
                        Target == t)
  plots[[t]] <- ggplot(plot_target,
                          aes(x = Treatment,
                              ymin = Expression + lower_error, 
                              ymax = Expression + upper_error)) +
    geom_bar(position = position_dodge(), 
             aes(y = Expression), 
             stat = "identity",
             fill = "white", color = "black") +
    geom_errorbar(position = position_dodge(), 
                  width = 0.3,
                  colour = "black") +
    geom_point(aes(y = 2^-(ddCT)),
               shape = 1, size = 2.5) +
    theme_classic() +
    labs(x = NULL, y = NULL, title = t) +
    theme(axis.text.x = element_blank(),
          title = element_text(size = 50)) +
    scale_y_log10(limits = c(plot_target$ylim_min[1], 
                             plot_target$ylim_max[1]))
}

ggarrange(plotlist = plots,
          ncol = 4,
          nrow = 4)
ggsave("qPCR_barplots.pdf", width = 7, height = 7)


#### compare to RNA-seq (analysis performed for reviewers)
plots2 <- list()
RNAseq<-readRDS("/Users/veronica/Documents/Leung/PARP13/CLIP analysis/all_DESeq_set2.rds")
WT3p_KO3p <- RNAseq["KO3P_WT3P"] %>% as.data.frame() %>%
  filter(KO3P_WT3P.symbol %in% c("TNFRSF10D", #TRAILR4
                                 "WRB", "ZIC2", "AGO2", "SERPING1", "PPP1R15A", 
                                 "ERAP2", "IFITM1", "IFIT1", "IFIT2", "ISG15", 
                                 "RSAD2","DDX58", #RIG-I
                                 "OAS2", "OASL", "CCL4")) %>%
  mutate(KO3P_WT3P.log2FoldChange = -KO3P_WT3P.log2FoldChange)
WT3p_KO3p_PCR <- summary_data_plot2 %>%
  select(Treatment, Target, Expression) %>%
  unique() %>%
  filter(Treatment %in% c("WT3p", "KO3p")) %>%
  spread(Treatment, Expression) %>%
  mutate(FC = log2(KO3p/WT3p))
WT3p_KO3p$KO3P_WT3P.symbol[which(WT3p_KO3p$KO3P_WT3P.symbol == "TNFRSF10D")] <- "TRAILR4"
WT3p_KO3p$KO3P_WT3P.symbol[which(WT3p_KO3p$KO3P_WT3P.symbol == "DDX58")] <- "RIG-I"

WT3p_KO3p_PCR <- merge(WT3p_KO3p_PCR, WT3p_KO3p, 
                       by.x = "Target", by.y = "KO3P_WT3P.symbol")
plots2[[1]] <- ggplot(WT3p_KO3p_PCR,
       aes(x = KO3P_WT3P.log2FoldChange,
           y = FC)) +
  geom_smooth(method = "lm") +
  geom_point(size = 2)  +
  geom_text_repel(aes(label = Target), 
                  max.overlaps = 5, 
                  color = "grey30", size = 2.5) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, 
           label = paste("~r^2==~", 
                         signif(cor(WT3p_KO3p_PCR$FC, 
                                    WT3p_KO3p_PCR$KO3P_WT3P.log2FoldChange,
                                    use = "pairwise.complete.obs")^2, 
                                digits = 2)), parse = T) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 2.5, 
           label = paste("~p==~", 
                         signif(summary(lm(KO3P_WT3P.log2FoldChange~FC, 
                                           data = WT3p_KO3p_PCR))$coefficients[2,4], 
                                digits = 2)), parse = T) +
  labs(x = NULL,
       y = "PCR log2FC",
       title = "WT3p vs KO3p")
summary(lm(KO3P_WT3P.log2FoldChange~FC, data = WT3p_KO3p_PCR))


WTss_KOss <- RNAseq["KOSS_WTSS"] %>% as.data.frame() %>%
  filter(KOSS_WTSS.symbol %in% c("TNFRSF10D", #TRAILR4
                                 "WRB", "ZIC2", "AGO2", "SERPING1", "PPP1R15A", 
                                 "ERAP2", "IFITM1", "IFIT1", "IFIT2", "ISG15", 
                                 "RSAD2","DDX58", #RIG-I
                                 "OAS2", "OASL", "CCL4")) %>%
  mutate(KOSS_WTSS.log2FoldChange = -KOSS_WTSS.log2FoldChange)
WTss_KOss_PCR <- summary_data_plot2 %>%
  select(Treatment, Target, Expression) %>%
  unique() %>%
  filter(Treatment %in% c("WTss", "KOss")) %>%
  spread(Treatment, Expression) %>%
  mutate(FC = log2(KOss/WTss))
WTss_KOss$KOSS_WTSS.symbol[which(WTss_KOss$KOSS_WTSS.symbol == "TNFRSF10D")] <- "TRAILR4"
WTss_KOss$KOSS_WTSS.symbol[which(WTss_KOss$KOSS_WTSS.symbol == "DDX58")] <- "RIG-I"

WTss_KOss_PCR <- merge(WTss_KOss_PCR, WTss_KOss, 
                       by.x = "Target", by.y = "KOSS_WTSS.symbol")
plots2[[2]] <- ggplot(WTss_KOss_PCR,
       aes(x = KOSS_WTSS.log2FoldChange,
           y = FC)) +
  geom_smooth(method = "lm") +
  geom_point(size = 2)  +
  geom_text_repel(aes(label = Target), 
                  max.overlaps = 5, 
                  color = "grey30", size = 2.5) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, 
           label = paste("~r^2==~", 
                         signif(cor(WTss_KOss_PCR$FC, 
                                    WTss_KOss_PCR$KOSS_WTSS.log2FoldChange,
                                    use = "pairwise.complete.obs")^2, 
                                digits = 2)), parse = T) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 2.5, 
           label = paste("~p==~", 
                         signif(summary(lm(KOSS_WTSS.log2FoldChange~FC, 
                                           data = WTss_KOss_PCR))$coefficients[2,4], 
                                digits = 2)), parse = T) +
  labs(x = NULL,
       y = NULL,
       title = "WTss vs KOss")
summary(lm(KOSS_WTSS.log2FoldChange~FC, data = WTss_KOss_PCR))

WTss_WT3p <- RNAseq["WT3P_WTSS"] %>% as.data.frame() %>%
  filter(WT3P_WTSS.symbol %in% c("TNFRSF10D", #TRAILR4
                                 "WRB", "ZIC2", "AGO2", "SERPING1", "PPP1R15A", 
                                 "ERAP2", "IFITM1", "IFIT1", "IFIT2", "ISG15", 
                                 "RSAD2","DDX58", #RIG-I
                                 "OAS2", "OASL", "CCL4")) %>%
  mutate(WT3P_WTSS.log2FoldChange = -WT3P_WTSS.log2FoldChange)
WTss_WT3p_PCR <- summary_data_plot2 %>%
  select(Treatment, Target, Expression) %>%
  unique() %>%
  filter(Treatment %in% c("WTss", "WT3p")) %>%
  spread(Treatment, Expression) %>%
  mutate(FC = log2(WT3p/WTss))
WTss_WT3p$WT3P_WTSS.symbol[which(WTss_WT3p$WT3P_WTSS.symbol == "TNFRSF10D")] <- "TRAILR4"
WTss_WT3p$WT3P_WTSS.symbol[which(WTss_WT3p$WT3P_WTSS.symbol == "DDX58")] <- "RIG-I"

WTss_WT3p_PCR <- merge(WTss_WT3p_PCR, WTss_WT3p, 
                       by.x = "Target", by.y = "WT3P_WTSS.symbol")
plots2[[3]] <- ggplot(WTss_WT3p_PCR,
       aes(x = WT3P_WTSS.log2FoldChange,
           y = FC)) +
  geom_smooth(method = "lm") +
  geom_point(size = 2)  +
  geom_text_repel(aes(label = Target), 
                  max.overlaps = 5, 
                  color = "grey30", size = 2.5) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, 
           label = paste("~r^2==~", 
                         signif(cor(WTss_WT3p_PCR$FC, 
                                    WTss_WT3p_PCR$WT3P_WTSS.log2FoldChange,
                                    use = "pairwise.complete.obs")^2, 
                                digits = 2)), parse = T) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 2.5, 
           label = paste("~p==~", 
                         signif(summary(lm(WT3P_WTSS.log2FoldChange~FC, 
                                           data = WTss_WT3p_PCR))$coefficients[2,4], 
                                digits = 2)), parse = T) +
  labs(x = "RNAseq log2FC",
       y = "PCR log2FC",
       title = "WTss vs WT3p")
summary(lm(WT3P_WTSS.log2FoldChange~FC, data = WTss_WT3p_PCR))


KOss_KO3p <- RNAseq["KO3P_KOSS"] %>% as.data.frame() %>%
  filter(KO3P_KOSS.symbol %in% c("TNFRSF10D", #TRAILR4
                                 "WRB", "ZIC2", "AGO2", "SERPING1", "PPP1R15A", 
                                 "ERAP2", "IFITM1", "IFIT1", "IFIT2", "ISG15", 
                                 "RSAD2","DDX58", #RIG-I
                                 "OAS2", "OASL", "CCL4")) %>%
  mutate(KO3P_KOSS.log2FoldChange = -KO3P_KOSS.log2FoldChange)
KOss_KO3p_PCR <- summary_data_plot2 %>%
  select(Treatment, Target, Expression) %>%
  unique() %>%
  filter(Treatment %in% c("KOss", "KO3p")) %>%
  spread(Treatment, Expression) %>%
  mutate(FC = log2(KO3p/KOss))
KOss_KO3p$KO3P_KOSS.symbol[which(KOss_KO3p$KO3P_KOSS.symbol == "TNFRSF10D")] <- "TRAILR4"
KOss_KO3p$KO3P_KOSS.symbol[which(KOss_KO3p$KO3P_KOSS.symbol == "DDX58")] <- "RIG-I"

KOss_KO3p_PCR <- merge(KOss_KO3p_PCR, KOss_KO3p, 
                       by.x = "Target", by.y = "KO3P_KOSS.symbol")
plots2[[4]] <- ggplot(KOss_KO3p_PCR,
       aes(x = KO3P_KOSS.log2FoldChange,
           y = FC)) +
  geom_smooth(method = "lm") +
  geom_point(size = 2)  +
  geom_text_repel(aes(label = Target), 
                  max.overlaps = 5, 
                  color = "grey30", size = 2.5) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, 
           label = paste("~r^2==~", 
                         signif(cor(KOss_KO3p_PCR$FC, 
                                    KOss_KO3p_PCR$KO3P_KOSS.log2FoldChange,
                                    use = "pairwise.complete.obs")^2, 
                                digits = 2)), parse = T) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 2.5, 
           label = paste("~p==~", 
                         signif(summary(lm(KO3P_KOSS.log2FoldChange~FC, 
                                           data = KOss_KO3p_PCR))$coefficients[2,4], 
                                digits = 2)), parse = T) +
  labs(x = "RNAseq log2FC",
       y = NULL,
       title = "KOss vs KO3p")
summary(lm(KO3P_KOSS.log2FoldChange~FC, data = KOss_KO3p_PCR))$coefficients[2,4]


ggarrange(plotlist = plots2,
          ncol = 2,
          nrow = 2)
ggsave("qPCR_RNAseq_FC_compare.pdf")

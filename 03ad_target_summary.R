library(tidyverse)
library(VennDiagram)


##### overlap of SS and 3P treatment PARP13 targets #####
WT3P_CLIP<-read.delim("results/WT3P_peaks.tsv", header = F) %>% 
    mutate(sample = "3p") %>% filter(!grepl(",", V8))
WTSS_CLIP<-read.delim("results/WTSS_peaks.tsv", header = F) %>% 
    mutate(sample = "ss") %>% filter(!grepl(",", V8))
jpeg("CLIP_venn.jpeg", 
     width = 2, height = 2, 
     quality = 100, res = 300, 
     units = "in")
grid.newpage(); draw.pairwise.venn(area1 = length(unique(WTSS_CLIP$V8)),
                                   area2 = length(unique(WT3P_CLIP$V8)),
                                   cross.area = sum(unique(WTSS_CLIP$V8) %in% unique(WT3P_CLIP$V8)),
                                   fill = c("blue", "red"),
                                   fontfamily = rep("sans",3),
                                   alpha = rep(0.3, 2))
dev.off()


###### localization of PARP13 targets #######
# APEX data derived from Kaewsapsak et al 2017
ER_APEX<-read.csv("Other_PARP13_datasets/Kaewsapsak_2017/ER_APEX.csv")
nuc_APEX<-read.csv("Other_PARP13_datasets/Kaewsapsak_2017/nuc_APEX.csv")
cyt_APEX<-read.csv("Other_PARP13_datasets/Kaewsapsak_2017/cyt_APEX.csv")

WTSS<-WTSS_CLIP[,7:8] %>%
    mutate(treat = "SS") %>% unique() %>%
    filter(!grepl(",", V8))
WT3P<-WT3P_CLIP[,7:8] %>%
    mutate(treat = "3P") %>% unique() %>%
    filter(!grepl(",", V8))
CLIP<-rbind(WTSS, WT3P) %>%
    group_by(V8) %>%
    mutate(combined = (n() > 1))
CLIP[which(CLIP$combined == "TRUE"), "treat"] <-"both"
CLIP<-CLIP[,1:3] %>% unique()

gather_CLIP_loc<-CLIP %>% mutate(ER= V8 %in% ER_APEX$Gene.name,
                                 cyt = V8 %in% cyt_APEX$Gene.name) %>%
    select(-V7) %>% 
    unique() %>%
    mutate(localization = "")
gather_CLIP_loc[which(gather_CLIP_loc$ER == T & 
                          gather_CLIP_loc$cyt == T), "localization"] <-"cyt & ER"
gather_CLIP_loc[which(gather_CLIP_loc$ER == T & 
                          gather_CLIP_loc$cyt == F), "localization"] <-"ER"
gather_CLIP_loc[which(gather_CLIP_loc$ER == F & 
                          gather_CLIP_loc$cyt == T), "localization"] <-"cytoplasm"
gather_CLIP_loc[which(gather_CLIP_loc$ER == F & 
                          gather_CLIP_loc$cyt == F), "localization"] <-"neither"

ggplot(gather_CLIP_loc, aes(x = (treat == "3P"))) +
    geom_bar(position = "fill", aes(fill = factor(localization, levels = c("neither", "cytoplasm", "cyt & ER", "ER")))) +
    theme_classic() +
    geom_text(stat='count', aes(label=..count..), vjust=-1, y = 1) +
    labs(x = NULL, y = "Proportion", fill = "RNA\nLocalization") +
    scale_y_continuous(lim= c(0,1.1), breaks = seq(0, 1, .25)) +
    scale_x_discrete(labels = c("+ssRNA", "+3p-RNA\nonly")) +
    theme(text = element_text(size = 16))
ggsave("localization_proportion.jpeg", width = 4.5, height = 3, units = "in")


### stats
see<-read.csv("PARP13 Table of Values.csv") # just to pull list of all genes
both<- CLIP %>% filter(treat == "both") %>% 
    select(V8) %>% unlist() %>% unique()

ER_loc<-matrix(c(sum(both %in% ER_APEX$Gene.name),
                 length(both) - sum(both %in% ER_APEX$Gene.name),
                 nrow(ER_APEX) - sum(both %in% ER_APEX$Gene.name),
                 sum(!(see$genename %in% ER_APEX$Gene.name | see$genename %in% both))),
               nrow = 2,
               dimnames = list(ER = c("loc", "not loc"),
                               PARP13 = c("bound", "not bound")))
fisher.test(ER_loc, alternative = "greater")$p.value #2.961706e-132


cyt_loc<-matrix(c(sum(WT3P$V8 %in% cyt_APEX$Gene.name),
                  length(WT3P$V8) - sum(WT3P$V8 %in% cyt_APEX$Gene.name),
                  nrow(cyt_APEX) - sum(WT3P$V8 %in% cyt_APEX$Gene.name),
                  sum(!(see$genename %in% cyt_APEX$Gene.name | see$genename %in% WT3P$V8))),
                nrow = 2,
                dimnames = list(cyt = c("loc", "not loc"),
                                PARP13 = c("bound", "not bound")))
fisher.test(cyt_loc, alternative = "greater")$p.value #1.036058e-69

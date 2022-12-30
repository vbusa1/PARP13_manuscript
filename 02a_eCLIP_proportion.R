library(tidyverse)

WTSS<-read.delim("results/WTSS_peaks.tsv", header = F, stringsAsFactors = F)
WT3P<-read.delim("results/WT3P_peaks.tsv", header = F, stringsAsFactors = F)

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
    theme(text = element_text(size = 20)) +
    scale_fill_discrete(labels = c("Intron", "Noncoding exon", "5' UTR", "CDS", "3' UTR"))
ggsave("PARP13_CLIP_gene_region_binding.png", units = "in", width = 5, height = 3)

chi<-regions %>% group_by(sample, V9) %>% 
    summarise(count = n()) %>% 
    spread(V9, count) %>%
    as.data.frame()
row.names(chi)<-chi[,1]
chi<-chi[,2:ncol(chi)]
chisq.test(chi)

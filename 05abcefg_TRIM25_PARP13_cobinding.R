library(tidyverse)
library(VennDiagram)
library(nearBynding)
library(GenomicRanges)
library(BiocIO)
library(rtracklayer)
theme_set(theme_classic())

WT3P_CLIP<-read.delim("results/WT3P_peaks.tsv", header = F) %>% 
    mutate(sample = "3p") %>% filter(!grepl(",", V8))
WTSS_CLIP<-read.delim("results/WTSS_peaks.tsv", header = F) %>% 
    mutate(sample = "ss") %>% filter(!grepl(",", V8))
CLIP<-rbind(WT3P_CLIP, WTSS_CLIP) %>% select("V8", "sample") %>% unique()

see<-read.csv("PARP13 Table of Values.csv") # to get all gene names

## Figure 5A
TRIM1<-import("GSE104949_RAW/GSM2810841_Trim25_exp1_significant_clusters.gtf")
TRIM1<-TRIM1@elementMetadata@listData[["gene_name"]]
TRIM2<-import("GSE104949_RAW/GSM2810842_Trim25_exp2_significant_clusters.gtf")
TRIM2<-TRIM2@elementMetadata@listData[["gene_name"]]
TRIM3<-import("GSE104949_RAW/GSM2810843_Trim25_exp3_significant_clusters.gtf")
TRIM3<-TRIM3@elementMetadata@listData[["gene_name"]]
TRIM25<-c(TRIM1, TRIM2, TRIM3)
TRIM25<-lapply(TRIM25, function(x){
    return(strsplit(x, ":::")[[1]][2])
}) %>% unlist() %>% unique()

jpeg("TRIM25_PARP13_venn.jpeg", width = 2, height = 2, quality = 100, res = 300, units = "in")
grid.newpage(); draw.pairwise.venn(area1 = length(unique(CLIP$V8)),
                                   area2 = length(TRIM25),
                                   cross.area = sum(TRIM25 %in% CLIP$V8),
                                   fill = c("blue", "red"),
                                   scale =  T,
                                   euler.d = T,
                                   fontfamily = rep("sans",3),
                                   alpha = rep(0.3, 2))
dev.off()

## stats assessment
gene_list<-unique(c(see$genename, CLIP$V8, TRIM25))
binding<-matrix(c(sum(TRIM25 %in% CLIP$V8),
                  sum(gene_list %in% CLIP$V8) - sum(TRIM25 %in% CLIP$V8),
                  sum(gene_list %in% TRIM25) - sum(TRIM25 %in% CLIP$V8),
                  sum(!(gene_list %in% CLIP$V8 | gene_list %in% TRIM25))),
                nrow = 2,
                dimnames = list(TRIM25 = c("bound", "not bound"),
                                PARP13 = c("bound", "not bound")))
fisher.test(binding, alternative = "greater")$p.value


## Figure 5B
CLIP<-rbind(WT3P_CLIP, WTSS_CLIP) %>% 
    select("V8", "sample") %>% 
    unique() %>% 
    group_by(V8) %>%
    mutate(combined = (dplyr::n() > 1))
CLIP[which(CLIP$combined == "TRUE"), "sample"] <-"both"
CLIP<-CLIP[,1:3] %>% unique()

ggplot(CLIP, aes(x = (sample == "3p"))) +
    geom_bar(position = "fill", aes(fill = V8 %in% TRIM25)) +
    theme_classic() +
    geom_text(stat='count', aes(label=..count..), vjust=-1, y = 1) +
    labs(x = "Treatment", y = "Proportion", fill = "bound by\nPARP13") +
    scale_y_continuous(lim= c(0,1.1), breaks = seq(0, 1, .25)) +
    scale_x_discrete(labels = c("+ssRNA", "+3p-ssRNA\nonly")) +
    theme(text = element_text(size = 16))+
    scale_fill_hue(labels = c("only", "and TRIM25"), l=55)
ggsave("TRIM25_PARP13_cobinding.jpeg", width = 4.5, height = 3, units = "in")





########### run nearBynding on PARP13 and TRIM25 ###########

# list of HEK293 expressed transcripts based on RNAseq from Sun et al 2018 GEO:GSE122425
transcript_list<-readLines("isoforms_HEK293.txt") 
GenomeMappingToChainFile(genome_gtf = "Homo_sapiens.GRCh38.98.gtf", # from Ensembl
                         out_chain_name = "HEK293_exon.chain",
                         RNA_fragment = "exon",
                         transcript_list = transcript_list,
                         alignment = "hg38")
getChainChrSize(chain = "HEK293_exon.chain",
                out_chr = "HEK293_exon.size")

###### prep TRIM25 data from Choudhury et al 2017 GEO:GSE104949
# separate by F and R strand
expt1<-read.table("GSE104949_RAW/TRIM25_expt1.bed", header = F, stringsAsFactors = F)
colnames(expt1)<-c('chr','start','end','value','strand')
expt1 <- with(expt1, GRanges(chr, IRanges(start, end), strand=strand, score = value))
expt1 %>% filter(strand == "+") %>% export("GSE104949_RAW/TRIM25_expt1_F.bedGraph", "bedGraph")
expt1 %>% filter(strand == "-") %>% export("GSE104949_RAW/TRIM25_expt1_R.bedGraph", "bedGraph")
expt2<-read.table("GSE104949_RAW/TRIM25_expt2.bed", header = F, stringsAsFactors = F)
colnames(expt2)<-c('chr','start','end','value','strand')
expt2 <- with(expt2, GRanges(chr, IRanges(start, end), strand=strand, score = value))
expt2 %>% filter(strand == "+") %>% export("GSE104949_RAW/TRIM25_expt2_F.bedGraph", "bedGraph")
expt2 %>% filter(strand == "-") %>% export("GSE104949_RAW/TRIM25_expt2_R.bedGraph", "bedGraph")
expt3<-read.table("GSE104949_RAW/TRIM25_expt3.bed", header = F, stringsAsFactors = F)
colnames(expt3)<-c('chr','start','end','value','strand')
expt3 <- with(expt3, GRanges(chr, IRanges(start, end), strand=strand, score = value))
expt3 %>% filter(strand == "+") %>% export("GSE104949_RAW/TRIM25_expt3_F.bedGraph", "bedGraph")
expt3 %>% filter(strand == "-") %>% export("GSE104949_RAW/TRIM25_expt3_R.bedGraph", "bedGraph")

liftOverToExomicBG(input = c("GSE104949_RAW/TRIM25_expt1_F.bedGraph",
                             "GSE104949_RAW/TRIM25_expt1_R.bedGraph"),
                   chain = "HEK293_exon.chain",
                   chrom_size = "HEK293_exon.size",
                   output_bg = "GSE104949_RAW/TRIM25_expt1.bedGraph")
liftOverToExomicBG(input = c("GSE104949_RAW/TRIM25_expt2_F.bedGraph",
                             "GSE104949_RAW/TRIM25_expt2_R.bedGraph"),
                   chain = "HEK293_exon.chain",
                   chrom_size = "HEK293_exon.size",
                   output_bg = "GSE104949_RAW/TRIM25_expt2.bedGraph")
liftOverToExomicBG(input = c("GSE104949_RAW/TRIM25_expt3_F.bedGraph",
                             "GSE104949_RAW/TRIM25_expt3_R.bedGraph"),
                   chain = "HEK293_exon.chain",
                   chrom_size = "HEK293_exon.size",
                   output_bg = "GSE104949_RAW/TRIM25_expt3.bedGraph")


#### prep and correlate PARP13
# must first convert to hg38
chain<-import("hg19ToHg38.over.chain") # publicly available on Ensembl

##### PARP13 +3p treatment vs TRIM25 #######
PARP13<-read.table("results/WT3P_peaks.bed", header = F, stringsAsFactors = F)
colnames(PARP13)<-c('chr','start','end','value','strand')
PARP13_GRanges <- with(PARP13, GRanges(chr, IRanges(start, end), strand=strand, score = value))
# separate by F and R strand
PARP13_GRanges %>% filter(strand == "+") %>% 
    liftOver(chain) %>% unlist() %>%
    export("PARP13_hg38_F.bedGraph", "bedGraph")
PARP13_GRanges %>% filter(strand == "-") %>% 
    liftOver(chain) %>% unlist() %>%
    export("PARP13_hg38_R.bedGraph", "bedGraph")

liftOverToExomicBG(input = c("PARP13_hg38_F.bedGraph",
                             "PARP13_hg38_R.bedGraph"),
                   chain = "HEK293_exon.chain",
                   chrom_size = "HEK293_exon.size",
                   output_bg = "PARP13_hg38_liftOver.bedGraph")

## correlate
runStereogene(track_files = c("PARP13_hg38_liftOver.bedGraph","GSE104949_RAW/TRIM25_expt1.bedGraph"),
              name_config = "HEK293_exon.cfg",
              nShuffle = 10000)
runStereogene(track_files = c("PARP13_hg38_liftOver.bedGraph","GSE104949_RAW/TRIM25_expt2.bedGraph"),
              name_config = "HEK293_exon.cfg",
              nShuffle = 10000)
runStereogene(track_files = c("PARP13_hg38_liftOver.bedGraph","GSE104949_RAW/TRIM25_expt3.bedGraph"),
              name_config = "HEK293_exon.cfg",
              nShuffle = 10000)

## visualize binding
visualizeStereogene(protein_file = c("TRIM25_expt1", "TRIM25_expt2", "TRIM25_expt3"),
                    context_file = "PARP13_hg38_liftOver",
                    x_lim = c(-800, 800),
                    y_lim = c(-2, 8),
                    out_file = "TRIM25_PARP13_3p",
                    legend = F)
# visualizeStereogene(protein_file = c("TRIM25_expt1", "TRIM25_expt2", "TRIM25_expt3"),
#                     context_file = "PARP13_hg38_liftOver",
#                     x_lim = c(-800, 800),
#                     y_lim = c(-8, 8),
#                     out_file = "TRIM25_PARP13_3p_heat",
#                     heatmap = T)




##### PARP13 +ss treatment vs TRIM25 #######
PARP13_ss<-read.table("results/WTSS_peaks.bed", header = F, stringsAsFactors = F)
colnames(PARP13_ss)<-c('chr','start','end','value','strand')
PARP13_ss_GRanges <- with(PARP13_ss, GRanges(chr, IRanges(start, end), strand=strand, score = value))
# separate by F and R strand
PARP13_ss_GRanges %>% filter(strand == "+") %>% 
    liftOver(chain) %>% unlist() %>% 
    export("PARP13_ss_hg38_F.bedGraph", "bedGraph")
PARP13_ss_GRanges %>% filter(strand == "-") %>% 
    liftOver(chain) %>% unlist() %>% 
    export("PARP13_ss_hg38_R.bedGraph", "bedGraph")

liftOverToExomicBG(input = c("PARP13_ss_hg38_F.bedGraph",
                             "PARP13_ss_hg38_R.bedGraph"),
                   chain = "HEK293_exon.chain",
                   chrom_size = "HEK293_exon.size",
                   output_bg = "PARP13_ss_hg38_liftOver.bedGraph")

## correlate
write_config(name_config = "HEK293_exon.cfg",
             chrom_size = "HEK293_exon.size")
runStereogene(track_files = c("PARP13_ss_hg38_liftOver.bedGraph","GSE104949_RAW/TRIM25_expt1.bedGraph"),
              name_config = "HEK293_exon.cfg",
              nShuffle = 10000)
runStereogene(track_files = c("PARP13_ss_hg38_liftOver.bedGraph","GSE104949_RAW/TRIM25_expt2.bedGraph"),
              name_config = "HEK293_exon.cfg",
              nShuffle = 10000)
runStereogene(track_files = c("PARP13_ss_hg38_liftOver.bedGraph","GSE104949_RAW/TRIM25_expt3.bedGraph"),
              name_config = "HEK293_exon.cfg",
              nShuffle = 10000)

## visualize binding
visualizeStereogene(protein_file = c("TRIM25_expt1", "TRIM25_expt2", "TRIM25_expt3"),
                    context_file = "PARP13_ss_hg38_liftOver",
                    x_lim = c(-800, 800),
                    y_lim = c(-2, 8),
                    out_file = "TRIM25_PARP13_ss",
                    legend = F)
# visualizeStereogene(protein_file = c("TRIM25_expt1", "TRIM25_expt2", "TRIM25_expt3"),
#                     context_file = "PARP13_ss_hg38_liftOver",
#                     x_lim = c(-800, 800),
#                     y_lim = c(-8, 8),
#                     out_file = "TRIM25_PARP13_ss_heat",
#                     heatmap = T)



########### AUTOCORRELATION ###############

# autocorrelation of PARP13
#created duplicate file so could run this (don't know why doesn't like to run same thing twice)
# 3p treatment
runStereogene(track_files = c("PARP13_hg38_liftOver_dup.bedGraph","PARP13_hg38_liftOver.bedGraph"),
              name_config = "HEK293_exon.cfg",
              nShuffle = 10000)
visualizeStereogene(protein_file = "PARP13_hg38_liftOver",
                    context_file = "PARP13_hg38_liftOver_dup",
                    x_lim = c(-800, 800),
                    y_lim = c(-1, 3),
                    out_file = "PARP13_3P_autocor",
                    legend = F)
visualizeStereogene(protein_file = "PARP13_hg38_liftOver",
                    context_file = "PARP13_hg38_liftOver_dup",
                    x_lim = c(-800, 800),
                    y_lim = c(-2,100),
                    out_file = "PARP13_3P_autocor",
                    legend = F)
# ss treatment
runStereogene(track_files = c("PARP13_ss_hg38_liftOver_dup.bedGraph","PARP13_ss_hg38_liftOver.bedGraph"),
              name_config = "HEK293_exon.cfg",
              nShuffle = 10000)
visualizeStereogene(protein_file = "PARP13_ss_hg38_liftOver",
                    context_file = "PARP13_ss_hg38_liftOver_dup",
                    x_lim = c(-800, 800),
                    y_lim = c(-1, 3),
                    out_file = "PARP13_SS_autocor",
                    legend = F)
visualizeStereogene(protein_file = "PARP13_ss_hg38_liftOver",
                    context_file = "PARP13_ss_hg38_liftOver_dup",
                    x_lim = c(-800, 800),
                    y_lim = c(-2, 100),
                    out_file = "PARP13_SS_autocor",
                    legend = F)

## autocorrelation of TRIM25
runStereogene(track_files = c("GSE104949_RAW/TRIM25_expt3_dup.bedGraph","GSE104949_RAW/TRIM25_expt3.bedGraph"),
              name_config = "HEK293_exon.cfg",
              nShuffle = 10000)
visualizeStereogene(protein_file = "TRIM25_expt3",
                    context_file = "TRIM25_expt3_dup",
                    x_lim = c(-800, 800),
                    y_lim = c(-2, 10),
                    out_file = "TRIM25_expt3_autocor",
                    legend = F)
visualizeStereogene(protein_file = "TRIM25_expt3",
                    context_file = "TRIM25_expt3_dup",
                    x_lim = c(-800, 800),
                    y_lim = c(-2, 100),
                    out_file = "TRIM25_expt3_autocor",
                    legend = F)
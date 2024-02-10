library(tidyverse)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)


######## luciferase ###########
data1<-read.csv("data/20211117_luciferase.csv") %>% 
    mutate(day = 1)
data2<-read.csv("data/20211120_luciferase.csv") %>% 
    mutate(day = 2)
data3<-read.csv("data/20211210_luciferase.csv") %>% 
    mutate(day = 3)

data<-rbind(data1, data2) %>% rbind(data3)

# remove background signal (no-plasmid control)
for(cell in unique(data$cell)){
    for(treatment in unique(data$treatment)){
        for(day in unique(data$day)){
            mean_ctrl<-mean(data[which(data$cell == cell &
                                           data$treatment == treatment &
                                           data$day == day &
                                           data$plasmid == "none"), "ctrl_luc"])
            mean_luc<-mean(data[which(data$cell == cell & 
                                          data$treatment == treatment &
                                          data$day == day &
                                          data$plasmid == "none"), "Rluc"])
            data[which(data$cell == cell &
                           data$day == day &
                           data$treatment == treatment), "ctrl_luc"]<-data[which(data$cell == cell &
                                                                                     data$day == day &
                                                                                     data$treatment == treatment), "ctrl_luc"] - mean_ctrl
            data[which(data$cell == cell &
                           data$day == day &
                           data$treatment == treatment), "Rluc"]<-data[which(data$cell == cell &
                                                                                 data$day == day &
                                                                                 data$treatment == treatment), "Rluc"] - mean_luc
        }
    }
}
data <- filter(data, plasmid !="none")

process<-gather(data, "luciferase", "value", -replicate, -plasmid, -cell, -treatment, -day) %>%
    filter(luciferase == "Rluc") %>% 
    mutate(normalized = 0, ctrl_normalized = 0)
for(i in 1:nrow(process)){
    process[i, "normalized"] <- process[i, "value"] / data[which(data$cell == process[i, "cell"] &
                                                                     data$day == process[i, "day"] &
                                                                     data$treatment == process[i, "treatment"] &
                                                                     data$plasmid == process[i, "plasmid"] &
                                                                     data$replicate == process[i, "replicate"]), "ctrl_luc"]
}
for(i in 1:nrow(process)){
    process[i, "ctrl_normalized"] <- process[i, "normalized"] / mean(process[which(process$cell == process[i, "cell"] &
                                                                                       process$day == process[i, "day"] &
                                                                                       process$treatment == process[i, "treatment"] &
                                                                                       process$plasmid == "no_insert"), "normalized"])
}

ReplicateAverages <- process %>%
    filter(plasmid != "no_insert") %>%
    select(-replicate, -value, -normalized, -luciferase) %>% 
    group_by(cell, plasmid, treatment, day) %>% 
    summarise_each(list(mean))

stat <- process %>%
    lm(ctrl_normalized ~ day, data = .)

stat.test <- process %>% mutate(resid = stat$residuals) %>%
    group_by(treatment, plasmid) %>%
    t_test(resid ~ cell, data = .) %>%
    add_significance() %>%
    filter(plasmid != "no_insert")
stat.test$y.position<-2.3

ggplot(filter(process, plasmid != "no_insert"), 
       aes(x=cell,y=ctrl_normalized, shape = factor(day))) + 
    geom_line(data=ReplicateAverages, aes(group = paste0(plasmid, treatment, day)), 
              alpha = .7, linetype = "dashed") +
    facet_wrap(~factor(treatment, levels = c("ssRNA", "3p-ssRNA")) + plasmid, 
               scales = "free_x", nrow = 1) +
    geom_beeswarm(aes(color=cell), cex=5, size = 3, dodge.width = 1, alpha = .3) + 
    geom_point(data=ReplicateAverages, aes(color=cell), size=3) + 
    stat_pvalue_manual(stat.test) +
    labs(x = NULL,
         y = "Renilla vs firefly luc expression\n normalized to no insert (luc)") +
    theme(legend.title = element_blank(), 
          text = element_text(size = 14),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.background = element_blank()) +
    ylim(0, 2.4)
ggsave("pooled_luciferase.png", height = 3.5, width = 5, units = "in")



######## qPCR ###########
data1<-read.csv("data/20211118_qPCR.csv") %>% 
    mutate(day = 1)
data2<-read.csv("data/20211206_qPCR.csv") %>% 
    mutate(day = 2)
qPCR<-rbind(data1, data2)
qPCR<-qPCR %>% separate(Sample.Name,
                        c("cell", "treatment", "UTR"),
                        "_") %>%
    filter(UTR != "none")
data3<-read.csv("data/20211210_qPCR.csv")
data3<-data3 %>% separate(Sample.Name,
                        c("cell", "treatment", "UTR", "day"),
                        "_")
data3$day<-as.integer(data3$day)+2
qPCR<-rbind(qPCR, data3)

ctrl_transcript<-"ctrl"
ctrl_UTR<-"no.ins"

summary_data<-qPCR %>% mutate(dCT=0, ddCT=0)

for(i in 1:nrow(summary_data)){
    summary_data[i,"dCT"]<-summary_data[i,"Ct"]-mean(summary_data[which(summary_data$cell == summary_data[i,"cell"] & 
                                                                            summary_data$UTR == ctrl_UTR & 
                                                                            summary_data$treatment == summary_data[i,"treatment"] & 
                                                                            summary_data$day == summary_data[i,"day"] &
                                                                            summary_data$Target.Name==summary_data[i,"Target.Name"]),"Ct"])
}

for(i in 1:nrow(summary_data)){
    summary_data[i,"ddCT"]<-summary_data[i,"dCT"]-mean(summary_data[which(summary_data$cell==summary_data[i,"cell"] & 
                                                                              summary_data$UTR==summary_data[i,"UTR"] &
                                                                              summary_data$day == summary_data[i,"day"] &
                                                                              summary_data$treatment==summary_data[i,"treatment"] &
                                                                              summary_data$Target.Name== ctrl_transcript), "dCT"])
}

summary_data<-summary_data %>% mutate(Expression=2^-ddCT) %>%
    filter(Target.Name != ctrl_transcript,
           UTR != ctrl_UTR)

ReplicateAverages <- summary_data[,c(1:3,6,9)] %>%
    group_by(cell, UTR, treatment, day) %>% 
    summarise_each(list(mean))

stat <- summary_data %>%
    lm(Expression ~ day, data = .)

stat.test <- summary_data %>% mutate(resid = stat$residuals) %>%
    group_by(treatment, UTR) %>%
    t_test(resid ~ cell, data = .) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance()
stat.test$y.position<-2.7

ggplot(summary_data, 
       aes(x=cell, y=Expression)) + 
    geom_line(data=ReplicateAverages, aes(group = paste0(UTR, treatment, day)), 
              alpha = .7, linetype = "dashed") +
    facet_wrap(~factor(treatment, levels = c("ssRNA", "3p-ssRNA")) + UTR, 
               scales = "free_x", nrow = 1) +
    geom_beeswarm(aes(color=cell), cex=5, size = 3, dodge.width = 1, alpha = .3) + 
    geom_point(data=ReplicateAverages, aes(color=cell), size=3) + 
    stat_pvalue_manual(stat.test) +
    theme_classic() +
    labs(x = NULL,
         y = "Renilla vs firefly luc expression\n normalized to no insert (qPCR)") +
    theme(legend.title = element_blank(), 
          text = element_text(size = 14),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.background = element_blank()) +
    ylim(0, 2.81) 
ggsave("pooled_qPCR.png", height = 3.5, width = 5, units = "in")



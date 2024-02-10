library(tidyverse)
library(ggpubr)
library(rstatix)

qPCR<-read.csv("data/20211217.csv")

ggplot(qPCR, aes(x=Target.Name, y=Ct))+geom_point() +
    facet_grid(~Sample.Name) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ctrl_transcript<-"GAPDH"
ctrl_cell<-"WT"

summary_data<-qPCR %>% mutate(dCT=0, ddCT=0)

for(i in 1:nrow(summary_data)){
    summary_data[i,"dCT"]<-summary_data[i,"Ct"] -
        mean(summary_data[which(summary_data$Sample.Name == ctrl_cell & 
                                    summary_data$Target.Name==summary_data[i,"Target.Name"]),"Ct"])
}

for(i in 1:nrow(summary_data)){
    summary_data[i,"ddCT"]<-summary_data[i,"dCT"] -
        mean(summary_data[which(summary_data$Sample.Name==summary_data[i,"Sample.Name"] &
                                    summary_data$Target.Name== ctrl_transcript), "dCT"])
}


summary_data<-summary_data %>% mutate(Expression=2^-ddCT)

plot<- summary_data %>% group_by(Target.Name, Sample.Name) %>%
    dplyr::summarize(upper_error = sd(Expression)/length(n) + mean(Expression),
                     lower_error = mean(Expression) - sd(Expression)/length(n),
                     Expression = mean(Expression)) %>%
    unique() %>%
    filter(Target.Name != "GAPDH")

# ### Create Figure 1E
# plot_TRAILR4 <- plot %>% filter(Target.Name == "TRAILR4")
# ggplot(plot_TRAILR4, aes(x = factor(Sample.Name, levels = c("WT", "KO")))) +
#     geom_bar(position=position_dodge(),
#              aes(y=Expression),
#              stat="identity",
#              fill = "white",
#              color = "black") +
#     geom_errorbar(aes(ymin = lower_error,
#                       ymax = upper_error),
#                   width = .5,
#                   color = "black") +
#     theme_classic() +
#     xlab(NULL) +
#     ylab("TRAILR4 expression\nnormalized to GAPDH") +
#     theme(legend.position = "none") +
#     scale_x_discrete(labels = c("WT", "PARP13\nKO"))
# ggsave("TRAILR4_expression_HEK_qPCR.jpeg", 
#        units = "in", height = 2, width = 2)

##### replicate data
qPCR2<-read.csv("data/20220106.csv")

ggplot(qPCR2, aes(x=Target.Name, y=Ct))+geom_point() +
    facet_grid(~Sample.Name) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
qPCR2<-filter(qPCR2, Target.Name != "UBA7") 
# need to get new primer for UBA7, has high Ct and two melting temps

summary_data2<-qPCR2 %>% mutate(dCT=0, ddCT=0)

for(i in 1:nrow(summary_data2)){
    summary_data2[i,"dCT"]<-summary_data2[i,"Ct"] -
        mean(summary_data2[which(summary_data2$Sample.Name == ctrl_cell & 
                                     summary_data2$Target.Name==summary_data2[i,"Target.Name"]),"Ct"])
}

for(i in 1:nrow(summary_data2)){
    summary_data2[i,"ddCT"]<-summary_data2[i,"dCT"] -
        mean(summary_data2[which(summary_data2$Sample.Name==summary_data2[i,"Sample.Name"] &
                                     summary_data2$Target.Name== ctrl_transcript), "dCT"])
}


summary_data2<-summary_data2 %>% mutate(Expression=2^-ddCT)

plot2<- summary_data2 %>% group_by(Target.Name, Sample.Name) %>%
    dplyr::summarize(upper_error = sd(Expression)/length(n) + mean(Expression),
                     lower_error = mean(Expression) - sd(Expression)/length(n),
                     Expression = mean(Expression)) %>%
    unique() %>%
    filter(Target.Name != "GAPDH")

### Create Figure 4B
plot_ISGylation<-rbind(plot, plot2) %>%
    filter(Target.Name != "TRAILR4")
ggplot(plot_ISGylation, aes(x = Target.Name,
                            group = Sample.Name,
                            fill = factor(Sample.Name, levels = c("KO", "WT")))) +
    geom_bar(position=position_dodge(),
             aes(y=Expression),
             stat="identity") +
    geom_errorbar(aes(ymin = lower_error,
                      ymax = upper_error),
                  width = .5,
                  color = "black",
                  position = position_dodge(width = 0.9)) +
    theme_classic() +
    xlab(NULL) +
    ylab("Expression\nnormalized to GAPDH") +
    theme(legend.title = element_blank(),
          text = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, h = 1)) +
    scale_fill_hue(labels = c("PARP13\nKO", "WT"), l=55)
ggsave("ISG_expression_HEK_qPCR.jpeg", 
       units = "in", height = 3.5, width = 6)


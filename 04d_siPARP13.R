library(tidyverse)
library(ggpubr)

###### 12 hr timepoints in HEK
qPCR1<-read.csv("data/20220308_HEK.csv")
qPCR2<-read.csv("data/20220328siPARP13.csv")
qPCR<-rbind(qPCR1, qPCR2)

# ggplot(qPCR, aes(x=Target.Name, y=Ct))+geom_point() +
#     facet_grid(~Sample.Name) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
# poor RNA quality in 48 hr timepoint
qPCR<-filter(qPCR, Sample.Name != 48 & Target.Name != "GAPDH")

ctrl_transcript<-"TBP"
ctrl_cell<-"0"

summary_data<-qPCR %>% mutate(dCT=0, ddCT=0)

for(i in 1:nrow(summary_data)){
    summary_data[i,"dCT"]<-summary_data[i,"Ct"]-mean(summary_data[which(summary_data$Sample.Name == ctrl_cell & 
                                                                            summary_data$Target.Name==summary_data[i,"Target.Name"]),"Ct"])
}

for(i in 1:nrow(summary_data)){
    summary_data[i,"ddCT"]<-summary_data[i,"dCT"]-mean(summary_data[which(summary_data$Sample.Name==summary_data[i,"Sample.Name"] &
                                                                              summary_data$Target.Name== ctrl_transcript), "dCT"])
}

summary_data<-summary_data %>% mutate(Expression=2^-ddCT)

plot<- summary_data %>% group_by(Target.Name, Sample.Name) %>%
    dplyr::summarize(upper_error = sd(Expression)/length(n) + mean(Expression),
                     lower_error = mean(Expression) - sd(Expression)/length(n),
                     Expression = mean(Expression)) %>%
    unique()

ggplot(plot, aes(x = Sample.Name, y = Expression, 
                 color = factor(Target.Name, 
                                levels = c("RIGI","TBP", "USP41", 
                                           "USP18", "TRIM25", "ISG15", 
                                           "PARP13.1", "UBE2L6")))) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = lower_error,
                      ymax = upper_error,
                      group = Sample.Name),
                  width = 0.7) +
    theme_classic() +
    scale_y_log10(limits = c(0.1, 2)) +
    xlab("Hours post-siPARP13 treatment") +
    ylab("Expression normalized\nto TBP and 0 hr") +
    theme(legend.title = element_blank(),
          text = element_text(size = 14),
          axis.text = element_text(size = 12)) +
    scale_color_hue(l=55)
ggsave("siPARP13_HEK_qPCR.jpeg", 
       units = "in", height = 3, width = 7)
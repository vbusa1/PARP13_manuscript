library(tidyverse)

# analyze qPCR
qPCR<-read.csv("data/20220405_P13IP.csv")
qPCR2<-read.csv("data/20220408_PARP13IP.csv")

ctrl_transcript<-"TBP"
ctrl_sample <- "input_1"

qPCR<-filter(qPCR, Sample.Name %in% c("IP_1", "input_1"))
qPCR$Ct<-as.numeric(qPCR$Ct)

summary_data<-qPCR %>% mutate(dCT=0, ddCT=0)
for(i in 1:nrow(summary_data)){
    summary_data[i,"dCT"]<-summary_data[i,"Ct"]-mean(summary_data[which(summary_data$Sample.Name == ctrl_sample & 
                                                                            summary_data$Target.Name==summary_data[i,"Target.Name"]),"Ct"])
}
for(i in 1:nrow(summary_data)){
    summary_data[i,"ddCT"]<-summary_data[i,"dCT"]-mean(summary_data[which(summary_data$Sample.Name==summary_data[i,"Sample.Name"] & 
                                                                              summary_data$Target.Name== ctrl_transcript), "dCT"])
}

summary_data<-summary_data %>% mutate(Expression=2^-ddCT) %>%
    filter(Target.Name != ctrl_transcript,
           Sample.Name != ctrl_sample)

plot<- summary_data %>% group_by(Target.Name) %>%
    dplyr::summarize(upper_error = sd(Expression)/length(n) + mean(Expression),
                     lower_error = mean(Expression) - sd(Expression)/length(n),
                     Expression = mean(Expression)) %>%
    unique()

qPCR2<-qPCR2 %>% filter(Target.Name != "ERI3") # high variability and high Ct

summary_data2<-qPCR2 %>% mutate(dCT=0, ddCT=0)
for(i in 1:nrow(summary_data2)){
    summary_data2[i,"dCT"]<-summary_data2[i,"Ct"]-mean(summary_data2[which(summary_data2$Sample.Name == ctrl_sample & 
                                                                            summary_data2$Target.Name==summary_data2[i,"Target.Name"]),"Ct"])
}
for(i in 1:nrow(summary_data2)){
    summary_data2[i,"ddCT"]<-summary_data2[i,"dCT"]-mean(summary_data2[which(summary_data2$Sample.Name==summary_data2[i,"Sample.Name"] & 
                                                                              summary_data2$Target.Name== ctrl_transcript), "dCT"])
}

summary_data2<-summary_data2 %>% mutate(Expression=2^-ddCT) %>%
    filter(Target.Name != ctrl_transcript,
           Sample.Name != ctrl_sample)

plot2<- summary_data2 %>% group_by(Target.Name) %>%
    dplyr::summarize(upper_error = sd(Expression)/length(n) + mean(Expression),
                     lower_error = mean(Expression) - sd(Expression)/length(n),
                     Expression = mean(Expression)) %>%
    unique()
plot<-rbind(plot, plot2)

ggplot(plot, aes(x = factor(Target.Name, 
                            levels = c("PARP13.1", "PARP13.2", "ATP6AP1", "EIF1", 
                                       "SLC3A2", "STT3A", "TMED3","XBP1", "HDAC3")),
                 ymin = lower_error,
                 ymax = upper_error)) +
    geom_bar(aes(y=Expression),
             stat="identity",
             color = "black",
             fill = "white",
             size = .7) +
    geom_errorbar(width=0.5,
                  size = .8) +
    theme_classic() +
    theme(text = element_text(size = 14),
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 65, h = 1)) +
    xlab(NULL) +
    ylab("IP enrichment relative to TBP") +
    scale_y_log10()
ggsave("PARP13_IP.jpeg", width = 3.5, height = 3.5, units = "in")



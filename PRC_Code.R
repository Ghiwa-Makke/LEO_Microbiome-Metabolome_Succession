#library:
library(vegan)
library(readr)

# RP
auc <- read.csv("Input_files/RP_PRC_input.csv", sep = ",")

env <- data.frame(Treatment = auc$Treatment, Time = auc$Time)
Treatment <- as.factor(env[ ,1])
Time <- as.factor(env[ ,2])
response.df <- auc[, -c(1,2,3,4)]
prc.obj <- prc(response.df, Treatment, Time)

summary(prc.obj)
plot(prc.obj)
##            0       1       2
### 1 -1.404e-16 -0.4836 -0.8452

#species_scores <- as.data.frame(scores(prc.obj, display = "species"))

#write.csv(species_scores, "PRC_feature_scores_RP.csv", row.names = TRUE)

# HILIC:
auc <- read.csv("Input_files/HILIC_PRC_input.csv", sep = ",")

env <- data.frame(Treatment = auc$Treatment, Time = auc$Time)
Treatment <- as.factor(env[ ,1])
Time <- as.factor(env[ ,2])
response.df <- auc[, -c(1,2,3,4)]
prc.obj <- prc(response.df, Treatment, Time)

summary(prc.obj)
plot(prc.obj)

#          0       1       2
# 1 -3.985e-17 -0.4061 -0.9085


#species_scores <- as.data.frame(scores(prc.obj, display = "species"))
#write.csv(species_scores, "PRC_feature_scores_RP.csv", row.names = TRUE)

# ASV:
auc <- read.csv("Input_files/OTU_PRC_input.csv", sep = ",")

env <- data.frame(Treatment = auc$Treatment, Time = auc$Time)
Treatment <- as.factor(env[ ,1])
Time <- as.factor(env[ ,2])
response.df <- auc[, -c(1,2,3,4)]
prc.obj <- prc(response.df, Treatment, Time)

summary(prc.obj)
plot(prc.obj)

#           0      1      2
# 1 -1.137e-16 0.4862 0.3853

#species_scores <- as.data.frame(scores(prc.obj, display = "species"))
#write.csv(species_scores, "PRC_feature_scores_ASV.csv", row.names = TRUE)


# Plotting 
Names <- c("Bare", "Biocrust", "Moss")
prc_values_RP <- c(0, 0.4836, 0.8452)
prc_values_HILIC <- c(0, 0.4061, 0.9085)
prc_values_ASV <- c(0, 0.4862, 0.3853)


prc_plot_data <- data.frame(Names, prc_values_HILIC, prc_values_RP, prc_values_ASV)
# Reshape the data using gather()
library(ggplot2)
library(tidyr)
prc_plot_gathered <- gather(prc_plot_data, key="Run", value="value", -Names)

# Create the plot
PRC_plot <- ggplot(prc_plot_gathered, aes(x=Names, y=value, color=Run, group=Run)) +
  geom_point(size = 4) +
  #geom_line() +
  #geom_smooth(method="loess") +
  #geom_line(curve_type="quadratic") +
  geom_smooth(linewidth=0.8)+
  xlab("Types") +
  ylab("Values") +
  scale_color_manual(values=c("blue", "red", "grey4"), labels=c("Microbial", "HILIC", "RP"))+
  ggtitle("PRC plot")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_minimal()
PRC_plot
ggsave("PRC_Plot.png", PRC_plot)





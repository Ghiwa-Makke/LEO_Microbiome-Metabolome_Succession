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
##            0       1       2
### 1 6.168e-16 1.876 3.278


# ASV:
auc <- read.csv("Input_files/ASV_PRC_input_.csv", sep = ",")

env <- data.frame(Treatment = auc$Treatment, Time = auc$Time)
Treatment <- as.factor(env[ ,1])
Time <- as.factor(env[ ,2])
response.df <- auc[, -c(1,2,3,4)]
prc.obj <- prc(response.df, Treatment, Time)

summary(prc.obj)

#           0      1      2
# 1 3.249e-16 1.431 0.9092


# ITS:
auc <- read.csv("Input_files/ITS_PRC_input.csv", sep = ",")

env <- data.frame(Treatment = auc$Treatment, Time = auc$Time)
Treatment <- as.factor(env[ ,1])
Time <- as.factor(env[ ,2])
response.df <- auc[, -c(1,2,3,4)]
prc.obj <- prc(response.df, Treatment, Time)

summary(prc.obj)
#0      1       2
#1 4.374e-16 0.19 0.1601

# Plotting 
Names <- c("Bare", "Biocrust", "Moss")
prc_values_RP <- c(0, 1.876, 3.278)
prc_values_ASV <- c(0, 1.431, 0.9092)
prc_values_fungal <- c(0, 0.19, 0.1601)

prc_plot_data <- data.frame(Names, prc_values_RP, prc_values_ASV, prc_values_fungal)
# Reshape the data using gather()
library(ggplot2)
library(tidyr)
prc_plot_gathered <- gather(prc_plot_data, key="Run", value="value", -Names)

# Create the plot
PRC_plot <- ggplot(prc_plot_gathered, aes(x=Names, y=value, color=Run, group=Run)) +
  geom_point(size = 3) +
  geom_smooth(linewidth=0.8)+
  xlab("Types") +
  ylab("Values") +
  #scale_color_manual(values=c("blue", "red", "grey4", "orange"), labels=c("Microbial", "HILIC", "RP", "MAGs"))+
  ggtitle("PRC plot")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_minimal()
PRC_plot
ggsave("PRC_Plot.png", PRC_plot)





#### Microbiome-Metabolome Interactions Reshape Ecosystem Properties During Biocrust-Moss Succession
### CCA RDA analysis

## 1. Libraries: 
library(tidyverse)
library("vegan")
library("ggplot2")
library("ggpmisc")

# set path variables
project_dir <- getwd()

## 2. import tables
#data table should be normalized with the samples in rows and features in columns
data <- read.csv("Input_files/LEO_RP_ML_Input_paretoLMM.csv", row.names = 1, sep = ",", check.names = F)

auc <- data %>%
  select(-Type, -Slope, -Location)

auc_shifted <- auc + abs(min(auc))


explan <- data %>%
  select(Type, Slope, Location)

#auc <- read.table(file = "clipboard", sep = "\t", header = TRUE, row.names = 1)
##3. Setting variables: 
Type = as.factor(explan[ , 1]) #this sets the type as factor and gets the info from expan table - 1st column
Slope = as.factor(explan[ , 2]) #this sets the slope as factor and gets the info from expan table - 2nd column
Location = as.factor(explan[ , 3]) #this sets the location as factor and gets the info from expan table - 3nd column

CCA_1 <- cca(auc_shifted ~ Type + Slope + Location, data = explan)
summary(CCA_1)
CCA_2 <- cca(auc_shifted ~ Type * Slope + Location  , data = explan)

##5. getting the statistical tables: 
CCA_stat_1 <- anova.cca(CCA_1, by="terms") # here by terms means look at the order of the factors # what is usually used
CCA_stat_2 <- anova.cca(CCA_2, by="terms")
write.csv(CCA_stat_2, "CCA_Stat_2.csv", row.names=TRUE)

##4. getting the CCA tables:
RDA_1 <- cca(auc ~ Type + Slope + Location, data=explan) 
RDA_2 <- cca(auc ~ Type * Slope * Location, data=explan)
summary(RDA_1)

##5. getting the statistical tables: 
RDA_stat_1 <- anova.cca(RDA_1, by="terms") # here by terms means look at the order of the factors # what is usually used
RDA_stat_2 <- anova.cca(RDA_2, by="terms")
write.csv(CCA_stat_2, "CCA_Stat_2.csv", row.names=TRUE)


# Numerical Variables: ----------------------------------------------------
metadata <- read.csv("env_metadata.csv")


## All together: -----------------------------------------------------------
explan <- metadata %>%
  column_to_rownames(var = "SampleID") %>%
  select(MossCoverage, Biocrust, Shannon_Diversity, WaterContent, pH, C_N)%>%
  rename(Moss_Coverage = MossCoverage) %>%
  rename(Biocrust_Coverage = Biocrust) %>%
  rename(Microbial_Diversity = Shannon_Diversity)%>%
  rename(Water_Content = WaterContent)

Moss_Coverage <- as.numeric(explan[ ,1])
Biocrust_Coverage <- as.numeric(explan[ ,2])
Microbial_Diversity <- as.numeric(explan[ ,3])
Water_Content <- as.numeric(explan[ ,4])
pH <- as.numeric(explan[ ,5])
C_N <- as.numeric(explan[ ,6])

explan_stand<-decostand(explan, method="standardize")

library(car)
auc_vector <- auc[, 1] 
vif(lm(auc_vector ~ Moss_Coverage + Biocrust_Coverage + Microbial_Diversity + Water_Content + pH + C_N, data = explan_stand))



RDA_1 <- rda(auc ~ Moss_Coverage + Biocrust_Coverage + Microbial_Diversity + Water_Content + pH + C_N, data = explan_stand)
Stat_1 <- anova.cca(RDA_1, by = "term")
RDA_2 <- rda(auc ~ Moss_Coverage * Biocrust_Coverage * Microbial_Diversity * Water_Content * pH * C_N, data = explan_stand)
stat_2 <- anova.cca(RDA_2, by = "term")

write.csv(Stat_1, "RDA_Stat_1_RP.csv", row.names=TRUE)


## Biotic Factors ----------------------------------------------------------
explan_biotic <- metadata %>%
  column_to_rownames(var = "SampleID") %>%
  select(MossCoverage, Biocrust, Shannon_Diversity)%>%
  rename(Moss_Coverage = MossCoverage) %>%
  rename(Biocrust_Coverage = Biocrust) %>%
  rename(Microbial_Diversity = Shannon_Diversity)

Moss_Coverage <- as.numeric(explan_biotic[ ,1])
Biocrust_Coverage <- as.numeric(explan_biotic[ ,2])
Microbial_Diversity <- as.numeric(explan_biotic[ ,3])

explan_biotic_stand<-decostand(explan_biotic, method="standardize")


biotic_RDA_1 <- rda(auc ~ Moss_Coverage + Biocrust_Coverage + Microbial_Diversity, data = explan_biotic_stand)
biotic_Stat_1 <- anova.cca(biotic_RDA_1, by = "term")
biotic_RDA_2 <- rda(auc ~ Moss_Coverage * Biocrust_Coverage * Microbial_Diversity, data = explan_biotic_stand)
biotic_stat_2 <- anova.cca(biotic_RDA_2, by = "term")

write.csv(biotic_stat_2, "RDA_Stat_2_RP_Biotic.csv", row.names=TRUE)



## Abiotic Factors ---------------------------------------------------------
explan_abiotic <- metadata %>%
  column_to_rownames(var = "SampleID") %>%
  select(WaterContent, pH, C_N) %>%
  rename(Water_Content = WaterContent)

Water_Content <- as.numeric(explan_abiotic[ ,1])
pH <- as.numeric(explan_abiotic[ ,2])
C_N <- as.numeric(explan_abiotic[ ,3])

explan_abiotic_stand<-decostand(explan_abiotic, method="standardize")

abiotic_RDA_1 <- rda(auc ~ Water_Content + pH + C_N, data = explan_abiotic_stand)
abiotic_Stat_1 <- anova.cca(abiotic_RDA_1, by = "term")
abiotic_RDA_2 <- rda(auc ~ Water_Content * pH * C_N, data = explan_abiotic_stand)
abiotic_stat_2 <- anova.cca(abiotic_RDA_2, by = "term")

write.csv(abiotic_stat_2, "RDA_Stat_2_RP_abiotic.csv", row.names=TRUE)



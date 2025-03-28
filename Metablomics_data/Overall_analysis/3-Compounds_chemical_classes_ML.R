### LEO 2022 Surface Samples Metabolomics - RP 
### 5. Differential Analysis
### LMM - classes 

# 1. Importing Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
}) 


# 2. Import data

project_dir <- getwd()
project_name <- 'LEO_RP'

# Load compounds table

Bare <- file.path('ML/Bare_ML_ft.csv')
Bare <- read_csv(Bare)
Bare<-Bare %>%
  rename(FeatureID = Bare_ML_ft)

Biocrust <- file.path('ML/Biocrust_ML_ft.csv')
Biocrust <- read_csv(Biocrust)
Biocrust<-Biocrust %>%
  rename(FeatureID = Biocrust_ML_ft)

Moss <- file.path('ML/Moss_ML_ft.csv')
Moss <- read_csv(Moss)
Moss<-Moss %>%
  rename(FeatureID = Moss_ML_ft)

classification <- read.csv("Input_Files/RP_Compound_Inchikey_classifications.csv") %>%
  rename(Superclass = Final_Superclass)
classification$Superclass[is.na(classification$Superclass)] <- "Unclassified"

Bare_class <- merge(Bare, classification, by= "FeatureID")
Biocrust_class <- merge(Biocrust, classification, by= "FeatureID")
Moss_class <- merge(Moss, classification, by= "FeatureID")

#so that I can use the facet wrap function to put them all in one figure



# 3. Counts ---------------------------------------------------------------

#Bare_Biocrust
Bare_.count<- Bare_class %>%
  select(FeatureID, Superclass) %>%
  count(Superclass)

# Biocrust
Biocrust_.count<- Biocrust_class %>%
  select(FeatureID, Superclass) %>%
  count(Superclass)
#Moss
Moss.count<- Moss_class %>%
  select(FeatureID, Superclass) %>%
  count(Superclass)


#Next I add labels for each isolate and combine them into one data file 
#so that I can use the facet wrap function to put them all in one figure
Bare_.count$Comp<- 'Bare'
Biocrust_.count$Comp<- 'Biocrust'
Moss.count$Comp<- 'Moss'

class.rank.all<- rbind(Bare_.count, Biocrust_.count, Moss.count)
class.rank.all$Comp <- factor(class.rank.all$Comp, levels = c('Bare', 'Biocrust', 'Moss'))


class.plot<- ggplot(class.rank.all)+
  geom_bar(aes(y=Superclass, x=n), stat="identity", alpha=0.7)+
  facet_wrap(~Comp)+
  #scale_fill_manual(values = c("lightsteelblue","tomato2"))+
  theme_classic()+xlab("Number of Features")+ylab("Classification")+
  labs(title="Class ranks")
class.plot

figure_file <- file.path(figures_dir, 'Class_rank_abundance_LMM.png')
ggsave(figure_file, class.plot, width=25, height=10, unit='cm')


# Combine data frames
Bare_class$Source<- 'Bare'
Biocrust_class$Source<- 'Biocrust'
Moss_class$Source<- 'Moss'

combined_df_metabo <- rbind(Bare_class, Biocrust_class, Moss_class)

# Count the number of compounds per classification for each sample type
count_df_metabo <- combined_df_metabo %>%
  group_by(Source, Superclass) %>%
  summarise(NumFeatures = n()) %>%
  ungroup()

# Load required packages



# Create the bubble plot with different colors for each sample type
Metabolite_class <- ggplot(count_df_metabo, aes(x = Source, y = Superclass)) +
  geom_point(aes(size = NumFeatures, fill = Source), shape = 21, color = "black", alpha = 0.7) +
  scale_size_continuous(range = c(3, 8), breaks = seq(1, max(count_df_metabo$NumFeatures), by = 2)) +
  scale_fill_manual(values = c("Bare" = '#3d1d0a', "Biocrust" = '#a8581b', "Moss" = '#116311')) +
  theme_minimal(base_size = 8) +
  labs(
    #title = "Number of Important Features by Sample Type and Classification",
    x = "",
    y = "Compound Classification",
    size = "Number of Features",
    fill = "Sample Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.box = "horizontal",
    legend.text = element_text(size = 8)
  )+ guides(fill = guide_legend(ncol = 1)) 

Metabolite_class

# Save the combined plot
ggsave("output_figures/ML_Bubble_Plot_Metabo.svg", Metabolite_class, dpi = 300, height = 2.5, width = 6)

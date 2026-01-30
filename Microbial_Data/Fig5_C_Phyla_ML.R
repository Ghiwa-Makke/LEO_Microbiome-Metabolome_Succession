### LEO 2022 Surface Samples Metabolomics - 16S
### Machine Learning discriminatiry Features 
### ML- Phyla 

# 1. Importing Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
}) 


# 2. Import data
Bare <- read_csv('16S_data/ML_Bare.csv')
Bare<-Bare %>%
  rename(ASV_ID = Bare)

Biocrust <- read_csv('16S_data/ML_Biocrust.csv')
Biocrust<-Biocrust %>%
  rename(ASV_ID = Biocrust)

Moss <- read_csv('16S_data/ML_Moss.csv')
Moss<-Moss %>%
  rename(ASV_ID = Moss)

classification <- read.csv("16S_data/Tax_data.csv")
classification$Phylum[is.na(classification$Phylum)] <- "Other"

Bare_phylum <- merge(Bare, classification, by= "ASV_ID")
Biocrust_phylum <- merge(Biocrust, classification, by= "ASV_ID")
Moss_phylum <- merge(Moss, classification, by= "ASV_ID")


# Combine data frames
Bare_phylum$Source<- 'Bare'
Biocrust_phylum$Source<- 'Biocrust'
Moss_phylum$Source<- 'Moss'

combined_df <- rbind(Bare_phylum, Biocrust_phylum, Moss_phylum)

# Count the number of compounds per classification for each sample type
combined_df <- combined_df %>%
  group_by(Source, Phylum) %>%
  summarise(NumFeatures = n()) %>%
  ungroup()

# Create the bubble plot with different colors for each sample type
Micro_Phyla <- ggplot(combined_df, aes(x = Source, y = Phylum)) +
  geom_point(aes(size = NumFeatures, fill = Source), shape = 21, color = "black", alpha = 0.7) +
  scale_size_continuous(range = c(3, 8), breaks = seq(1, max(combined_df$NumFeatures), by = 2)) +
  scale_fill_manual(values = c("Bare" = '#3d1d0a', "Biocrust" = '#a8581b', "Moss" = '#116311')) +
  theme_minimal(base_size = 8) +
  labs(
    #title = "Number of Important Features by Sample Type and Classification",
    x = "",
    y = "Phyla",
    size = "Number of ASVs",
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

Micro_Phyla

# Save the combined plot
ggsave("output_figures/ML_Bubble_Plot_Micro.svg", Metabolite_class, dpi = 300, height = 2.5, width = 6)

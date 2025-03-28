# Assuming your matrix is called 'feature_matrix', Compound_Classification dataframe is called 'compound_classification', and metadata dataframe is called 'metadata'
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(readxl)
# 2. Import data

project_dir <- getwd()

RP_norm <- read.csv("Input_files/RP_norm_pareto.csv", row.names = 1, check.names = F)%>%
  t()%>%
  as.data.frame()
  
RP_stand <- scale(RP_norm) %>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = "SampleID")

HILIC_norm <- read.csv("Input_files/HILIC_norm_pareto.csv", row.names = 1, check.names = F)%>%
  t()%>%
  as.data.frame()

HILIC_stand <- scale(HILIC_norm) %>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = "SampleID")

Norm_full_matrix <- merge(RP_stand,HILIC_stand, by = "SampleID")%>%
  column_to_rownames(var = "SampleID")

# Load compounds table
norm_matrix_stand <- as.data.frame(t(Norm_full_matrix))%>% 
  rownames_to_column(var = "FeatureID")

# compounds_table <- file.path(tables_dir, 'gap_filled_compounds_table.csv')
# compounds_table <- read_csv(compounds_table)

# Import metadata and fix names
metadata <- read.csv("Input_files/fixed_metadata.csv")

metadata <- metadata %>%
  arrange(Type)

RP_classification <- read.csv("Input_Files/RP_Compound_Inchikey_classifications.csv")%>%
  select(FeatureID, Final_Superclass)%>%
  rename(Superclass = Final_Superclass)
HILIC_classification <- read.csv("Input_Files/HILIC_Compound_Inchikey_classifications.csv") %>%
  select(FeatureID, Final_Superclass)%>%
  rename(Superclass = Final_Superclass)

classification <- rbind(RP_classification, HILIC_classification) %>%
  distinct()

classification_count2 <- classification %>%
  select(FeatureID, Superclass)%>%
  distinct_all() %>%
  filter(!is.na(Superclass))


# Replace NA values in 'Class' column with "Other"
classification$Superclass[is.na(classification$Superclass)] <- "Other" 

# norm_matrix <- as.data.frame(t(norm_matrix))%>%
#   rownames_to_column(var = "FeatureID")


joined_data_classes <- left_join(norm_matrix_stand, classification, by = "FeatureID")%>%
  pivot_longer(cols = 2:55, names_to = 'SampleName', values_to = 'intensities')

Class_annotation_count <- joined_data_classes %>%
  select(FeatureID, Superclass)%>%
  distinct_all() %>%
  filter(Superclass != "Other")

Class_annotation_count2 <- joined_data_classes %>%
  select(FeatureID, Superclass)%>%
  distinct_all() %>%
  filter(Superclass == "Other")

# Group by sample and class, and calculate the average intensity for each class per sample
class_intensity <- joined_data_classes %>%
  group_by(SampleName, Superclass) %>%
  summarise(avg_intensity = mean(intensities))

class_matrix <- class_intensity %>%
  pivot_wider(names_from = SampleName, values_from = avg_intensity)%>%
  column_to_rownames(var = "Superclass")

# Step 3: Make a heatmap showing the change in intensities of the classes for the different sample types


# Heatmap of Classes: 
# Initialize graphical device
dev.off()

# Set sampleID as row.names to annotate heatmap
col_annot <- metadata %>% 
  select(SampleName, Type)

mapcolor <- colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)[100:1]

annot_colors <- list(
  Type = c(Bare = "#3d1d0a", Biocrust = "#a8581b", Moss = "#116311"))

row_annotations <- rownames(class_matrix)
matrix <- as.data.frame(t(class_matrix))%>%
  rownames_to_column("SampleName")%>%
  right_join(col_annot, by="SampleName")%>%
  arrange(Type)%>%
  select(-Type)%>%
  #select(-SampleID)%>%
  column_to_rownames(var = "SampleName")
new_matrix <- as.data.frame(t(matrix))

col_annot <- metadata %>% 
  select(SampleName, Type) %>%
  column_to_rownames(var = "SampleName")
figure_file <- file.path('Output_figures', 'Compound_Superclass_Heatmap_stand2.png')
png(figure_file, width = 2200, height = 1000, res = 300)
p <- pheatmap(new_matrix,
              clustering_distance_rows = 'correlation',
              cluster_cols = F,
              scale = 'row',
              border_color = NA,
              annotation_col = col_annot,
              annotation_colors = annot_colors,
              color = mapcolor,
              gaps_col = c(18,36),
              #cutree_rows = 3,
              show_colnames = F,
              main = 'Compound Superclass - Heatmap',
              fontsize_row = 8,  # Adjust fontsize for row annotation
              fontsize_col = 8,  # Adjust fontsize for column annotation
              fontsize = 8
)

gridExtra::grid.arrange(p$gtable, vp=viewport(width=0.9, height=1))
dev.off()




## Pathways ####
RP_Kegg_pathways <- read.csv("Input_files/RP_Compound_Kegg_pathways.csv")
HILIC_Kegg_pathways <- read.csv("Input_files/HILIC_Compound_Kegg_pathways.csv")
KEGG_db <- read_xlsx("Input_files/KEGG_patheway_Classes_all.xlsx")

RP_Kegg_pathways_all <- left_join(RP_Kegg_pathways, KEGG_db, by = "KEGG_pathway")%>%
  select(FeatureID, KEGG_pathway, Pathway_Subclass, KEGG_id)
HILIC_Kegg_pathways_all <- left_join(HILIC_Kegg_pathways, KEGG_db, by = "KEGG_pathway")%>%
  select(FeatureID, KEGG_pathway, Pathway_Subclass, KEGG_id)

Kegg_pathways_all <- rbind(RP_Kegg_pathways_all, HILIC_Kegg_pathways_all) 
# Step 1: Group features based on Class categories
# Replace NA values in 'Class' column with "Other"
Kegg_pathways_all$KEGG_pathway[is.na(Kegg_pathways_all$KEGG_pathway)] <- "No Kegg annotation"


# Step 2: Make an average of the intensities of each class categories per sample
# Join the feature matrix with Compound_Classification dataframe to get the class information for each feature

joined_data <- left_join(norm_matrix_stand, Kegg_pathways_all, by = "FeatureID")%>%
  pivot_longer(cols = 2:55, names_to = 'SampleName', values_to = 'intensities')

Kegg_count <- joined_data %>%
  select(FeatureID, KEGG_id)%>%
  distinct_all() %>%
  filter(!is.na(KEGG_id))

Kegg_count2 <- joined_data %>%
  select(FeatureID, KEGG_id)%>%
  distinct_all() %>%
  filter(is.na(KEGG_id))

# Group by sample and class, and calculate the average intensity for each class per sample
pathway_intensity <- joined_data %>%
  group_by(SampleName, Pathway_Subclass) %>%
  summarise(avg_intensity = mean(intensities))
pathway_intensity$Pathway_Subclass[is.na(pathway_intensity$Pathway_Subclass)] <- "No Kegg annotation"  

pathway_matrix <- pathway_intensity %>%
  pivot_wider(names_from = SampleName, values_from = avg_intensity)%>%
  column_to_rownames(var = "Pathway_Subclass")


# Step 3: Make a heatmap showing the change in intensities of the classes for the different sample types


# Heatmap of Classes: 
# Initialize graphical device
dev.off()

# Set sampleID as row.names to annotate heatmap
col_annot <- metadata %>% 
  select(SampleName, Type)

mapcolor <- colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)[100:1]

annot_colors <- list(
  Type = c(Bare = "#3d1d0a", Biocrust = "#a8581b", Moss = "#116311"))

row_annotations <- rownames(pathway_matrix)
matrix <- as.data.frame(t(pathway_matrix))%>%
  rownames_to_column("SampleName")%>%
  right_join(col_annot, by="SampleName")%>%
  arrange(Type)%>%
  select(-Type)%>%
  column_to_rownames(var = "SampleName")
new_matrix <- as.data.frame(t(matrix))

col_annot <- metadata %>% 
  select(SampleName, Type) %>%
  column_to_rownames(var = "SampleName")
figure_file <- file.path('Output_figures', 'Pathway_Classes_heatmap_allFeatures.png')
png(figure_file, width = 2200, height = 1000, res = 300)
p <- pheatmap(new_matrix,
         clustering_distance_rows = 'correlation',
         cluster_cols = F,
         scale = 'row',
         annotation_col = col_annot,
         annotation_colors = annot_colors,
         color = mapcolor,
         gaps_col = c(18,36),
         #cutree_rows = 5,
         main = 'KEGG Pathway Classes - Heatmap',
         fontsize_row = 8,  # Adjust fontsize for row annotation
         fontsize_col = 8,  # Adjust fontsize for column annotation
         fontsize = 8,
         border_color = NA,
         show_colnames = F
)
gridExtra::grid.arrange(p$gtable, vp=viewport(width=0.9, height=1))
dev.off()


#gridExtra::grid.arrange(p$gtable, vp=viewport(width=0.9, height=1))


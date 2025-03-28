#### LEO - Surface samples 2022 - 16S data - Picrust visualization
#### Microbiome-Metabolome Interactions Reshape Ecosystem Properties During Biocrust-Moss Succession 

# Load libraries: 
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)


# Load your metadata and dataset
metadata <- read_delim("sample_metadata_2.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 

metadata_filtered <- metadata %>%
  filter(Type != "blank")

# data(kegg_abundance)
kegg_abundance <- ko2kegg_abundance("output/KO_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv") 

# Keep only the columns corresponding to the filtered sample IDs
# Assuming the sample IDs are in column names
kegg_abundance_filtered <- kegg_abundance %>%
  select(all_of(metadata$sample_name))

# pathway pca -------------------------------------------------------------
source(file = 'Pathway_PCA_function.R')
Plot1 <- pathway_pca(abundance = kegg_abundance_filtered , metadata = metadata_filtered, group = "Type", colors = c("Bare" = "#3d1d0a", "Biocrust" = "#a8581b", "Moss" = "#116311"))
Plot1


# Pathway daa -------------------------------------------------------------


## Bare vs Biocrust --------------------------------------------------------
# Define the pair of groups
group1 <- "Bare"
group2 <- "Biocrust"

# Subset metadata and dataset for these two groups
subset_metadata <- metadata_filtered[metadata_filtered$Type %in% c(group1, group2), ]
subset_dataset <- kegg_abundance_filtered %>%
  select(all_of(subset_metadata$sample_name))


# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df_BvC <- pathway_daa(abundance = subset_dataset, metadata = subset_metadata, group = "Type", daa_method = "ALDEx2", select = NULL, reference = NULL) 

# Filter results for ALDEx2_Wilcoxon rank test method
# Please check the unique(daa_results_df$method) and choose one
unique(daa_results_df_BvC$method)
daa_sub_method_results_df_BvC <- daa_results_df_BvC[daa_results_df_BvC$method == "ALDEx2_Wilcoxon rank test", ]

## splitting into chuncks
rows_per_chunk <- 100

# Create a sequence of indices for splitting
indices <- split(seq_len(nrow(daa_sub_method_results_df_BvC)), ceiling(seq_len(nrow(daa_sub_method_results_df_BvC))/rows_per_chunk))

# Split the dataframe using the indices
chunked_list <- lapply(indices, function(idx) daa_sub_method_results_df_BvC[idx, , drop = FALSE])

# Check the number of chunks
length(chunked_list)

# Access the chunks
daa_sub_method_results_df_BvC_1 <- chunked_list[[1]]
daa_sub_method_results_df_BvC_2 <- chunked_list[[2]]
daa_sub_method_results_df_BvC_3 <- chunked_list[[3]]

# Annotate pathway results using KO to KEGG conversion
daa_annotated_sub_method_results_df_BvC_1 <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df_BvC_1, ko_to_kegg = TRUE)

# Annotate pathway results using KO to KEGG conversion
daa_annotated_sub_method_results_df_BvC1 <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df_BvC_1, ko_to_kegg = TRUE)

daa_annotated_sub_method_results_df_BvC2 <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df_BvC_2, ko_to_kegg = TRUE)

daa_annotated_sub_method_results_df_BvC3 <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df_BvC_3, ko_to_kegg = TRUE)

daa_annotated_sub_method_results_df_BvC_full <- rbind(daa_annotated_sub_method_results_df_BvC1, daa_annotated_sub_method_results_df_BvC2, daa_annotated_sub_method_results_df_BvC3)


### Keep only Metabolism pathways
daa_annotated_BvC_soil <- daa_annotated_sub_method_results_df_BvC_full %>%
  filter(grepl("^Metabolism", pathway_class))

## Remove NA pathways
daa_annotated_BvC_soil <- daa_annotated_BvC_soil[!is.na(daa_annotated_BvC_soil$pathway_name),]

## round the pvalue adjusted to 5 decimal places
daa_annotated_BvC_soil$p_adjust <- round(daa_annotated_BvC_soil$p_adjust,5)

# identifies the top 20 features in daa_annotated_BvC_soil with the lowest adjusted p-values (p_adjust)
low_p_feature <- daa_annotated_BvC_soil[order(daa_annotated_BvC_soil$p_adjust), ]$feature[1:20]

source('pathway_heatmap_functon.R')

subset_dataset <- subset_dataset %>%
  rownames_to_column(var = "feature")

# Filter features with p < 0.05
feature_with_p_0.05 <- daa_annotated_BvC_soil %>% 
  filter(p_adjust < 0.05)

library(ggplot2)
library(ggh4x)

daa_annotated_BvC_soil <- daa_annotated_BvC_soil %>%
  mutate(pathway_class = sub("^Metabolism;", "", pathway_class))

heatmap_BC <- pathway_heatmap(
  abundance = subset_dataset %>% 
    right_join(
      daa_annotated_BvC_soil %>%
        select(all_of(c("feature","pathway_class"))),
      by = "feature"
    ) %>% 
    filter(feature %in% feature_with_p_0.05$feature) %>% 
    select(-"feature") %>% 
    column_to_rownames("pathway_class"),
  metadata = subset_metadata2, 
  group = "Type", 
  colors = c("Bare" = "#3d1d0a", "Biocrust" = "#a8581b95")
)
heatmap_BC



## Biocrust_Moss -----------------------------------------------------------
# Define the pair of groups
group1 <- "Biocrust"
group2 <- "Moss"

# Subset metadata and dataset for these two groups
subset_metadata2 <- metadata_filtered[metadata_filtered$Type %in% c(group1, group2), ]
subset_dataset2 <- kegg_abundance_filtered %>%
  select(all_of(subset_metadata2$sample_name))

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df_CvM <- pathway_daa(abundance = subset_dataset2, metadata = subset_metadata2, group = "Type", daa_method = "ALDEx2", select = NULL, reference = NULL) 

# Filter results for ALDEx2_Welch's t test method
# Please check the unique(daa_results_df$method) and choose one
unique(daa_results_df_CvM$method)
daa_sub_method_results_df_CvM <- daa_results_df_CvM[daa_results_df_CvM$method == "ALDEx2_Wilcoxon rank test", ]


## splitting into chuncks
rows_per_chunk <- 100

# Create a sequence of indices for splitting
indices <- split(seq_len(nrow(daa_sub_method_results_df_CvM)), ceiling(seq_len(nrow(daa_sub_method_results_df_CvM))/rows_per_chunk))

# Split the dataframe using the indices
chunked_list <- lapply(indices, function(idx) daa_sub_method_results_df_CvM[idx, , drop = FALSE])

# Check the number of chunks
length(chunked_list)

# Access the chunks
daa_sub_method_results_df_CvM_1 <- chunked_list[[1]]
daa_sub_method_results_df_CvM_2 <- chunked_list[[2]]
daa_sub_method_results_df_CvM_3 <- chunked_list[[3]]

# Annotate pathway results using KO to KEGG conversion
daa_annotated_sub_method_results_df_CvM_1 <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df_CvM_1, ko_to_kegg = TRUE)

daa_annotated_sub_method_results_df_CvM_2 <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df_CvM_2, ko_to_kegg = TRUE)

daa_annotated_sub_method_results_df_CvM_3 <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df_CvM_3, ko_to_kegg = TRUE)

daa_annotated_sub_method_results_df_CvM_full <- rbind( daa_annotated_sub_method_results_df_CvM_2, daa_annotated_sub_method_results_df_CvM_3)

## keep only metabolism pathways
daa_annotated_CvM_soil <- daa_annotated_sub_method_results_df_CvM_full %>%
  filter(grepl("^Metabolism", pathway_class))


## remove NA pathways
daa_annotated_CvM_soil <- daa_annotated_CvM_soil[!is.na(daa_annotated_CvM_soil$pathway_name),]

## round pvalue to 5 digits
daa_annotated_CvM_soil$p_adjust <- round(daa_annotated_CvM_soil$p_adjust,5)

low_p_feature <- daa_annotated_CvM_soil[order(daa_annotated_CvM_soil$p_adjust), ]$feature[1:20]


subset_dataset2 <- subset_dataset2 %>%
  rownames_to_column(var = "feature")

# Filter features with p < 0.05
feature_with_p_0.05_2 <- daa_annotated_CvM_soil %>% 
  filter(p_adjust < 0.05)
library(ggplot2)
library(ggh4x)

# Create the heatmap pathway_name

daa_annotated_CvM_soil <- daa_annotated_CvM_soil %>%
  mutate(pathway_class = sub("^Metabolism;", "", pathway_class))

# Create the heatmap pathway_class
source('pathway_heatmap_functon.R')
heatmap_CM <- pathway_heatmap(
  abundance = subset_dataset2 %>% 
    right_join(
      daa_annotated_CvM_soil %>%
        select(all_of(c("feature","pathway_class"))),
      by = "feature"
    ) %>% 
    filter(feature %in% feature_with_p_0.05_2$feature) %>% 
    select(-"feature") %>% 
    column_to_rownames("pathway_class"),
  metadata = subset_metadata2, 
  group = "Type", 
  colors = c("Biocrust" = "#a8581b95", "Moss" = "#11631195")
)
heatmap_CM


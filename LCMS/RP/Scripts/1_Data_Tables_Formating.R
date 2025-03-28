### LEO 2022 Surface Samples Metabolomics - RP
### 1. Data Editing 

# Load Libraries ----------------------------------------------------------
suppressPackageStartupMessages({  
  library(tidyverse)
  library(readxl)
  source('Scripts/functions_cdis_exploration_1.R')
}) 

# Function to add spaces between chemical elements and quantities, but keep two-letter elements intact
add_spaces_to_formula <- function(formula) {
  # Add spaces between numbers and the next element
  spaced_formula <- gsub("([0-9])([A-Z])", "\\1 \\2", formula)
  # Add spaces between elements but avoid separating two-letter elements
  spaced_formula <- gsub("([A-Z][a-z]*)(?=[A-Z])", "\\1 ", spaced_formula, perl = TRUE)
  return(spaced_formula)
}


# Import Data Tables ------------------------------------------------------
#set path variables
project_dir <- getwd()

# Give a name for your project so the results are all grouped together in a directory with the same name
project_name <- 'LEO_RP'
figures_dir <- file.path(project_dir, paste0(project_name, '_output_figures'))
tables_dir <- file.path(project_dir,  paste0(project_name, '_output_tables'))

# Create output directories
dir.create(figures_dir, showWarnings = FALSE)
dir.create(tables_dir, showWarnings = FALSE)

# The data for in this script is the **Compounds table** table exported from **Compound Discoverer**. For this script to work export and excel table with the columns: **Name**, **Formula**, **Annot. deltaMass [ppm]**, **Calc. MW**, **RT [min]**, **MS2**, **Area**, **Gap Status** and **Gap Fill Status**.
# the **metadata table** can be exported directly from the **Samples tab** in **Compound Discoverer** by copying to a blank text file.

# Import data tables
cd_results_file <- file.path(project_dir, 'Input_files/LEO_RP_Compound_Annotation_Table.xlsx')
cd_results_table <- read_xlsx(cd_results_file)


# Select columns needed for downstream analysis
cd_results_table <- cd_results_table %>%
  arrange(desc(`Calc. MW`)) %>% 
  mutate(Feature_ID = paste0('Feature',formatC(n():0001, width = 4, flag = '0'))) %>%
  select(FeatureID, Feature_ID, Final_name, Final_formula, `Calc. MW`, contains('Annotation source'), contains('Results'), contains('Area:'), contains('Gap Status:'), contains('Gap Fill Status:'), Annotation_Confidence) %>% 
  
  # Differentiate between features that share the same name using "peak#" at the end of the name
  group_by(Final_name) %>% 
  add_count(Final_name) %>% 
  
  # Create variable with names for plotting (useful in following scripts)
  mutate(name4plot = ifelse(is.na(Final_name), FeatureID, ifelse(n == 1, Final_name, paste0(Final_name, '-isomer', n():1)))) %>% 
  select(-n) %>% 
  ungroup() %>%
  #rename(Formula = Final_formula) %>%
  rename(Name = Final_name)

# Apply the function to the Formula column
cd_results_table$Formula <- sapply(cd_results_table$Final_formula, add_spaces_to_formula)

table_file <- file.path(tables_dir, 'RP_cd_features_table.csv')
write_csv(cd_results_table, table_file)

# Import metadata and fix names
metadata_file <- file.path(project_dir, 'Input_files/LEO_RP_metadata.csv')
metadata <- read_csv(metadata_file)

# Select only the useful columns and fix column names
metadata$SampleID <- gsub("^RP_Neg_|\\.raw \\(F\\d{1,3}\\)$", "", metadata$SampleID)
metadata$SampleName <- gsub("^RP_Neg_|\\.raw \\(F\\d{1,3}\\)$", "", metadata$SampleName)


metadata <- metadata %>%  
  mutate(SampleID = str_remove(SampleID, 'Area: '),
         SampleID = str_remove(SampleID, '.raw'), 
         SampleID = str_remove(SampleID, 'RP_')) %>%
  mutate(SampleName = str_remove(SampleName, 'Area: '),
         SampleName = str_remove(SampleName, '.raw'), 
         SampleName = str_remove(SampleName, 'RP_')) %>%
filter(!grepl("Blank|Pooled|Not", SampleName, ignore.case = TRUE))

table_file <- file.path(tables_dir, 'fixed_metadata.csv')
write_csv(metadata, table_file)

# Thermodynamics ----------------------------------------------------------
# Split formula column into elemental counts
cd_results_table <- separate_formula(cd_results_table)
# Calculate ratios and thermodynamic indices 
cd_results_table <- calc_ratios_n_idxs(cd_results_table)
# Calculate classes
cd_results_table <- calc_classes(cd_results_table)


# Gap-Filled Table - Needed for Statistical Analysis ----------------------
# Gather area under the curve (AUC) values per sample
compounds_table <- cd_results_table %>% 
  select(-contains('Gap Status:'), -contains('Gap Fill Status:')) %>% 
  pivot_longer(contains('Area:'), names_to = 'SampleID', values_to = 'AUC') %>% 
  filter(AUC > 0) %>% 
  mutate(SampleID = str_remove(SampleID, 'Area: '),
         SampleID = str_remove(SampleID, '.raw.*'),
         SampleID = str_remove(SampleID, 'RP_'))%>%
  filter(!grepl("Blank|Pooled|S55", SampleID, ignore.case = TRUE)) %>%
  left_join(metadata, by = "SampleID")

# Save gap-filled table to be used in Statistical Analysis
table_file <- file.path(tables_dir, 'gap_filled_compounds_table.csv')
write_csv(compounds_table, table_file)


# Non-gap filled table ----------------------------------------------------
label = FALSE

# Gather "Gap Status" for filtering 
if(label == FALSE){
  gap_status <- cd_results_table %>% 
    select(FeatureID, contains('Gap Status:')) %>% 
    gather(contains('Gap Status:'), key = 'SampleID', value = 'gap_status')
  gap_status$SampleID <- str_remove(gap_status$SampleID, 'Gap Status: ')
  gap_status$SampleID <- str_remove(gap_status$SampleID, '.raw.*')
  gap_status$SampleID <- str_remove(gap_status$SampleID, 'RP_') 
  gap_status <- gap_status %>%
  filter(!grepl("Blank|Pooled|S55", SampleID, ignore.case = TRUE))
  
  ## Filtering
  
  compounds_table <- left_join(compounds_table, gap_status, by = c('FeatureID', 'SampleID')) %>% 
    filter(gap_status != 'Full gap')
  
  # Plotting types of gaps detected
  
  gap_fill_status <- cd_results_table%>% 
    select(FeatureID, Feature_ID, Name, contains('Gap Fill Status:')) %>% 
    pivot_longer(contains('Gap Fill Status:'), names_to = 'SampleID', values_to = 'gap_fill') %>% 
    mutate(SampleID = str_remove(SampleID, 'Gap Fill Status: '),
           SampleID = str_remove(SampleID, '.raw.*'),
           SampleID = str_remove(SampleID, 'RP_')) %>% 
    mutate(gap_fill_type = case_when(gap_fill == 32 ~ 'Filled by spectrum noise',
                                     gap_fill == 128 ~ 'Filled by re-detected peak',
                                     gap_fill == 0 ~ 'No gap to fill',
                                     gap_fill == 16 ~ 'Filled by simulated peak',
                                     gap_fill == 8 ~ 'Filled by trace area',
                                     gap_fill == 64 ~ 'Filled by matching ion')) %>% 
    mutate(gap_fill_type = factor(gap_fill_type, levels = c('No gap to fill', 'Filled by re-detected peak', 
                                                            'Filled by simulated peak', 'Filled by matching ion',
                                                            'Filled by trace area', 'Filled by spectrum noise'))) %>%
    filter(!grepl("Blank|Pooled|S55", SampleID, ignore.case = TRUE))
}

write.csv(gap_fill_status, file = "LEO_RP_output_tables/gap_fill_status.csv")
non_gap_cmpds <- merge(gap_fill_status, gap_status, by = c("FeatureID", "SampleID"))

remove_full_gaps <- non_gap_cmpds %>%
  filter(gap_status != 'Full gap')
remove_real_gaps <- non_gap_cmpds%>%
  filter(gap_fill_type !='Filled by spectrum noise')%>%
  filter(gap_fill_type !='Filled by trace area')
  
## Filtering

compounds_table_non_gap <- left_join(compounds_table, non_gap_cmpds, by = c('FeatureID', 'Feature_ID', 'SampleID', "Name", "gap_status")) %>% 
  filter(gap_fill_type !='Filled by spectrum noise')%>%
  filter(gap_fill_type !='Filled by trace area')


# Add metadata information
#compounds_table_non_gap <- left_join(compounds_table_non_gap, metadata, by = 'SampleID')
table_file <- file.path(tables_dir, 'compounds_table_non_real_gap_filled.csv')
write_csv(compounds_table_non_gap, table_file)


## make zero matrix 
zero_auc_table <- compounds_table_non_gap %>% 
  select(FeatureID, SampleName, AUC) %>% 
  pivot_wider(names_from = 'SampleName', values_from = 'AUC', values_fill = 0) %>%  
  column_to_rownames(var = 'FeatureID')

write.csv(zero_auc_table, "LEO_RP_output_tables/RP_zero_auc_table.csv")
library(compositions)
metabolite_data_clr <- as.data.frame(clr(t(zero_auc_table)))
write.csv(metabolite_data_clr, "LEO_RP_output_tables/RP_CLR_gap.csv")



## separate by type
Bare_matrix <- compounds_table_non_gap %>%
  select(Feature_ID, SampleName, AUC, Type) %>%
  filter(Type == "Bare") %>%
  select(-Type) %>%
  pivot_wider(names_from = 'SampleName', values_from = 'AUC', values_fill = 0) %>%  
  column_to_rownames(var = 'Feature_ID')

Bare_matrix_CLR <- as.data.frame(clr(t(Bare_matrix))) %>%
  t() %>%
  as.data.frame()

# Step 1: Find the minimum value of the CLR data
min_clr <- min(Bare_matrix_CLR)

# Step 2: Shift the data to ensure all values are >= 1
shifted_data <- (Bare_matrix_CLR - min_clr) + 1

write.csv(shifted_data, "LEO_RP_output_tables/tables_per_type/Bare_clr_rescaled.csv")


Biocrust_matrix <- compounds_table_non_gap %>%
  select(Feature_ID, SampleName, AUC, Type) %>%
  filter(Type == "Biocrust") %>%
  select(-Type) %>%
  pivot_wider(names_from = 'SampleName', values_from = 'AUC', values_fill = 0) %>%  
  column_to_rownames(var = 'Feature_ID')

Biocrust_matrix_CLR <- as.data.frame(clr(t(Biocrust_matrix))) %>%
  t() %>%
  as.data.frame()

# Step 1: Find the minimum value of the CLR data
min_clr <- min(Biocrust_matrix_CLR)

# Step 2: Shift the data to ensure all values are >= 1
shifted_data <- (Biocrust_matrix_CLR - min_clr) + 1

write.csv(shifted_data, "LEO_RP_output_tables/tables_per_type/Biocrust_clr_rescaled.csv")

Moss_matrix <- compounds_table_non_gap %>%
  select(Feature_ID, SampleName, AUC, Type) %>%
  filter(Type == "Moss") %>%
  select(-Type) %>%
  pivot_wider(names_from = 'SampleName', values_from = 'AUC', values_fill = 0) %>%  
  column_to_rownames(var = 'Feature_ID')

Moss_matrix_CLR <- as.data.frame(clr(t(Moss_matrix))) %>%
  t() %>%
  as.data.frame()

# Step 1: Find the minimum value of the CLR data
min_clr <- min(Moss_matrix_CLR)

# Step 2: Shift the data to ensure all values are >= 1
shifted_data <- (Moss_matrix_CLR - min_clr) + 1

write.csv(shifted_data, "LEO_RP_output_tables/tables_per_type/Moss_clr_rescaled.csv")

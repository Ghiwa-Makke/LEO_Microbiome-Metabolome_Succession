##### Microbiome-Metabolome Interactions Reshape Ecosystem Properties During Biocrust-Moss Succession
##### LEO - surface samples - LC - HILIC


# Data Reformatting -------------------------------------------------------

# Load Libraries --
suppressPackageStartupMessages({  
  library(tidyverse)
  library(readxl)
  library(vegan)
  library(rstatix)
  library(ggpubr)
  library(cowplot)
  source('Custom_Functions.R')
}) 


#set path variables
project_dir <- getwd()

# Give a name for your project so the results are all grouped together in a directory with the same name
project_name <- 'HILIC_LEO'

# Import data tables
cd_results_table <- read_xlsx('Input_files/HILIC_Compounds_Table_final.xlsx')

# Select columns needed for downstream analysis
cd_results_table <- cd_results_table %>%
  arrange(desc(`Calc. MW`)) %>% 
  select(FeatureID, Final_name, Final_formula, `Calc. MW`, contains('Annotation source'), contains('Results'), contains('Area:'), contains('Gap Status:'), contains('Gap Fill Status:'), Annotation_Confidence) %>% 
  
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

write_csv(cd_results_table, 'output_tables/HILIC_cd_features_table.csv')

# Import metadata and fix names
metadata <- read_csv('Input_files/HILIC_Matadata.csv')

# Select only the useful columns and fix column names
metadata$SampleID <- gsub("^HILIC_Neg_|\\.raw \\(F\\d{1,3}\\)$", "", metadata$SampleID)
metadata$SampleName <- gsub("^HILIC_Neg_|\\.raw \\(F\\d{1,3}\\)$", "", metadata$SampleName)


metadata <- metadata %>%  
  mutate(SampleID = str_remove(SampleID, 'Area: '),
         SampleID = str_remove(SampleID, '.raw'), 
         SampleID = str_remove(SampleID, 'HILIC_')) %>%
  mutate(SampleName = str_remove(SampleName, 'Area: '),
         SampleName = str_remove(SampleName, '.raw'), 
         SampleName = str_remove(SampleName, 'HILIC_')) %>%
  filter(!grepl("Blank|Pooled|Not", SampleName, ignore.case = TRUE))

write_csv(metadata, 'output_tables/fixed_metadata.csv')


# Split formula column into elemental counts
cd_results_table <- separate_formula(cd_results_table)
# Calculate ratios and thermodynamic indices 
cd_results_table <- calc_ratios_n_idxs(cd_results_table)
# Calculate classes
cd_results_table <- calc_classes(cd_results_table)

# Gap-Filled Table - Needed for Statistical Analysis -
# Gather area under the curve (AUC) values per sample
compounds_table <- cd_results_table %>% 
  select(-contains('Gap Status:'), -contains('Gap Fill Status:')) %>% 
  pivot_longer(contains('Area:'), names_to = 'SampleID', values_to = 'AUC') %>% 
  filter(AUC > 0) %>% 
  mutate(SampleID = str_remove(SampleID, 'Area: '),
         SampleID = str_remove(SampleID, '.raw.*'),
         SampleID = str_remove(SampleID, 'HILIC_'))%>%
  filter(!grepl("Blank|Pooled|S55", SampleID, ignore.case = TRUE)) %>%
  left_join(metadata, by = "SampleID")

# Save gap-filled table to be used in Statistical Analysis
write_csv(compounds_table, 'output_tables/gap_filled_compounds_table.csv')


# Non-gap filled table ---
label = FALSE

# Gather "Gap Status" for filtering 
if(label == FALSE){
  gap_status <- cd_results_table %>% 
    select(FeatureID, contains('Gap Status:')) %>% 
    gather(contains('Gap Status:'), key = 'SampleID', value = 'gap_status')
  gap_status$SampleID <- str_remove(gap_status$SampleID, 'Gap Status: ')
  gap_status$SampleID <- str_remove(gap_status$SampleID, '.raw.*')
  gap_status$SampleID <- str_remove(gap_status$SampleID, 'HILIC_') 
  gap_status <- gap_status %>%
    filter(!grepl("Blank|Pooled|S55", SampleID, ignore.case = TRUE))
  
  ## Filtering
  
  compounds_table <- left_join(compounds_table, gap_status, by = c('FeatureID', 'SampleID')) %>% 
    filter(gap_status != 'Full gap')
  
  # Plotting types of gaps detected
  
  gap_fill_status <- cd_results_table%>% 
    select(FeatureID, Name, contains('Gap Fill Status:')) %>% 
    pivot_longer(contains('Gap Fill Status:'), names_to = 'SampleID', values_to = 'gap_fill') %>% 
    mutate(SampleID = str_remove(SampleID, 'Gap Fill Status: '),
           SampleID = str_remove(SampleID, '.raw.*'),
           SampleID = str_remove(SampleID, 'HILIC_')) %>% 
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

write.csv(gap_fill_status, file = "output_tables/gap_fill_status.csv")
non_gap_cmpds <- merge(gap_fill_status, gap_status, by = c("FeatureID", "SampleID"))

remove_full_gaps <- non_gap_cmpds %>%
  filter(gap_status != 'Full gap')
remove_real_gaps <- non_gap_cmpds%>%
  filter(gap_fill_type !='Filled by spectrum noise')%>%
  filter(gap_fill_type !='Filled by trace area')

## Filtering

compounds_table_non_gap <- left_join(compounds_table, non_gap_cmpds, by = c('FeatureID', 'SampleID', "Name", "gap_status")) %>% 
  filter(gap_fill_type !='Filled by spectrum noise')%>%
  filter(gap_fill_type !='Filled by trace area')

# Add metadata information
#compounds_table_non_gap <- left_join(compounds_table_non_gap, metadata, by = 'SampleID')
write_csv(compounds_table_non_gap, 'output_tables/compounds_table_non_real_gap_filled.csv')

# 2. Normalization --------------------------------------------------------

# Load compounds table
compounds_table <- read_csv('output_tables/gap_filled_compounds_table.csv')

# Import metadata and fix names
metadata <- read_csv('output_tables/fixed_metadata.csv')

# Create a new tibble with the AUC per each mass from each sample
auc_table <- compounds_table %>% 
  select(FeatureID, SampleName, AUC) %>% 
  pivot_wider(names_from = 'SampleName', values_from = 'AUC') %>%  
  column_to_rownames(var = 'FeatureID')

write.csv(auc_table, 'output_tables/raw_auc_table.csv', row.names = TRUE)

## Pareto
HILIC_raw <- auc_table%>%
  t()%>%
  as.data.frame()

HILIC_clean <- HILIC_raw %>%
  log_transform() %>%
  pareto_scale() 
# Remove the "X" in front of column names if they exist
colnames(HILIC_clean) <- gsub("^X", "", colnames(HILIC_clean))

write.csv(HILIC_clean, file = "output_tables/HILIC_norm_pareto.csv", row.names = TRUE)

HILIC_full <- HILIC_clean %>%
  rownames_to_column('SampleName') %>%
  left_join(metadata)%>%
  select(-SampleID)

write.csv(HILIC_full, "output_tables/HILIC_LMM_norm_auc.csv", row.names = F)


# 3. Ordination -----------------------------------------------------------

norm_matrix <- read_csv('output_tables/HILIC_LMM_norm_auc.csv')%>% 
  column_to_rownames(var = 'SampleName')%>%
  select(-Type, -Slope, - Location)

# Import metadata and fix names
metadata <- read_csv('output_tables/fixed_metadata.csv')

# adjust the dataframe name based on which normalization matrix chosen 

# Calculating distance matrix
dm <- (norm_matrix) %>% 
  vegdist(., method = 'euclidean')

# Performing NMDS analysis
nmds_res <- metaMDS(dm,
                    k = 2,
                    maxit = 999,
                    trymax = 500,
                    wascores = TRUE)

nmds_res$stress
stressplot(nmds_res) 


# Extracting NMDS scores and plotting
nmds_scores <- as.data.frame(scores(nmds_res, display = 'sites')) %>%
  rownames_to_column(var = 'SampleName') %>% 
  left_join(metadata, by = 'SampleName')

# Save the nmds_scores data frame to a CSV file

table_file <- file.path(tables_dir, 'nmds_scores.csv')
write_csv(nmds_scores, table_file)

my_colors <- c("#3d1d0a", "#a8581b", "#116311")

nmds_plot <- nmds_scores %>% 
  ggplot() +
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 color = Type),
             size = 1.5) +
  scale_color_manual(values = setNames(my_colors, unique(nmds_scores$Type))) +
  labs(title = 'NMDS plot') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size = 8),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Type, fill = Type), 
               geom = "polygon", level = 0.95, alpha = 0.1, 
               linetype = 1, linewidth = 1) +
  
  # Adjust the legend position
  guides(color = guide_legend(override.aes = list(shape = NA))) + # Remove shape legend
  # Set fill colors for the ellipses
  scale_fill_manual(values = setNames(my_colors, unique(nmds_scores$Type)))

nmds_plot


ggsave(filename = "output_figures/HILIC_nmds_metabolites.svg", plot = nmds_plot, width = 4, height =3, dpi = 300)


# 4. Diversity ------------------------------------------------------------

# Load compounds table
compounds_table <- read_csv('output_tables/compounds_table_non_real_gap_filled.csv')

compounds_table_HILIC <- compounds_table %>% 
  select(- gap_status, - Type, - Slope, - Location, - SampleID, -gap_fill, -gap_fill_type) %>%
  pivot_wider(names_from = 'SampleName', values_from = 'AUC')

write_csv(compounds_table_HILIC, 'output_tables/cmpd_HILIC.csv')  

# Import metadata and fix names
metadata <- read_csv('output_tables/fixed_metadata.csv') 

raw_intensity <- read_csv('output_tables/raw_auc_table.csv')%>%
  column_to_rownames(var = '...1')

samples <- intersect(metadata$SampleName, names(compounds_table_HILIC))

#### Get intensity matrix for diversity analysis 
intensity_matrix <- compounds_table_HILIC %>%
  select(name4plot, all_of(samples)) %>%
  pivot_longer(!name4plot, names_to = 'SampleName', values_to = 'intensity') %>%
  pivot_wider(names_from = 'name4plot', values_from = 'intensity') %>%
  column_to_rownames(var = 'SampleName')

intensity_matrix[is.na(intensity_matrix)] <- 0

# Sum normalize intensities
sample_sum = rowSums(intensity_matrix)
norm_intensity_matrix <- intensity_matrix /sample_sum

# Get color palettes
my_colors <- get_palette('Dark2', length(unique(metadata$'Type')))
names(my_colors) <- unique(metadata$'Type') 


# Diversity Index

Shannon_table <- tibble(SampleName = rownames(norm_intensity_matrix),
                        Shannon = diversity(norm_intensity_matrix, index = 'shannon'))

write_csv(Shannon_table, 'output_tables/Shannon_diversity_HILIC.csv')

#Shannon_table <- read.csv("output_tables/Shannon_diversity_HILIC.csv")

Shannon_plot_table <- Shannon_table %>%
  filter(SampleName!="Mass")%>%
  pivot_longer(!SampleName, names_to = 'index', values_to = 'values') %>% 
  left_join(metadata, by = 'SampleName')

stat_table <- Shannon_plot_table %>% 
  filter(!is.na(Type)) %>%
  select(values, Type) %>% 
  wilcox_test(values ~ Type) %>% 
  add_xy_position(step.increase = 0.3) %>% 
  add_significance()

Shannon_plot <-  Shannon_plot_table %>%
  ggplot() +
  geom_boxplot(aes(x = Type,
                   y = values,
                   fill = Type)) +
  scale_fill_manual(values = c("#3d1d0a", "#a8581b", "#116311")) +
  stat_pvalue_manual(data = stat_table, inherit.aes = FALSE, hide.ns = TRUE) +
  labs(title = 'Metabolite Diversity',
       y = 'Shannon diversity') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size = 8),
        axis.text.x = element_text(face = 'bold', size = 8, colour = "black"),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  xlab(NULL) +
  guides(fill = "none")


Shannon_plot

ggsave(filename = "output_figures/HILIC_shannon_metabolites.svg", plot = Shannon_plot, width = 3, height = 3, dpi = 300)


# 5. Differentail analysis ------------------------------------------------

## Reshape the data for Linear Model 
lm_data_HILIC <- HILIC_full %>%
  pivot_longer(cols = ends_with('HILIC'), names_to = 'ft', values_to = 'area') %>%
  modify_at('Type', factor, levels = c("Bare", 'Biocrust', "Moss")) %>%
  modify_at('Slope', factor, levels = c('West', 'Center', 'East')) %>%
  modify_at('Location', factor) %>%
  group_by(ft) %>%
  nest() 

## Fit Linear Mixed effect model "Area" is response variable and "Type is fixed effect while "Slope and Location are random effect" 
# this is applied for each of the features in the nested table above
lmer_HILIC <- lm_data_HILIC %>%
  mutate(mod = purrr::map(data, function(x) lmerTest::lmer(area ~ Type + (1|Slope/Location), data = x))) 


# defining contrasts of categorical variables
Bare_Biocrust <- c(0,1,0) 
Bare_Moss <- c(0,0,1)
Biocrust_Moss <- c(0,-1,1)

## get the p-values for each comparison for each feature
lmer_HILIC_contra <- lmer_HILIC %>%
  mutate(p_Bare_Biocrust = map_dbl(mod, function(x) lmerTest::contest(x, L = Bare_Biocrust)$`Pr(>F`)) %>%
  mutate(p_Bare_Moss = map_dbl(mod, function(x) lmerTest::contest(x, L = Bare_Moss)$`Pr(>F`)) %>%
  mutate(p_Biocrust_Moss = map_dbl(mod, function(x) lmerTest::contest(x, L = Biocrust_Moss)$`Pr(>F`))

#get the bumber of significant features based on Linear Mixed effect model for different comaprisons 
lmer_HILIC_contra %>% #check the significant ones
  dplyr::select(-data, - mod) %>%
  column_to_rownames('ft')  %>%
  apply(., 2, function(x) sum(x < 0.05))

# get the adj. p-values for multiple hypothesis testing using the Bonferroni method
lmer_HILIC_contra_adj <- lmer_HILIC_contra %>% # pvalue corrections
  pivot_longer(cols = starts_with('p_'), names_to = 'test') %>%
  ungroup() %>%
  mutate(padj = p.adjust(value, method = 'bonferroni')) %>%
  mutate(new_test = paste0(gsub('p_', '', test), '_padj')) %>%
  dplyr::select(-data, -mod, -test, -value) %>%
  pivot_wider(names_from = 'new_test', values_from = 'padj')

# get the bumber of significant features based on adj. pvalue of Linear Mixed effect model for different comaprisons
lmer_HILIC_contra_adj %>%
  column_to_rownames('ft')  %>%
  apply(., 2, function(x) sum(x < 0.05))

Compound_Classification <- read.csv("Input_files/HILIC_Compound_Inchikey_classifications.csv")


HILIC_fold <- HILIC_raw %>%
  rownames_to_column('SampleName') %>%
  left_join(metadata) %>%
  pivot_longer(cols = ends_with('HILIC'), names_to = 'ft', values_to = 'area') %>%
  group_by(ft) %>%
  nest() %>%
  mutate(B_C_l2fc = map_dbl(data, function(x) calc_log2fc(x, factor = 'Type', groups = c('Bare', 'Biocrust')))) %>%
  mutate(B_M_l2fc = map_dbl(data, function(x) calc_log2fc(x, factor = 'Type', groups = c('Bare', 'Moss')))) %>%
  mutate(C_M_l2fc = map_dbl(data, function(x) calc_log2fc(x, factor = 'Type', groups = c('Biocrust', 'Moss'))))

## l2FC positive value (increasing) means Biocrust is higher than Bare in BvsC comparison

### Bare vs Biocrust
BvC <- lmer_HILIC_contra_adj %>%
  left_join(HILIC_fold) %>%
  dplyr::select(ft, B_C_l2fc, Bare_Biocrust_padj) %>%
  rename('padj' = Bare_Biocrust_padj, 'L2FC' = B_C_l2fc) %>%
  mutate(Comment = case_when(L2FC < 0 ~ 'Downregulated',                                                                    L2FC > 0 ~ 'Upregulated'))

dim(BvC)  

BvC_df <- left_join(BvC, Compound_Classification, by = c("ft" = "FeatureID"))
BvC_sig <- BvC_df %>%
  filter(padj < 0.05) %>%
  filter(abs(L2FC) > 2)
dim(BvC_sig)

write_csv(BvC_sig, file = 'output_tables/BvC_sig_HILIC.csv')

### Biocrust vs Moss
CvM <- lmer_HILIC_contra_adj %>%
  left_join(HILIC_fold) %>%
  dplyr::select(ft, C_M_l2fc, Biocrust_Moss_padj) %>%
  rename('padj' = Biocrust_Moss_padj, 'L2FC' = C_M_l2fc) %>%
  mutate(Comment = case_when(L2FC < 0 ~ 'Downregulated',                                                                    
                             L2FC > 0 ~ 'Upregulated'))

dim(CvM)
CvM_df <- left_join(CvM, Compound_Classification, by = c("ft" = "FeatureID"))

CvM_sig <- CvM_df %>%
  filter(padj < 0.05) %>%
  filter(abs(L2FC) > 2)
dim(CvM_sig)

write_csv(CvM_sig, file = 'output_tables/CvM_sig_HILIC.csv')


## Bare vs Moss
BvM <- lmer_HILIC_contra_adj %>%
  left_join(HILIC_fold) %>%
  dplyr::select(ft, B_M_l2fc, Bare_Moss_padj) %>%
  rename('padj' = Bare_Moss_padj, 'L2FC' = B_M_l2fc)%>%
  mutate(Comment = case_when(L2FC < 0 ~ 'Downregulated',                                                                    
                             L2FC > 0 ~ 'Upregulated'))

dim(BvM)  
BvM_df <- left_join(BvM, Compound_Classification, by = c("ft" = "FeatureID"))

BvM_sig <- BvM_df %>%
  filter(padj < 0.05) %>%
  filter(abs(L2FC) > 2)
dim(BvM_sig)

write_csv(BvM_sig, file = 'output_tables/BvM_sig_HILIC.csv')


## Significant compounds  --------------------------------------------------

Bare_Biocrust_LMM <- read_csv('output_tables/BvC_sig_HILIC.csv')
Bare_Biocrust_LMM$Final_Superclass[is.na(Bare_Biocrust_LMM$Final_Superclass)] <- "Unclassified"

Bare_Moss_LMM <- read_csv('output_tables/BvM_sig_HILIC.csv')
Bare_Moss_LMM$Final_Superclass[is.na(Bare_Moss_LMM$Final_Superclass)] <- "Unclassified"

Biocrust_Moss_LMM <- read_csv('output_tables/CvM_sig_HILIC.csv')
Biocrust_Moss_LMM$Final_Superclass[is.na(Biocrust_Moss_LMM$Final_Superclass)] <- "Unclassified"

## count 
#Bare_Biocrust
Bare_Biocrust.count<- Bare_Biocrust_LMM %>%
  filter(padj < 0.05) %>%
  mutate(Abundance = recode(Comment, 
                            "Downregulated" = "Decreasing", 
                            "Upregulated" = "Increasing")) %>%
  select(ft, Final_Superclass, Abundance, padj) %>%
  rename(Superclass = Final_Superclass)%>%
  count(Superclass, Abundance)%>%
  mutate(n = ifelse(Abundance == "Decreasing", -n, n))

#Bare_Moss
Bare_Moss.count<- Bare_Moss_LMM %>%
  filter(padj < 0.05) %>%
  mutate(Abundance = recode(Comment, 
                            "Downregulated" = "Decreasing", 
                            "Upregulated" = "Increasing")) %>%
  select(ft, Final_Superclass, Abundance, padj) %>%
  rename(Superclass = Final_Superclass)%>%
  count(Superclass, Abundance)%>%
  mutate(n = ifelse(Abundance == "Decreasing", -n, n))

#Biocrust_Moss
Biocrust_Moss.count<- Biocrust_Moss_LMM %>%
  filter(padj < 0.05) %>%
  mutate(Abundance = recode(Comment, 
                            "Downregulated" = "Decreasing", 
                            "Upregulated" = "Increasing")) %>%
  select(ft, Final_Superclass, Abundance, padj) %>%
  rename(Superclass = Final_Superclass)%>%
  count(Superclass, Abundance)%>%
  mutate(n = ifelse(Abundance == "Decreasing", -n, n))

#Next I add labels for each isolate and combine them into one data file 
#so that I can use the facet wrap function to put them all in one figure
Bare_Biocrust.count$Comp<- 'Biocrust compared to Bare'
Biocrust_Moss.count$Comp<- 'Moss compared to Biocrust'
Bare_Moss.count$Comp<- 'Moss compared to Bare'

class.rank.all<- rbind(Bare_Biocrust.count, Biocrust_Moss.count, Bare_Moss.count)
class.rank.all$Comp <- factor(class.rank.all$Comp, levels = c('Biocrust compared to Bare', 'Moss compared to Biocrust', 'Moss compared to Bare'))


class.plot<- ggplot(class.rank.all)+
  geom_bar(aes(y=Superclass, x=n,fill= Abundance), stat="identity", alpha=0.7)+
  facet_wrap(~Comp)+
  scale_fill_manual(values = c("lightsteelblue","tomato2"))+
  theme_classic()+xlab("Number of Features")+ylab("Superclass")+
  labs(title="Superclass ranks")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(color = "black", size = 6.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8))
class.plot

ggsave('output_figures/Class_rank_abundance_LMM.svg', class.plot, width = 10, height = 4, dpi = 300)


# Supplementary Figure S3 ----------------------------------------------------------------

# First, create top row (plots A and B)
top_row <- plot_grid(nmds_plot, Shannon_plot, labels = c("A", "B"), 
                     label_size = 9, label_fontface = "bold", align = 'hv', ncol = 2)

# Combine top row with bottom plot (C)
final_plot <- plot_grid(top_row, class.plot, labels = c("", "C"), 
                        label_size = 9, label_fontface = "bold", ncol = 1,
                        rel_heights = c(1, 1), align = 'hv')

# Display the combined plot
print(final_plot)
ggsave("output_figures/SupFigS3.png", final_plot, width = 7.25, height = 6, dpi = 300)



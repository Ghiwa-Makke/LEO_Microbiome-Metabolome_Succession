##### Microbiome-Metabolome Interactions Reshape Ecosystem Properties During Biocrust-Moss Succession
##### LEO - surface samples - LC -RP 


# 1. Data Formatting  --------------------------------------------------------

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

#set path variables
project_dir <- getwd()

# Give a name for your project so the results are all grouped together in a directory with the same name
project_name <- 'LEO_RP'

# Import data tables
cd_results_table <- read_xlsx('Input_files/LEO_RP_Compound_Annotation_Table.xlsx')

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

write_csv(cd_results_table, "output_tables/RP_cd_features_table.csv")

# Import metadata and fix names
metadata <- read_csv('Input_files/LEO_RP_metadata.csv')


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

write_csv(metadata, 'output_tables/fixed_metadata.csv')


# Split formula column into elemental counts
cd_results_table <- separate_formula(cd_results_table)
# Calculate ratios and thermodynamic indices 
cd_results_table <- calc_ratios_n_idxs(cd_results_table)
# Calculate classes
cd_results_table <- calc_classes(cd_results_table)

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

write.csv(gap_fill_status, file = "output_tables/gap_fill_status.csv")
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
write_csv(compounds_table_non_gap, 'output_tables/compounds_table_non_real_gap_filled.csv')


# 2. Normalization --------------------------------------------------------
## functions:
#Pareto Scaling:
PS_helper <- function(x) {
  (x - mean(x)) / sqrt(sd(x, na.rm = T))
}	

pareto_scale <- function(x){
  mtb_scaled <- data.frame(apply(x, 2, PS_helper))
  return(mtb_scaled)
}

#Auto Scaling:
AS_helper <- function(x) {
  (x - mean(x)) / sd(x, na.rm = T)
} 

auto_scale <- function(x){
  mtb_scaled <- apply(x, 2, AS_helper) 
  return(mtb_scaled)
}

#Log Transformation Functions:
log_helper <- function(x, min.val) {
  log2((x + sqrt(x ^ 2 + min.val ^ 2)) / 2)
}

#Log Scaling:
log_transform <- function(x){
  x_nz <- x[ ,which(apply(x, 2, sum) != 0)]
  min.val <- min(abs(x_nz[x_nz!=0]))/10
  x_log_trans <- data.frame(apply(x_nz, 2, log_helper, min.val))
  return(x_log_trans)
}


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
RP_raw <- auc_table%>%
  t()%>%
  as.data.frame()

RP_clean <- RP_raw %>%
  log_transform() %>%
  pareto_scale() 
# Remove the "X" in front of column names if they exist
colnames(RP_clean) <- gsub("^X", "", colnames(RP_clean))

write.csv(RP_clean, file = "output_tables/RP_norm_pareto.csv", row.names = TRUE)

RP_full <- RP_clean %>%
  rownames_to_column('SampleName') %>%
  left_join(metadata)%>%
  select(-SampleID)

write.csv(RP_full, "output_tables/LMM_norm_auc.csv", row.names = F)


# 3. Ordination -----------------------------------------------------------
suppressPackageStartupMessages({  
  library(tidyverse)
  library(vegan)
}) 
## Import data
norm_matrix <- read_csv('output_tables/LMM_norm_auc.csv')%>% 
  column_to_rownames(var = 'SampleName')

# Import metadata and fix names
metadata <- read_csv('output_tables/fixed_metadata.csv')

# adjust the dataframe name based on which normalization matrix chosen 
norm_matrix <- norm_matrix%>% 
  #column_to_rownames(var = 'SampleName')%>%
  select(-Type, - Location, -Slope)
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
#nmds_scores <- read.csv('output_tables/nmds_scores_pareto.csv')
write_csv(nmds_scores, 'output_tables/nmds_scores_pareto.csv')

my_colors <- c("#3d1d0a", "#a8581b", "#116311")

nmds_plot <- nmds_scores %>% 
  ggplot() +
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 color = Type),
             #shape = Slope),
             size = 1.5) +
  scale_color_manual(values = setNames(my_colors, unique(nmds_scores$Type))) +
  labs(title = 'NMDS plot') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size = 8),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8), #) +
        legend.position = "right") + # This line removes the legend
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Type, fill = Type), 
               geom = "polygon", level = 0.95, alpha = 0.1, 
               linetype = 1, linewidth = 1) +
  
  # Adjust the legend position
  #guides(color = guide_legend(override.aes = list(shape = NA))) + # Remove shape legend
  # Set fill colors for the ellipses
  scale_fill_manual(values = setNames(my_colors, unique(nmds_scores$Type)))

nmds_plot

ggsave(filename = "output_figures/RP_nmds_metabolites.svg", plot = nmds_plot, width = 4, height =3, dpi = 300)


# 4. Diversity ---------------------------------------------------------------

# Importing Libraries
suppressPackageStartupMessages({  
  library(tidyverse)
  library(vegan)
  library(rstatix)
  library(ggpubr)
  library(SYNCSA)
  source('Scripts/functions_cdis_exploration_1.R')
}) 

compounds_table <- read_csv('output_tables/compounds_table_non_real_gap_filled.csv')

compounds_table_RP <- compounds_table %>% 
  select(- gap_status, - Type, - Slope, - Location, - SampleID, -gap_fill, -gap_fill_type) %>%
  pivot_wider(names_from = 'SampleName', values_from = 'AUC')

write_csv(compounds_table_RP, 'output_tables/cmpd_RP.csv')  

# Import metadata and fix names
metadata <- read_csv('output_tables/fixed_metadata.csv')

raw_intensity <- read_csv('output_tables/raw_auc_table.csv')%>%
  column_to_rownames(var = '...1')

#### Get intensity matrix for diversity analysis 
intensity_matrix <- compounds_table_RP %>%
  select(name4plot, all_of(metadata$SampleName)) %>%
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

# Shannon Diversity Index

Shannon_table <- tibble(SampleName = rownames(norm_intensity_matrix),
                        Shannon = diversity(norm_intensity_matrix, index = 'shannon'))

write_csv(Shannon_table, 'output_tables/Shannon_diversity_RP.csv')

#Shannon_table <- read.csv("output_tables/Shannon_diversity_RP.csv")

Shannon_plot_table <- Shannon_table %>%
  filter(SampleName!="Mass")%>%
  pivot_longer(!SampleName, names_to = 'index', values_to = 'values') %>% 
  left_join(metadata, by = 'SampleName')

## get significance
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
  guides(fill = FALSE)


Shannon_plot

ggsave(filename = "output_figures/RP_shannon_metabolites.svg", plot = Shannon_plot, width = 3, height = 3, dpi = 300)


# Differential Analysis ---------------------------------------------------

library(tidyverse)

RP_full <- read.csv("output_tables/LMM_norm_auc.csv")

## linear mixed model ------------------------------------------------------

## Reshape the data for Linear Model 
lm_data_RP <- RP_full %>%
  pivot_longer(cols = ends_with('RP'), names_to = 'ft', values_to = 'area') %>%
  modify_at('Type', factor, levels = c("Bare", 'Biocrust', "Moss")) %>%
  modify_at('Slope', factor, levels = c('West', 'Center', 'East')) %>%
  modify_at('Location', factor) %>%
  group_by(ft) %>%
  nest() 

## Fit Linear Mixed effect model "Area" is response variable and "Type is fixed effect while "Slope and Location are random effect" 
# this is applied for each of the features in the nested table above
lmer_RP <- lm_data_RP %>%
  mutate(mod = purrr::map(data, function(x) lmerTest::lmer(area ~ Type + (1|Slope/Location), data = x))) 


# defining contrasts of categorical variables
Bare_Biocrust <- c(0,1,0) 
Bare_Moss <- c(0,0,1)
Biocrust_Moss <- c(0,-1,1)

## get the p-values for each comparison for each feature
lmer_RP_contra <- lmer_RP %>%
  mutate(p_Bare_Biocrust = map_dbl(mod, function(x) lmerTest::contest(x, L = Bare_Biocrust)$`Pr(>F`)) %>%
  mutate(p_Bare_Moss = map_dbl(mod, function(x) lmerTest::contest(x, L = Bare_Moss)$`Pr(>F`)) %>%
  mutate(p_Biocrust_Moss = map_dbl(mod, function(x) lmerTest::contest(x, L = Biocrust_Moss)$`Pr(>F`))

#get the number of significant features based on Linear Mixed effect model for different comaprisons 
lmer_RP_contra %>% #check the significant ones
  dplyr::select(-data, - mod) %>%
  column_to_rownames('ft')  %>%
  apply(., 2, function(x) sum(x < 0.05))

# get the adj. p-values for multiple hypothesis testing using the Bonferroni method
lmer_RP_contra_adj <- lmer_RP_contra %>% # pvalue corrections
  pivot_longer(cols = starts_with('p_'), names_to = 'test') %>%
  ungroup() %>%
  mutate(padj = p.adjust(value, method = 'bonferroni')) %>%
  mutate(new_test = paste0(gsub('p_', '', test), '_padj')) %>%
  dplyr::select(-data, -mod, -test, -value) %>%
  pivot_wider(names_from = 'new_test', values_from = 'padj')

# get the bumber of significant features based on adj. pvalue of Linear Mixed effect model for different comaprisons
lmer_RP_contra_adj %>%
  column_to_rownames('ft')  %>%
  apply(., 2, function(x) sum(x < 0.05))


## Log 2 FC ----------------------------------------------------------------
Compound_Classification <- read.csv("Input_files/RP_Compound_Inchikey_classifications.csv")

calc_log2fc <- function(x, factor, groups){
  if(length(groups) != 2){
    stop('Must Be Exactly 2 Groups')
  }
  group1 <- x[x[factor] == groups[1], ][['area']]
  group2 <- x[x[factor] == groups[2], ][['area']]
  l2fc <- log2(mean(group2)/mean(group1))
  return(l2fc)
}

RP_fold <- RP_raw %>%
  rownames_to_column('SampleName') %>%
  left_join(meta_RP) %>%
  pivot_longer(cols = ends_with('RP'), names_to = 'ft', values_to = 'area') %>%
  group_by(ft) %>%
  nest() %>%
  mutate(B_C_l2fc = map_dbl(data, function(x) calc_log2fc(x, factor = 'Type', groups = c('Bare', 'Biocrust')))) %>%
  mutate(B_M_l2fc = map_dbl(data, function(x) calc_log2fc(x, factor = 'Type', groups = c('Bare', 'Moss')))) %>%
  mutate(C_M_l2fc = map_dbl(data, function(x) calc_log2fc(x, factor = 'Type', groups = c('Biocrust', 'Moss'))))

## l2FC positive value (increasing) means Biocrust is higher than Bare in BvsC comparison
### Bare vs Biocrust
BvC <- lmer_RP_contra_adj %>%
  left_join(RP_fold) %>%
  dplyr::select(ft, B_C_l2fc, Bare_Biocrust_padj) %>%
  rename('padj' = Bare_Biocrust_padj, 'L2FC' = B_C_l2fc) %>%
  mutate(Comment = case_when(L2FC < 0 ~ 'Downregulated',                                                                    L2FC > 0 ~ 'Upregulated'))

BvC_df <- left_join(BvC, Compound_Classification, by = c("ft" = "FeatureID"))

BvC_sig <- BvC_df %>%
  filter(padj < 0.05) %>%
  filter(abs(L2FC) > 2)
dim(BvC_sig)

write_csv(BvC_sig, file = 'output_tables/Pareto_LMM_BvC_sig_RP.csv')

### Biocrust vs Moss
CvM <- lmer_RP_contra_adj %>%
  left_join(RP_fold) %>%
  dplyr::select(ft, C_M_l2fc, Biocrust_Moss_padj) %>%
  rename('padj' = Biocrust_Moss_padj, 'L2FC' = C_M_l2fc) %>%
  mutate(Comment = case_when(L2FC < 0 ~ 'Downregulated',                                                                    L2FC > 0 ~ 'Upregulated'))
CvM_df <- left_join(CvM, Compound_Classification, by = c("ft" = "FeatureID"))

CvM_sig <- CvM_df %>%
  filter(padj < 0.05) %>%
  filter(abs(L2FC) > 2)
dim(CvM_sig)

write_csv(CvM_sig, file = 'output_tables/Pareto_LMM_CvM_sig_RP.csv')

## Bare vs Moss
BvM <- lmer_RP_contra_adj %>%
  left_join(RP_fold) %>%
  dplyr::select(ft, B_M_l2fc, Bare_Moss_padj) %>%
  rename('padj' = Bare_Moss_padj, 'L2FC' = B_M_l2fc)%>%
  mutate(Comment = case_when(L2FC < 0 ~ 'Downregulated',                                                                    L2FC > 0 ~ 'Upregulated'))

dim(BvM)  
BvM_df <- left_join(BvM, Compound_Classification, by = c("ft" = "FeatureID"))

BvM_sig <- BvM_df %>%
  filter(padj < 0.05) %>%
  filter(abs(L2FC) > 2)
dim(BvM_sig)

write_csv(BvM_sig, file = 'output_tables/Pareto_LMM_BvM_sig_RP.csv')



## Significant compound classs ----------------------------------------------

Bare_Biocrust_LMM <- read_csv('output_tables/Pareto_LMM_BvC_sig_RP.csv')
Bare_Biocrust_LMM$Final_Superclass[is.na(Bare_Biocrust_LMM$Final_Superclass)] <- "Unclassified"

Bare_Moss_LMM <- read_csv('output_tables/Pareto_LMM_BvM_sig_RP.csv')
Bare_Moss_LMM$Final_Superclass[is.na(Bare_Moss_LMM$Final_Superclass)] <- "Unclassified"

Biocrust_Moss_LMM <- read_csv('output_tables/Pareto_LMM_CvM_sig_RP.csv')
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


# Figure 3 ----------------------------------------------------------------

library(cowplot)
# Add labels to each plot
plot_A <- nmds_plot + ggtitle("A") + theme(plot.title = element_text(hjust = 0, size = 9, face = "bold"))
plot_B <- Shannon_plot + ggtitle("B") + theme(plot.title = element_text(hjust = 0, size = 9, face = "bold"))
plot_C <- class.plot + ggtitle("C") + theme(plot.title = element_text(hjust = -1.5, size = 9, face = "bold"))

# Combine A and B horizontally first
top_row <- plot_grid(plot_A, plot_B, labels = NULL, ncol = 2, align = 'hv')
# Combine the top row with plot C below, spanning full width
combined_plot <- plot_grid(top_row, plot_C, labels = NULL, ncol = 1, rel_heights = c(1, 1))
combined_plot
ggsave("combined_figure_final.png", combined_plot, width = 7.25, height = 5, dpi = 300)



# First, create top row (plots A and B)
top_row <- plot_grid(nmds_plot, Shannon_plot, labels = c("A", "B"), 
                     label_size = 9, label_fontface = "bold", align = 'hv', ncol = 2)

# Combine top row with bottom plot (C)
final_plot <- plot_grid(top_row, class.plot, labels = c("", "C"), 
                        label_size = 9, label_fontface = "bold", ncol = 1,
                        rel_heights = c(1, 1), align = 'hv')

# Display the combined plot
print(final_plot)
ggsave("combined_figure_final.png", final_plot, width = 7.25, height = 6, dpi = 300)

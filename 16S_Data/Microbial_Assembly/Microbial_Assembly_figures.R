


# functions ---------------------------------------------------------------

library(reshape2)
require(vegan)
library(ggsci)
library(rstatix)
library(ggpubr)
library(ggrepel)
library(grid)
library(gridExtra)
library(phyloseq)
library(ggridges)
library(patchwork)
library(corrplot)
library(ggh4x)
library(tidyverse)

distance_func <- function(bnti_matrix, geo_matrix, geo_col){
  
  geo_selected <- geo_matrix %>% 
    select(geo_col, Sample_name) %>% 
    rename(Selected = geo_col) %>% 
    filter(!is.na(Selected))
  
  bnti_filt <- bnti_matrix[rownames(bnti_matrix) %in% geo_selected$Sample_name, 
                           colnames(bnti_matrix) %in% geo_selected$Sample_name]
  
  geo_selected <- geo_selected %>% 
    filter(Sample_name %in% rownames(bnti_filt))
  
  bnti_filt <- bnti_filt[geo_selected$Sample_name, geo_selected$Sample_name]
  
  n <- length(rownames(bnti_filt))
  dist_mat <- as.data.frame(matrix(nrow = n, ncol = n))
  rownames(dist_mat) <- rownames(bnti_filt)
  colnames(dist_mat) <- colnames(bnti_filt)
  
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      dist_mat[j, i] <- abs(geo_selected$Selected[j] - geo_selected$Selected[i])
    }
  }
  dist_mat[upper.tri(dist_mat)] = t(dist_mat)[upper.tri(dist_mat)]
  return(dist_mat)
  
}


mantel_fun <- function(bnti_matrix, geo_matrix){
  
  bnti_filt <- bnti_matrix[rownames(bnti_matrix) %in% rownames(geo_matrix), 
                           colnames(bnti_matrix) %in% colnames(geo_matrix)]
  
  mantel_test <- mantel(bnti_filt, geo_matrix, 
                        method = "pearson", permutations = 9999)
  
  return(mantel_test)
}


# set colors --------------------------------------------------------------

list_colors <- c('Bare' = '#3d1d0a', 'Biocrust' = '#a8581b','Moss' = '#116311')

# Microbial ---------------------------------------------------------------
## Load in bNTI ----
weig <- read_csv("output/leo/weighted_BNTI_ASV_leo_rel_abundance.csv") %>% 
  column_to_rownames(var = "...1")

# Load in RCBC
weig.rcbc <- read_csv("output/leo/RCBC_ASV_leo_rel_abundance.csv")%>% 
  column_to_rownames(var = "sampleid")

# metadata
metadata <- read_csv("input/leo/metadata_samples.csv")
metadata$Type <- factor(metadata$Type, levels = c("Bare","Biocrust","Moss"))

#Create a named vector to map sampleid to samplename
name_map <- setNames(metadata$sample_name, metadata$sampleid)
# Replace column and row names in the correlation matrix
colnames(weig) <- name_map[colnames(weig)]
rownames(weig) <- name_map[rownames(weig)]

colnames(weig.rcbc) <- name_map[colnames(weig.rcbc)]
rownames(weig.rcbc) <- name_map[rownames(weig.rcbc)]


## Data Reprocessing ----

# Matching the order of RCBC to bNTI results
weig.rcbc = weig.rcbc[row.names(weig), colnames(weig), drop = F]

# Setting the RCBC diagonals to NA

diag(weig.rcbc) = NA

# Reflecting null matrices

weig[upper.tri(weig)] = t(weig)[upper.tri(weig)]
weig.rcbc[upper.tri(weig.rcbc)] = t(weig.rcbc)[upper.tri(weig.rcbc)]

# Removing significant bNTI results from the RCBC results

weig.rcbc[abs(weig) > 2] = NA


## Preprocessing Data for plotting ----

# Melting data
weig = melt(as.matrix(weig))

weig.rcbc = melt(as.matrix(weig.rcbc))

#Combine data and remove null values

bnti_data <- weig %>% 
  filter(!is.na(value))

rcbc_data <- weig.rcbc %>% 
  filter(!is.na(value))

# Adding Habitat Information
bnti_data$Habitat = "Bare"
bnti_data$Habitat[grep("C", bnti_data$Var2)] = "Biocrust"
bnti_data$Habitat[grep("M", bnti_data$Var2)] = "Moss"

bnti_data$Habitat <- factor(bnti_data$Habitat, levels = c("Bare", "Biocrust", "Moss"))

rcbc_data$Habitat = "Bare"
rcbc_data$Habitat[grep("C", rcbc_data$Var2)] = "Biocrust"
rcbc_data$Habitat[grep("M", rcbc_data$Var2)] = "Moss"

rcbc_data$Habitat <- factor(rcbc_data$Habitat, levels = c("Bare", "Biocrust", "Moss"))


### BNTI plot all microbes ----- 

# Boxplots by like samples (example: Palsa vs Palsa)

bnti_data_within <- bnti_data %>% 
  filter(
    (str_detect(Var1, 'B') & str_detect(Var2, 'B')) |
      (str_detect(Var1, 'C') & str_detect(Var2, 'C')) |
      (str_detect(Var1, 'C') & str_detect(Var2, 'C'))
  ) 


stat_df <- bnti_data_within %>% 
  wilcox_test(value~Habitat) %>% 
  adjust_pvalue(method = 'bonferroni') %>% 
  add_significance() %>%
  add_xy_position() %>% 
  mutate(p.adj = paste0("P = ", p.adj))


#Violin plot

bp_group_vio_so <- bnti_data_within %>% 
  #mutate(Type = "Metabolome") %>% 
  ggplot(aes(x = Habitat, y = value, fill = Habitat))+
  geom_boxplot(width = 0.02, show.legend = FALSE) +
  geom_violin(alpha = 0.5, show.legend = FALSE) +
  geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
  labs(title = expression(bold(paste(beta, "NTI ", microbes))),
    y = expression(bold(paste(beta, 'NTI'))),
    fill = 'Habitat')+
  scale_y_continuous(limits = c(-3.5, 25)) +
  scale_fill_manual(values = list_colors)+
  stat_pvalue_manual(stat_df,
                     label = 'p.adj.signif', inherit.aes = FALSE, hide.ns = TRUE,
                     size = 4)+
  #facet_wrap(~Type)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
        
bp_group_vio_so

ggsave("violin_plot.png", plot = bp_group_vio_so, width = 5, height = 4, dpi = 300)

library(cowplot)




## Determining ecological processes fractionation
bnti_data <- bnti_data %>% 
  rename(BNTI = value)
  

rcbc_data <- rcbc_data %>% 
  rename(RCBC = value)
## Joining BNTI and RCBC datasets
data_full_microbe <- full_join(bnti_data, rcbc_data)

## removing duplicated comparisons

temp_data <- data_full_microbe

for(i in 1:nrow(temp_data)){
  temp_data$comb[i] <- paste0(sort(c(as.character(temp_data$Var1[i]), 
                                     as.character(temp_data$Var2[i]))), 
                              collapse = ';')
}

data_full_filt_micro <- temp_data %>% 
  group_by(comb) %>% 
  summarise(BNTI = mean(BNTI),
            RCBC = mean(RCBC, na.rm = TRUE)) %>% 
  separate(comb, into = c('Var1', 'Var2'), sep = ';')


eco_microbe <- data_full_filt_micro %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100)


eco_microbe_full <- data_full_microbe %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated"))


##### Plot by sampling - type ####


Bare <- metadata %>% 
  filter(Type == "Bare") %>% 
  pull(sample_name)

Biocrust <- metadata %>% 
  filter(Type == "Biocrust") %>% 
  pull(sample_name)

Moss <- metadata %>% 
  filter(Type == "Moss") %>% 
  pull(sample_name)


# Boxplots by like samples (example: Palsa vs Palsa)

bnti_data_bare <- data_full_filt_micro %>% 
  filter(Var1 %in% Bare & Var2 %in% Bare) %>% 
  mutate(Type = "Bare")

bnti_data_biocrust <- data_full_filt_micro  %>% 
  filter(Var1 %in% Biocrust & Var2 %in% Biocrust) %>% 
  mutate(Type = "Biocrust")

bnti_data_moss <- data_full_filt_micro %>% 
  filter(Var1 %in% Moss & Var2 %in% Moss) %>% 
  mutate(Type = "Moss")


# Joining data

bnti_final_microbe <- rbind(bnti_data_bare, bnti_data_biocrust, bnti_data_moss)


bnti_final_microbe$Type <- factor(bnti_final_microbe$Type, levels = c("Bare", "Biocrust", "Moss"))


new_joined <- bnti_final_microbe %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated"),
         #Type = case_when(Var1 %in% Bare & Bare %in% Bare ~ 'Bare',
          #                Var1 %in% Biocrust & Var2 %in% Biocrust ~ 'Biocrust',
           #               Var1 %in% Moss & Var2 %in% Moss ~ 'Moss'),
         #Type = factor(Type, levels = c('Bare', 'Biocrust', 'Moss')))
)


########## Ecological processes by Type

eco_type <- bnti_final_microbe %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  group_by(Type) %>% 
   filter(!is.na(ecological)) %>% ## removed NA 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100) %>% 
  mutate(Analysis = "Type") %>% 
  rename(set = Type)


## Bar plot -ecological processees - Month

ecobar_type <- ggplot(data = eco_type, aes(x = set, y = percentage, fill = ecological))+
  geom_col()+#stat = "identity", position = "dodge")+
  labs(title = expression(bold("Ecological processes - Microbiome")),
       y = expression(bold("Assembly patterns (%)")),
       fill = "Assembly process")+
  xlab(NULL)+
  ylim(0,100)+
  scale_fill_manual(values = c('Homogeneous selection' = '#ED553B', 
                               'Variable selection' = '#F9A95A', 
                               'Undominated' = '#3CAEA3', 
                               'Homogenizing dispersal' = '#5894C1' ,
                               'Dispersal limitation' = '#89C6F5'))+
  theme_bw()+
  # Facet by Type with a strip for each Type
  #facet_wrap(~Type, ncol = 1, strip.position = "top") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', face = 'bold', size = 8))

ecobar_type

ggsave("ecobar_type.png", plot = ecobar_type, width = 5, height = 4, dpi = 300)


## Mantel Heatmap

### Loading biogeochemical data


envs_data <- read_csv("input/leo/env_data.csv")

geo <- envs_data %>% 
  select(Sample_name, sampleid, water_content, MossCoverage, Biocrust_Coverage, C_N, OC, TN, pH)
## load bnti

bacteria_bnti <- read_csv("output/leo/weighted_BNTI_ASV_leo_rel_abundance.csv") %>% 
  column_to_rownames("...1")

bacteria_bnti <- bacteria_bnti[colnames(bacteria_bnti),] %>% 
  rename(!!(set_names(geo$sampleid, geo$Sample_name))) 

rownames(bacteria_bnti) <- colnames(bacteria_bnti)

# Reflecting  matrices

bacteria_r <- bacteria_bnti
bacteria_r[upper.tri(bacteria_bnti)] <-  t(bacteria_bnti)[upper.tri(bacteria_bnti)]



### Analysis by Habitat
#### Bare - ALL 


#Filtering Bare data

geo_bare <- geo %>% 
  filter(str_detect(Sample_name, "B"))

sel_columns <- c('water_content', 'MossCoverage',
                 'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')


matrix_list <- map(sel_columns, function(x){
  if(x == 'bacteria'){
    mat <- bacteria_r[rownames(bacteria_r) %in% geo_bare$Sample_name, 
                      colnames(bacteria_r) %in% geo_bare$Sample_name,]
  } else {
    mat <- distance_func(bacteria_r, geo_bare, x)
  }
  
  return(mat)
})

names(matrix_list) <- c('water_content', 'MossCoverage',
                        'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')


mantel_stats <- tibble(Geochemistry = character(), Mantel_stat = numeric(), p_value = numeric(), pval.adj = numeric(), .rows = 0)

for(i in 1:length(matrix_list)){
  
  mat <- as.data.frame(matrix_list[i])
  colnames(mat) <- str_remove(colnames(mat), '.*\\.')
  
  mantel_result <- mantel_fun(bacteria_r, mat)
  temp <- tibble(Geochemistry = names(matrix_list)[i], Mantel_stat = mantel_result$statistic, 
                 p_value = mantel_result$signif)
  
  mantel_stats <- rbind(mantel_stats, temp)
}

mantel_stats$pval.adj <- p.adjust(mantel_stats$p_value, method = 'bonferroni')


data_plot_Bare_bacteria <- mantel_stats %>% 
  mutate(Habitat = 'Bare',
         Type = 'Microbial \u03B2NTI')



#### Bacteria - Biocrust


#Filtering Biocrust
geo_crust <- geo %>% 
  filter(str_detect(Sample_name, "C"))

sel_columns <- c('water_content', 'MossCoverage',
                 'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')

matrix_list <- map(sel_columns, function(x){
  if(x == 'bacteria'){
    mat <- bacteria_r[rownames(bacteria_r) %in% geo_crust$Sample_name, 
                      colnames(bacteria_r) %in% geo_crust$Sample_name,]
  } else {
    mat <- distance_func(bacteria_r, geo_crust, x)
  }
  
  return(mat)
})

names(matrix_list) <- c('water_content', 'MossCoverage',
                        'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')


mantel_stats <- tibble(Geochemistry = character(), Mantel_stat = numeric(), p_value = numeric(), pval.adj = numeric(), .rows = 0)

for(i in 1:length(matrix_list)){
  
  mat <- as.data.frame(matrix_list[i])
  colnames(mat) <- str_remove(colnames(mat), '.*\\.')
  
  mantel_result <- mantel_fun(bacteria_r, mat)
  temp <- tibble(Geochemistry = names(matrix_list)[i], Mantel_stat = mantel_result$statistic, 
                 p_value = mantel_result$signif)
  
  mantel_stats <- rbind(mantel_stats, temp)
}

mantel_stats$pval.adj <- p.adjust(mantel_stats$p_value, method = 'bonferroni')


data_plot_Crust_bacteria <- mantel_stats %>% 
  mutate(Habitat = 'Biocrust',
         Type = 'Microbial \u03B2NTI')


#### Bacteria - Fen

#Filtering fen data 
geo_moss <- geo %>% 
  filter(str_detect(Sample_name, "M"))

sel_columns <- c('water_content', 'MossCoverage',
                 'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')

matrix_list <- map(sel_columns, function(x){
  if(x == 'bacteria'){
    mat <- bacteria_r[rownames(bacteria_r) %in% geo_moss$Sample_name, 
                      colnames(bacteria_r) %in% geo_moss$Sample_name,]
  } else {
    mat <- distance_func(bacteria_r, geo_moss, x)
  }
  
  return(mat)
})

names(matrix_list) <- c('water_content', 'MossCoverage',
                        'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')

mantel_stats <- tibble(Geochemistry = character(), Mantel_stat = numeric(), p_value = numeric(), pval.adj = numeric(), .rows = 0)

for(i in 1:length(matrix_list)){
  
  mat <- as.data.frame(matrix_list[i])
  colnames(mat) <- str_remove(colnames(mat), '.*\\.')
  
  mantel_result <- mantel_fun(bacteria_r, mat)
  temp <- tibble(Geochemistry = names(matrix_list)[i], Mantel_stat = mantel_result$statistic, 
                 p_value = mantel_result$signif)
  
  mantel_stats <- rbind(mantel_stats, temp)
}

mantel_stats$pval.adj <- p.adjust(mantel_stats$p_value, method = 'bonferroni')


data_plot_Moss_bacteria <- mantel_stats %>% 
  mutate(Habitat = 'Moss',
         Type = 'Microbial \u03B2NTI')




data_plot <- rbind(data_plot_Bare_bacteria, data_plot_Crust_bacteria, data_plot_Moss_bacteria) %>% 
  mutate(Signif = case_when(pval.adj <= 0.0001 ~ '****',
                            pval.adj <= 0.001 ~  '***',
                            pval.adj <= 0.01 ~    '**',
                            pval.adj <= 0.05 ~ '*')) %>% 
  mutate(Habitat = factor(Habitat, levels = c('Bare', 'Biocrust', 'Moss'))) #%>% 
#rbind(list('microbial_bNTI', NA, NA, NA, 'Microbial \u03B2NTI', 'Bare', NA)) %>% 
#rbind(list('microbial_bNTI', NA, NA, NA, 'Microbial \u03B2NTI', 'Crust', NA)) %>% 
#rbind(list('microbial_bNTI', NA, NA, NA, 'Microbial \u03B2NTI', 'Moss', NA))

color_strips <- strip_themed(background_x = elem_list_rect(fill = c("#3d1d0a", "#a8581b", "#116311")))


mantel_heatmap <- data_plot %>% 
  # mutate(Geochemistry = case_when(Geochemistry == 'soil_t' ~ 'Soil temperature',
  #                                 Geochemistry == 'depth' ~ 'Depth',
  #                                 Geochemistry == 'microbial_bNTI' ~ 'Microbial \u03B2NTI',
  #                                 Geochemistry == 'C:N_ratio' ~ 'C:N ratio',
  #                                 Geochemistry == 'precipitation' ~ 'Precipitation',
  #                                 TRUE ~ Geochemistry),
  mutate(Geochemistry = factor(Geochemistry, levels = c('water_content', 'MossCoverage',
                                                        'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')),
         Type = str_replace(Type, " ", "\n")) %>%  
  ggplot(aes(x = Type, y = Geochemistry, fill = Mantel_stat))+
  geom_tile(color = "white")+
  scale_fill_distiller(palette = 'RdBu', 
                       limits = c(-0.55, 0.55), na.value = 'gray80')+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_text(aes(label = Signif))+
  labs(title = expression(bold("Environmental Data Correlation")),
       fill = 'Mantel \nStatistic')+
  theme(plot.title = element_text(face="bold", hjust = 0.5),
        strip.text.x = element_text(face = "bold", color = "white", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank()) +
  facet_wrap2(~Habitat, strip = color_strips)

mantel_heatmap


# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(patchwork)

# Combine the plots
grid.arrange(bp_group_vio_so, ecobar_type, mantel_heatmap, ncol=3)


# Stack bp_group_vio_so and mantel_heatmap vertically
vertical_stack <- bp_group_vio_so / mantel_heatmap

# Combine the stacked plots with ecobar_type horizontally
combined_plot <- vertical_stack | ecobar_type

# Save the plot
ggsave("combined_plot.png", plot = combined_plot, width = 12, height = 6, dpi = 300)
ggsave("combined_plot.svg", plot = combined_plot, width = 12, height = 6, dpi = 300)

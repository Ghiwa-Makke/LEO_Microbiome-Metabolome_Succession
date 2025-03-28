#### Averaging Samples across Slopes: 


# Libraries and Functions -------------------------------------------------

library(tidyverse)
library(reshape2)
library(rstatix)
library(ggpubr)
library(vegan)
library(ggplot2)
library(ggh4x)

# Function to calculate standard deviations
compute_sd <- function(x, y, z) {
  sapply(1:nrow(x), function(i) {
    sapply(1:ncol(x), function(j) {
      sd(c(x[i, j], y[i, j], z[i, j]), na.rm = TRUE)
    })
  })
}

distance_func <- function(bnti_matrix, geo_matrix, geo_col){
  
  geo_selected <- geo_matrix %>% 
    select(geo_col, SampleID) %>% 
    rename(Selected = geo_col) %>% 
    filter(!is.na(Selected))
  
  bnti_filt <- bnti_matrix[rownames(bnti_matrix) %in% geo_selected$SampleID, 
                           colnames(bnti_matrix) %in% geo_selected$SampleID]
  
  geo_selected <- geo_selected %>% 
    filter(SampleID %in% rownames(bnti_filt))
  
  bnti_filt <- bnti_filt[geo_selected$SampleID, geo_selected$SampleID]
  
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

eco_colors <- set_names(c('#e76f51', '#f4a261', '#264653', '#2a9d8f', 'red'),
                        nm = c('Variable selection', 'Homogeneous selection',
                               'Dispersal limitation', 'Homogenizing dispersal', 
                               'Undominated'))

type_color <- c('Bare' = '#3d1d0a', 'Biocrust' = '#a8581b','Moss' = '#116311')


# BNTI ---------------------------------------------------------------------


# Load data from all three slopes
char_west = read_csv("West/output/LEO_West_TWCD_bNTI_999.csv") %>% 
  column_to_rownames(var = '...1')
char_center = read_csv("Center/output/LEO_Center_TWCD_bNTI_999.csv") %>% 
  column_to_rownames(var = '...1')
char_east = read_csv("East/output/LEO_East_TWCD_bNTI_999.csv") %>% 
  column_to_rownames(var = '...1')

# Calculate average and standard deviation matrices
average_matrix <- (char_west + char_center + char_east) / 3
sd_matrix <- compute_sd(char_west, char_center, char_east)

# Reflecting upper triangle to lower triangle to make matrices symmetric
average_matrix[upper.tri(average_matrix)] <- t(average_matrix)[upper.tri(average_matrix)]
sd_matrix[lower.tri(sd_matrix)] <- t(sd_matrix)[lower.tri(sd_matrix)]

# Set the diagonal to NA if necessary
diag(average_matrix) <- NA
diag(sd_matrix) <- NA

# Save or print your results
write.csv(average_matrix, "Average_slopes/average_bNTI_matrix.csv")
write.csv(sd_matrix, "Average_slopes/sd_bNTI_matrix.csv")


char_metabo <- average_matrix 

rm(char_west, char_center, char_east)

# RCBC --------------------------------------------------------------------


# Load data from all three slopes
rcbc_west = read_csv("West/output/LEO_West_TWCD_RCBC_9999.csv") %>% 
  column_to_rownames(var = '...1')
rcbc_center = read_csv("Center/output/LEO_Center_TWCD_RCBC_9999.csv") %>% 
  column_to_rownames(var = '...1')
rcbc_east = read_csv("East/output/LEO_East_TWCD_RCBC_9999.csv") %>% 
  column_to_rownames(var = '...1')

# Calculate average and standard deviation matrices
average_matrix <- (rcbc_west + rcbc_center + rcbc_east) / 3
sd_matrix <- compute_sd(rcbc_west, rcbc_center, rcbc_east)

# Reflecting upper triangle to lower triangle to make matrices symmetric
average_matrix[upper.tri(average_matrix)] <- t(average_matrix)[upper.tri(average_matrix)]
sd_matrix[lower.tri(sd_matrix)] <- t(sd_matrix)[lower.tri(sd_matrix)]

# Set the diagonal to NA if necessary
diag(average_matrix) <- NA
diag(sd_matrix) <- NA

# Save or print your results
write.csv(average_matrix, "Average_slopes/average_rcbc_matrix.csv")
write.csv(sd_matrix, "Average_slopes/sd_rcbc_matrix.csv")


char.rcbc_metabo <- average_matrix 

rm(rcbc_west, rcbc_center, rcbc_east)
rm(average_matrix, sd_matrix)

# Removing significant bNTI results from the RCBC results
char.rcbc_metabo[abs(char_metabo) > 2] = NA


#### Preprocessing Data for plotting ####
# ##################################### #

# Melting data
char_metabo = melt(as.matrix(char_metabo))

char.rcbc_metabo = melt(as.matrix(char.rcbc_metabo))


# Removing null values
char_metabo = char_metabo[!is.na(char_metabo$value),]

char.rcbc_metabo = char.rcbc_metabo[!is.na(char.rcbc_metabo$value),]

### pLOTTING
metadata <- read_csv("West/input/metadata_West.csv")
metadata$Type[metadata$Type == "Crust"] <- "Biocrust"

metadata$Type <- factor(metadata$Type, levels = c("Bare","Biocrust","Moss"))


###### Ecological processes per whole system ###############

# Determining ecological processes fractionation 
## This is for the whole system

## Renaming columns

bnti_data <- char_metabo %>% 
  rename(BNTI = value)

rcbc_data <- char.rcbc_metabo %>% 
  rename(RCBC = value)

## Joining BNTI and RCBC datasets

data_full_metabo <- full_join(bnti_data, rcbc_data)

## removing duplicated comparisons

temp_data <- data_full_metabo

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


## Determining ecological processes fractionation

eco_metabo <- data_full_filt_micro %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100)

################################## eco processes per sample


eco_metabo_full <- data_full_metabo %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated"))


##### Plot by sampling - type ####


Bare <- metadata %>% 
  filter(Type == "Bare") %>% 
  pull(SampleID)

Biocrust <- metadata %>% 
  filter(Type == "Biocrust") %>% 
  pull(SampleID)

Moss <- metadata %>% 
  filter(Type == "Moss") %>% 
  pull(SampleID)



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

bnti_final_metabo <- rbind(bnti_data_bare, bnti_data_biocrust, bnti_data_moss)


bnti_final_metabo$Type <- factor(bnti_final_metabo$Type, levels = c("Bare", "Biocrust", "Moss"))

# New plot

# New figure ---- within samples


new_joined <- bnti_final_metabo %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated"),
         Type = case_when(Var1 %in% Bare & Bare %in% Bare ~ 'Bare',
                          Var1 %in% Biocrust & Var2 %in% Biocrust ~ 'Biocrust',
                          Var1 %in% Moss & Var2 %in% Moss ~ 'Moss'),
         Type = factor(Type, levels = c('Bare', 'Biocrust', 'Moss')))

violin_plot_avg <- new_joined %>% 
  mutate(ecological = factor(ecological, levels = c('Variable selection', 
                                                    'Homogeneous selection',
                                                    'Dispersal limitation', 
                                                    'Homogenizing dispersal',
                                                    'Undominated'))) %>% 
  ggplot() +
  geom_violin(aes(x = Type,
                  y = BNTI),
              alpha = 0.3,
              show.legend = FALSE) +
  geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
  geom_point(aes(x = Type,
                 y = BNTI,
                 color = ecological),
             size = 2) +
  #scale_fill_manual(values = month_color)+
  scale_color_manual(values = eco_colors) +
  labs(title = 'Averaged Slopes',
       y = expression(bold(paste(beta,"-nearest taxon index"))),
       color = "Ecological processes") +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right") 

violin_plot_avg  
ggsave("Average_slopes/Violin_processes.png", plot = violin_plot_avg, width = 7, height = 4, dpi = 300)


## Violin_box plot ----

# Adding Habitat Information
bnti_data$Habitat = "Bare"
bnti_data$Habitat[grep("C$", bnti_data$Var2)] = "Biocrust"
bnti_data$Habitat[grep("M", bnti_data$Var2)] = "Moss"

bnti_data$Habitat <- factor(bnti_data$Habitat, levels = c("Bare", "Biocrust", "Moss"))

rcbc_data$Habitat = "Bare"
rcbc_data$Habitat[grep("C$", rcbc_data$Var2)] = "Biocrust"
rcbc_data$Habitat[grep("M", rcbc_data$Var2)] = "Moss"

rcbc_data$Habitat <- factor(rcbc_data$Habitat, levels = c("Bare", "Biocrust", "Moss"))

# Boxplots by like samples (example: Palsa vs Palsa)

bnti_data_within <- bnti_data %>% 
  filter(
    (str_detect(Var1, 'B') & str_detect(Var2, 'B')) |
      (str_detect(Var1, 'C') & str_detect(Var2, 'C')) |
      (str_detect(Var1, 'M') & str_detect(Var2, 'M'))
  ) 


stat_df <- bnti_data_within %>% 
  wilcox_test(BNTI~Habitat) %>% 
  adjust_pvalue(method = 'bonferroni') %>% 
  add_significance() %>%
  add_xy_position(step.increase = 0.11) %>% 
  mutate(p.adj = paste0("P = ", p.adj))

#Violin plot

violin_box_plot <- bnti_data_within %>% 
  #mutate(Type = "Metabolome") %>% 
  ggplot(aes(x = Habitat, y = BNTI, fill = Habitat))+
  geom_boxplot(width = 0.02, show.legend = FALSE) +
  geom_violin(alpha = 0.5, show.legend = FALSE) +
  geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
  labs(title = expression(bold(paste(beta, "NTI ", Metabolome))),
       y = expression(bold(paste(beta, 'NTI'))),
       fill = 'Habitat')+
  #scale_y_continuous(limits = c(-10, 10)) +
  scale_fill_manual(values = list_colors)+
  stat_pvalue_manual(stat_df,
                     label = 'p.adj.signif', inherit.aes = FALSE, hide.ns = TRUE,
                     size = 4)+
  #facet_wrap(~Type)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

violin_box_plot


ggsave("Average_slopes/Violin_boxPlot.png", plot = violin_box_plot, width = 5, height = 4, dpi = 300)


library(patchwork)

violin_combined_plot <- violin_plot_avg | violin_box_plot
violin_combined_plot

ggsave("Average_slopes/violin_combined_plot.png", plot = violin_combined_plot, width = 10, height = 4, dpi = 300)

### Density_plot ----
density_plot_data <- data_full_metabo %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  left_join(metadata, by = c('Var1' = 'SampleID')) %>% 
  mutate(Type = factor(Type, levels = c('Bare', 'Biocrust', 'Moss')))

density_plot <- density_plot_data %>% 
  mutate(ecological = factor(ecological, levels = c('Variable selection', 
                                                    'Homogeneous selection',
                                                    'Dispersal limitation', 
                                                    'Homogenizing dispersal',
                                                    'Undominated'))) %>% 
  ggplot() +
  # geom_boxplot(aes(x = month,
  #                 y = BNTI,
  #                 fill = month),
  #             alpha = 0.3,
  #             show.legend = FALSE) +
  # ggbeeswarm::geom_beeswarm(aes(x = month,
  #                               y = BNTI,
  #                               color = ecological)) +
  ggridges::geom_density_ridges(aes(x= BNTI,
                                    y = Type,
                                    fill = Type),
                                alpha = 0.5) +
  geom_vline(xintercept = c(-2,2), color = "red", lty = 2)+
  # geom_point(aes(x = month,
  #                y = BNTI,
  #                color = ecological),
  #            size = 2) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = type_color)+
  #scale_color_manual(values = assembly_color) +
  labs(title = 'Averaged Slopes',
       y = expression(bold(paste(beta,"-nearest taxon index"))),
       color = "Ecological processes") +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none') 

density_plot  

ggsave("Average_slopes/density_plot.png", plot = density_plot, width = 5, height = 4, dpi = 300)


### Overlayed density plot ----
overlayed_density_plot <- density_plot_data %>% 
  mutate(ecological = factor(ecological, levels = c('Variable selection', 
                                                    'Homogeneous selection',
                                                    'Dispersal limitation', 
                                                    'Homogenizing dispersal',
                                                    'Undominated'))) %>% 
  ggplot() +
  # Overlay densities of `Type`
  geom_density(aes(x = BNTI, 
                   fill = Type, 
                   color = Type),
               alpha = 0.3, 
               adjust = 1) +
  geom_vline(xintercept = c(-2, 2), color = "red", lty = 2) +
  scale_fill_manual(values = type_color) +
  scale_color_manual(values = type_color) +
  labs(title = 'Averaged Slopes',
       x = expression(bold(paste(beta,"-nearest taxon index"))),
       y = NULL,
       color = "Type",
       fill = "Type") +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        panel.grid = element_blank(),
        legend.position = 'right') 

overlayed_density_plot

ggsave("Average_slopes/overlayed_density_plot.png", plot = overlayed_density_plot, width = 5, height = 4, dpi = 300)


density_combined_plot <- density_plot | overlayed_density_plot
density_combined_plot

ggsave("Average_slopes/density_combined_plot.png", plot = density_combined_plot, width = 10, height = 4, dpi = 300)

## Barplot Ecological processes by Type -----

eco_type <- bnti_final_metabo %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  group_by(Type) %>% 
  # filter(!is.na(ecological)) %>% ## removed NA 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100) %>% 
  mutate(Analysis = "Type") %>% 
  rename(set = Type)

## Bar plot -ecological processees - Month
ecobar_type <- ggplot(data = eco_type, aes(x = set, y = percentage, fill = ecological))+
  geom_col()+#stat = "identity", position = "dodge")+
  labs(title = expression(bold("Ecological processes - Metabolome")),
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
  theme(plot.title = element_text(hjust = 0.5, size = 8),
        axis.title.y = element_text(hjust = 0.5, size = 8),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', face = 'bold', size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

ecobar_type

ggsave("Average_slopes/ecobar_type.png", plot = ecobar_type, width = 5, height = 4, dpi = 300)


# Mantel Heatmap ----
### Loading biogeochemical data
envs_data <- read_csv("West/input/env_data.csv")

geo <- envs_data %>% 
  filter(grepl("^S_\\d+W[a-zA-Z]+$", SampleID))%>%
  select(SampleID, water_content, MossCoverage, Biocrust_Coverage, C_N, OC, TN, pH)

## load bnti

bnti <- read_csv("Average_slopes/average_bNTI_matrix.csv") %>% 
  column_to_rownames("...1")

### Analysis by Type
#### Bare   ----
#Filtering Bare data
geo_bare <- geo %>% 
  filter(str_detect(SampleID, "B"))

sel_columns <- c('water_content', 'MossCoverage',
                 'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')

matrix_list <- map(sel_columns, function(x){
  if(x == 'metabolite'){
    mat <- bnti[rownames(bnti) %in% geo_bare$SampleID, 
                     colnames(bnti) %in% geo_bare$SampleID,]
  } else {
    mat <- distance_func(bnti, geo_bare, x)
  }
  
  return(mat)
})

names(matrix_list) <- c('water_content', 'MossCoverage',
                        'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')


mantel_stats <- tibble(Geochemistry = character(), Mantel_stat = numeric(), p_value = numeric(), pval.adj = numeric(), .rows = 0)

for(i in 1:length(matrix_list)){
  
  mat <- as.data.frame(matrix_list[i])
  colnames(mat) <- str_remove(colnames(mat), '.*\\.')
  
  mantel_result <- mantel_fun(bnti, mat)
  temp <- tibble(Geochemistry = names(matrix_list)[i], Mantel_stat = mantel_result$statistic, 
                 p_value = mantel_result$signif)
  
  mantel_stats <- rbind(mantel_stats, temp)
}

mantel_stats$pval.adj <- p.adjust(mantel_stats$p_value, method = 'bonferroni')

data_plot_Bare <- mantel_stats %>% 
  mutate(Habitat = 'Bare',
         Type = 'Metabolite \u03B2NTI')

#### Crust   ----
#Filtering Crust data
geo_Crust <- geo %>% 
  filter(str_detect(SampleID, "C$"))

sel_columns <- c('water_content', 'MossCoverage',
                 'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')

matrix_list <- map(sel_columns, function(x){
  if(x == 'metabolite'){
    mat <- bnti[rownames(bnti) %in% geo_Crust$SampleID, 
                     colnames(bnti) %in% geo_Crust$SampleID,]
  } else {
    mat <- distance_func(bnti, geo_Crust, x)
  }
  
  return(mat)
})

names(matrix_list) <- c('water_content', 'MossCoverage',
                        'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')


mantel_stats <- tibble(Geochemistry = character(), Mantel_stat = numeric(), p_value = numeric(), pval.adj = numeric(), .rows = 0)

for(i in 1:length(matrix_list)){
  
  mat <- as.data.frame(matrix_list[i])
  colnames(mat) <- str_remove(colnames(mat), '.*\\.')
  
  mantel_result <- mantel_fun(bnti, mat)
  temp <- tibble(Geochemistry = names(matrix_list)[i], Mantel_stat = mantel_result$statistic, 
                 p_value = mantel_result$signif)
  
  mantel_stats <- rbind(mantel_stats, temp)
}

mantel_stats$pval.adj <- p.adjust(mantel_stats$p_value, method = 'bonferroni')

data_plot_Crust <- mantel_stats %>% 
  mutate(Habitat = 'Crust',
         Type = 'Metabolite \u03B2NTI')

#### Moss   ----
#Filtering Moss data
geo_Moss <- geo %>% 
  filter(str_detect(SampleID, "M$"))

sel_columns <- c('water_content', 'MossCoverage',
                 'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')

matrix_list <- map(sel_columns, function(x){
  if(x == 'metabolite'){
    mat <- bnti[rownames(bnti) %in% geo_Moss$SampleID, 
                     colnames(bnti) %in% geo_Moss$SampleID,]
  } else {
    mat <- distance_func(bnti, geo_Moss, x)
  }
  
  return(mat)
})

names(matrix_list) <- c('water_content', 'MossCoverage',
                        'Biocrust_Coverage', 'C_N' , 'TN', 'pH', 'OC')


mantel_stats <- tibble(Geochemistry = character(), Mantel_stat = numeric(), p_value = numeric(), pval.adj = numeric(), .rows = 0)

for(i in 1:length(matrix_list)){
  
  mat <- as.data.frame(matrix_list[i])
  colnames(mat) <- str_remove(colnames(mat), '.*\\.')
  
  mantel_result <- mantel_fun(bnti, mat)
  temp <- tibble(Geochemistry = names(matrix_list)[i], Mantel_stat = mantel_result$statistic, 
                 p_value = mantel_result$signif)
  
  mantel_stats <- rbind(mantel_stats, temp)
}

mantel_stats$pval.adj <- p.adjust(mantel_stats$p_value, method = 'bonferroni')

data_plot_Moss <- mantel_stats %>% 
  mutate(Habitat = 'Moss',
         Type = 'Metabolite \u03B2NTI')


data <- rbind(data_plot_Bare, data_plot_Crust, data_plot_Moss) %>% 
  mutate(Signif = case_when(pval.adj <= 0.0001 ~ '****',
                            pval.adj <= 0.001 ~  '***',
                            pval.adj <= 0.01 ~    '**',
                            pval.adj <= 0.05 ~ '*')) %>% 
  mutate(Habitat = factor(Habitat, levels = c('Bare', 'Crust', 'Moss')))


color_strips <- strip_themed(background_x = elem_list_rect(fill = c("#3d1d0a", "#a8581b", "#116311")))


mantel_heatmap <- data %>% 
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
                       #limits = c(-0.58, 0.79), 
                       na.value = 'gray80')+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_text(aes(label = Signif))+
  labs(title = expression(bold("Environmental Data Correlation")),
       fill = 'Mantel \nStatistic')+
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 8),
        strip.text.x = element_text(face = "bold", color = "white", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  facet_wrap2(~Habitat, strip = color_strips)

mantel_heatmap

ggsave("Average_slopes/mantel_heatmap_avg.png", plot = mantel_heatmap, width = 6, height = 4, dpi = 300)


library(patchwork)
# Stack bp_group_vio_so and mantel_heatmap vertically
vertical_stack <- violin_box_plot / mantel_heatmap
# Combine the stacked plots with ecobar_type horizontally
combined_plot <- (vertical_stack | ecobar_type) +
  plot_annotation(
    #title = "Averaged Slope", 
    #theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
    tag_levels = 'A' ) &
  theme(plot.tag = element_text(size = 9, face = "bold"))

combined_plot
# Save the plot
ggsave("Average_slopes/combined_plot_avg.png", plot = combined_plot, width = 9, height =5, dpi = 300)
ggsave("Average_slopes/combined_plot_avg.svg", plot = combined_plot, width = 12, height = 6, dpi = 300)

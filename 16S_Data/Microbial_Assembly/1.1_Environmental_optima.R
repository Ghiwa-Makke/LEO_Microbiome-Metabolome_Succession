require(reshape2)
require(ggplot2)
require(vegan)
library(dbplyr)
library(ggsci)
library(rstatix)
library(ggpubr)
library(ggrepel)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(SYNCSA)
library(ggnewscale)
library(analogue)
library(picante)
library(phyloseq)
library(adephylo)
#library(ggtree)
library(treeio)
library(future)
library(future.apply)
library(furrr)
library(tidyverse)
library(mixOmics)
library(patchwork)



## Loading data

tree <- read_tree('input/phylogenetic_tree.tree')

## Working with rarefied matrix

abundance_matrix <- read_csv("input/ASV_matrix_rarefied_LEO.csv")

matrix_ready <- abundance_matrix %>% 
  dplyr::select(sampleid, abundance, ASV) %>% 
  pivot_wider(names_from = 'ASV', values_from = 'abundance') %>%
  #mutate(sampleid = str_replace(sampleid, "-", "_")) %>% 
  column_to_rownames(var = 'sampleid')

matrix_ready <- matrix_ready[,colSums(matrix_ready > 0) >= 4]

env <- read_csv('input/env_data.csv') %>% 
  # mutate(C_N = total_C/total_N,
  #        sample = str_replace_all(sample, '-','_')) %>%
  column_to_rownames(var = 'sampleid') 

metadata <- read_csv("input/metadata_samples.csv") 


## Checking if tree is rooted

is.rooted(tree)

## Rooting tree

tree_root <- phytools::midpoint.root(tree)

is.rooted(tree_root)


### Sanity check

all(rownames(matrix_ready) %in% metadata$sampleid)

## TRUE



### Testing phylogenetic signal of all data before inferring ecological processes


# Matching OTUs and tree

matrix_ready2 <- t(matrix_ready) %>% 
  as.data.frame()

phylo_match <- picante::match.phylo.data(tree_root, matrix_ready2)

match_otu <- as.data.frame(t(phylo_match$data))

match_tree <- phylo_match$phy

match_env <- env[rownames(match_otu),]


## Changing 

### Calculate the abundance-weighted mean environmental optima for each MAG


# Calculate the weighted average optima for each species for each environmental variable

optima_C <- as.data.frame(optima(match_otu, match_env$OC))

optima_N <- as.data.frame(optima(match_otu, match_env$TN))

optima_CN <- as.data.frame(optima(match_otu, match_env$C_N))

optima_ph <- as.data.frame(optima(match_otu, match_env$pH))

optima_moisture <- as.data.frame(optima(match_otu, match_env$water_content))

optima_Biocrust <- as.data.frame(optima(match_otu, match_env$Biocrust_Coverage))

optima_Moss <- as.data.frame(optima(match_otu, match_env$MossCoverage))

optima_metabo <- as.data.frame(optima(match_otu, match_env$Metabolite_Diversity))


# Merging optima

all_optima <- list('Carbon' = optima_C,
                   'Nitrogen' = optima_N, 
                   'C_N' = optima_CN, 
                   'pH' = optima_ph,
                   'Water_Content' = optima_moisture,
                   'Biocrust_Coverage' = optima_Biocrust,
                   'Moss_Coverage' = optima_Moss,
                   'Metabolite_Diversity' = optima_metabo)


niche_matrix <- purrr::map(all_optima, function(x){
  dist(x, method = 'euclidean')
})


## Testing Doherty et al. code

OTU_distance <- distTips(match_tree, tips = "all", "patristic") #  time consuming

print('finished distTips step')

## Normalize distance so that they range from 0-1

OTU_distance_mx <- as.matrix(OTU_distance)

xmax = max(OTU_distance_mx)
xmin = min(OTU_distance_mx)


## Note: The distance matrix is super large and 'apply' is a function that DOESN'T allow 
## paralelization

## Changing the code to use 'future'

plan(multisession, workers = 94)

# Increase the maximum allowed size for future.globals.maxSize
options(future.globals.maxSize = 2 * 1024^3)  # Set to 2 GB

# Use future_apply for parallel processing
ASV_distance_normalized <- future_apply(
  OTU_distance_mx,
  MARGIN = c(1, 2),
  FUN = function(x) {
    ((x - xmin) / (xmax - xmin))
  }
)


save.image(file = "Envs_optima_ASV_9_2_24.RData")



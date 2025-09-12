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


args <- commandArgs(trailingOnly = TRUE)

mat_idx <- as.numeric(args[1])


## Loading data


load('Envs_optima_ASV_9_2_24.RData')


### Mantel Correlog

## Changing the code to use 'future'

# plan(multisession, workers = 94)
# 
# # Increase the maximum allowed size for future.globals.maxSize
# options(future.globals.maxSize = 2 * 1024^3)  # Set to 2 GB

sel_matrix <- niche_matrix[[mat_idx]]

names_matrix <- names(niche_matrix)[mat_idx]

correlogs <- vegan::mantel.correlog(sel_matrix,
                                    ASV_distance_normalized,
                                    r.type = 'pearson',
                                    nperm = 999,
                                    n.class = 50,
                                    cutoff = FALSE)


filename <- paste0('output/Envs_optima_result_ASV_', names_matrix ,'.png')

png(filename, 
    res = 300, width = 2400, height = 1800)

plot(correlogs)
title(main = names_matrix)
dev.off()






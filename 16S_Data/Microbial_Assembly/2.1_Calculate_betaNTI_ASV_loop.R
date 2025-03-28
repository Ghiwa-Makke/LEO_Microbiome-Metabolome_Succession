
# Loading libraries


library(picante); packageVersion("picante") # handles tree manipulations
library(vegan); packageVersion("vegan") # for ecological applications
#library(ggtree)
library(treeio)
library(phyloseq)
library(adephylo)
library(tidyverse); packageVersion("tidyverse") # for dataframe processing

args = commandArgs(trailingOnly = TRUE)
start = as.numeric(args[1])
end = as.numeric(args[2])

## Working directory

setwd("/xdisk/tfaily/ghiwamakke/LEO-assembly/")

## Load data 

tree <- read_tree('input/phylogenetic_tree.tree')

## Working with rarefied matrix

abundance_matrix <- read_csv("input/ASV_matrix_rarefied_LEO.csv")

matrix_ready <- abundance_matrix %>% 
  dplyr::select(sampleid, abundance, ASV) %>% 
  pivot_wider(names_from = 'ASV', values_from = 'abundance') %>%
  #mutate(sampleid = str_replace(sampleid, "-", "_")) %>% 
  column_to_rownames(var = 'sampleid')


matrix_ready <- matrix_ready[,colSums(matrix_ready > 0) >= 4]


metadata <- read_csv("input/metadata_samples.csv") 


## Rerooting the tree between Archaea and Bacteria 


## Checking if tree is rooted

is.rooted(tree)

## Rooting tree

tree_root <- phytools::midpoint.root(tree)

is.rooted(tree_root)


### Sanity check

all(rownames(matrix_ready) %in% metadata$sampleid)

## TRUE


## Calculating relative abundance of ASV
## Based on what Hannah did in my paper for MAGs coverage


rel_abundance <- as.data.frame(t(matrix_ready)) %>% 
  rownames_to_column(var = 'ASV') %>% 
  pivot_longer(!ASV, names_to = 'sampleid', values_to = 'counts') %>% 
  filter(counts > 0) %>% 
  group_by(sampleid) %>% 
  mutate(relative_abun = counts / sum(counts)) %>% 
  select(ASV, sampleid, relative_abun) %>% 
  pivot_wider(names_from = sampleid, values_from = relative_abun, values_fill = 0) %>% 
  column_to_rownames(var = 'ASV')

# Matching OTUs and tree


## Community data matrix = MAGs as rownames and Samples as colnames

phylo_match <- match.phylo.data(tree_root, rel_abundance)

#"Dropping taxa from the data because they are not present in the phylogeny:"
## Total 1801 ASV dropped out - keeping 15061


## Step 1: First we calculate the empirical betaMNTD
# This is the observed betaMNTD which we will compare to our null distribution that we will generate below


beta.mntd.weighted <- as.matrix(comdistnt(t(phylo_match$data),
                                          cophenetic(phylo_match$phy),
                                          abundance.weighted=T))
dim(beta.mntd.weighted)

## 32 x 32


## Step 2: Calculate null expectation (randomized) betaMNTD
# This is our null distribution. Calculating this will take a while - best to run on server

# Running cophenetic outside of the for-loop
coph = cophenetic(phylo_match$phy)


#null_bMNTD <- NULL # make a null object to hold list of results

for(i in start:end){
  matrix_i <- as.matrix(comdistnt(t(phylo_match$data),
                              taxaShuffle(coph),
                              abundance.weighted=T,
                              exclude.conspecifics = F))
  
  filename <- paste0('output/beta_null_reps/beta_rep_',
                     str_pad(i, 4, pad = '0'), '.csv')
  
  matrix_i %>% 
    as.data.frame() %>% 
    rownames_to_column(var = ' ') %>% 
    select(` `, everything()) %>% 
    write_csv(filename)
  
  print(paste0(Sys.time(),
               ' Finished null bMNTD calculation for rep ' , i))
}

print(paste0('Matrices from ', start, ' to ', end, ' finished'))


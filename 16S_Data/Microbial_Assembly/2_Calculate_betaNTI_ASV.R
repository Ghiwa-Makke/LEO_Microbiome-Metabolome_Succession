
# Loading libraries


library(picante); packageVersion("picante") # handles tree manipulations
library(vegan); packageVersion("vegan") # for ecological applications
#library(ggtree)
library(treeio)
library(picante)
library(phyloseq)
library(adephylo)
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(doParallel)

## Working directory

#setwd("/xdisk/tfaily/vfreirezapata/Assembly_MAGs/Assembly_AVS_updated_2024/")

## Load data 

tree <- read_tree('input/phylogenetic_tree.tree')

## Working with rarefied matrix

abundance_matrix <- read_csv("input/ASV_matrix_rarefied_LEO.csv")

matrix_ready <- abundance_matrix %>% 
  dplyr::select(sampleid, abundance, ASV) %>% 
  pivot_wider(names_from = 'ASV', values_from = 'abundance') %>%
  #mutate(sampleid = str_replace(sampleid, "-", "_")) %>% 
  column_to_rownames(var = 'sampleid')

matrix_ready <- matrix_ready[,colSums(matrix_ready) > 1]


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

## Setting parallel execution

cl <- parallel::makeCluster(92)
registerDoParallel(cl)

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


beta.reps <- 1000

# Running cophenetic outside of the for-loop
coph = cophenetic(phylo_match$phy)


null_bMNTD <- NULL # make a null object to hold list of results
null_bMNTD <- foreach (i=1:beta.reps, .packages = c("picante")) %dopar% {
  as.matrix(comdistnt(t(phylo_match$data),
                      taxaShuffle(coph),
                      abundance.weighted=T,
                      exclude.conspecifics = F))
} # end of parallel computation

print('finished null bMNTD calculation')

# reformat data output from foreach
# Unlisting the list of null beta.reps; This creates an 3d array
# of size: # samples x # samples x # beta.reps

rand.weighted.bMNTD.comp <- array(as.numeric(unlist(null_bMNTD)), 
                                  dim=c(ncol(phylo_match$data),
                                        ncol(phylo_match$data), 
                                        beta.reps))

dim(rand.weighted.bMNTD.comp) # should be # samples x # samples x # beta.reps


## Step 3: Now calculate betaNTI by normalizing observed against the null expectation
# Reminder: 
# -2 < bNTI < 2 indicates stochastic process which can be further refined using values from RCbray in next step.
# bNTI < -2 indicates homogeneous selection
# bNTI > 2 indicates heterogeneous selection


weighted.bNTI <- matrix(c(NA),nrow=ncol(phylo_match$data),ncol=ncol(phylo_match$data));

dim(weighted.bNTI); # should be # samples x # samples

for (columns in 1:(ncol(phylo_match$data)-1)) {
  for (rows in (columns+1):ncol(phylo_match$data)) { # columns + 1 so samples are not compared to themselves
    
    # Pull out the randomized betaMNTD values
    rand.vals <- rand.weighted.bMNTD.comp[rows,columns,]; # length of this object = beta.reps
    # Normalize the observed against the randomized values
    weighted.bNTI[rows,columns] <- (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals); # normalize
    rm("rand.vals");
    
  };
};

# Rename columns and rows to get a samplexsample table
rownames(weighted.bNTI) = colnames(phylo_match$data);
colnames(weighted.bNTI) = colnames(phylo_match$data)
# weighted.bNTI;

## Save Data


write.csv(weighted.bNTI, "output/weighted_BNTI_ASV_leo_rel_abundance.csv")




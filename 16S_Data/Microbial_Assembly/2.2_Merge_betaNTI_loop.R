# Loading libraries


library(picante); packageVersion("picante") # handles tree manipulations
library(vegan); packageVersion("vegan") # for ecological applications
#library(ggtree)
library(treeio)
library(phyloseq)
library(adephylo)
library(tidyverse); packageVersion("tidyverse") # for dataframe processing



beta.mntd.weighted <- as.matrix(comdistnt(t(phylo_match$data),
                                          cophenetic(phylo_match$phy),
                                          abundance.weighted=T))
dim(beta.mntd.weighted)

# reformat data output from foreach
# Unlisting the list of null beta.reps; This creates an 3d array
# of size: # samples x # samples x # beta.reps

## Working directory

setwd("/xdisk/tfaily/vfreirezapata/Assembly_MAGs/Assembly_AVS_updated_2024/")


beta_null_files <- list.files('output/beta_null_reps',
                              full.names = TRUE)

null_bMNTD <- purrr::map(beta_null_files, function(x){
  df <- read_csv(x) %>% 
    column_to_rownames(var = '...1')
})

beta.reps <- 1000

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


write.csv(weighted.bNTI, "output/weighted_BNTI_ASV_rel_abundance.csv")




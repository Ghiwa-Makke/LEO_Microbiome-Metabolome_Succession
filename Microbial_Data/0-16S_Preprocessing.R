#### ##### LEO - Surface samples 2022 - 16S data
## Preprocessing 

## Making Phyloseq obj -----------------------------------------------------


### Load Libraries: 
library(phyloseq)
library(ape)

#Import Biom file from QIIME2. 
biom_file = "16S_data/feature-table_w_tax.biom"
biomot = import_biom(biom_file)

#Import metadata file with blanks and controls removed
map_file = ("16S_data/sample_metadata_minusblanks_controls.txt")
map_file
bmsd = import_qiime_sample_data(map_file)
sample_names(biomot)
class(bmsd)
bmsd # metadatafile

# Editing Column Names
bmsd$Type <- bmsd$covering
bmsd$Slope <- bmsd$slope
bmsd$Location <- bmsd$length
# Remove the old column
bmsd$covering <- NULL
bmsd$slope <- NULL
bmsd$length <- NULL
# Change values in the 'Type' column from 'Crust' to 'Biocrust'
bmsd$Type[bmsd$Type == "Crust"] <- "Biocrust"
bmsd$Slope[bmsd$Slope == "Central"] <- "Center"

biomot # phyloseq object
nsamples(biomot) # number of samples 62 (includes all)
ntaxa(biomot) #taxa we have 3045
sample_sums(biomot)
sample_names(biomot)
sample_names(bmsd)
nsamples(bmsd) # samples we have in metadata 54

#Merge biomot and bmsd
Merge = merge_phyloseq(biomot, bmsd) # merging taxa with metadata
Merge
sample_variables(Merge)
head(sample_data(Merge))
sample_data(Merge)
rank_names(Merge)
TaxaNames<-taxa_names(Merge)
Tax<-tax_table(Merge)
Tax

## The sample variables are also including sample ID as a variable.
#import tree
phy_tree<-read.tree("16S_data/tree.nwk")

#Merge Physeq
Merged_Phyloseq_Object<-merge_phyloseq(Merge,phy_tree,Tax,TaxaNames)
Merged_Phyloseq_Object

ntaxa(Merged_Phyloseq_Object) # 3045 taxa
nsamples(Merged_Phyloseq_Object) # 53 samples (removed outlier and other blanks)
sample_variables(Merged_Phyloseq_Object)
sample_data(Merged_Phyloseq_Object)
rank_names(Merged_Phyloseq_Object)
head(sample_data(Merged_Phyloseq_Object))
colnames(tax_table(Merged_Phyloseq_Object))<-c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")

#check col names
colnames(tax_table(Merged_Phyloseq_Object))

#checking sample variables
sample_data(Merged_Phyloseq_Object)$Type
sample_data(Merged_Phyloseq_Object)$Slope
sample_data(Merged_Phyloseq_Object)$Location
sample_variables(Merged_Phyloseq_Object)

#save this merged file to .RData filebecause it is computationally intensive to create
save(Merged_Phyloseq_Object, file= "LEO_16S.RData", compress = "bzip2")

## PS cleaning -------------------------------------------------------------
library(phyloseq)
library(tidyverse)

#Load in the data:
load( "LEO_16S.RData")

ps_raw <- Merged_Phyloseq_Object

taxa_names(ps_raw) <- paste0('ASV', 1:nrow(tax_table(ps_raw))) #Assign ASV numbers instead of codes for the rows in taxa table

#Deconstruct phyloseq object:
otu_new <- as.data.frame(otu_table(ps_raw)) # ASV in rows - Samples in columns
tax_new <- as.data.frame(tax_table(ps_raw)) # taxa classification (ASV in rows)
meta_new <- data.frame(sample_data(ps_raw)) # Metadata of the samples (Samples in rows)
tree_new <- phy_tree(ps_raw) 

tree_tips <- phy_tree(ps_raw)$tip.label

write.tree(tree_new, file = "phylogenetic_tree.tree")

# removing prefix of classification
tax_fixed <- tax_new %>%
  rownames_to_column('ASV') %>%
  modify_at(colnames(tax_new), ~gsub('k__|p__|c__|o__|f__|g__|s__', '', .x)) %>%
  mutate_all(~na_if(., "")) %>%
  column_to_rownames('ASV')

## Changing NA's into previous classification  
renames <- rownames(tax_fixed[is.na(tax_fixed[, 'Class']),]) 
taxdf <- tax_fixed[renames,]
renamed_class <- unname(sapply(rownames(taxdf), function(x) paste0('p_', taxdf[x, 'Phylum'], '_', x)))
tax_fixed[renames, 'Class'] <- renamed_class

renames <- rownames(tax_fixed[is.na(tax_fixed[, 'Order']),]) 
taxdf <- tax_fixed[renames,]
renamed_Order <- unname(sapply(rownames(taxdf), function(x) paste0('c_', taxdf[x, 'Class'], '_', x)))
tax_fixed[renames, 'Order'] <- renamed_Order

renames <- rownames(tax_fixed[is.na(tax_fixed[, 'Family']),]) 
taxdf <- tax_fixed[renames,]
renamed_family <- unname(sapply(rownames(taxdf), function(x) paste0('o_', taxdf[x, 'Order'], '_', x)))
tax_fixed[renames, 'Family'] <- renamed_family

renames <- rownames(tax_fixed[is.na(tax_fixed[, 'Genus']),]) 
taxdf <- tax_fixed[renames,]
renamed_Genus <- unname(sapply(rownames(taxdf), function(x) paste0('f_', taxdf[x, 'Family'], '_', x)))
tax_fixed[renames, 'Genus'] <- renamed_Genus

renames <- rownames(tax_fixed[is.na(tax_fixed[, 'Species']),]) 
taxdf <- tax_fixed[renames,]
renamed_Species <- unname(sapply(rownames(taxdf), function(x) paste0('g_', taxdf[x, 'Genus'], '_', x)))
tax_fixed[renames, 'Species'] <- renamed_Species


tax_final <- tax_fixed

names_taxa <- tax_table(tax_final)
taxa_names(names_taxa) <- rownames(tax_final)
# Get the column names of tax_final
original_colnames <- colnames(tax_final)

# Rename the columns of names_taxa to match the original column names
colnames(names_taxa) <- original_colnames

ps_new <- phyloseq(otu_table(otu_new, taxa_are_rows = T),
                   sample_data(meta_new),
                   names_taxa,
                   phy_tree(tree_new))

save(ps_new, file = 'PhyloSeq_edited2.RDATA')
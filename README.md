# LEO - Microbiome-Metabolome Succession

This is a repository containing the code for replicating the analysis and the figures for the paper: 

**Microbiome-Metabolome Interactions Reshape Ecosystem Properties During Biocrust-Moss Succession**

## 16S - Microbiome analysis
1. **`Analysis_1`** includes R code and input data used for Analysis of 16S and plotting Figure 2 of the paper
2. **`Picrust_Analysis`** includes R code, R functions and input data used for Picrust analysis visualization
3.  **`Microbial_Assembly`** includes R codes and input data used for community assembly analysis of the microbial data, this is based on Freire-Zapata et. al (https://www.nature.com/articles/s41564-024-01800-z)


## LC-MS/MS - Metabolome Analysis
1. **`R scripts from 1 to 2`** were used for pre-processing of LC-MS/MS data 
to calculate the ecological processes contributing to microbiome assembly. The input of these scripts include a phylogenetic tree derived from metagenome assembled geneomes and their abundances. Author Hannah Holland Moritz

2. Normalization
3. NMDS
4. Metabolite diversity
5. Differential Analysis (Linear Mixed effect Model and L2FC)
6. CCA and RDA analysis
7. PRC analysis
8. Heatmap (compound classes)
9. Heatmap (KEGG pathways) 
10. Van Krevelen Diagram

## FT-ICR-MS Metabolome analysis

### Assembly Metabolome
1. **`R scripts from 1 to 5`** were adapted from Danczak et al. 2020, https://github.com/danczakre/Meta-Metabolome_Ecology to calculate the assembly processes contributing to the peat metabolome. The input of the first script is a report of the relative abundances of masses detected through FT-ICR MS. The outputs of each script serve as input of the following ones.


## co-occurance Network: 
2. Calculating Correlations


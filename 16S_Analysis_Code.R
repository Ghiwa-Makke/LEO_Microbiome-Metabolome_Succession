#### ##### LEO - Surface samples 2022 - 16S data
#### Microbiome-Metabolome Interactions Reshape Ecosystem Properties During Biocrust-Moss Succession 


# 0 - Data Prep ---------------------------------------------------------------

## Making Phyloseq obj -----------------------------------------------------


### Load Libraries: 
library(phyloseq)
library(ape)

#Import Biom file from QIIME2. 
biom_file = "feature-table_w_tax.biom"
biomot = import_biom(biom_file)

#Import metadata file with blanks and controls removed
map_file = ("sample_metadata_minusblanks_controls.txt")
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
phy_tree<-read.tree("tree.nwk")

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
save(Merged_Phyloseq_Object, file= "WholeLEO2022Jan.RData", compress = "bzip2")





## PS cleaning -------------------------------------------------------------
library(phyloseq)
library(tidyverse)

#Load in the data:
load( "WholeLEO2022Jan.RData")

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



# 1 - Diversity ---------------------------------------------------------------
### Load Libraries
library(phyloseq)
library(magrittr)
library(tidyverse)
library(vegan)
library(microViz)

### load data
load("PhyloSeq_edited2.RDATA")
ps_raw = ps_new

# grouping at different levels 
ps_phylum <- ps_raw %>%
  tax_glom(taxrank = 'Phylum') # group the tax information at the phylum level # new variable 
ps_class <- ps_raw %>%
  tax_glom(taxrank = 'Class')
ps_order <- ps_raw %>%
  tax_glom(taxrank = 'Order')

# filtering out the singlton - doubletons ... at each level 
ps_counts_ph <- ps_phylum %>%
  filter_taxa(., function(x) sum(x > 3) > 0, TRUE)
ps_counts_c <- ps_class %>%
  filter_taxa(., function(x) sum(x > 3) > 0, TRUE)
ps_counts_o <- ps_order %>%
  filter_taxa(., function(x) sum(x > 3) > 0, TRUE)

# calculating the relative abundance (normalization)
ps_relab_ph <- ps_counts_ph %>%
  transform_sample_counts(function(x) x/sum(x)) #36
ps_relab_c <- ps_counts_c %>%
  transform_sample_counts(function(x) x/sum(x)) #188
ps_relab_o <- ps_counts_o %>%
  transform_sample_counts(function(x) x/sum(x)) #535

# filtering out the ones that have mean relative abundance below (0.0001) across all samples
# keeps the ones that are more abundant in our samples 
ps_ph <- ps_relab_ph %>%
  filter_taxa(function(x) mean(x)> 1e-4, TRUE) #31
ps_c <- ps_relab_c %>%
  filter_taxa(function(x) mean(x)> 1e-4, TRUE) #124
ps_o <- ps_relab_o %>%
  filter_taxa(function(x) mean(x)> 1e-4, TRUE) #269


## Alpha Diversity  --------------------------------------------------------
adiv <- ps_raw %>% # using raw data 
  rarefy_even_depth(rngseed = 12) %>% # subsamples (randomized) the data to an even depth 
  estimate_richness(measures = c('Observed', 'Shannon', 'Simpson')) %>% 
  rownames_to_column('sample') %>% # making rowname into a column called sample 
  left_join(., as.data.frame(sample_data(ps_raw)), by = c('sample' = 'sample.id')) # adds the alpha diversity to the table 

#27 OTUs were removed because they are no longer present in any sample after random subsampling

write.csv(adiv, file = "output_tables/diversity.microbial.csv", row.names = FALSE)
 adiv <- read.csv("output_tables/diversity.microbial.csv")
## Plotting diversity
adiv_plot_data <- adiv %>% # reshaping the data format to be easier for plotting 
  pivot_longer(cols = c('Observed', 'Shannon', 'Simpson'), names_to = 'measure')

Slope_order <- c("West", "Center", "East")
adiv_plot_data$Slope <- factor(adiv_plot_data$Slope, levels = Slope_order)

# Define custom colors for the 'Type' categories
type_colors <- c("Bare" = "#3d1d0a", "Biocrust" = "#a8581b", "Moss" = "#0b7b4c")

df <- adiv_plot_data%>%
  filter(measure == "Shannon") %>%
  select(sample_name, Type, value)

df$Type <- as.factor(df$Type)
sum(is.na(df$value))
df_clean <- df %>%
  filter(!is.na(Type) & !is.na(value))

stat.test <- df %>% 
  #group_by(Type) %>% 
  rstatix::tukey_hsd(value ~ Type) %>% 
  rstatix::add_y_position(step.increase = 1)
stat.test <- as.data.frame(stat.test)

write.csv(df_clean, "microbial_shannon_diversity.csv", row.names = F)
write.csv(stat.test, "microbial_shannon_diversity_stats.csv", row.names = F)


## Diversity Plot ----------------------------------------------------------

microb_shannon_diversity_plot <- ggplot(adiv, aes(x = Type, y = Shannon, fill = Type)) +
  geom_boxplot() +
  scale_fill_manual(values = type_colors) +
  labs(title = "Microbial Diversity",
       x = "",
       y = "Shannon Diversity") +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size = 8),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab(NULL) +
  guides(fill = FALSE) +
  ggpubr::stat_pvalue_manual(stat.test,
                             label = 'p.adj.signif', inherit.aes = FALSE,
                             hide.ns = TRUE) #+
#ggtitle("Microbial Shannon Diversity")
microb_shannon_diversity_plot
ggsave("output_figures/Microbial_shannon_diversity.svg", plot = shannon_diversity_plot, width = 3, height = 3, dpi = 300)


# 2 - Phyla realtive abundance barplot ------------------------------------------------------------
## Phylum level
# tranform the data using psmelt 
ph_bar <- ps_relab_ph %>% # using the dataset that was relative abundance and filtered for most abundant 
  psmelt() %>% # convert to long format 
  mutate(Type_abbr = case_when(
    Type == "Bare" ~ "B",
    Type == "Biocrust" ~ "C",
    Type == "Moss" ~ "M",
    TRUE ~ Type
  )) %>%
  mutate(Slope_abbr = case_when(
    Slope == "West" ~ "W",
    Slope == "East" ~ "E",
    Slope == "Center" ~ "C",
    TRUE ~ Slope
  )) %>%
  mutate(int = paste0(Type_abbr, Slope_abbr, Location)) %>%
  mutate(int2 = paste0(Slope_abbr, Type_abbr, Location))

ph_top <- ph_bar %>%
  group_by(Phylum, Type) %>%
  summarise(mean = mean(Abundance)) %>% # calculate the mean abundnace 
  ungroup() %>%
  group_by(Type) %>%
  nest() %>%
  mutate(top = map(data, function(x) x %>% arrange(desc(mean)) %>% slice_head(n = 20) %>% pull(Phylum))) # selects the top 20 by mean abundance within each type 

final_top <- c(unlist(ph_top$top)[!duplicated(unlist(ph_top$top))], 'Other')

clowns <- c(c(RColorBrewer::brewer.pal(6, 'Dark2')),
            c(RColorBrewer::brewer.pal(6, 'Accent')),
            c(RColorBrewer::brewer.pal(12, 'Set3')), 
            "#808080")

names(clowns) <- final_top

# adding a column based on the final top classes and the rest are labeled as others
ph_bar_edit <- ph_bar %>%
  mutate(Phylum_glom = ifelse(Phylum %in% final_top, Phylum, 'Other'))


## Bar Plot ----------------------------------------------------------------
h1_ph <- ggplot(ph_bar_edit, aes(x = int, y = Abundance, fill = Phylum_glom)) +
  geom_col() +
  cowplot::theme_cowplot() +
  labs(x = "")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(angle = 0, size = 8),
        axis.title.y = element_text(angle = 90, size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.background = element_rect(fill = "white"),
        legend.key.height = unit(1, "lines")) +
  scale_fill_manual(values = clowns,
                    name = 'Phylum')  +
  scale_y_continuous(expand = c(0,0),
                     labels = scales::label_percent()) +
  geom_vline(xintercept = 18.5, linewidth = 1)+
  geom_vline(xintercept = 36.5, linewidth = 1)+
  geom_vline(xintercept = 0.45, linewidth = 1)+
  geom_vline(xintercept = 53.5, linewidth = 1)+
  geom_hline(yintercept = 1, linewidth = 1)+
  geom_hline(yintercept = 0, linewidth = 1) +
  guides(fill = guide_legend(ncol = 2))

h2_ph <- ggplot(ph_bar_edit, aes(x = int, y = 1, fill = Type)) +
  geom_bar(stat = 'identity', width = 1) +
  theme_void() +
  theme(legend.position = 'none') +
  annotate('text', x = 9, y = 23, label = 'Bare', size = 4, color = "white") +
  annotate('text', x = 27, y = 23, label = 'Crust', size = 4, color = "white") +  
  annotate('text', x = 45, y = 23, label = 'Moss', size = 4, color= "white") +
  scale_fill_manual(values = c('#3d1d0a', '#a8581b', '#116311'))


leg_ph <- cowplot::get_legend(h1_ph)  

h1n_ph <- h1_ph + theme(legend.position = 'none')

h3_ph <- cowplot::plot_grid(h2_ph, NULL, h1n_ph, align = 'v', ncol = 1, axis = 'lr', rel_heights = c(1,-0.2, 12)) 

Phylum_byType <- cowplot::plot_grid(h3_ph, leg_ph, ncol = 2, rel_widths = c(10,6))
Phylum_byType
ggsave("output_figures/barplot_Phylum_byType2.svg", plot = Phylum_byType, width = 7, height = 3, dpi = 300)


# 3 - Ordination --------------------------------------------------------------
my_colors <- c('#3d1d0a', '#a8581b', '#116311')
ordC <- ordinate(ps_c, method = 'PCoA')
ord_c_plot <- plot_ordination(ps_c, ordC, color = 'Type') +
  stat_ellipse(aes(fill = Type), alpha = 0.2, geom = 'polygon') +
  geom_point(size = 1) +
  scale_fill_manual(values = my_colors)+
  scale_color_manual(values = my_colors) +
  labs(#title = "PCoA Plot from Class Composition",
    x= "PC1 (30.7%)",
    y = "PC2 (13.5%)")+
  #theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  # Centered title
ord_c_plot
ggsave("output_figures/PCoA_Type_Classes1.svg", plot = ord_c_plot, dpi = 300, width = 4, height = 2.5)


# 4 - Differential Analysis ---------------------------------------------------

# Differential Abundance Testing ---

finalASVs <- rownames(tax_table(ps_c)) # extracts the ASV identifiers of interest from the taxonomic table

# Subsetting the data: 
ps_counts_c_wn <- ps_counts_c %>%
  tax_names2rank('ASV') %>% # filters data to keep only ASVs
  subset_taxa(ASV %in% finalASVs) # Subsets the data to include only ASVs identified in the previous step.

# data preparation for Linear Mixed Effects Models (LMM) 
lmer_data <- ps_counts_c_wn %>%
  psmelt() %>%
  modify_at('Type', factor, levels = c('Bare', 'Biocrust', 'Moss')) %>% # type as factor
  modify_at('Slope', factor) %>% # slope as factor
  modify_at('Location', factor) %>% # location as factor
  group_by(ASV) %>%
  nest() # make tables for each ASV 

# hist(lmer_data[1]) 
# boxplot(lmer_data[1])

# Fitting Linear Mixed-Effects Models (LMMs):

## Type + (1|Slope/Location) -
lmer_mod <- lmer_data %>%
  mutate(mod = map(data, function(x) try(lme4::glmer.nb(Abundance ~ Type + (1|Slope/Location), data = x))))

BarevsCrust <- matrix(c(0,1,0), nrow = 1)
BarevsMoss <- matrix(c(0,0,1), nrow = 1)
CrustvsMoss <- matrix(c(0,-1,1), nrow = 1)


lmer_test <- lmer_mod %>%
  mutate(works = map_lgl(mod, function(x) class(x) != 'try-error')) %>%
  filter(works) %>%
  # Coefficients for Contrasts
  #IF NEGATIVE: Higher in Bare (1st Group)
  #IF POSITIVE: Higher in Crust (2nd Group)
  mutate(p_BarevsCrust = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = BarevsCrust))$test$pvalue)) %>%
  mutate(coef_BarevsCrust = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = BarevsCrust))$test$coefficients)) %>%
  
  #IF NEGATIVE: Higher in Bare (1st Group)
  #IF POSITIVE: Higher in Moss (2nd Group)
  mutate(p_BarevsMoss = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = BarevsMoss))$test$pvalue)) %>%
  mutate(coef_BarevsMoss = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = BarevsMoss))$test$coefficients)) %>%
  #IF NEGATIVE: Higher in Crust (1st Group)
  #IF POSITIVE: Higher in Moss (2nd Group)
  mutate(p_CrustvsMoss = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = CrustvsMoss))$test$pvalue)) %>%
  mutate(coef_CrustvsMoss = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = CrustvsMoss))$test$coefficients))

#First checking the coefficient direction of 1 feature
lmer_test[90,c(1,9:10)]

ggplot(lmer_test$data[[90]], aes(x = Type, y = Abundance, color = Type)) +
  geom_point() +
  stat_summary(geom = 'point', stat = 'mean', shape = 17, size = 5, color = 'black')

# Data Transformation and Adjustment:
lmer_test_adj <- lmer_test %>%
  dplyr::select(-works, -data, -mod) %>%
  pivot_longer(starts_with('p_'), names_to = 'test', values_to = 'pval') %>%
  ungroup() %>%
  mutate(padj = p.adjust(pval, method = 'BH')) %>%
  dplyr::select(-pval) %>%
  pivot_wider(names_from = 'test', values_from = 'padj')

# Counting Significant ASVs (adjusted p-values less than 0.05):
lmer_test_adj %>%
  column_to_rownames('ASV') %>%
  apply(., 2, function(x) sum(x<0.05))

#Taxonomic Data and Significant ASVs: (extracts taxonomic information from the ps_c phyloseq object and combines it with the ASV identifiers)
alltax <- data.frame(tax_table(ps_c)) %>%
  rownames_to_column('ASV')%>%
  select(ASV, Kingdom, Phylum, Class)

BarevsCrust_sig <- lmer_test_adj %>%
  select(ASV, coef_BarevsCrust, p_BarevsCrust)%>%
  filter(p_BarevsCrust < 0.05) %>%
  mutate(Comment = ifelse(coef_BarevsCrust>0, "Increasing", "Decreasing"))%>%
  left_join(alltax)

write.csv(BarevsCrust_sig, file = "output_tables/Microbial_BarevCrust_sig.csv", row.names = FALSE)

BarevsMoss_sig <- lmer_test_adj %>%
  select(ASV, coef_BarevsMoss, p_BarevsMoss)%>%
  filter(p_BarevsMoss < 0.05) %>%
  mutate(Comment = ifelse(coef_BarevsMoss>0, "Increasing", "Decreasing"))%>%
  left_join(alltax)

write.csv(BarevsMoss_sig, file = "output_tables/Microbial_BarevMoss_sig.csv", row.names = FALSE)

CrustvsMoss_sig <- lmer_test_adj %>%
  select(ASV, coef_CrustvsMoss, p_CrustvsMoss)%>%
  filter(p_CrustvsMoss < 0.05) %>%
  mutate(Comment = ifelse(coef_CrustvsMoss>0, "Increasing", "Decreasing"))%>%
  left_join(alltax)

write.csv(CrustvsMoss_sig, file = "output_tables/Microbial_CrustvMoss_sig.csv", row.names = FALSE)


## Significantly Changing Phyla --------------------------------------------
library(tidyverse)

Bare_Biocrust_LMM <- read.csv('output_tables/Microbial_BarevCrust_sig.csv')
Bare_Moss_LMM <- read.csv('output_tables/Microbial_BarevMoss_sig.csv')
Biocrust_Moss_LMM <- read.csv('output_tables/Microbial_CrustvMoss_sig.csv')

Bare_Biocrust_LMM$Phylum[is.na(Bare_Biocrust_LMM$Phylum)] <- "Unclassified"
Bare_Moss_LMM$Phylum[is.na(Bare_Moss_LMM$Phylum)] <- "Unclassified"
Biocrust_Moss_LMM$Phylum[is.na(Biocrust_Moss_LMM$Phylum)] <- "Unclassified"

#### Counts
#Bare_Biocrust
Bare_Biocrust.count<- Bare_Biocrust_LMM %>%
  filter(p_BarevsCrust < 0.05) %>%
  select(ASV, Phylum, X, p_BarevsCrust) %>%
  rename(Abundance = X)%>%
  count(Phylum, Abundance)%>%
  mutate(n = ifelse(Abundance == "Decreasing", -n, n))

#Bare_Moss
Bare_Moss.count<- Bare_Moss_LMM %>%
  filter(p_BarevsMoss < 0.05) %>%
  select(ASV, Phylum, Comment, p_BarevsMoss) %>%
  rename(Abundance = Comment)%>%
  count(Phylum, Abundance)%>%
  mutate(n = ifelse(Abundance == "Decreasing", -n, n))


#Biocrust_Moss
Biocrust_Moss.count<- Biocrust_Moss_LMM %>%
  filter(p_CrustvsMoss < 0.05) %>%
  select(ASV, Phylum, Comment, p_CrustvsMoss) %>%
  rename(Abundance = Comment)%>%
  count(Phylum, Abundance)%>%
  mutate(n = ifelse(Abundance == "Decreasing", -n, n))

#Next I add labels for each isolate and combine them into one data file 
#so that I can use the facet wrap function to put them all in one figure
Bare_Biocrust.count$Comp<- 'Biocrust compared to Bare'
Biocrust_Moss.count$Comp<- 'Moss compared to Biocrust'
Bare_Moss.count$Comp<- 'Moss compared to Bare'

class.rank.all<- rbind(Bare_Biocrust.count, Biocrust_Moss.count, Bare_Moss.count)
class.rank.all$Comp <- factor(class.rank.all$Comp, levels = c('Biocrust compared to Bare', 'Moss compared to Biocrust', 'Moss compared to Bare'))


## Sig_Phyla_plot ----------------------------------------------------------

class.plot<- ggplot(class.rank.all)+
  geom_bar(aes(y=fct_rev(Phylum), x=n,fill= Abundance), stat="identity", alpha=0.7)+
  facet_wrap(~Comp)+
  scale_fill_manual(values = c("lightsteelblue","tomato2"))+
  theme_classic()+xlab("Number of ASVs")+ylab("Phylum")+
  #labs(title="Superclass ranks")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(color = "black", size = 8))
class.plot

figure_file <- file.path('output_figures/phyla_abundance_LMM_l2fc.png')
ggsave(figure_file, class.plot, width=25, height=10, unit='cm')

figure_file <- file.path('output_figures/phyla_abundance_LMM_l2fc2.svg')
ggsave(figure_file, class.plot, width = 7, height = 3.5, dpi = 300)


# Figure 2 ----------------------------------------------------------------
library(patchwork)

# Assuming you have already created the three plots: ord_c_plot, shannon_diversity_plot, and h3_plot

# Arrange the plots using patchwork
combined_plot <- (ord_c_plot + microb_shannon_diversity_plot) /
  Phylum_byType /
  class.plot+
  plot_layout(heights = c(1, 2, 1.8)) +  # Set the height ratio for the rows
  plot_annotation(tag_levels = 'A') #+  # Automatically labels panels as A, B, C
#theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))  # Add margins to the plot


# Print the combined plot
print(combined_plot)
ggsave(filename = "output_figures/Microbial_Figure2.png", plot = combined_plot, width = 7, height = 8, dpi = 300)




#### ##### LEO - Surface samples 2022 - 16S data
## 16S - Differential Analysis 

# libraries:
library(phyloseq)
library(tidyverse)
library(microViz)


# 16S - Differential Analysis ---------------------------------------------
### load data
load("16S_data/PhyloSeq_edited2.RDATA")
ps_raw = ps_new

# grouping at class level
ps_class <- ps_raw %>%
  tax_glom(taxrank = 'Class')
# filtering out the singlton - doubletons 
ps_counts_c <- ps_class %>%
  filter_taxa(., function(x) sum(x > 3) > 0, TRUE)
# calculating the relative abundance (normalization)
ps_relab_c <- ps_counts_c %>%
  transform_sample_counts(function(x) x/sum(x)) #188

# filtering out the ones that have mean relative abundance below (0.0001) across all samples
# keeps the ones that are more abundant in our samples 
ps_c <- ps_relab_c %>%
  filter_taxa(function(x) mean(x)> 1e-4, TRUE) #124

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
  select(ASV, Phylum, Comment , p_BarevsCrust) %>%
  rename(Abundance = Comment)%>%
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







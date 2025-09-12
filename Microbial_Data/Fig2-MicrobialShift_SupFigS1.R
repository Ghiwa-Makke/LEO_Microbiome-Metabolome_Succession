# Load Libraries -------------
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)

# Ordinations ---------------------------------------------------------------
## 16S --------------
load("16S_data/PhyloSeq_edited2.RDATA")
ps_raw = ps_new

# grouping at different levels 
ps_phylum <- ps_raw %>%
  tax_glom(taxrank = 'Phylum') # group the tax information at the phylum level # new variable 
ps_class <- ps_raw %>%
  tax_glom(taxrank = 'Class')

# filtering out the singlton - doubletons ... at each level 
ps_counts_ph <- ps_phylum %>%
  filter_taxa(., function(x) sum(x > 3) > 0, TRUE)
ps_counts_c <- ps_class %>%
  filter_taxa(., function(x) sum(x > 3) > 0, TRUE)

# calculating the relative abundance (normalization)
ps_relab_ph <- ps_counts_ph %>%
  transform_sample_counts(function(x) x/sum(x)) #36
ps_relab_c <- ps_counts_c %>%
  transform_sample_counts(function(x) x/sum(x)) #188


# filtering out the ones that have mean relative abundance below (0.0001) across all samples
# keeps the ones that are more abundant in our samples 
ps_ph <- ps_relab_ph %>%
  filter_taxa(function(x) mean(x)> 1e-4, TRUE) #31
ps_c <- ps_relab_c %>%
  filter_taxa(function(x) mean(x)> 1e-4, TRUE) #124


# Beta Diversity --

### PCoA ----------
my_colors <- c('#3d1d0a', '#a8581b', '#116311')
ordC <- ordinate(ps_c, method = "MDS", distance = "bray")

# Extract percent variance explained (handles correction if present)
vals <- ordC$values
pct  <- if (!is.null(vals$Rel_corr_eig)) 100 * vals$Rel_corr_eig else 100 * vals$Relative_eig
pc1  <- round(pct[1], 1)
pc2  <- round(pct[2], 1)

Ord_16S <- plot_ordination(ps_c, ordC, color = 'Type') +
  geom_point(shape = 20, size = 1) +
  stat_ellipse(aes(fill = Type), alpha = 0.2, geom = 'polygon') +
  scale_fill_manual(values = my_colors)+
  scale_color_manual(values = my_colors) +
  labs(title = "PCoA Ordination - 16S rRNA",
       x = paste0("PC1 (", pc1, "%)"),
       y = paste0("PC2 (", pc2, "%)"))+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, face = "bold"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7.5),
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 7.5),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  # Centered title
Ord_16S


## MetaG -Singlem ------------
# Read SingleM full output
singlem_all <- read_tsv("Metagenomic_data/singlem.otu_table.tsv")

genes <- singlem_all %>%
  select(gene)%>%
  distinct()

# Filter for gene of interest
gene_of_interest <- "S3.40.ribosomal_protein_L11_rplK"
singlem_gene <- singlem_all %>% 
  filter(gene == gene_of_interest, !str_detect(taxonomy, "Eukaryota"))

# Generate OTU IDs
singlem_id <- singlem_gene %>%
  select(sequence) %>% 
  distinct() %>% 
  mutate(otu_id = paste0("singlem_rplK_OTU_", str_pad(n():1, 6, pad = "0")))

singlem_gene_id <- singlem_gene %>% 
  inner_join(singlem_id, by = "sequence") %>% 
  select(otu_id, everything())

# Save filtered table (optional)
#write_csv(singlem_gene_id, "singlem_rplK.csv")

# Extract sample list and assign Type and Slope
metadata <- singlem_gene_id %>%
  select(sample) %>%
  mutate(sample = substr(sample, 1, 3)) %>%
  distinct() %>%
  mutate(
    Type = case_when(
      str_ends(sample, "B") ~ "Bare",
      str_ends(sample, "C") ~ "Crust",
      str_ends(sample, "M") ~ "Moss",
      TRUE ~ "Other"
    ),
    Slope = case_when(
      str_sub(sample, 2, 2) == "C" ~ "Center",
      str_sub(sample, 2, 2) == "W" ~ "West",
      str_sub(sample, 2, 2) == "E" ~ "East",
      TRUE ~ "Unknown"
    )
  )

# Prepare OTU table
otu_table_long <- singlem_gene_id %>%
  mutate(sample = substr(sample, 1, 3)) %>%
  group_by(sample, otu_id) %>%
  summarise(value = sum(num_hits), .groups = "drop")

# Check for missing samples
missing_samples <- setdiff(metadata$sample, unique(otu_table_long$sample))
if (length(missing_samples) > 0) {
  warning("Gene not represented in these samples: ", paste(missing_samples, collapse = ", "))
}

# Pivot to wide format
otu_table_wide <- otu_table_long %>%
  pivot_wider(names_from = otu_id, values_from = value, values_fill = 0) %>%
  column_to_rownames("sample")

# Rarefy to even depth
set.seed(123)
singlem_rarefied <- rrarefy(otu_table_wide, 200) %>% as.data.frame()

# Filter OTUs in at least 4 samples
singlem_filt <- singlem_rarefied[, colSums(singlem_rarefied > 0) >= 4]


# 2. Create phyloseq object

ps_singlem <- phyloseq(
  otu_table(singlem_filt, taxa_are_rows = FALSE),
  sample_data(metadata %>% filter(sample %in% rownames(singlem_filt)) %>% column_to_rownames("sample"))
)


### PCoA: ------


# Colors
facet_colors <- c("Bare" = "#3d1d0a", "Crust" = "#a8581b", "Moss" = "#116311")

# --- PCoA: All samples 
ord_all <- ordinate(ps_singlem, method = "MDS", distance = "bray")

ordination_df_all <- plot_ordination(ps_singlem, ord_all, justDF = TRUE) %>%
  rownames_to_column("sample") 

ord_metag <- ggplot(ordination_df_all, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = Type), shape = 19, size = 1.25) +
  stat_ellipse(aes(fill = Type, colour = Type), alpha = 0.2, geom = 'polygon') +
  scale_fill_manual(values = facet_colors) +
  scale_color_manual(values = facet_colors) +
  labs(title = "PCoA Ordination - SingleM",
       x = paste0("PC1 (", round(ord_all$values$Relative_eig[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(ord_all$values$Relative_eig[2]*100, 1), "%)")) +
  coord_cartesian(xlim = range(ordination_df_all$Axis.1) * 1.2,
                  ylim = range(ordination_df_all$Axis.2) * 1.2)+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, face = "bold"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7.5),
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 7.5),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  # Centered title

ord_metag

## Fungal ITS ---------

## Load data:
vic = read.table('ITS_data/Clean_DADA2_ITS.csv', sep = ',', header = 1, row.names = 1)
vim = read.csv('ITS_data/metadata_its.csv', sep = ',', header = 1,fileEncoding = 'UTF-8-BOM')
taxa <- read.csv("ITS_data/Clean_DADA2_ITS_taxa.csv")

### standerdize : 
vic.s <- decostand(t(vic), method = "hellinger")

### PCoA: ------

# Convert to data frame 
vic_df <- as.data.frame(t(vic.s))

# Add OTU_IDs: OTU1, OTU2, ...
vic_df$OTU_ID <- paste0("OTU", seq_len(nrow(vic_df)))

## ID matching:
fungal_ID <- vic_df%>%
  select(OTU_ID)%>%
  rownames_to_column(var = "OTU")

taxa_names1 <- taxa %>%
  left_join(fungal_ID, by="OTU")%>%
  select(-OTU)

# use OTU_ID as rownames

rownames(vic_df) <- vic_df$OTU_ID

vic_df <- vic_df[, -ncol(vic_df)]  # remove OTU_ID column (already used as rowname)


vic_t <- as.data.frame(t(vic_df))
common_samples <- intersect(rownames(vic_t), vim$Sample)# Ensure sample names match

# Filter to shared samples
vic_t <- vic_t[common_samples, ]
vim_filt <- vim %>%
  filter(Sample %in% common_samples) %>%
  column_to_rownames("Sample")

# Create phyloseq object
ps_vic <- phyloseq(
  otu_table(vic_t, taxa_are_rows = FALSE),
  sample_data(vim_filt)
)

# Run PCoA and plot

# Define colors
facet_colors <- c("Bare" = "#3d1d0a", "Crust" = "#a8581b", "Moss" = "#116311")

# Run PCoA
ord_vic <- ordinate(ps_vic, method = "MDS", distance = "bray")

# Extract coordinates
ordination_df <- plot_ordination(ps_vic, ord_vic, justDF = TRUE) %>%
  rownames_to_column("Sample") 

# Variance explained
eig_vals <- ord_vic$values$Eigenvalues
var_exp <- eig_vals / sum(eig_vals)
percent1 <- round(var_exp[1] * 100, 1)
percent2 <- round(var_exp[2] * 100, 1)

# Plot
ord_its <- ggplot(ordination_df, aes(x = Axis.1, y = Axis.2, fill = Type, colour = Type)) +
  geom_point( shape = 19, size = 1.25) +
  stat_ellipse(aes(fill = Type), alpha = 0.2, geom = 'polygon') +
  scale_fill_manual(values = facet_colors) +
  scale_color_manual(values = facet_colors) +
  labs(title = "PCoA Ordination - Fungal ITS",
       x = paste0("PC1 (", percent1, "%)"),
       y = paste0("PC2 (", percent2, "%)")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, face = "bold"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7.5),
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 7.5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  # Centered title
ord_its


# Ordination - combined:-----------
library(patchwork)

# Arrange the plots using patchwork
combined_plot <- (Ord_16S + ord_metag + ord_its)+
  plot_annotation(tag_levels = 'A') 

# Print the combined plot
print(combined_plot)

ggsave(filename = "output_figures/Fig2-ABC.png", plot = combined_plot, width = 15, height = 4, dpi = 300)



# Taxonomy BarPlots -------
#define color scheme

phyl_colors <- c(
  "Acidobacteria"   = "#1b9e77",  "Acidobacteriota" = "#1b9e77",
  "Actinobacteria" = "#D95F02",   "Actinomycetota" = "#d95f02",
  "Armatimonadetes" = "#7570B3",  "Armatimonadota" = "#7570B3",
  "Bacteroidetes"   = "#E7298A",   "Bacteroidota"    = "#E7298A",
  "Chloroflexi" = "#E6AB02", "Chloroflexota" = "#E6AB02",
  "Cyanobacteria" = "#66A61E", "Cyanobacteriota" = "#66A61E", 
  "Patescibacteriota" = "#FFED6F",
  "Planctomycetes" = "#BEAED4", "Planctomycetota" = "#BEAED4", 
  "Proteobacteria"  = "#386CB0", "Pseudomonadota" = "#386CB0", 
  "Thermoproteota" = "#FDC086",  "[Thermi]" = "#FDC086", 
  "Verrucomicrobia" = "#CCEBC5", 
  "Euryarchaeota" = "#FCCDE5", 
  "Gemmatimonadetes" = "#FFFF99", 
  "Armatimonadetes" = "bisque3", 
  "Other"           = "#808080"
)

### 16S  --------------------------
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
  mutate(top = map(data, function(x) x %>% arrange(desc(mean)) %>% slice_head(n = 10) %>% pull(Phylum))) # selects the top 20 by mean abundance within each type 

final_top <- c(unlist(ph_top$top)[!duplicated(unlist(ph_top$top))], 'Other')

# adding a column based on the final top classes and the rest are labeled as others
ph_bar_edit <- ph_bar %>%
  mutate(Phylum_glom = ifelse(Phylum %in% final_top, Phylum, 'Other'))

# Start from per-sample relative abundance (ph_bar_edit), then average within Type
ph_type_ra <- ph_bar_edit %>%
  group_by(Type, Phylum_glom) %>%
  summarise(r_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = factor(Type, levels = c("Bare","Biocrust","Moss")))

# Plot: 3 stacked bars, one per Type
Phylum_byType_3bars <- ggplot(ph_type_ra,
                              aes(x = Type, y = r_abundance, fill = Phylum_glom)) +
  geom_col(color = "black", linewidth = 0.2) +
  scale_fill_manual(values = phyl_colors, name = "Phylum") +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(title = "Taxonomy profile - 16S rRNA",
    x = NULL, y = "Relative abundance") +
  cowplot::theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.justification = c(0.5, 1.5),  
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.box = "vertical",
        legend.box.just = "left",
        legend.key.spacing.y = unit(.05, 'lines'),
        legend.key.spacing.x = unit(.5, 'lines'),
        legend.key.size = unit(0.3, 'lines'),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8))+
  guides(fill = guide_legend(ncol = 2))

Phylum_byType_3bars
ggsave("output_figures/barplot_Phylum_byType_3bars.svg",
       plot = Phylum_byType_3bars, width = 4.5, height = 6, dpi = 300)


## Singlem- MetaG -----
singlem_taxonomy <- singlem_all %>% 
  mutate(sample = substr(sample, 1, 3)) %>% 
  group_by(sample, taxonomy) %>%
  summarise(coverage = sum(coverage), .groups = 'drop') %>%
  separate(taxonomy, 
           into = c('Root', 'Domain', 'Phylum', 'Class', 'Order', 
                    'Family', 'Genus', 'Species'),
           sep = '; ') %>% 
  filter(!is.na(Phylum),
         Domain != 'd__Eukaryota') %>% 
  inner_join(metadata, by = 'sample') %>% 
  select(sample, coverage, Domain, Phylum, Class, Order, Family, Genus, Species) %>% 
  mutate(across(Domain:Species, ~str_remove(.x, '.__')))




singlem_taxonomy <- singlem_taxonomy %>% 
  mutate(Phylum_lump = fct_lump_n(Phylum, 10, w = coverage))

singlem_ra <- singlem_taxonomy %>% 
  group_by(sample) %>% 
  mutate(r_abundance = coverage / sum (coverage)) %>% 
  group_by(sample, Phylum_lump) %>% 
  summarise(r_abundance = sum(r_abundance)) %>% 
  inner_join(metadata, by ="sample")

singlem_ra <- singlem_ra %>%
  mutate(sample = factor(sample, levels = unique(sample[order(Type)])))

facet_colors <- c("Bare" = "#3d1d0a", 
                  "Crust" = "#a8581b", 
                  "Moss" = "#116311")


### type_plot
singlem_tax_plot2 <- singlem_ra %>%  
  group_by(Phylum_lump, Type) %>%
  summarise(r_abundance = mean(r_abundance)) %>%
  ggplot(aes(x = Type, y = r_abundance, fill = Phylum_lump)) +
  geom_col(color = 'black',
           linewidth = 0.2) +
  scale_fill_manual(values = phyl_colors) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(title = 'Taxonomy profile - SingleM',
    y = 'Relative Abundance',
    fill = 'Phylum') +
  guides(fill = guide_legend(ncol = 2)) +
  cowplot::theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.justification = c(0.2, 1),  
        legend.key.spacing.y = unit(.05, 'lines'),
        legend.key.spacing.x = unit(.5, 'lines'),
        legend.key.size = unit(0.3, 'lines'),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.box = "vertical",
        legend.box.just = "left",
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8))

singlem_tax_plot2


## Fungal ITS -----
counts_long <- vic %>%
  rownames_to_column("OTU") %>%
  left_join(fungal_ID, by="OTU")%>%
  select(-OTU)%>%
  pivot_longer(-OTU_ID, names_to = "Sample", values_to = "count")

## get taxa classifications
its_tax_long <- counts_long %>%
  left_join(taxa_names1, by = "OTU_ID") %>%
  left_join(vim %>% select(Sample, Type), by = "Sample")

## calculate relative abundance
its_ra <- its_tax_long %>%
  mutate(Phyla = ifelse(is.na(Phyla) | Phyla == "", "Unclassified", Phyla)) %>%
  group_by(Sample) %>%
  mutate(sample_total = sum(count, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(sample_total > 0) %>%
  group_by(Sample, Type, Phyla) %>%
  summarise(r_abundance = sum(count, na.rm = TRUE) / first(sample_total), .groups = "drop")


## Collapsed by Type
its_tax_plot_type <- its_ra %>%
  group_by(Type, Phyla) %>%
  summarise(r_abundance = mean(r_abundance), .groups = "drop") %>%
  ggplot(aes(x = Type, y = r_abundance, fill = Phyla)) +
  geom_col(color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(title = "Taxonomy profile - ITS", 
    y = "Relative Abundance", 
    fill = "Phylum") +
  cowplot::theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7.5),
        legend.justification = c(0.5, 1),
        legend.key.spacing.y = unit(.05, 'lines'),
        legend.key.spacing.x = unit(.5, 'lines'),
        legend.key.size = unit(0.3, 'lines'),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8))+
  guides(fill = guide_legend(ncol = 1))

its_tax_plot_type


# Taxonomy Combined ------

# Arrange the plots using patchwork
combined_plot <- (Phylum_byType_3bars + singlem_tax_plot2 + its_tax_plot_type)+
  plot_annotation(tag_levels = 'A') 

# Print the combined plot
print(combined_plot)

ggsave(filename = "output_figures/Fig2-DEF.png", plot = combined_plot, width = 15, height = 8, dpi = 300)



# figure2 ----
figure2 <- (ord_c_plot2 + pcoa_all + pcoa_plot) /
  (Phylum_byType_3bars + singlem_tax_plot2 + its_tax_plot_type) +
  plot_layout(height = c(1.5,2)) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(size = 9, face = "bold")
  )
figure2

ggsave(filename = "output_figures/figure2-Full_updated2.png", plot = figure2, width = 7.25, height = 6, dpi = 300, units = "in")




# supplamentary figure 1 ----
## 16S ----

adiv <- ps_raw %>% # using raw data 
  rarefy_even_depth(rngseed = 12) %>% # subsamples (randomized) the data to an even depth 
  estimate_richness(measures = c('Observed', 'Shannon', 'Simpson')) %>% 
  rownames_to_column('sample') %>% 
  left_join(., as.data.frame(sample_data(ps_raw)), by = c('sample' = 'sample.id'))

adiv_plot_data <- adiv %>% # reshaping the data format to be easier for plotting 
  pivot_longer(cols = c('Observed', 'Shannon', 'Simpson'), names_to = 'measure')


# Define custom colors for the 'Type' categories
type_colors <- c("Bare" = "#3d1d0a", "Biocrust" = "#a8581b", "Moss" = "#0b7b4c")

df <- adiv_plot_data%>%
  filter(measure == "Shannon") %>%
  select(sample_name, Type, value)

df$Type <- as.factor(df$Type)
sum(is.na(df$value))
df_clean <- df %>%
  filter(!is.na(Type) & !is.na(value))

stat.test_16 <- df %>% 
  rstatix::tukey_hsd(value ~ Type) %>% 
  rstatix::add_y_position(step.increase = 1)
stat.test_16 <- as.data.frame(stat.test_16)

#write.csv(df_clean, "microbial_shannon_diversity.csv", row.names = F)
#write.csv(stat.test, "microbial_shannon_diversity_stats.csv", row.names = F)


microb_shannon_diversity_plot <- ggplot(adiv, aes(x = Type, y = Shannon, fill = Type)) +
  geom_boxplot() +
  scale_fill_manual(values = type_colors) +
  labs(title = "Shannon Diversity - 16S",
       x = "",
       y = "Shannon Diversity") +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size = 8),
        axis.text.x = element_text(size = 7, colour = "black"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab(NULL) +
  guides(fill = FALSE) +
  ggpubr::stat_pvalue_manual(stat.test_16,
                             label = 'p.adj.signif', inherit.aes = FALSE,
                             hide.ns = TRUE) #+
#ggtitle("Microbial Shannon Diversity")
microb_shannon_diversity_plot


## MetaG - Singlem ----------
# Split SingleM data by gene
singlem_split <- singlem_all %>%
  split(.$gene)

# Calculate alpha diversity per gene, across all and by slope
# Loop to calculate alpha diversity per gene
alpha_by_gene <- imap_dfr(singlem_split, function(df, gene) {
  tryCatch({
    # Format abundance table
    df <- df %>%
      mutate(sample = substr(sample, 1, 3))
    
    otu_table <- df %>%
      group_by(sample, sequence) %>%
      summarise(value = sum(num_hits), .groups = "drop") %>%
      pivot_wider(names_from = sequence, values_from = value, values_fill = 0) %>%
      column_to_rownames("sample")
    
    # Skip genes with too few samples
    if (nrow(otu_table) < 2) return(NULL)
    
    rare_depth <- min(rowSums(otu_table))
    if (rare_depth < 10) return(NULL)
    
    # Rarefy
    rarefied <- rrarefy(otu_table, rare_depth)
    
    # Fix sample names that start with digits (X1CB)
    sample_ids <- str_remove(rownames(rarefied), "^X")
    rownames(rarefied) <- sample_ids
    
    # Match metadata
    meta <- metadata %>%
      mutate(sample = str_remove(sample, "^X")) %>%
      filter(sample %in% sample_ids) %>%
      column_to_rownames("sample")
    
    if (nrow(meta) < 2) return(NULL)
    
    # Build phyloseq
    physeq <- phyloseq(
      sample_data(meta),
      otu_table(rarefied, taxa_are_rows = FALSE)
    )
    
    # Compute richness and merge
    estimate_richness(physeq, measures = c("Observed", "Shannon", "Chao1")) %>%
      rownames_to_column("sample") %>%
      pivot_longer(cols = c("Observed", "Shannon", "Chao1"),
                   names_to = "Metric", values_to = "Value") %>%
      mutate(sample = str_remove(sample, "^X")) %>%
      left_join(metadata, by = "sample") %>%
      mutate(gene = gene)
  },
  error = function(e) {
    message("Skipping gene: ", gene, " (", conditionMessage(e), ")")
    return(NULL)
  })
})

alpha_avg_per_sample <- alpha_by_gene %>%
  group_by(sample, Metric, Type) %>%
  summarise(mean_diversity = mean(Value, na.rm = TRUE), .groups = "drop")

shannon_avg_type <- alpha_avg_per_sample %>%
  filter(Metric == "Shannon")

shannon_all <- ggplot(shannon_avg_type, aes(x = Type, y = mean_diversity, fill = Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Bare" = "#3d1d0a", "Crust" = "#a8581b", "Moss" = "#116311")) +
  labs(title = "Shannon Diversity - SingleM",
       x = "", 
       y = "Mean Shannon Diversity", fill = "Type"
  ) +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size = 8),
        axis.text.x = element_text(size = 7, colour = "black"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab(NULL) +
  guides(fill = FALSE) 

shannon_all

## Fungal -----
shannon_df <- estimate_richness(ps_vic, measures = "Shannon") %>%
  rownames_to_column("Sample") %>%
  left_join(vim, by = "Sample")  
stat.test_fungal <- shannon_df %>% 
  rstatix::tukey_hsd(Shannon ~ Type) %>% 
  rstatix::add_y_position(step.increase = 1)
stat.test_fungal <- as.data.frame(stat.test_fungal)

# --- Plot Shannon diversity ---

fungal_shannon_diversity_plot_cleaned <- ggplot(shannon_df, aes(x = Type, y = Shannon, fill = Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Bare" = "#3d1d0a", "Crust" = "#a8581b", "Moss" = "#116311")) +
  labs(title = "Shannon Diversity - ITS",
       x = "",
       y = "Shannon Diversity") +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size = 8),
        axis.text.x = element_text(size = 7, colour = "black"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab(NULL) +
  guides(fill = FALSE) +
  ggpubr::stat_pvalue_manual(stat.test_fungal,
                             label = 'p.adj.signif', inherit.aes = FALSE,
                             hide.ns = TRUE) 

fungal_shannon_diversity_plot_cleaned


# Combined-shannon ----------
sup_fig1 <- (microb_shannon_diversity_plot + shannon_all + fungal_shannon_diversity_plot_cleaned)+
  plot_annotation(tag_levels = 'A')&
  theme(
    plot.tag = element_text(size = 9, face = "bold")
  )

sup_fig1
ggsave(filename = "output_figures/sup_fig1.png", plot = sup_fig1, width = 7.25, height = 3 , dpi = 300, units = "in")

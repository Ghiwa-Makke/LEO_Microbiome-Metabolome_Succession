# Packages
library(tidyverse)
library(vegan)
library(ggrepel)


# PCA Functional Plot -----
cols <- c(Bare = "#3d1d0a", Crust = "#a8581b", Moss = "#116311")

# 1) metadata 
meta <- read.csv("Metagenomic_data/Metadata.csv")                 # Sample, Type

# 2) Trait completeness scaled to [0,1] # derived from DRAM output
cnps_df <- read_csv('Metagenomic_data/traits_from_annotations_presence_by_bin.csv') %>% 
  filter(!str_detect(m_id, 'Trait'),
         def != 'K08094 K13831')

## C_N_S_P definition based on Ayala-Ortiz et al 2025 https://doi.org/10.1038/s43247-025-02637-y 
cnps_ids <- read_csv('Metagenomic_data/definitions_C_N_S_P_11_14_24.csv') %>% 
  arrange(Module_id) %>% 
  group_by(Category_1) %>% 
  mutate(order = n():1,
         rx_per_path = n())


cpns_presence <- cnps_df %>% 
  #inner_join(mag_phyl, by = 'bin') %>% 
  inner_join(cnps_ids, by = c('m_id' = 'Module_id')) %>% 
  distinct() %>% 
  filter((Category_2 %in% c('Methane metabolism', 'Carbon fixation', 'Sulfur metabolism', 'Nitrogen metabolism')))

cpns_binplot <- cpns_presence %>% 
  filter(present) %>% 
  group_by(bin, Category_2, Category_1, rx_per_path) %>% 
  summarise(rx_present = n()) %>% 
  ungroup() %>% 
  mutate(perc_present = rx_present / rx_per_path)

traits01 <- cpns_binplot %>%
  select(bin, Category_2, Category_1, perc_present)


# Trimmed_means -
mags_coverm <- read_tsv("Metagenomic_data/trimmed_means_coverm_aln.tsv") %>%
  pivot_longer(cols = -Genome, names_to = "Sample", values_to = "abund") %>%
  mutate(Sample = str_remove(Sample, "_raw.fastq.gz Trimmed Mean*")) %>%   # Keep first 3 chars and prefix "S"
  mutate(Sample = paste0("S", Sample))%>%
  group_by(Sample, Genome) %>%
  summarise(abund= mean(abund), .groups = "drop") %>%
  pivot_wider(names_from = Genome, values_from = abund) %>%
  column_to_rownames("Sample")


coverm <- mags_coverm %>%
  rownames_to_column("Sample")%>%
  pivot_longer(cols = -Sample, names_to = "bin", values_to = "abund")%>%
  left_join(meta, by ="Sample")


samp_func_weighted <- coverm %>%
  filter(abund > 0) %>%
  inner_join(traits01, by = "bin") %>%
  group_by(Sample, Category_2) %>%
  summarise(score = weighted.mean(perc_present, w = abund, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = Category_2, values_from = score, values_fill = 0)


X_fun_w <- samp_func_weighted %>%
  column_to_rownames("Sample") %>%
  as.matrix()

# Merge plotting metadata
meta2 <- meta %>% filter(Sample %in% rownames(X_fun_w))

# standardization then PCA via rda:
X_fun_w_std <- decostand(X_fun_w, method = "standardize")  # column-wise (functions) z-score
pca_fun   <- rda(X_fun_w_std)                            # PCA

# Extract scores for plotting (scaling=2: good for interpreting correlations among arrows)
site_scr <- scores(pca_fun, display = "sites", scaling = 2) %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% left_join(meta2, by = "Sample")

sp_scr <- scores(pca_fun, display = "species", scaling = 2) %>% as.data.frame() %>%
  rownames_to_column("Function") %>%
  mutate(len = sqrt(PC1^2 + PC2^2)) %>%
  slice_max(len, n = 15)   # plot top N arrows to reduce clutter

# Percent variance
eig   <- eigenvals(pca_fun)
varPC <- round(100 * eig / sum(eig), 1)

## Plot -----
library(ggplot2)
biplot <- ggplot(site_scr, aes(PC1, PC2, color = Type)) +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  geom_point(size = 2, alpha = 0.9) +
  stat_ellipse(aes(fill = Type), geom = "polygon", alpha = 0.15,
               type = "norm", level = 0.68, colour = NA) +
  stat_ellipse(linewidth = 0.6, type = "norm", level = 0.68, show.legend = FALSE) +
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_segment(data = sp_scr,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.02, "npc")), inherit.aes = FALSE) +
  geom_text_repel(data = sp_scr, aes(x = PC1, y = PC2, label = Function),
                  min.segment.length = 0, size = 3, inherit.aes = FALSE) +
  labs(x = paste0("PC1 (", varPC[1], "%)"),
       y = paste0("PC2 (", varPC[2], "%)")
  ) +
  theme_bw(base_size = 8, base_family = "Helvetica")+
  theme(legend.position = "right", 
        legend.text = element_text(size = 8),
        plot.margin = margin(10, 0, 1, 0, "mm"))

biplot

#ggsave("Output_figures/Functional_pathways_PCA_biplot.svg", plot, dpi = 300, width= 4.5, height = 3.5)

# Phylogenetic Tree -----

## Load required libraries
library(ggtree)
library(tidyverse)
library(ggnewscale)
library(ggtreeExtra)
library(ggpubr)
library(ape)

library(patchwork)
library(treedataverse)
library(ggh4x)


checkm2 <- read_tsv('Metagenomic_data/quality_report.tsv') %>% 
  mutate(bin_type =case_when(
    Completeness >= 70 & Contamination <= 10 ~ 'High quality',
    Completeness >= 50 & Contamination <= 10 ~ 'Medium quality',
    TRUE ~ 'Low quality'
  ))
load_mag_taxonomy <- function(){
  gtdb_bac <- read_tsv('Metagenomic_data/gtdbtk.bac120.decorated.tree-taxonomy',
                       col_names = FALSE) %>% 
    rename(bin = X1) %>% 
    filter(str_detect(bin, 'bin')) %>% 
    separate_wider_delim(X2,
                         delim = '; ',
                         names = c('Domain', 'Phylum', 'Class', 'Order',
                                   'Family', 'Genus', 'Species'),
                         too_few = "align_start")
  
  
  gtdb_ar <- read_tsv('Metagenomic_data/gtdbtk.ar53.decorated.tree-taxonomy',
                      col_names = FALSE)  %>% 
    rename(bin = X1) %>% 
    filter(str_detect(bin, 'bin')) %>% 
    separate_wider_delim(X2,
                         delim = '; ',
                         names = c('Domain', 'Phylum', 'Class', 'Order',
                                   'Family', 'Genus', 'Species'),
                         too_few = "align_start")
  
  gtdb_all <- rbind(gtdb_bac, gtdb_ar) %>% 
    mutate(across(Domain:Species, ~str_remove(.x, '.__')))
  
  return(gtdb_all)
}
mag_taxonomy <- load_mag_taxonomy() %>% 
  mutate(Phylum_lump = fct_lump_n(Phylum, 10)) %>% 
  group_by(Phylum) %>% 
  mutate(n = n()) %>% 
  mutate(Phylum_n = str_c(Phylum, ' (', n, ')'))

bac_tree <- read.tree('Metagenomic_data/gtdbtk.bac120.decorated.tree')
ar_tree <- read.tree('Metagenomic_data/gtdbtk.ar53.decorated.tree')


tree_tips <- tip.label(bac_tree)
drop_tips <- tree_tips[which(!str_detect(tree_tips, 'bin'))]

mybac_tree <- drop.tip(bac_tree, drop_tips)

ar_tree_tips <- tip.label(ar_tree)
ar_drop_tips <- ar_tree_tips[which(!str_detect(ar_tree_tips, 'bin'))]

myar_tree <- drop.tip(ar_tree, ar_drop_tips)

phyla_colors <- set_names(ggpubr::get_palette("Paired", length(unique(na.omit(mag_taxonomy$Phylum)))), 
                          sort(unique(na.omit(mag_taxonomy$Phylum))))

phyla_colors[3] = 'blue'
phyla_colors[11] = "green"

bactree_df <- as_tibble(mybac_tree) %>% 
  left_join(mag_taxonomy, by = c('label' = 'bin')) %>%
  as.treedata()

artree_df <- as_tibble(myar_tree) %>% 
  left_join(mag_taxonomy, by = c('label' = 'bin')) %>%
  as.treedata()

all_tree <- bind.tree(myar_tree, mybac_tree) %>% 
  as_tibble() %>% 
  left_join(mag_taxonomy, by = c('label' = 'bin')) %>%
  as.treedata()

gtdb_checkm <- checkm2 %>% 
  select(Name, Completeness, Contamination, Genome_Size, GC_Content) %>% 
  inner_join(mag_taxonomy, by = c('Name' = 'bin'))

## plot ----
tree_annotated <- ggtree(all_tree, layout = 'fan', open.angle = 20) +
  geom_tippoint(aes(fill = Phylum), size = 2.5, shape = 21) +
  scale_fill_manual(values = phyla_colors) +
  guides(fill = guide_legend(ncol = 5), ) +
  ggnewscale::new_scale_fill() +
  geom_fruit(data = gtdb_checkm,
             geom = geom_col,
             mapping = aes(y = Name,
                           x = Completeness,
                           color = Completeness),
             #color = 'white',
             offset = 0.05,
             pwidth = 0.2) +
  scale_color_distiller(palette = 'PuBuGn', labels = ~scales::percent(.x, scale = 1),
                        limits = c(50, 100), direction = 1)+
  ggnewscale::new_scale_fill() +
  geom_fruit(data = gtdb_checkm,
             geom = geom_col,
             mapping = aes(y = Name,
                           x = Contamination,
                           fill = Contamination),
             color = 'black',
             offset = 0.05,
             pwidth = 0.2) +
  scale_fill_distiller(palette = 'Greys', 
                       labels = ~scales::percent(.x, scale = 1), direction = 1,
                       limits = c(0, 10)) +
  ggnewscale::new_scale_fill() +
  geom_fruit(data = gtdb_checkm,
             geom = geom_col,
             mapping = aes(y = Name,
                           x = Genome_Size,
                           fill = Genome_Size),
             color = 'purple4',
             offset = 0.05,
             pwidth = 0.5) +
  scale_fill_distiller(palette = 'Purples', direction = 1,
                       labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.direction = 'vertical',
        legend.position = 'bottom',
        legend.key.width = unit(2, 'lines'),
        legend.title.position = 'top',
        legend.spacing = unit(.5, 'lines'),
        legend.key.spacing.y = unit(0.05, 'lines'),
        legend.key.spacing.x = unit(.5, 'lines'),
        legend.key.size = unit(0.3, 'lines')
  )


tree_annotated

# Show and save

#ggsave('output_figures/Ext_figure_2_version2.svg', tree_annotated, dpi = 300, width = 200, height = 185, units = 'mm')
#ggsave('output_figures/Ext_figure_2_version2.png', tree_annotated, dpi = 300, width = 200, height = 185, units = 'mm')

tree_annotated <- ggtree(all_tree, layout = 'fan', open.angle = 20) +
  geom_tippoint(aes(fill = Phylum), size = 3, shape = 21) +
  scale_fill_manual(values = phyla_colors) +
  guides(fill = guide_legend(ncol = 3, order = 1)) +
  ggnewscale::new_scale_fill() +
  geom_fruit(data = gtdb_checkm,
             geom = geom_col,
             mapping = aes(y = Name,
                           x = Completeness,
                           color = Completeness),
             #color = 'white',
             offset = 0.05,
             pwidth = 0.2) +
  scale_color_distiller(palette = 'PuBuGn', labels = ~scales::percent(.x, scale = 1),
                        limits = c(50, 100), , direction = 1) +
  guides(fill = guide_legend(order = 2)) +
  ggnewscale::new_scale_color() +
  geom_fruit(data = gtdb_checkm,
             geom = geom_col,
             mapping = aes(y = Name,
                           x = Contamination,
                           color = Contamination),
             #color = 'grey',
             offset = 0.05,
             pwidth = 0.1) +
  scale_color_distiller(palette = 'Greys', 
                        labels = ~scales::percent(.x, scale = 1), direction = 1,
                        limits = c(0, 10)) +
  guides(fill = guide_legend(order = 3)) +
  ggnewscale::new_scale_color() +
  geom_fruit(data = gtdb_checkm,
             geom = geom_col,
             mapping = aes(y = Name,
                           x = Genome_Size,
                           color = Genome_Size),
             #color = 'white',
             offset = 0.05,
             pwidth = 0.3) +
  scale_color_distiller(palette = 'Purples', direction = 1,
                        labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  guides(fill = guide_legend(order = 4)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        #legend.key.size = unit(1, 'lines'),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        #legend.direction = 'horizontal',
        #legend.position = 'bottom',
        legend.direction = "horizontal",         # <- key
        legend.position = "right",
        legend.box = "vertical",             # lay two groups side-by-side
        legend.key.width = unit(2, 'lines'),
        legend.title.position = 'top',
        legend.spacing = unit(.05, 'lines'),
        legend.key.spacing.y = unit(0.05, 'lines'),
        legend.key.spacing.x = unit(.5, 'lines'),
        legend.box.spacing = unit(0, 'lines'),
        plot.margin = margin(0, 0, 0, 0, "mm"),
        legend.key.size = unit(0.3, 'lines'))

tree_annotated

#ggsave('output_figures/Ext_figure_2_version2.svg', tree_annotated, dpi = 300, width = 150, height = 150, units = 'mm', limitsize = FALSE)

tree_core <- tree_annotated + theme(legend.position = "none",
                                    plot.margin = margin(0, 0, 0, 0, "mm"))

tree_core
## 2) Make a version that shows ONLY the Phylum (fill) legend
tree_phyla <- ggtree(all_tree, layout = 'fan', open.angle = 20) +
  geom_tippoint(aes(fill = Phylum), size = 3, shape = 21) +
  scale_fill_manual(values = phyla_colors) +
  guides(fill = guide_legend(ncol = 3, order = 1))+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal",         # <- key
        legend.position = "bottom",
        legend.box = "vertical",             # lay two groups side-by-side
        legend.key.width = unit(2, 'lines'),
        legend.title.position = 'top',
        legend.spacing = unit(.05, 'lines'),
        legend.key.spacing.y = unit(0.05, 'lines'),
        legend.key.spacing.x = unit(.5, 'lines'),
        legend.box.spacing = unit(0, 'lines'),
        plot.margin = margin(0, 0, 0, 0, "mm"),
        legend.key.size = unit(0.3, 'lines'))

phylum_leg <- cowplot::get_legend(tree_phyla)

## 3) Make a version that shows ONLY the three colourbars
tree_circles <- ggtree(all_tree, layout = 'fan', open.angle = 20) +
  geom_tippoint(aes(fill = Phylum), size = 3, shape = 21) +
  scale_fill_manual(values = phyla_colors) +
  guides(fill = "none") +
  ggnewscale::new_scale_fill() +
  geom_fruit(data = gtdb_checkm,
             geom = geom_col,
             mapping = aes(y = Name,
                           x = Completeness,
                           color = Completeness),
             #color = 'white',
             offset = 0.05,
             pwidth = 0.2) +
  scale_color_distiller(palette = 'PuBuGn', labels = ~scales::percent(.x, scale = 1),
                        limits = c(50, 100), , direction = 1) +
  guides(fill = guide_legend(order = 2)) +
  ggnewscale::new_scale_color() +
  geom_fruit(data = gtdb_checkm,
             geom = geom_col,
             mapping = aes(y = Name,
                           x = Contamination,
                           color = Contamination),
             #color = 'grey',
             offset = 0.05,
             pwidth = 0.1) +
  scale_color_distiller(palette = 'Greys', 
                        labels = ~scales::percent(.x, scale = 1), direction = 1,
                        limits = c(0, 10)) +
  guides(fill = guide_legend(order = 3)) +
  ggnewscale::new_scale_color() +
  geom_fruit(data = gtdb_checkm,
             geom = geom_col,
             mapping = aes(y = Name,
                           x = Genome_Size,
                           color = Genome_Size),
             #color = 'white',
             offset = 0.05,
             pwidth = 0.3) +
  scale_color_distiller(palette = 'Purples', direction = 1,
                        labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  guides(fill = guide_legend(order = 4)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        #legend.key.size = unit(1, 'lines'),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        #legend.direction = 'horizontal',
        #legend.position = 'bottom',
        legend.direction = "horizontal",         # <- key
        legend.position = "right",
        legend.box = "vertical",             # lay two groups side-by-side
        legend.key.width = unit(2, 'lines'),
        legend.title.position = 'top',
        legend.spacing = unit(.1, 'lines'),
        legend.key.spacing.y = unit(0.1, 'lines'),
        legend.key.spacing.x = unit(.5, 'lines'),
        legend.box.spacing = unit(0.1, 'lines'),
        plot.margin = margin(0, 0, 0, 0, "mm"),
        legend.key.size = unit(0.2, 'lines'))

other_leg <- cowplot::get_legend(tree_circles)

## 4) Assemble: top row A|B, bottom row (Phylum | other)
top_row <- cowplot::plot_grid(
  tree_core,                                    # panel A
  plot,       # panel B
  ncol = 2, 
  rel_heights = c(1.5, 1),
  rel_widths = c(1, 1.3),
  labels = c("A", "B"), label_size = 9, label_fontface = "bold"
)

legend_row <- cowplot::plot_grid(phylum_leg, other_leg, ncol = 2, rel_widths = c(2, 1))

final <- cowplot::plot_grid(top_row, legend_row, ncol = 1, rel_heights = c(1, 0.3))


legL <- cowplot::ggdraw(phylum_leg) +
  theme(plot.margin = margin(t = 4, r = 0, b = 0, l = 0, unit = "mm"))

legR <- cowplot::ggdraw(other_leg) +
  theme(plot.margin = margin(t = 4, r = 0, b = 0, l = 0, unit = "mm"))

legend_row <- cowplot::plot_grid(legL, legR, ncol = 2, 
                                 rel_widths = c(1.5, 0.75), 
                                 align = "h")

# Add extra space *under B* using a spacer stacked below the PCA:
right_col <- biplot / patchwork::plot_spacer() +
  plot_layout(heights = c(1, 0.10))  # <- 10% extra gap under B (tune as needed)

## --- 3) make sure the *top* row doesn’t carry a huge bottom margin ---
top_row_clean <- cowplot::ggdraw() +
  cowplot::draw_plot(
    cowplot::plot_grid(
      tree_core,                                     # panel A
      right_col,        # panel B
      ncol = 2, 
      rel_widths = c(1, 1),
      labels = c("A", "B"), label_size = 9, label_fontface = "bold"
    ),
    x = 0, y = 0, width = 1, height = 1
  ) +
  theme(plot.margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "mm"))  # tiny gap


## --- 4) stack rows with vertical alignment and enough space for legends ---
final <- cowplot::plot_grid(
  top_row_clean,
  legend_row,
  ncol = 1, 
  rel_heights = c(1.2, 0.5), 
  align = "v"
)
final

ggsave("output_figures/Fig3_3.png", final,
       dpi = 300, width = 7.25, height = 4.5, units = "in", limitsize = FALSE, bg = "white")

# Stress traits ------
library(tidyverse)
library(ggplot2)
library(ggh4x)

# Load data
rule_matrix <- read_tsv("Metagenomic_data/traits_hmm_hits.tsv")       
trait_defs <- read_csv("Metagenomic_data/microtrait_trait_definitions.csv")  

# Merge trait rules with definitions and taxonomy
trait_rules_long <- rule_matrix %>%
  pivot_longer(-id, names_to = "rule", values_to = "present")%>%
  rename("bin" = "id")%>%
  mutate(present = as.logical(present))  # Convert 0/1 to FALSE/TRUE

trait_data <- trait_rules_long %>%
  rename(trait = rule) %>%
  left_join(trait_defs, by = c("trait" = "microtrait_rule-name")) %>%
  left_join(mag_taxonomy, by = "bin") %>%
  filter(present == TRUE)

# Filter MAGs based on abundance thresholds
prev_min_counts <- 0
prev_min_samples <- 2
mags_coverm <- mags_coverm[, colSums(mags_coverm > prev_min_counts) >= prev_min_samples]

abund_long <- as.data.frame(t(mags_coverm)) %>%
  rownames_to_column("bin") %>%
  pivot_longer(-bin, names_to = "Sample", values_to = "count") %>%
  left_join(meta, by = "Sample") %>%
  filter(count > 0)

# Join trait presence with abundance
trait_abund_data <- trait_data %>%
  inner_join(abund_long, by = "bin")  # adds count and Type
trait_abund_data <- trait_abund_data %>%
  filter(!is.na(`microtrait_trait-name3`)) %>%
  separate(`microtrait_trait-name3`, into = c("Category1", "Category2", "Category3", "Trait"),
           sep = ":", fill = "right", remove = FALSE)

trait_abund_stress <- trait_abund_data %>%
  filter(Category1 == "Stress Tolerance")
trait_abund_summary <- trait_abund_stress %>%
  group_by(Type, Phylum_lump, trait, Category3) %>%
  summarise(total_abundance = sum(count, na.rm = TRUE), .groups = "drop") %>%
  mutate(prop_abundance = total_abundance / sum(total_abundance)) %>%
  ungroup()

trait_abund_summary %>%
  arrange(Type, desc(prop_abundance)) %>%
  head(20)  # or View() in RStudio

stress_heatmap <- ggplot(trait_abund_summary, 
                         aes(x = Phylum_lump, y = fct_rev(trait), fill = prop_abundance)) +
  geom_tile(color = "white") +
  facet_grid2(rows = vars(Category3), cols = vars(Type),
              scales = "free_y", space = "free_y",
              switch = "y", labeller = label_wrap_gen(width = 20)) +
  scale_fill_gradient(
    low = "pink", 
    high = "blue",
    na.value = "white",
    limits = c(0, max(trait_abund_summary$prop_abundance, na.rm = TRUE)),
    labels = scales::percent_format(accuracy = 0.5)
  ) +
  labs(x = "Phylum", y = "Stress Trait", fill = "Abundance %") +
  theme_bw(base_size = 8, base_family = "Halvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(face = "italic"),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
    strip.text.x = element_text(size = 8),
    strip.placement = "outside",
    legend.position = "right",
    legend.key.width = unit(1, 'lines'),
    panel.spacing.x = unit(0.6, "lines"),
    panel.spacing.y = unit(0.2, "lines"),panel.grid.major = element_blank()
  )

stress_heatmap
ggsave("output_figures/microtrait_stress_tolerance.svg", stress_heatmap, dpi = 300, width = 8.5, height = 5)


# Figure 4 -------------------

## 4) Assemble: top row A|B, bottom row (Phylum | other)

legL <- cowplot::ggdraw(phylum_leg) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))

legR <- cowplot::ggdraw(other_leg) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))

legend_row <- cowplot::plot_grid(legL, legR, ncol = 2, 
                                 rel_widths = c(1.5, 0.75), 
                                 align = "h")

# Add extra space *under B* using a spacer stacked below the PCA:
right_col <- biplot / patchwork::plot_spacer() +
  plot_layout(heights = c(1, 0.02)) 

## --- 3) make sure the *top* row doesn’t carry a huge bottom margin ---
top_row_clean <- cowplot::ggdraw() +
  cowplot::draw_plot(
    cowplot::plot_grid(
      tree_core,                                     # panel A
      right_col,        # panel B
      ncol = 2, 
      rel_widths = c(1, 1),
      labels = c("A", "B"), label_size = 9, label_fontface = "bold"
    ),
    x = 0, y = 0, width = 1, height = 1
  ) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))  # tiny gap

## --- 3) Bottom row: (C) stress heatmap ---------------------------------------
bottom_row <- cowplot::plot_grid(
  stress_heatmap,
  labels = "C",
  label_size = 9,
  label_fontface = "bold"
)

# --- 4) Stack all rows --------------------------------------------------------
final <- cowplot::plot_grid(
  top_row_clean,
  legend_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 0.4, 1),  # tune vertical proportions
  align = "v"
)

final

ggsave("output_figures/Fig3_final.png", final,
       dpi = 300, width = 7.25, height = 9, units = "in", bg = "white", limitsize = FALSE)

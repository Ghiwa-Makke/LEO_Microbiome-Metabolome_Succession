# 1. Loading libraries

library(patchwork)
library(ggpubr)
library(igraph)
library(ggfx)
library(ggh4x)
library(ggnetwork)
library(readxl)
library(MetaNet)
library(tidyverse)
library(ggClusterNet)

## Load data

colors <- c('Bare' = '#3d1d0acc', 'Crust' = '#a8581bcc', 'Moss' = '#116311cc')

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
mag_taxonomy <- load_mag_taxonomy()

bin_code <- read_csv('Metagenomic_data/mag_network_names_updated.csv')

mag_phyl <- mag_taxonomy %>%
  mutate(Phylum = str_remove(Phylum, 'p__'),
         Phylum_lump = fct_lump_n(Phylum, n = 12))
mag_phyl$bin <- gsub("[_\\.]", "", mag_phyl$bin)


mag_phyl <- bin_code %>%
  left_join(mag_phyl, by = c('Genome' = 'bin'))

phyla_colors <- set_names(get_palette('Paired', 14),
                          nm =sort(unique(mag_phyl$Phylum_lump)))

# MAG-MAG co-occurance --------------

# Default_No_iDirect

## Nodes ----

nodes_bare <- read_delim('Metagenomic_data/Bare_node_attribute.txt') %>% mutate(net = 'Bare')
nodes_crust <- read_delim('Metagenomic_data/Crust_node_attribute.txt') %>% mutate(net = 'Crust')
nodes_moss <- read_delim('Metagenomic_data/Moss_node_attribute.txt') %>% mutate(net = 'Moss')

nodes_all <- bind_rows(nodes_bare, nodes_crust, nodes_moss) %>% 
  inner_join(bin_code, by = c('Name' = 'ID'))

## Edges ------

edges_bare <- read_delim('Metagenomic_data/Bare_edge_attribute.txt') %>% 
  select(from, to, interaction, value) %>% mutate(net = 'Bare')
edges_crust <- read_delim('Metagenomic_data/Crust_edge_attribute.txt') %>% 
  select(from, to, interaction, value) %>% mutate(net = 'Crust')
edges_moss <- read_delim('Metagenomic_data/Moss_edge_attribute.txt') %>% 
  select(from, to, interaction, value) %>% mutate(net = 'Moss')

edges_all <- bind_rows(edges_bare, edges_crust, edges_moss) %>% 
  mutate(interaction = ifelse(interaction == 'pp', 'Positive', 'Negative'))

## 2.3 Networks
# Define roles of interest
core_roles <- c("Module hubs", "Connectors", "Network hubs")

nodes_ready <- nodes_all %>% 
  group_by(`No. module`) %>% 
  mutate(nodes_per_module = n(),
         module = paste0('Module ', `No. module`),
         module = factor(module, levels = paste0('Module ', 0:16)),
         node_topological_role = case_when(Zi <= 2.5 & Pi <= 0.62 ~ "Peripherals",
                                           Zi <= 2.5 & Pi > 0.62 ~ "Connectors",
                                           Zi > 2.5 & Pi <= 0.62 ~ "Module hubs",
                                           Zi > 2.5 & Pi > 0.62 ~ "Network hubs")) %>% 
  inner_join(mag_phyl, by = 'Genome') %>% 
  ungroup()

summary_def_no_id <- nodes_ready %>%
  count(net, module, node_topological_role) %>%
  group_by(net) %>%
  mutate(n_modules = n_distinct(module)) %>%
  ungroup() %>%
  group_by(net, node_topological_role, n_modules) %>%
  summarise(n_nodes = sum(n), .groups = "drop") %>%
  mutate(dataset = 'default_no_iDirect')  # Replace with correct label in each section

# Filter from your full node dataset (adjust name as needed)
core_nodes_def_no_id  <- nodes_ready %>%
  filter(node_topological_role %in% core_roles)%>%
  select(Name, ID, Domain, Phylum, Class, Phylum_lump, node_topological_role) %>%
  mutate(dataset = 'default_no_iDirect')

module_sizes_def_no_id <- nodes_ready %>%
  count(net, module, name = "n_nodes_in_module")%>%
  mutate(dataset = 'default_no_iDirect')

i <- 1
module_colors <- map_chr(levels(nodes_ready$module), function(mod){
  n <- nodes_ready %>% 
    filter(module == mod) %>% 
    pull(nodes_per_module) %>% 
    unique
  
  if(n > 10){
    col <- get_palette('Dark2', 6)[i]
    i <<- i +1
  } else {
    col = 'transparent'
  }
  
  return(col)
})

names(module_colors) <- levels(nodes_ready$module)

phyl_colors <- set_names(get_palette('Paired', 18),
                         nm = sort(unique(nodes_ready$Phylum)))

nodes_split <- nodes_ready %>% 
  split(.$net)

edges_split <- edges_all %>% 
  split(.$net)


# Filter nodes and edges for Bare
nod_bare <- nodes_split$Bare %>%
  select(Name, Phylum_lump, node.degree, node_topological_role, module) %>%
  distinct(Name, .keep_all = TRUE)

edg_bare <- edges_split$Bare %>%
  select(from, to, interaction, value)%>%
  filter(from %in% nod_bare$Name, to %in% nod_bare$Name)

# Create igraph object
graph_bare <- graph_from_data_frame(edg_bare, directed = FALSE, vertices = nod_bare)
# Filter nodes and edges for Crust
nod_Crust <- nodes_split$Crust %>%
  select(Name, Phylum_lump, node.degree, node_topological_role, module) %>%
  distinct(Name, .keep_all = TRUE)

edg_Crust <- edges_split$Crust %>%
  select(from, to, interaction, value)%>%
  filter(from %in% nod_Crust$Name, to %in% nod_Crust$Name)

# Create igraph object
graph_Crust <- graph_from_data_frame(edg_Crust, directed = FALSE, vertices = nod_Crust)

# Filter nodes and edges for Moss
nod_Moss <- nodes_split$Moss %>%
  select(Name, Phylum_lump, node.degree, node_topological_role, module) %>%
  distinct(Name, .keep_all = TRUE)

edg_Moss <- edges_split$Moss %>%
  select(from, to, interaction, value)%>%
  filter(from %in% nod_Moss$Name, to %in% nod_Moss$Name)

# Create igraph object
graph_Moss <- graph_from_data_frame(edg_Moss, directed = FALSE, vertices = nod_Moss)

# Get ggnetwork layout for each graph and tag with cover type
layout_bare <- ggnetwork(graph_bare) %>% mutate(net = "Bare")
layout_crust <- ggnetwork(graph_Crust) %>% mutate(net = "Biocrust")
layout_moss <- ggnetwork(graph_Moss) %>% mutate(net = "Moss")

# Combine into one layout
layout_all <- bind_rows(layout_bare, layout_crust, layout_moss)


### Network Plot -----
single_net_plot <- layout_all %>% 
  mutate(node_size = node.degree/max(node.degree) * .5) %>% 
  ggplot(aes(x = x,
             y = y,
             xend = xend,
             yend = yend)) +
  ggforce::geom_mark_ellipse(aes(group = module,
                                 fill = module),
                             color = 'transparent',
                             show.legend = TRUE,
                             alpha = .1) +
  scale_fill_manual(values = module_colors) +
  with_outer_glow(
    geom_nodes(data = . %>% filter(node_topological_role == 'Connectors'),
               aes(x = x,
                   y = y,
                   shape = node_topological_role,
                   size = node_size),
               show.legend = TRUE),
    colour = 'red',
    sigma = 5,
    expand = 7
  ) +
  #labs(title = paste0(name, ' months')) +
  facet_wrap(~net) +
  geom_edges(aes(color = interaction)) +
  scale_color_manual(values = c('steelblue', 'indianred2')) +
  scale_size_continuous(range = c(.5, 2)) +
  ggnewscale::new_scale_fill() +
  geom_nodes(aes(fill = Phylum_lump,
                 shape = node_topological_role,
                 size = node_size),
  ) +
  scale_fill_manual(values = phyla_colors) +
  scale_shape_manual(values = c('Connectors' = 23, 
                                'Network hubs' = 24, 
                                'Peripherals' = 21)) +
  labs(title = 'Co-ocurrence networks') +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  guides(size = 'none') +
  theme_blank() +
  theme(legend.position = 'Bottom',
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))

single_net_plot


# Create dummy plot for legend extraction using layout_all
leg_plot <- layout_all %>%
  distinct(name, Phylum_lump, node_topological_role) %>%  # Ensure uniqueness
  mutate(
    TopologicalRole = node_topological_role,
    Phylum = Phylum_lump
  ) %>%
  ggplot() +
  geom_point(
    aes(x = Phylum, y = TopologicalRole, fill = Phylum, shape = TopologicalRole),
    size = 3
  ) +
  scale_shape_manual(
    values = c('Connectors' = 23, 
               'Network hubs' = 24, 
               'Peripherals' = 21),
    name = "Topological Role"
  ) +
  scale_fill_manual(
    values = phyla_colors,
    name = "Phylum"
  ) +
  # Fake correlation types for legend only
  geom_hline(aes(yintercept = 1, color = 'Negative')) +
  geom_hline(aes(yintercept = 2, color = 'Positive')) +
  scale_color_manual(
    values = c('Negative' = 'steelblue', 'Positive' = 'indianred2'),
    name = "Correlation\n(Co-occurrence network)",
    guide = guide_legend(nrow = 4)
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21), nrow = 4),
    shape = guide_legend(nrow = 4),
    color = guide_legend(nrow = 4)
  ) +
  theme_blank() +
  theme(
    legend.position = 'bottom',
    legend.title.position = 'top',
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.spacing.y = unit(.01, 'lines'),
    legend.key.spacing.x = unit(.5, 'lines'),
    legend.key.size = unit(0.3, 'lines')
  )
leg_plot
library(ggplot2)
library(grid)
library(gtable)

# Convert plot to grob
g <- ggplotGrob(leg_plot)

# Find which grob is the legend
legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")

# Extract it
leg_ready <- g$grobs[[legend_index]]


library(patchwork)
library(cowplot)

###Panel_A ----
final_plot_cooccurance <- plot_grid(
  single_net_plot + theme(legend.position = "none"),  # Main plot
  leg_ready,                                          # Extracted legend
  ncol = 1,                                           # Stack vertically
  rel_heights = c(1, 0.3)                             # Custom height ratio
)

# Display
print(final_plot_cooccurance)


ggsave('output_figs/Panel_A.png', final_plot_cooccurance, height = 4.5, width = 7.25, units = "in", dpi = 300, bg = "white")

# MIP -----
create_network_zipi <- function(edges, nodes){
  
  graph <- graph_from_data_frame(edges, 
                                 directed = FALSE, 
                                 vertices = nodes)
  
  node_par <- net_par(graph, mode = 'all')$v_index
  
  greedy_modules <- module_detect(graph, method = 'cluster_fast_greedy')
  
  nodes_zipi <- zp_analyse(greedy_modules) %>% 
    MetaNet::get_v() %>% 
    mutate(module = paste0('Module_', module)) %>% 
    left_join(node_par)
  
  hubs <- nodes_zipi %>% 
    filter(roles == 2) %>% 
    pull(name)
  
  updated_edges <- edges %>% 
    mutate(connection_with = ifelse(to %in% hubs | from %in% hubs, 
                                    'Link with hubs', 'Link with other MAG'))
  
  updated_graph <- graph_from_data_frame(updated_edges, 
                                         directed = FALSE, 
                                         vertices = nodes_zipi)
  
  return(updated_graph)
  
}

network_stats <- function(graph){
  
  cluster <- cluster_fast_greedy(graph)
  
  df <- net_properties(graph) %>%
    as.data.frame %>% 
    rownames_to_column(var = 'index') %>% 
    add_row(index = 'transitivity',
            value = transitivity(graph)) %>% 
    add_row(index = 'num.modules',
            value = length(cluster)) %>%
    add_row(index = 'modularity',
            value = modularity(cluster)) 
  
  return(df)
  
}
## Metabolic interactions network
global_results <- read_tsv('Metagenomic_data/smetana_pairwise_summary.tsv') %>% 
  mutate(comparison = str_remove(comparison, '_results.*\\.tsv$'),
         # Split on underscore, reconstruct MAG_1 and MAG_2
         MAG_1 = str_c(word(comparison, 1, sep = "_"), "_", word(comparison, 2, sep = "_")),
         MAG_2 = str_c(word(comparison, 3, sep = "_"), "_", word(comparison, 4, sep = "_")),
         mip = as.numeric(mip),
         mro = as.numeric(mro))

detailed_results <- read_tsv('Metagenomic_data/smetana_pairwise_detailed_summary.tsv')%>%
  mutate(receiver = str_remove(receiver, '_model*'),
         donor = str_remove(donor, '_model*')) 

receiver_counts <- detailed_results %>% 
  group_by(receiver) %>% 
  count(name = 'receiver_counts')

donor_counts <- detailed_results %>% 
  group_by(donor) %>% 
  count(name = 'donor_counts')

donor_or_receiver <- full_join(receiver_counts, donor_counts, by = c('receiver' = 'donor')) %>% 
  rename(bin = receiver) %>% 
  mutate(across(c(receiver_counts, donor_counts), ~ifelse(is.na(.x), 0, .x)),
         main_role = ifelse(receiver_counts > donor_counts, 'Receiver', 'Donor'),
         bin = str_remove(bin, '_genes.*'))

global_filtered <- global_results %>% 
  filter(mip >= 5)

edges <- global_filtered %>% 
  filter(mip >= 5) %>% 
  select(from = MAG_1, to = MAG_2, mip)

nodes <- tibble(Name = unique(c(edges$from, edges$to))) %>% 
  inner_join(mag_taxonomy, by = c('Name' = 'bin')) %>% 
  inner_join(donor_or_receiver, by = c('Name' = 'bin'))

nodes$Name <- gsub("[_\\.]", "", nodes$Name)
edges$from <- gsub("[_\\.]", "", edges$from)
edges$to <- gsub("[_\\.]", "", edges$to)

# Filtering for those in the Type networks

## default-no idirect
## co-occurance data: 
nodes_bare <- read_delim('Metagenomic_data/Bare_node_attribute.txt') %>% mutate(net = 'Bare')%>%
  inner_join(bin_code, by=c('Name'='ID'))
nodes_crust <- read_delim('Metagenomic_data/Crust_node_attribute.txt') %>% mutate(net = 'Crust')%>%
  inner_join(bin_code, by=c('Name'='ID'))
nodes_moss <- read_delim('Metagenomic_data/Moss_node_attribute.txt') %>% mutate(net = 'Moss')%>%
  inner_join(bin_code, by=c('Name'='ID'))

nodes_met_bare <- nodes %>% 
  filter(Name %in% nodes_bare$Genome)

edges_met_bare <- edges %>% 
  filter((from %in% nodes_met_bare$Name & to %in% nodes_met_bare$Name))
edges_met_bare_clean <- edges_met_bare %>%
  mutate(edge_min = pmin(from, to),
         edge_max = pmax(from, to)) %>%
  group_by(edge_min, edge_max) %>%
  summarise(mip = max(mip), .groups = "drop") %>%
  rename(from = edge_min, to = edge_max)

nodes_met_crust <- nodes %>% 
  filter(Name %in% nodes_crust$Genome)

edges_met_crust <- edges %>% 
  filter((from %in% nodes_met_crust$Name & to %in% nodes_met_crust$Name))
edges_met_crust_clean <- edges_met_crust %>%
  mutate(edge_min = pmin(from, to),
         edge_max = pmax(from, to)) %>%
  group_by(edge_min, edge_max) %>%
  summarise(mip = max(mip), .groups = "drop") %>%
  rename(from = edge_min, to = edge_max)


nodes_met_moss <- nodes %>% 
  filter(Name %in% nodes_moss$Genome)

edges_met_moss <- edges %>% 
  filter((from %in% nodes_met_moss$Name & to %in% nodes_met_moss$Name))
edges_met_moss_clean <- edges_met_moss %>%
  mutate(edge_min = pmin(from, to),
         edge_max = pmax(from, to)) %>%
  group_by(edge_min, edge_max) %>%
  summarise(mip = max(mip), .groups = "drop") %>%
  rename(from = edge_min, to = edge_max)

## 2.3 Create graph objects
bare_met_network <- create_network_zipi(edges_met_bare_clean, nodes_met_bare)
crust_met_network <- create_network_zipi(edges_met_crust_clean, nodes_met_crust)
moss_met_network <- create_network_zipi(edges_met_moss_clean, nodes_met_moss)
## 2.4 Plotting as facets
met_list <- list('Bare' = bare_met_network,
                 'Biocrust' = crust_met_network,
                 'Moss' = moss_met_network)

roles_df <- tibble(role_name = c("Peripherals", "Network hubs", "Module hubs", "Connectors"),
                   roles = 1:4)

# ----------------- BARE ------------------
nodes_bare_df <- igraph::as_data_frame(bare_met_network, what = "vertices") %>%
  mutate(Name = name,
         role_name = roles,   # Use roles directly as labels
         origin = "Bare") %>%
  filter(!is.na(role_name))

lay_bare <- layout_on_sphere(bare_met_network)
coords_bare <- as.data.frame(lay_bare)
colnames(coords_bare) <- c("x", "y")
coords_bare$Name <- V(bare_met_network)$name

df_bare <- left_join(coords_bare, nodes_bare_df, by = "Name")

# ----------------- CRUST ------------------
nodes_crust_df <- igraph::as_data_frame(crust_met_network, what = "vertices") %>%
  mutate(Name = name,
         role_name = roles,
         origin = "Biocrust") %>%
  filter(!is.na(role_name))

lay_crust <- layout_on_sphere(crust_met_network)
coords_crust <- as.data.frame(lay_crust)
colnames(coords_crust) <- c("x", "y")
coords_crust$Name <- V(crust_met_network)$name

df_crust <- left_join(coords_crust, nodes_crust_df, by = "Name")

# ----------------- MOSS ------------------
nodes_moss_df <- igraph::as_data_frame(moss_met_network, what = "vertices") %>%
  mutate(Name = name,
         role_name = roles,
         origin = "Moss") %>%
  filter(!is.na(role_name))

lay_moss <- layout_on_sphere(moss_met_network)
coords_moss <- as.data.frame(lay_moss)
colnames(coords_moss) <- c("x", "y")
coords_moss$Name <- V(moss_met_network)$name

df_moss <- left_join(coords_moss, nodes_moss_df, by = "Name")

# ----------------- Combine all networks ------------------
met_net_df <- bind_rows(df_bare, df_crust, df_moss)



get_edge_df <- function(graph, origin_label) {
  edgelist <- igraph::as_data_frame(graph, what = "edges")
  layout_coords <- layout_on_sphere(graph)
  node_names <- igraph::V(graph)$name
  
  coords_df <- tibble(Name = node_names, x = layout_coords[, 1], y = layout_coords[, 2])
  
  edgelist %>%
    left_join(coords_df, by = c("from" = "Name")) %>%
    rename(x = x, y = y) %>%
    left_join(coords_df, by = c("to" = "Name"), suffix = c("", ".to")) %>%
    rename(xend = x.to, yend = y.to) %>%
    mutate(origin = origin_label)
}


edges_bare_df <- get_edge_df(bare_met_network, "Bare")
edges_crust_df <- get_edge_df(crust_met_network, "Biocrust")
edges_moss_df <- get_edge_df(moss_met_network, "Moss")

edges_all_df <- bind_rows(edges_bare_df, edges_crust_df, edges_moss_df)

# Normalize node size
met_net_df <- met_net_df %>%
  mutate(node_size = Degree / max(Degree, na.rm = TRUE) * 0.5)

### plot -----
met_network_plot <- ggplot() +
  # Plot edges
  geom_segment(data = edges_all_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey70", size = 0.2, alpha = 0.6) +
  
  # Highlight Network hubs with glow (optional)
  with_outer_glow(
    geom_point(data = met_net_df %>% filter(role_name == "Network hubs"),
               aes(x = x, y = y, size = node_size),
               shape = 24),
    colour = "red", sigma = 5, expand = 5
  ) +
  
  # Plot all nodes with fill by phylum and shape by role
  geom_point(data = met_net_df,
             aes(x = x, y = y, fill = Phylum, shape = role_name, size = node_size)) +
  
  facet_wrap(~origin) +
  
  # Manual shape scale
  scale_shape_manual(values = c(
    "Connectors" = 23,
    "Network hubs" = 24,
    "Peripherals" = 21
  )) +
  
  # Phylum colors
  scale_fill_manual(values = phyla_colors) +
  
  scale_size_continuous(range = c(1, 3)) +
  
  labs(title = "Metabolic Interaction Potential Network",
       x = NULL, y = NULL, fill = "Phylum", shape = "Topological Role") +
  
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(title = "Topological Role", direction = "vertical"),
         size = "none") +
  
  theme_void(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    legend.position = "bottom",
    strip.text = element_text(size = 8),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.spacing.y = unit(.01, 'lines'),
    legend.key.spacing.x = unit(.5, 'lines'),
    legend.key.size = unit(0.3, 'lines'),
    panel.spacing = unit(1, "cm"),
    
    
  )+
  coord_fixed(ratio = 1.1)

# Display the plot
met_network_plot
ggsave('output_figs/Panle_B.png', met_network_plot, height = 4.5, width = 7.25, units = "in", dpi = 300, bg= "white")


# Figure 5 -----

library(patchwork)
library(cowplot)

figure3 <- (final_plot_cooccurance /  met_network_plot)+  
  plot_layout(height = c(1,0.8)) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(size = 9, face = "bold")
  )
figure3

  
ggsave('Fig5.png', figure3, height = 7, width = 7.25, units = "in", dpi = 300, bg = "white")



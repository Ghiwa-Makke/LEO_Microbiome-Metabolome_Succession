##### Microbiome-Metabolome Interactions Reshape Ecosystem Properties During Biocrust-Moss Succession
##### LEO - surface samples - LC -RP 
### Van Krevelen Diagrams from Machine learning Outputs

# 1. Importing Libraries
suppressPackageStartupMessages({  
  library(tidyverse)
  library(ggsci)
}) 

# 2. Import data

project_dir <- getwd()
project_name <- 'LEO_RP'

## Imprt ML lists
compounds_table <- read.csv("output_tables/cmpd_RP.csv") 
compounds_table <- compounds_table %>%
  select(FeatureID, Name, Final_formula, Formula, Annotation_Confidence, H_to_C, O_to_C, Class, AI_mod, AI)
## Bare
ML_Bare <- read.csv("Input_files/ML/Bare_ML_ft.csv")%>%
  rename(FeatureID = Bare_ML_ft)
Bare <- left_join(ML_Bare, compounds_table, by = "FeatureID") %>%
  filter(!is.na(H_to_C))

Bare$Type <- "Bare"
## Biocrust
ML_Biocrust <- read.csv("Input_files/ML/Biocrust_ML_ft.csv")%>%
  rename(FeatureID = Biocrust_ML_ft)
Biocrust <- left_join(ML_Biocrust, compounds_table, by = "FeatureID")%>%
  filter(!is.na(H_to_C))

Biocrust$Type <- "Biocrust"
## Moss
ML_Moss <- read.csv("Input_files/ML/Moss_ML_ft.csv")%>%
  rename(FeatureID = Moss_ML_ft)
Moss <- left_join(ML_Moss, compounds_table, by = "FeatureID")%>%
  filter(!is.na(H_to_C))

Moss$Type <- "Moss"

# Define the slope and intercept for the first oblique line
slope_1 <- -0.47
intercept_1 <- 1.09

# Define the slope and intercept for the second oblique line
slope_2 <- -0.6
intercept_2 <- 1.5


# Contour plot ------------------------------------------------------------
## Determine shared axis limits, starting from 0
x_limits <- c(0, max(c(Bare$O_to_C, Biocrust$O_to_C, Moss$O_to_C), na.rm = TRUE))
y_limits <- c(0, max(c(Bare$H_to_C, Biocrust$H_to_C, Moss$H_to_C), na.rm = TRUE))


library(ggplot2)

# Define the color palette for the compound classes
class_colors <- c("Carbohydrate" = "blue3", "Lignin" = "red3", "Lipid" = "green4", 
                  "Protein" = "purple", "Tannin" = "orange3", "Unsaturated HC" = "skyblue", "Other" = "grey")  # Replace with your actual Class names and desired colors

# Plot for Bare
p_bare <- ggplot(data = Bare, aes(x = O_to_C, y = H_to_C)) +
  geom_point(aes(color = Class), size = 2) +
  coord_cartesian(xlim = x_limits, ylim = y_limits) + # Apply shared limits
  geom_density2d(aes(fill = ..level..))+
  scale_color_manual(values = class_colors) +  # Apply fixed color scale
  theme_bw(base_family = "helvatica", base_size = 8) +
  labs(title = "Bare",
       x = "O:C ratio",
       y = "H:C ratio")+
  #geom_text_repel(aes(label = FeatureID), size = 2, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10)+ # Add labels with repelling # Add labels
  geom_abline(slope = slope_1, intercept = intercept_1, color = "blue", linetype = "dashed") + # Add first oblique line
  geom_abline(slope = slope_2, intercept = intercept_2, color = "red", linetype = "dashed")+
  annotate("text", x = 1.6, y = 0.1, label = "AI > 0.5", color = "blue", size = 3.5) +  # Add text for blue line
  annotate("text", x = 1.9, y = 0.6, label = "AI < 0.3", color = "red", size = 3.5)+    # Add text for red line
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.direction = 'vertical',
        legend.position = 'right',
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
  )

p_bare
# Plot for Biocrust
p_crust <- ggplot(data = Biocrust, aes(x = O_to_C, y = H_to_C)) +
  geom_point(aes(color = Class), size = 1.5) +
  geom_density2d(aes(fill = ..level..))+
  coord_cartesian(xlim = x_limits, ylim = y_limits) + # Apply shared limits
  scale_color_manual(values = class_colors) +  # Apply fixed color scale
  theme_bw(base_family = "helvatica", base_size = 8) +
  labs(title = "Biocrust",
       x = "O:C ratio",
       y = "H:C ratio")+
  geom_abline(slope = slope_1, intercept = intercept_1, color = "blue", linetype = "dashed") + # Add first oblique line
  geom_abline(slope = slope_2, intercept = intercept_2, color = "red", linetype = "dashed")+
  annotate("text", x = 1.6, y = 0.1, label = "AI > 0.5", color = "blue", size = 3.5) +  # Add text for blue line
  annotate("text", x = 1.9, y = 0.6, label = "AI < 0.3", color = "red", size = 3.5)+    # Add text for red line
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.direction = 'vertical',
        legend.position = 'right',
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
  )

# Plot for Moss
p_moss <- ggplot(data = Moss, aes(x = O_to_C, y = H_to_C)) +
  geom_point(aes(color = Class), size = 2) +
  geom_density2d(aes(fill = ..level..))+
  coord_cartesian(xlim = x_limits, ylim = y_limits) + # Apply shared limits
  scale_color_manual(values = class_colors) +  # Apply fixed color scale
  theme_bw(base_family = "helvatica", base_size = 8) +
  labs(title = "Moss",
       x = "O:C ratio",
       y = "H:C ratio")+
  geom_abline(slope = slope_1, intercept = intercept_1, color = "blue", linetype = "dashed") + # Add first oblique line
  geom_abline(slope = slope_2, intercept = intercept_2, color = "red", linetype = "dashed")+
  annotate("text", x = 1.6, y = 0.1, label = "AI > 0.5", color = "blue", size = 3.5) +  # Add text for blue line
  annotate("text", x = 1.9, y = 0.6, label = "AI < 0.3", color = "red", size = 3.5)+    # Add text for red line
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.direction = 'vertical',
        legend.position = 'right',
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
  )
p_moss
# Display plots
p_bare
p_crust
p_moss
library(patchwork) # For combining plots
# Combine the plots into a single row
combined_plot <- p_bare + p_crust + p_moss + 
  plot_layout(nrow = 3) # Set the number of columns

# Display the combined plot
combined_plot

# Save the combined plot as a PNG file
ggsave("output_figures/SupFigS7.png", plot = combined_plot, width = 6, height = 7, dpi = 300)



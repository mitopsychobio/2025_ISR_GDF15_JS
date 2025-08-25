# script to get the average of the spearman rhos of each Factor# vs AGE. Then to compare against GDF15 vs Age, CHOP, ATF4, and others. 

#This code currently gets the log2 TMM expression of GDF15, then centers and scales it, and then compared the pc2 score (which also gets centerd and scaled) for each tissue

# Uncommented out many of the "print" plots

rm(list = ls())

# BiocManager::install("edgeR")

library(tidyverse)
library(edgeR)
library(corrr)
library(dplyr)
library(tibble)
library(ggpubr)
library(readr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(beepr)
library(ggplot2)
library(ggpubr)
library(FSA)      # For Dunn's test
library(dplyr)

# Helper to get the *current script* path (works in RStudio, Rscript, or knitr)
get_script_path <- function() {
  # RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  # knitr / Quarto
  if (requireNamespace("knitr", quietly = TRUE)) {
    p <- tryCatch(knitr::current_input(), error = function(e) "")
    if (nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  # Rscript
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  # Fallback: unknown
  ""
}

# ---- 1) Use script name as Folder_name (fallback to date if unknown) ----
script_path <- get_script_path()
Folder_name <- if (nzchar(script_path)) {
  sub("\\.[^.]*$", "", basename(script_path))
} else {
  format(Sys.Date(), "%Y_%m_%d")
}

############################################################################################

folder_path <- here("Results", "gtex", Folder_name)
if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)


script_path <- here("Scripts", "gtex", "Attibutes_Phenos_Merged_plus_COD.R")
source(script_path)


# Set the path to the main folder
CSV_folder_path <- here("Data", "gtex", "Age_comparisons", "Comparing_AGE_Other_Factors_BH")

# List all subdirectories (12 folders)
subfolders <- list.dirs(CSV_folder_path, full.names = TRUE, recursive = FALSE)

# Initialize an empty list to store the dataframes
data_list <- list()

# Loop through each subfolder to read the specific CSV file
for (subfolder in subfolders) {
  # Find the CSV file in the current subfolder
  csv_file <- list.files(subfolder, pattern = "\\.csv$", full.names = TRUE)
  
  # Check if there’s exactly one CSV file, then read it
  if (length(csv_file) == 1) {
    # Extract the file name without the extension to use as the dataframe name
    file_name <- tools::file_path_sans_ext(basename(csv_file))
    
    # Read the CSV file and store it in the list with the file name as key
    data_list[[file_name]] <- read.csv(csv_file)
  } else {
    message(paste("Skipping folder", subfolder, "due to multiple or no CSV files."))
  }
}

# Initialize an empty dataframe for the result, starting with the "Tissue" column
All_Spearmans <- data.frame(Tissue = unique(data_list[[1]]$Tissue))

# Loop through each dataframe in data_list
for (name in names(data_list)) {
  # Get the current dataframe
  df <- data_list[[name]]
  
  # Find the column with "Spearman_Rho_chosen" in its name
  p_value_col <- grep("Spearman_Rho_chosen", colnames(df), value = TRUE)
  
  # Check if we found the exact column and "Tissue" column exists
  if (length(p_value_col) == 1 && "Tissue" %in% colnames(df)) {
    # Select the relevant columns (Tissue and Spearman_Rho_chosen column)
    df_selected <- df[, c("Tissue", p_value_col)]
    
    # Rename the Spearman_Rho_chosen column to the CSV file name
    colnames(df_selected)[2] <- name
    
    # Merge this dataframe with All_Spearmans by "Tissue"
    All_Spearmans <- merge(All_Spearmans, df_selected, by = "Tissue", all = TRUE)
  } else {
    message(paste("Skipping file", name, "due to missing or multiple 'Spearman_Rho_chosen' columns."))
  }
}


# Set the path to the main folder
ATF4_csv <- here("Data", "gtex", "Age_comparisons", "Comparing_AGE_Other_Factors_BH_ATF4","2024_10_30_All_Factors_AGE_ATF4", "Factor1_vs_AGE_SpearmanRhos.csv")
ATF4_df <- read.csv(ATF4_csv)
# Merge the Spearman_Rho_ATF4 column from ATF4_df into All_Spearmans by Tissue
All_Spearmans <- merge(All_Spearmans, ATF4_df[, c("Tissue", "Spearman_Rho_ATF4")], by = "Tissue", all.x = TRUE)



# Set the path to the main folder
DDIT3_csv <- here("Data", "gtex", "Age_comparisons", "Comparing_AGE_Other_Factors_BH_DDIT3", "2024_10_30_All_Factors_AGE_DDIT3", "Factor1_vs_AGE_SpearmanRhos.csv")
DDIT3_df <- read.csv(DDIT3_csv)
# Merge the Spearman_Rho_DDIT3 column from DDIT3_df into All_Spearmans by Tissue
All_Spearmans <- merge(All_Spearmans, DDIT3_df[, c("Tissue", "Spearman_Rho_DDIT3")], by = "Tissue", all.x = TRUE)


# Set the path to the main folder
ATF5_csv <-  here("Data", "gtex", "Age_comparisons", "Comparing_AGE_Other_Factors_BH_ATF5", "2024_10_30_All_Factors_AGE_ATF5", "Factor1_vs_AGE_SpearmanRhos.csv")
ATF5_df <- read.csv(ATF5_csv)
# Merge the Spearman_Rho_ATF5 column from ATF5_df into All_Spearmans by Tissue
All_Spearmans <- merge(All_Spearmans, ATF5_df[, c("Tissue", "Spearman_Rho_ATF5")], by = "Tissue", all.x = TRUE)

GDF15_csv <-  here("Data", "gtex", "Age_comparisons", "Comparing_AGE_Other_Factors_BH", "2024_10_28_All_Factors_AGE_Factor1", "Factor1_vs_AGE_SpearmanRhos.csv")
GDF15_df <- read.csv(GDF15_csv)
All_Spearmans <- merge(All_Spearmans, GDF15_df[, c("Tissue", "Spearman_Rho_GDF15")], by = "Tissue", all.x = TRUE)

Proliferation_csv <-  here("Data", "gtex", "Age_comparisons", "2024_12_02_All_Factors_AGE_AvgOf3ProliferationGenes", "Factor1_vs_AGE_SpearmanRhos.csv")
Proliferation_df <- read.csv(Proliferation_csv)
Proliferation_df <- Proliferation_df %>%
  rename(Spearman_Rho_Proliferation = Spearman_Rho)
All_Spearmans <- merge(All_Spearmans, Proliferation_df[, c("Tissue", "Spearman_Rho_Proliferation")], by = "Tissue", all.x = TRUE)

Senescence_csv <-  here("Data", "gtex", "Age_comparisons", "2024_12_02_All_Factors_AGE_AvgOf3SenescentGenes", "Factor1_vs_AGE_SpearmanRhos.csv")
Senescence_df <- read.csv(Senescence_csv)
Senescence_df <- Senescence_df %>%
  rename(Spearman_Rho_Senescence = Spearman_Rho)
All_Spearmans <- merge(All_Spearmans, Senescence_df[, c("Tissue", "Spearman_Rho_Senescence")], by = "Tissue", all.x = TRUE)


# # Calculate mean and SEM for each specified column, assuming All_Spearmans is your dataframe
# plot_columns <- c('Spearman_Rho_GDF15', 'Factor1_vs_AGE_SpearmanRhos', 'Factor2_vs_AGE_SpearmanRhos',
#                   'Factor3_vs_AGE_SpearmanRhos', 'Factor4_vs_AGE_SpearmanRhos', 'Factor5_vs_AGE_SpearmanRhos',
#                   'Factor6_vs_AGE_SpearmanRhos', 'Factor7_vs_AGE_SpearmanRhos', 'Factor8_vs_AGE_SpearmanRhos',
#                   'Factor9_vs_AGE_SpearmanRhos', 'Factor10_vs_AGE_SpearmanRhos', 'Factor11_vs_AGE_SpearmanRhos',
#                   'Factor12_vs_AGE_SpearmanRhos', 'Spearman_Rho_ATF4', 'Spearman_Rho_DDIT3', 'Spearman_Rho_ATF5')

# Calculate mean and SEM for each specified column, assuming All_Spearmans is your dataframe
plot_columns <- c('Spearman_Rho_GDF15', 'Spearman_Rho_Senescence', 'Spearman_Rho_Proliferation', 
                  'Spearman_Rho_ATF4', 'Spearman_Rho_ATF5', 'Spearman_Rho_DDIT3', 'Factor1_vs_AGE_SpearmanRhos',
                  'Factor12_vs_AGE_SpearmanRhos', 'Factor6_vs_AGE_SpearmanRhos', 'Factor2_vs_AGE_SpearmanRhos',
                  'Factor11_vs_AGE_SpearmanRhos', 'Factor4_vs_AGE_SpearmanRhos', 'Factor9_vs_AGE_SpearmanRhos',
                  'Factor8_vs_AGE_SpearmanRhos', 'Factor3_vs_AGE_SpearmanRhos', 'Factor10_vs_AGE_SpearmanRhos',
                  'Factor7_vs_AGE_SpearmanRhos', 'Factor5_vs_AGE_SpearmanRhos'
                  )

# Compute means and SEM for each column
means <- sapply(All_Spearmans[plot_columns], mean, na.rm = TRUE)
sems <- sapply(All_Spearmans[plot_columns], function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))

# Convert to a dataframe for plotting
plot_data <- data.frame(
  Factor = plot_columns,
  Mean = means,
  SEM = sems
)




# Define custom colors for each specific column
color_mapping <- c(
  "Spearman_Rho_GDF15" = "pink",
  "Factor1_vs_AGE_SpearmanRhos" = "red",
  "Factor2_vs_AGE_SpearmanRhos" = "#808080",
  "Factor3_vs_AGE_SpearmanRhos" = "#808080",
  "Factor4_vs_AGE_SpearmanRhos" = "#CCCCCC",
  "Factor5_vs_AGE_SpearmanRhos" = "#CCCCCC",
  "Factor6_vs_AGE_SpearmanRhos" = "#CCCCCC",
  "Factor7_vs_AGE_SpearmanRhos" = "#808080",
  "Factor8_vs_AGE_SpearmanRhos" = "#CCCCCC",
  "Factor9_vs_AGE_SpearmanRhos" = "#808080",
  "Factor10_vs_AGE_SpearmanRhos" = "#CCCCCC",
  "Factor11_vs_AGE_SpearmanRhos" = "#808080",
  "Factor12_vs_AGE_SpearmanRhos" = "#808080",
  "Spearman_Rho_ATF4" = "orange",  # tomato red
  "Spearman_Rho_DDIT3" = "lightblue",  # firebrick red
  "Spearman_Rho_ATF5" = "brown",    # dark red
  "Spearman_Rho_Senescence" = "blue",   
  "Spearman_Rho_Proliferation" = "purple"
)

# plot_data$Factor <- factor(plot_data$Factor, levels = plot_columns)

# Ensure the levels of Factor match the keys in color_mapping
plot_data$Factor <- factor(plot_data$Factor, levels = names(color_mapping))


plot_columns_csv <- c('Tissue', 'Spearman_Rho_GDF15', 'Spearman_Rho_Senescence', 'Spearman_Rho_Proliferation', 
                  'Spearman_Rho_ATF4', 'Spearman_Rho_ATF5', 'Spearman_Rho_DDIT3', 'Factor1_vs_AGE_SpearmanRhos',
                  'Factor12_vs_AGE_SpearmanRhos', 'Factor6_vs_AGE_SpearmanRhos', 'Factor2_vs_AGE_SpearmanRhos',
                 'Factor4_vs_AGE_SpearmanRhos', 'Factor9_vs_AGE_SpearmanRhos',
                  'Factor8_vs_AGE_SpearmanRhos', 'Factor3_vs_AGE_SpearmanRhos', 'Factor10_vs_AGE_SpearmanRhos',
                  'Factor7_vs_AGE_SpearmanRhos', 'Factor5_vs_AGE_SpearmanRhos',  'Factor11_vs_AGE_SpearmanRhos'
)



# Ensure Factor column has the desired order
# Reverse the order
plot_columns_csv_reversed <- rev(plot_columns_csv)

# Ensure Factor column has the reversed order
plot_data$Factor <- factor(plot_data$Factor, levels = plot_columns_csv_reversed)

# Plotting with custom colors
p2 <- ggplot(plot_data, aes(x = Mean, y = Factor, fill = Factor, color = Factor)) +
  geom_point(
    aes(shape = factor(
      case_when(
        Factor == "Factor1_vs_AGE_SpearmanRhos" ~ 8,    # Star for Factor1_vs_AGE_SpearmanRhos
        Factor == "Spearman_Rho_GDF15" ~ 22,            # Square for Spearman_Rho_GDF15
        TRUE ~ 21                                       # Circle for all other factors
      )
    )), 
    size = 6
  ) + # shape = 21 for filled points
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_errorbar(aes(xmin = Mean - SEM, xmax = Mean + SEM), width = 0.2) +
  scale_fill_manual(values = color_mapping) +  # Use color_mapping for fill colors
  scale_color_manual(values = color_mapping) +  # Use color_mapping for point outlines
  scale_x_continuous(
    breaks = seq(-1, 1, by = 0.2),  # Major grid lines at 0.2 intervals
    minor_breaks = seq(-1, 1, by = 0.1)  # Minor grid lines at 0.1 intervals
  ) +
  labs(
    x = "Spearman Rho Mean ± SEM", 
    y = "Factors and Genes", 
    title = "Average Spearman Rho Values with SEM"
  ) +
  theme_minimal() +
  # theme(
  #   axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 8),
  #   legend.position = "none"  # Remove legend if colors are self-explanatory
  # )
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),  # Black axis titles
    axis.text = element_text(size = 12, color = "black"),  # Black axis text
    axis.ticks = element_line(size = 0.25, color = "black"),  # Thin black axis ticks
    panel.grid.major.y = element_blank(),  # Remove major x gridlines
    panel.grid.minor.y = element_blank(),  # Remove minor x gridlines
    panel.grid.major.x = element_line(size = 0.1, color = "black"),  # Thin black major y gridlines
    panel.grid.minor.x = element_line(size = 0.08, color = "black"),  # Thin black minor y gridlines
    panel.background = element_blank(),  # Remove panel background
    panel.border = element_blank(),  # Remove panel border
    legend.position = "none"  # Hide legend
  )

# Print the plot
print(p2)




ggsave(filename = paste0(folder_path, "/", "SpearmanRho_Averages_with_SEM_Circle_Square_Triangle.png"), plot = p2, width = 7, height = 8, dpi = 300)

write.csv(plot_data, file = paste0(folder_path, "/", "SpearmanRho_Averages_with_SEM.csv"))







# Plotting with custom colors
p <- ggplot(plot_data, aes(x = Mean, y = Factor, fill = Factor, color = Factor)) +
  geom_point(shape = 21, size = 6) +  # shape = 21 for filled points
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_errorbar(aes(xmin = Mean - SEM, xmax = Mean + SEM), width = 0.2) +
  scale_fill_manual(values = color_mapping) +  # Use color_mapping for fill colors
  scale_color_manual(values = color_mapping) +  # Use color_mapping for point outlines
  scale_x_continuous(
    breaks = seq(-1, 1, by = 0.2),  # Major grid lines at 0.2 intervals
    minor_breaks = seq(-1, 1, by = 0.1)  # Minor grid lines at 0.1 intervals
  ) +
  labs(
    x = "Spearman Rho Mean ± SEM", 
    y = "Factors and Genes", 
    title = "Average Spearman Rho Values with SEM"
  ) +
  theme_minimal() +
  # theme(
  #   axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 8),
  #   legend.position = "none"  # Remove legend if colors are self-explanatory
  # )
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),  # Black axis titles
    axis.text = element_text(size = 12, color = "black"),  # Black axis text
    axis.ticks = element_line(size = 0.25, color = "black"),  # Thin black axis ticks
    panel.grid.major.y = element_blank(),  # Remove major x gridlines
    panel.grid.minor.y = element_blank(),  # Remove minor x gridlines
    panel.grid.major.x = element_line(size = 0.1, color = "black"),  # Thin black major y gridlines
    panel.grid.minor.x = element_line(size = 0.08, color = "black"),  # Thin black minor y gridlines
    panel.background = element_blank(),  # Remove panel background
    panel.border = element_blank(),  # Remove panel border
    legend.position = "none"  # Hide legend
  )

# Print the plot
print(p)

ggsave(filename = paste0(folder_path, "/", "SpearmanRho_Averages_with_SEM_CIRCLES.png"), plot = p, width = 7, height = 8, dpi = 300)






# Plotting with custom colors
p3 <- ggplot(plot_data, aes(x = Mean, y = Factor, fill = Factor, color = Factor)) +
  geom_point(
    aes(shape = factor(
      case_when(
        Factor == "Factor1_vs_AGE_SpearmanRhos" ~ 21,    # Star for Factor1_vs_AGE_SpearmanRhos
        # Factor == "Spearman_Rho_GDF15" ~ 22,            # Square for Spearman_Rho_GDF15
        TRUE ~ 8                                       # Circle for all other factors
      )
    )), 
    size = 6
  ) + # shape = 21 for filled points
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_errorbar(aes(xmin = Mean - SEM, xmax = Mean + SEM), width = 0.2) +
  scale_fill_manual(values = color_mapping) +  # Use color_mapping for fill colors
  scale_color_manual(values = color_mapping) +  # Use color_mapping for point outlines
  scale_x_continuous(
    breaks = seq(-1, 1, by = 0.2),  # Major grid lines at 0.2 intervals
    minor_breaks = seq(-1, 1, by = 0.1)  # Minor grid lines at 0.1 intervals
  ) +
  labs(
    x = "Spearman Rho Mean ± SEM", 
    y = "Factors and Genes", 
    title = "Average Spearman Rho Values with SEM"
  ) +
  theme_minimal() +
  # theme(
  #   axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 8),
  #   legend.position = "none"  # Remove legend if colors are self-explanatory
  # )
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),  # Black axis titles
    axis.text = element_text(size = 12, color = "black"),  # Black axis text
    axis.ticks = element_line(size = 0.25, color = "black"),  # Thin black axis ticks
    panel.grid.major.y = element_blank(),  # Remove major x gridlines
    panel.grid.minor.y = element_blank(),  # Remove minor x gridlines
    panel.grid.major.x = element_line(size = 0.1, color = "black"),  # Thin black major y gridlines
    panel.grid.minor.x = element_line(size = 0.08, color = "black"),  # Thin black minor y gridlines
    panel.background = element_blank(),  # Remove panel background
    panel.border = element_blank(),  # Remove panel border
    legend.position = "none"  # Hide legend
  )

# Print the plot
print(p3)




ggsave(filename = paste0(folder_path, "/", "SpearmanRho_Averages_with_SEM_Triangle.png"), plot = p3, width = 7, height = 8, dpi = 300)
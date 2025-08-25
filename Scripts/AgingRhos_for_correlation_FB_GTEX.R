


#This code currently gets the log2 TMM expression of GDF15, then centers and scales it, and then compared the pc2 score (which also gets centerd and scaled) for each tissue

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
library(boot)
library(here)

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


csv_tissue_data_read_in <- here("Data", "AgingRhos_correlate_Fibroblast_vs_GTEX.csv")
boxplot_data <- read.csv(csv_tissue_data_read_in)

# Calculate Spearman correlation
spearman_result <- cor.test(boxplot_data$Fibroblasts, boxplot_data$GTEx, method = "spearman")
spearman_rho <- spearman_result$estimate
spearman_p <- spearman_result$p.value

# Fit a linear model
lm_model <- lm(GTEx ~ Fibroblasts, data = boxplot_data)
r_squared <- summary(lm_model)$r.squared


# Create the color mapping
color_mapping <- c(
  "GDF15" = "pink",
  "Factor1" = "red",
  "Factor2" = "#808080",
  "Factor3" = "#CCCCCC",
  "Factor4" = "#8C8C8C",
  "Factor5" = "#808080",
  "Factor6" = "#CCCCCC",
  "Factor7" = "#CCCCCC",
  "Factor8" = "#808080",
  "Factor9" = "#CCCCCC",
  "Factor10" = "#808080",
  "Factor11" = "#CCCCCC",
  "Factor12" = "#808080",
  "ATF4" = "orange",
  "DDIT3" = "lightblue",
  "ATF5" = "brown",
  "sen_mean" = "blue",
  "prolif_mean" = "purple"
)

# Create the scatterplot with 95% confidence interval and colored points
p <- ggplot(boxplot_data, aes(x = Fibroblasts, y = GTEx, color = factor_gene)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  scale_color_manual(values = color_mapping) +
  labs(title = "Fibroblasts vs GTEx Spearman Correlation",
       x = "Fibroblasts",
       y = "GTEx") +
  annotate("text", x = min(boxplot_data$Fibroblasts), y = max(boxplot_data$GTEx), 
           label = paste0("Spearman rho: ", round(spearman_rho, 4), "\n",
                          "p-value: ", formatC(spearman_p, format = "e", digits = 2), "\n",
                          "R-squared: ", round(r_squared, 4)),
           hjust = 0, size = 4.5, color = "black", fontface = "italic") +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_blank())

# Display the plot
print(p)

ggsave(filename = file.path(paste0(folder_path,"/FB_vs_GTEx_Aging_Correlation.png")), plot = p, width = 14, height = 10)


# Filter the dataframe to include only rows up to and including 'Factor1'
boxplot_data2 <- boxplot_data[1:which(boxplot_data$factor_gene == "Factor1"), ]

# Calculate Spearman correlation
spearman_result <- cor.test(boxplot_data2$Fibroblasts, boxplot_data2$GTEx, method = "spearman")
spearman_rho <- spearman_result$estimate
spearman_p <- spearman_result$p.value

# Fit a linear model
lm_model <- lm(GTEx ~ Fibroblasts, data = boxplot_data2)
r_squared <- summary(lm_model)$r.squared

# Create the color mapping
color_mapping <- c(
  "GDF15" = "pink",
  "Factor1" = "red",
  # "Factor2" = "#808080",
  # "Factor3" = "#CCCCCC",
  # "Factor4" = "#8C8C8C",
  # "Factor5" = "#808080",
  # "Factor6" = "#CCCCCC",
  # "Factor7" = "#CCCCCC",
  # "Factor8" = "#808080",
  # "Factor9" = "#CCCCCC",
  # "Factor10" = "#808080",
  # "Factor11" = "#CCCCCC",
  # "Factor12" = "#808080",
  "ATF4" = "orange",
  "DDIT3" = "lightblue",
  "ATF5" = "brown",
  "sen_mean" = "blue",
  "prolif_mean" = "purple"
)

# Create the scatterplot with 95% confidence interval and colored points
p2 <- ggplot(boxplot_data2, aes(x = Fibroblasts, y = GTEx, color = factor_gene)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  scale_color_manual(values = color_mapping) +
  labs(title = "Fibroblasts vs GTEx Spearman Correlation",
       x = "Fibroblasts",
       y = "GTEx") +
  annotate("text", x = min(boxplot_data2$Fibroblasts), y = max(boxplot_data2$GTEx), 
           label = paste0("Spearman rho: ", round(spearman_rho, 4), "\n",
                          "p-value: ", formatC(spearman_p, format = "e", digits = 2), "\n",
                          "R-squared: ", round(r_squared, 4)),
           hjust = 0, size = 4.5, color = "black", fontface = "italic") +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_blank())

# Display the plot
print(p2)

ggsave(filename = file.path(paste0(folder_path,"/FB_vs_GTEx_Aging_Correlation_No_Other_Factors.png")), plot = p2, width = 14, height = 10)




# 
# Control F and then replace "Factor1" with Factor_
rm(list = ls())

# BiocManager::install("edgeR")
?effsize::cohen.d
library(tidyverse)
library(edgeR)
library(corrr)
library(tibble)
library(ggpubr)
library(readr)
library(circlize)
library(beepr)
library(ggplot2)
library(dplyr)
library(effsize)


condition_column_names <- c("DTHPLCE")
control <- "Hospital inpatient"
condition <- "Emergency room"
condition_for_titles <- "deathplace"
method_multiple_comparisons <- "BH"
date_of_csv_folder <- "2024_12_13"
CohenFalse_HedgesTrue <- "TRUE"


#------------------------------------------------------
#------------------------------------------------------ NAMING FOLDER AFTER THE NAME OF THE SCRIPT
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
  file_path_sans_ext(basename(script_path))
} else {
  format(Sys.Date(), "%Y_%m_%d")
}

#------------------------------------------------------

Plot_Save = "ON"
#------------------------------------------------------


Deathplace_option <- "DTHPLCE_csv_generator_for_all_tissues" # does not include Other 
# Deathplace_option <- "DTHPLCE_all_dthplce_csv_generator_for_all_tissues" # does include Other


folder_path_1 <- here("Results", "gtex")
if (!dir.exists(folder_path_1)) dir.create(folder_path_1, recursive = TRUE)

folder_path_2 <- here("Results", "gtex", Folder_name)
if (!dir.exists(folder_path_2)) dir.create(folder_path_2, recursive = TRUE)

condition_path_title <- here("Results", "gtex", Folder_name, condition_for_titles)
if (!dir.exists(condition_path_title)) dir.create(condition_path_title, recursive = TRUE)

folder_path <- condition_path_title


csv_readin <- here("Results", "gtex", "CSV_generator", Deathplace_option, "All_Tissue_data_DTHPLCE.csv")
data <- read.csv(csv_readin)

#----------------------------------------------------------------------------------------------------------

# Assuming your dataframe is named 'data'
data <- data %>%
  group_by(Tissue, group) %>%
  mutate(N = n())

data <- data %>%
  group_by(Tissue) %>%
  filter(all(N >= 20))


# Initialize result storage
results_list <- list()
plots_list <- list()
cohen_d_results <- data.frame()


# Define the correct group names for comparisons
comparisons <- list(
  c("Hospital inpatient_Factor1", "Emergency room_Factor1"),
  c("Hospital inpatient_GDF15", "Emergency room_GDF15")
)


# Iterate through each tissue
for (tissue in unique(data$Tissue)) {
  tissue_data <- data %>% filter(Tissue == tissue)
  
  # Initialize a vector to store unadjusted p-values
  tissue_p_values <- c()
  comparison_names <- c()
  
  for (comp in comparisons) {
    group1 <- comp[1]
    group2 <- comp[2]
    
    # Subset data for the comparison
    subset_data <- tissue_data %>% filter(group %in% c(group1, group2))
    
    if (nrow(subset_data) > 0) {
      # Set group levels explicitly to ensure order
      subset_data$group <- factor(subset_data$group, levels = c(group1, group2))
      
      # Perform Wilcoxon test
      wilcox_test <- wilcox.test(value ~ group, data = subset_data, exact = FALSE)
      
      # Store the unadjusted p-value and corresponding comparison
      tissue_p_values <- c(tissue_p_values, wilcox_test$p.value)
      comparison_names <- c(comparison_names, paste(group1, "vs", group2))
      
      
      
    } else {
      print(paste("No data for comparison:", group1, "vs", group2, "in tissue", tissue))
    }
  }
  
  # Apply p-value adjustment for multiple comparisons within this tissue
  if (length(tissue_p_values) > 0) {
    adjusted_p_values <- p.adjust(tissue_p_values, method = method_multiple_comparisons)
    
    # Store results for each comparison with adjusted p-values
    for (i in seq_along(comparison_names)) {
      comp <- strsplit(comparison_names[i], " vs ")[[1]]
      group1 <- comp[1]
      group2 <- comp[2]
      
      subset_data <- tissue_data %>% filter(group %in% c(group1, group2))
      
      if (nrow(subset_data) > 0) {
        # Set group levels explicitly to ensure order
        subset_data$group <- factor(subset_data$group, levels = c(group1, group2))
        
        # Calculate Cohen's D
        cohen_d <- cohen.d(subset_data$value, subset_data$group, hedges.correction = CohenFalse_HedgesTrue)
        
        # Save results
        p_values <- data.frame(
          Tissue = tissue,
          Comparison = comparison_names[i],
          Unadjusted_P = tissue_p_values[i],
          Adjusted_P = adjusted_p_values[i]
        )
        results_list[[paste(tissue, comparison_names[i], sep = "_")]] <- p_values
        
        cohen_d_results <- rbind(cohen_d_results, data.frame(
          Tissue = tissue,
          Comparison = comparison_names[i],
          Cohens_D = cohen_d$estimate,
          Adjusted_P = adjusted_p_values[i]
        ))
        
        # Create plot
        plot <- ggplot(subset_data, aes(x = group, y = value)) +
          geom_boxplot() +
          geom_jitter(width = 0.2, alpha = 0.5) +
          ggtitle(paste("Tissue:", tissue, "-", group1, "vs", group2)) +
          annotate("text", x = 1.5, y = max(subset_data$value, na.rm = TRUE), 
                   label = paste("Adj P:", signif(adjusted_p_values[i], 4), "\nMethod: ", method_multiple_comparisons), 
                   hjust = 0.5) +
          theme_minimal()
        
        plots_list[[paste(tissue, group1, group2, sep = "_")]] <- plot
      }
    }
  }
}


# Combine all p-values and adjusted p-values into a single dataframe
final_p_values <- bind_rows(results_list)

if (length(final_p_values) > 0) {
  final_p_values$Adjusted_P <- p.adjust(final_p_values$Unadjusted_P, method = method_multiple_comparisons)
}


cohen_d_results <- cohen_d_results %>%
  mutate(Group = paste(Tissue, Comparison, sep = "_"))

cohen_d_results <- cohen_d_results %>%
  select(-Adjusted_P, -Comparison)

final_p_values <- final_p_values %>%
  mutate(Group = paste(Tissue, Comparison, sep = "_"))

final_p_values <- final_p_values %>%
  select(-Comparison, -Tissue)

cohen_d_results <- cohen_d_results %>%
  left_join(final_p_values %>% select(Group, Adjusted_P), by = "Group")

cohen_d_results <- cohen_d_results %>%
  mutate(Comparison = str_extract(Group, "Hospital.*"))


# Apply p-value adjustment for multiple comparisons within this tissue
if (length(tissue_p_values) > 0) {
  adjusted_p_values <- p.adjust(tissue_p_values, method = method_multiple_comparisons)
  
  # Store results for each comparison with adjusted p-values
  for (i in seq_along(comparison_names)) {
    comp <- strsplit(comparison_names[i], " vs ")[[1]]
    group1 <- comp[1]
    group2 <- comp[2]
    
    subset_data <- tissue_data %>% filter(group %in% c(group1, group2))
    
    if (nrow(subset_data) > 0) {

     
      # Create plot
      plot <- ggplot(subset_data, aes(x = group, y = value)) +
        geom_boxplot() +
        geom_jitter(width = 0.2, alpha = 0.5) +
        ggtitle(paste("Tissue:", tissue, "-", group1, "vs", group2)) +
        annotate("text", x = 1.5, y = max(subset_data$value, na.rm = TRUE), 
                 label = paste("Adj P:", signif(adjusted_p_values[i], 4), "\nMethod: ", method_multiple_comparisons), 
                 hjust = 0.5) +
        theme_minimal()
      
      plots_list[[paste(tissue, group1, group2, sep = "_")]] <- plot
    }
  }
}







# Save each plot in the specified folder
for (name in names(plots_list)) {
  ggsave(
    filename = file.path(folder_path, paste0(name, ".png")),
    plot = plots_list[[name]],
    width = 6,
    height = 4
  )
}

# Create Cohen's D plots with significance coloring
cohen_d_results <- cohen_d_results %>%
  mutate(Significant = ifelse(Adjusted_P < 0.05, "Significant", "Not Significant"))

# Plot for Factor1 comparisons
cohen_d_factor1 <- cohen_d_results %>%
  filter(grepl("Factor1", Group)) %>%
  arrange(Cohens_D)

plot_factor1 <- ggplot(cohen_d_factor1, aes(x = Cohens_D, y = reorder(Tissue, Cohens_D), color = Significant)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  ggtitle("Hedge's G for Factor1 Comparisons") +
  xlab("Hedge's G") + ylab("Tissue") +
  theme_minimal()

# Plot for GDF15 comparisons
cohen_d_gdf15 <- cohen_d_results %>%
  filter(grepl("GDF15", Group)) %>%
  arrange(Cohens_D)

plot_gdf15 <- ggplot(cohen_d_gdf15, aes(x = Cohens_D, y = reorder(Tissue, Cohens_D), color = Significant)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  ggtitle("Hedge's G for GDF15 Comparisons") +
  xlab("Hedge's G") + ylab("Tissue") +
  theme_minimal()

# Display Cohen's D plots
print(plot_factor1)
print(plot_gdf15)

# Save Cohen's D plots
ggsave(filename = file.path(paste0(folder_path, "/1_Hedges_G_Factor1_", method_multiple_comparisons, ".png")), plot = plot_factor1, width = 6, height = 8)
ggsave(filename = file.path(paste0(folder_path, "/1_Hedges_G_GDF15_", method_multiple_comparisons, ".png")), plot = plot_gdf15, width = 6, height = 8)


# Display or save final results
write.csv(final_p_values, file = file.path(folder_path, "1_p_values_results.csv"), row.names = FALSE)
write.csv(cohen_d_results, file = file.path(folder_path, "1_Hedges_G_results.csv"), row.names = FALSE)

# Combine Factor1 and GDF15 data for combined plot
combined_data <- cohen_d_results %>%
  mutate(
    Type = ifelse(grepl("Factor1", Group), "Factor1", "GDF15"),
    Shape = ifelse(Type == "Factor1", "Circle", "Square"),
    Color = ifelse(Significant == "Significant", "Red", "Gray")
  ) %>%
  arrange(Cohens_D) %>%
  mutate(Tissue = factor(Tissue, levels = unique(cohen_d_factor1$Tissue)))

# Create combined plot
combined_plot <- ggplot(combined_data, aes(x = Cohens_D, y = Tissue)) +
  geom_point(aes(
    shape = Shape,
    color = Color,   # Outline color
    fill = Color     # Fill color
  ),
  size = 6,          # Size of the points
  stroke = 1,        # Bold outline thickness
  alpha = 0.6        # Transparency for the fill
  ) +
  scale_shape_manual(values = c("Circle" = 21, "Square" = 22)) +  # Circle = filled circle, Square = filled square
  scale_color_manual(values = c("Red" = "red", "Gray" = "gray")) +  # Outline color
  scale_fill_manual(values = c("Red" = "red", "Gray" = "gray")) +   # Fill color
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  ggtitle("Hedge's G for Factor1 and GDF15 Comparisons") +
  xlab("Hedge's G") +
  ylab("Tissue") +
  theme_minimal() +
  theme(
    axis.line = element_line(size = 0.25),  # Thinner axis lines
    panel.grid.major = element_line(size = 0.25),  # Thinner major grid lines
    panel.grid.minor = element_line(size = 0.1)  # Thinner minor grid lines
  )

# Display the plot
print(combined_plot)

# Save the combined plot
ggsave(filename = file.path(paste0(folder_path, "/1_Hedges_G_Combined_", method_multiple_comparisons, ".png")), plot = combined_plot, width = 14, height = 10)



# CHI SQUARE
# Step 1: Create the contingency table
contingency_table <- cohen_d_results %>%
  group_by(Comparison, Significant) %>%
  summarise(Count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Significant, values_from = Count, values_fill = 0)

# Rename columns for clarity
contingency_table <- contingency_table %>%
  rename(
    Not_Significant = `Not Significant`,
    Significant = Significant
  )

# Step 2: Perform the chi-square test
chisq_result <- chisq.test(contingency_table[, c("Not_Significant", "Significant")])

# Print the results
print("Contingency Table:")
print(contingency_table)

print("Chi-Square Test Result:")
print(chisq_result)


data_transformed <- data %>%
  select(SAMPID, AGE, SEX, Tissue, fa_vs_gdf15, value) %>%
  pivot_wider(names_from = fa_vs_gdf15, values_from = value)

Averages <- data_transformed %>%
  group_by(Tissue) %>%
  summarise(
    Avg_GDF15 = mean(GDF15, na.rm = TRUE),
    SD_GDF15 = sd(GDF15, na.rm = TRUE),
    SEM_GDF15 = SD_GDF15 / sqrt(n()),
    N_GDF15 = sum(!is.na(GDF15)),
    
    Avg_Factor1 = mean(Factor1, na.rm = TRUE),
    SD_Factor1 = sd(Factor1, na.rm = TRUE),
    SEM_Factor1 = SD_Factor1 / sqrt(n()),
    N_Factor1 = sum(!is.na(Factor1)),
    
    Avg_AGE = mean(AGE, na.rm = TRUE),
    SD_AGE = sd(AGE, na.rm = TRUE),
    SEM_AGE = SD_AGE / sqrt(n()),
    N_AGE = sum(!is.na(AGE)),
    
    Avg_SEX = mean(SEX, na.rm = TRUE),
    SD_SEX = sd(SEX, na.rm = TRUE),
    SEM_SEX = SD_SEX / sqrt(n()),
    N_SEX = sum(!is.na(SEX))
  )

write.csv(Averages, file = file.path(folder_path, "Averages_of_ISR_and_GDF15.csv"), row.names = FALSE)





# Count number of SEX values for all SAMPIDs
sex_counts_overall <- data %>%
  group_by(SEX) %>%
  summarise(count = n(), .groups = "drop")

# Count number of SEX values within each DTHPLCE category
sex_counts_by_dthplce <- data %>%
  group_by(DTHPLCE, SEX) %>%
  summarise(count = n(), .groups = "drop")

# View results
sex_counts_overall
sex_counts_by_dthplce


df_unique <- data %>% distinct(SAMPID, SEX, DTHPLCE)

# Overall counts
sex_counts_overall <- df_unique %>%
  count(SEX)

# Counts within each DTHPLCE
sex_counts_by_dthplce <- df_unique %>%
  count(DTHPLCE, SEX)

sex_counts_overall
sex_counts_by_dthplce

df_unique$Individual <- sub("^((?:[^.]+\\.){2}).*", "\\1", df_unique$SAMPID)

df_unique_Individual <- df_unique %>% distinct(Individual, SEX, DTHPLCE)

# Overall counts
sex_counts_overall <- df_unique_Individual %>%
  count(SEX)

# Counts within each DTHPLCE
sex_counts_by_dthplce <- df_unique_Individual %>%
  count(DTHPLCE, SEX)

sex_counts_overall
sex_counts_by_dthplce


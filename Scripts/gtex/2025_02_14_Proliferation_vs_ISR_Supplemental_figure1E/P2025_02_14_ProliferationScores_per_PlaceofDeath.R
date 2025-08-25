# This is a script to take the proliferation scores that were generated via Jacks method and compare with ISR scores vs Tissues in GTEx


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
library(effsize)

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


csv_prolif_of_tissues <- file.path(paste0(path, "/1_Scripts/gtex/proliferation_index_Jack_Devine_namesChanged.csv"))
proliferation_scores <- read.csv(csv_prolif_of_tissues)

# DO WE WANT LOG TRANSFORMED PROLIFERATION SCORES?!?!?!  Yes

#If yes:
proliferation_scores <- proliferation_scores %>%
  mutate(proliferation_score = log2(proliferation_score))

#If NO: COMMENT OUT THE 2 LINES ABOVE!

csv_ISR_of_tissues <- file.path(paste0(path, "/1_Scripts/gtex/deathplace/2025_01_14_BH_CohenFalse_HedgesTrue_False/Averages_of_ISR_and_GDF15.csv"))
ISR_scores <- read.csv(csv_ISR_of_tissues)

csv_sig <- file.path(paste0(path, "/1_Scripts/gtex/deathplace/2025_01_14_BH_CohenFalse_HedgesTrue_False/1_cohens_d_results.csv"))
sig_columns <- read.csv(csv_sig)

# Need to match the different Tissue names - Jacks are different - did by hand. 

# Perform a left join to add the proliferation_score column to ISR_score
merged_df <- merge(ISR_scores, proliferation_scores[, c("Tissue", "proliferation_score")], 
                   by = "Tissue", 
                   all.x = TRUE)

# Split into Factor1_sig and GDF15_sig
Factor1_sig <- sig_columns %>%
  filter(Comparison == "Hospital inpatient_Factor1 vs Emergency room_Factor1") %>%
  distinct(Tissue, .keep_all = TRUE)

GDF15_sig <- sig_columns %>%
  filter(Comparison == "Hospital inpatient_GDF15 vs Emergency room_GDF15") %>%
  distinct(Tissue, .keep_all = TRUE)

# Add the Significant column from Factor1_sig to merged_df
merged_df <- merged_df %>%
  left_join(Factor1_sig %>% select(Tissue, Significant), by = "Tissue")

merged_df <- merged_df %>%
  left_join(Factor1_sig %>% select(Tissue, Cohens_D), by = "Tissue") %>%
  mutate(deathplace_comparison = ifelse(Cohens_D > 0, "HI>ER", "ER>HI"))

# Drop the Cohens_D column if not needed
merged_df <- merged_df %>% select(-Cohens_D)


# Define tissues to remove
tissues_to_remove <- c("pancreas", "kidney_cortex", "uterus", "vagina", "spleen", "small_intestines_terminal_ileum")

# Filter out the unwanted tissues
filtered_df <- merged_df %>% 
  filter(!Tissue %in% tissues_to_remove)

merged_df <- filtered_df

spearman_result <- cor.test(merged_df$proliferation_score, merged_df$Avg_Factor1, method = "spearman")
spearman_rho <- spearman_result$estimate
spearman_p <- spearman_result$p.value

# This is the spearman rho of the ISR score of each tissue vs the proliferation score including significant AND NON significant tissues (in terms of their HI>ER or ER>HI)
p1 <- ggplot(merged_df, aes(x = proliferation_score, y = Avg_Factor1)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "prolif_score vs AvgISR Spearman All Tissues (Sig plus Non Sig)",
       x = "proliferation_score(log2)",
       y = "Avg_Factor1") +
  annotate("text", x = min(merged_df$proliferation_score), y = max(merged_df$Avg_Factor1), 
           label = paste0("Spearman rho: ", round(spearman_rho, 4), "\n",
                          "p-value: ", formatC(spearman_p, format = "e", digits = 2)),
           hjust = 0, size = 4.5, color = "black", fontface = "italic") +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_blank())

# Display the plot
print(p1)
ggsave(filename = file.path(paste0(folder_path_1stfolder,"/spearman_rho_ISRscore_vs_proliferation_both_non_and_sig.png")), plot = p1, width = 14, height = 10)


# Here we are Filtering by Significance= Significant: 

data <- merged_df %>%
  select(Tissue, Avg_Factor1, proliferation_score, Significant, deathplace_comparison)

# Scale the proliferation_score column
# data <- data %>%
#   mutate(proliferation_score = scale(proliferation_score))

data_sig <- data %>%
  filter(Significant != "Not Significant")

# Calculate Spearman correlation
spearman_result <- cor.test(data_sig$proliferation_score, data_sig$Avg_Factor1, method = "spearman")
spearman_rho <- spearman_result$estimate
spearman_p <- spearman_result$p.value


# Create the scatterplot with 95% confidence interval and colored points
p2 <- ggplot(data_sig, aes(x = proliferation_score, y = Avg_Factor1)) +
  geom_point(alpha = 0.7, size = 6) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "proliferation_score vs Avg_Factor1 Spearman Correlation",
       x = "proliferation_score(log2)",
       y = "Avg_Factor1") +
  annotate("text", x = min(data_sig$proliferation_score), y = max(data_sig$Avg_Factor1), 
           label = paste0("Spearman rho: ", round(spearman_rho, 4), "\n",
                          "p-value: ", formatC(spearman_p, format = "e", digits = 2)),
           hjust = 0, size = 4.5, color = "black", fontface = "italic") +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_blank())

# Display the plot
print(p2)

ggsave(filename = file.path(paste0(folder_path_1stfolder,"/significant_Tissues_spearman_rho_ISRscore_vs_proliferation.png")), plot = p2, width = 14, height = 10)


#####################################################################################################################################################################
# Now we want to split the data in two separate ways :
# 1. Have a cutoff of 3 for proliferation score to generate two groupings of tissues: Proliferative vs Non-Proliferative And compare the ISR scores between these two groups
    # Need to figure out how to split these into such groups, options include:
          #                                                                   - splitting between suspected tissues switching from prolif to non prolif, such as between breast mammary tissue and thyroid (non prolif)
          #                                                                   - splitting at the median

# 2. Split the tissues into the groups from plot 2G, where some tissues have significant HI>ER and some have ER>HI. Then, compare the proliferation scores of the two groups. Get Hedges G. 



# Calculate the median of the proliferation_score column
median_proliferation_score <- median(data_sig$proliferation_score, na.rm = TRUE)

# Print the result
print(median_proliferation_score)


#####
# 1 #
#####
data_sig <- data_sig %>%
  mutate(Proliferative_Group = ifelse(proliferation_score > median_proliferation_score, "Proliferative", "Non-Proliferative"))

# Check normality for Avg_Factor1 in both groups
normality_results <- data_sig %>%
  group_by(Proliferative_Group) %>%
  summarise(shapiro_p_value = shapiro.test(Avg_Factor1)$p.value)

# View the normality results
print(normality_results)


  # Perform Wilcoxon rank-sum test
test_result <- wilcox.test(Avg_Factor1 ~ Proliferative_Group, data = data_sig)

# View test result
print(test_result)

# Filter data for Proliferative and Non-Proliferative groups
subset_data <- data_sig %>% 
  filter(Proliferative_Group %in% c("Proliferative", "Non-Proliferative"))

# # Calculate Cohen's d with Hedges' g correction
hedges_g <- cohen.d(subset_data$Avg_Factor1,
                    as.factor(subset_data$Proliferative_Group),
                    hedges.correction = TRUE)
g_value <- hedges_g$estimate

# #Calculate Cohen's d with Hedges' g correction It is already true, and I get an error if I state = True like in the uncommented out above version
# hedges_g <- cohen.d(subset_data$Avg_Factor1,
#                     as.factor(subset_data$Proliferative_Group))
# 
# str(hedges_g)

# # Extract the estimate for Hedges' g
# g_value <- hedges_g$hedges.g[1, "effect"]

# Print the result
print(g_value)

# View the result
print(hedges_g)


# Define custom colors
custom_colors <- c("Non-Proliferative" = "blue", "Proliferative" = "purple")


# Create the plot
plot <- ggplot(data_sig, aes(x = Proliferative_Group, y = Avg_Factor1, color = Proliferative_Group)) +
  geom_jitter(width = 0.2, size = 6, alpha = 0.7) +  # Larger jittered points
  stat_summary(fun = median, geom = "crossbar", width = 0.4, 
               aes(ymin = ..y.., ymax = ..y.., color = Proliferative_Group), 
               size = 0.5, alpha = 0.5) +  # More transparent median lines
  labs(title = "Comparison of Avg_Factor1 by Proliferative Group",
       x = "Proliferative Group",
       y = "Avg_Factor1") +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove legend
  scale_color_manual(values = custom_colors)  # Apply custom colors


# Display the plot
print(plot)

ggsave(filename = file.path(paste0(folder_path_1stfolder,"/ISRscore_per_prolif_group_", median_proliferation_score, "_Hedges_G_", g_value,".png")), plot = plot, width = 14, height = 10)






#####
# 2 #
#####

# Check normality for Avg_Factor1 in both groups
normality_results2 <- data_sig %>%
  group_by(deathplace_comparison) %>%
  summarise(shapiro_p_value = shapiro.test(proliferation_score)$p.value)

print(normality_results2)

  # Perform Wilcoxon rank-sum test
test_result2 <- wilcox.test(proliferation_score ~ deathplace_comparison, data = data_sig)


# View test result
print(test_result2)


# Filter data for Proliferative and Non-Proliferative groups
subset_data2 <- data_sig %>% 
  filter(deathplace_comparison %in% c("HI>ER", "ER>HI"))

# # Calculate Cohen's d with Hedges' g correction
hedges_g2 <- cohen.d(subset_data$proliferation_score,
                    as.factor(subset_data2$deathplace_comparison),
                    hedges.correction = TRUE)

g2_value <- hedges_g2$estimate


# #Calculate Cohen's d with Hedges' g correction It is already true, and I get an error if I state = True like in the uncommented out above version
# hedges_g2 <- cohen.d(subset_data2$proliferation_score,
#                     as.factor(subset_data2$deathplace_comparison))
# 
# str(hedges_g2)

# # Extract the estimate for Hedges' g
# g2_value <- hedges_g2$hedges.g[1, "effect"]

# Print the result
print(g2_value)

# View the result
print(hedges_g)





# View the result
print(hedges_g2)
# g2_value <- hedges_g2$estimate

# Define custom colors
custom_colors <- c("HI>ER" = "orange", "ER>HI" = "maroon")

# Create the plot
plot <- ggplot(data_sig, aes(x = deathplace_comparison, y = proliferation_score, color = deathplace_comparison)) +
  geom_jitter(width = 0.2, size = 6, alpha = 0.7) +  # Larger jittered points
  stat_summary(fun = median, geom = "crossbar", width = 0.4, 
               aes(ymin = ..y.., ymax = ..y.., color = deathplace_comparison), 
               size = 0.5, alpha = 0.5) +  # More transparent median lines
  labs(title = "Comparison of proliferation_score by deathplace_comparison Group",
       x = "deathplace_comparison Group",
       y = "proliferation_score") +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove legend
  scale_color_manual(values = custom_colors)  # Apply custom colors


# Display the plot
print(plot)

ggsave(filename = file.path(paste0(folder_path_1stfolder,"/prolifScore_per_dthplce_comparison_Hedges_G_", g2_value,".png")), plot = plot, width = 14, height = 10)



#################################################################################
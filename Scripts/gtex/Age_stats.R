


#This code currently gets the log2 TMM expression of GDF15, then centers and scales it, and then compared the pc2 score (which also gets centerd and scaled) for each tissue

# Control F and then replace "Factor1" with Factor_
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

folder_path1 <- here("Results", "gtex", "Age_stats")
if (!dir.exists(folder_path1)) dir.create(folder_path1, recursive = TRUE)

folder_path <- here("Results", "gtex", "Age_stats", Folder_name)
if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

# Read necessary input files

script_path_cod <- here("Scripts", "gtex", "Attibutes_Phenos_Merged_plus_COD.R")
source(script_path_cod)
#------------------------------------------------------


condition_column_names <- c("DTHPLCE")

folder_path_1stfolder <- here("Results", "gtex", "Age_stats", Folder_name, condition_column_names, "_COD")
if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)



DTHPLCE <- annotations_merged %>%
  select(SAMPID, DTHPLCE, AGE) %>%
  mutate(SAMPID = sub("^([^.]*\\.[^.]*).*", "\\1", SAMPID)) %>%  # Extract up to second decimal
  distinct()

# Count rows for "Hospital inpatient" and "Emergency room"
place_counts <- annotations_merged %>%
  filter(DTHPLCE %in% c("Hospital inpatient", "Emergency room")) %>%
  count(DTHPLCE)

# Print the counts
print(place_counts)


# Count rows for "Hospital inpatient" and "Emergency room"
place_counts <- DTHPLCE %>%
  filter(DTHPLCE %in% c("Hospital inpatient", "Emergency room")) %>%
  count(DTHPLCE)

# Print the counts
print(place_counts)

# Filter and calculate statistics for 'Hospital inpatient' and 'Emergency room'
age_stats <- DTHPLCE %>%
  filter(DTHPLCE %in% c("Hospital inpatient", "Emergency room")) %>%
  group_by(DTHPLCE) %>%
  summarise(
    N = n(),  # Sample size
    Mean_AGE = mean(AGE, na.rm = TRUE),
    SD_AGE = sd(AGE, na.rm = TRUE),
    SEM_AGE = SD_AGE / sqrt(n())
  )

write.csv(age_stats, paste0(folder_path, "/_Age_States_DTHPLCE.csv"))


# Filter and calculate statistics for 'Hospital inpatient' and 'Emergency room'
age_stats2 <- DTHPLCE %>%
  group_by(DTHPLCE) %>%
  summarise(
    N = n(),  # Sample size
    Mean_AGE = mean(AGE, na.rm = TRUE),
    SD_AGE = sd(AGE, na.rm = TRUE),
    SEM_AGE = SD_AGE / sqrt(n())
  )

# Filter and calculate statistics for 'Hospital inpatient' and 'Emergency room'
age_stats3 <- DTHPLCE %>%
  # group_by(DTHPLCE) %>%
  summarise(
    N = n(),  # Sample size
    Mean_AGE = mean(AGE, na.rm = TRUE),
    SD_AGE = sd(AGE, na.rm = TRUE),
    SEM_AGE = SD_AGE / sqrt(n())
  )

print(age_stats3)
  
DTHPLCE2 <- annotations_merged %>%
  select(SAMPID, DTHPLCE, AGE, SEX) %>%
  mutate(SAMPID = sub("^([^.]*\\.[^.]*).*", "\\1", SAMPID)) %>%  # Extract up to second decimal
  distinct()
write.csv(DTHPLCE2, paste0(folder_path, "/Data_for_comparing_number_individuals.csv"))


# Count number of SEX values for all SAMPIDs
sex_counts_overall <- DTHPLCE2 %>%
  group_by(SEX) %>%
  summarise(count = n(), .groups = "drop")



rm(list = ls())

#This code reads in the different sample attributes and subject phenotypes, and adds in cause of death column to the merged attributes and phenotpes. 


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

# Set the base path
# Set the base path
path <- "/Users/janell/Library/CloudStorage/GoogleDrive-janells29@gmail.com/.shortcut-targets-by-id/1pUjuYVHpRn5G2HeHNW2X0DBAVzvL-vMa/MitoLab - General/ Members Folders/Janell/R_Code/ISRGDF15/2024_08_29_fb_gtex"
setwd(path)

Title <- "Hospital_Timing"

date <- Sys.Date()
date <- gsub("-", "_", as.character(date))

Plot_Save = "ON"
#------------------------------------------------------
Folder_name <- paste0(date, "_Time_In_Hospital")
#------------------------------------------------------

folder_path_1stfolder <- paste0(path, "/1_Scripts/gtex/", Title)

# Check if the folder exists
if (!dir.exists(folder_path_1stfolder)) {
  # If it doesn't exist, create the folder
  dir.create(folder_path_1stfolder, recursive = TRUE)
  cat("Folder created at:", folder_path_1stfolder, "\n")
} else {
  cat("Folder already exists at:", folder_path_1stfolder, "\n")
}

folder_path <- paste0(folder_path_1stfolder, "/", Folder_name)

# Check if the folder exists
if (!dir.exists(folder_path)) {
  # If it doesn't exist, create the folder
  dir.create(folder_path, recursive = TRUE)
  cat("Folder created at:", folder_path, "\n")
} else {
  cat("Folder already exists at:", folder_path, "\n")
}







# Read necessary input files
gtex_input <- paste0(path, "/2_InputFiles/gtex")
ISR_read_in <- file.path(paste0(gtex_input, "/Total_ISR_Gene_List_plus_gdf15.csv"))
Total_ISR_List <- read.csv(ISR_read_in) %>%
  select(-X)

setwd(gtex_input)
getwd()

annotations_subject_phenotypes <- read_tsv("GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.tsv")
annotations_sample_attributes <- read_tsv("GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.tsv") %>%
  filter(SMAFRZE == "RNASEQ") %>%
  mutate(SUBJID = substring(SAMPID,1,10)) %>%
  filter(SMRIN >= 6) %>%
  mutate(SUBJID_tmp = case_when(
    substr(SUBJID, nchar(SUBJID), nchar(SUBJID)) == "-" ~ substr(SUBJID, 1, nchar(SUBJID)-1),
    TRUE ~ SUBJID
  )) %>%
  select(-SUBJID) %>%
  rename(SUBJID = SUBJID_tmp)

annotations_subject_phenotypes$SUBJID <- gsub("\\-", ".", annotations_subject_phenotypes$SUBJID)
annotations_sample_attributes$SAMPID <- gsub("\\-", ".", annotations_sample_attributes$SAMPID)

annotations_sample_attributes <- annotations_sample_attributes %>%
  mutate(SUBJID = sub("^([^.]+\\.[^.]+)\\..*$", "\\1", SAMPID))

annotations_merged <- annotations_sample_attributes %>%
  inner_join(annotations_subject_phenotypes, by = "SUBJID")

# Initialize a dataframe to store Spearman's rho and p-values for all tissues
all_results <- data.frame(Tissue = character(), PC = character(), Spearman_Rho = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

setwd(path)
getwd()
# Path to the folder containing individual tissue folders
individual_tissues_folder <- paste0(path, "/2_InputFiles/gtex/individual_tissues_folders")


Timing_and_DTHPLCE <- annotations_subject_phenotypes %>%
  select(SUBJID, AGE, WGHT, HGHT, DTHPLCE, DTHFUCODD)

# Calculate mean DTHFUCODD and sample size for each DTHPLCE
summary_stats <- Timing_and_DTHPLCE %>%
  group_by(DTHPLCE) %>%
  summarise(
    avg_DTHFUCODD = mean(DTHFUCODD, na.rm = TRUE), # Average of DTHFUCODD
    n = sum(!is.na(DTHFUCODD)) # Count of non-NA values
  )


# Filter only Hospital inpatient and Emergency room groups
filtered_data <- Timing_and_DTHPLCE %>%
  filter(DTHPLCE %in% c("Hospital inpatient", "Emergency room"))

# Calculate summary statistics
summary_stats <- filtered_data %>%
  group_by(DTHPLCE) %>%
  summarise(
    avg_DTHFUCODD = mean(DTHFUCODD, na.rm = TRUE),  # Mean
    sd_DTHFUCODD = sd(DTHFUCODD, na.rm = TRUE),      # Standard Deviation
    sem_DTHFUCODD = sd_DTHFUCODD / sqrt(n()),        # Standard Error of the Mean (SEM)
    n = sum(!is.na(DTHFUCODD))                      # Sample size
  )

# Print summary stats
print(summary_stats)

# Create a plot with individual points and summary statistics
ggplot(filtered_data, aes(x = DTHPLCE, y = DTHFUCODD, color = DTHPLCE)) +
  geom_jitter(width = 0.2, alpha = 0.6) +  # Jittered points for visibility
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") + # Mean point
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") + # Error bars for SEM
  theme_minimal() +
  labs(
    title = "DTHFUCODD by Death Place",
    x = "Death Place",
    y = "DTHFUCODD"
  ) +
  theme(legend.position = "none")

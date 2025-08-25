


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
path <- "/Users/janell/Library/CloudStorage/GoogleDrive-janells29@gmail.com/.shortcut-targets-by-id/1pUjuYVHpRn5G2HeHNW2X0DBAVzvL-vMa/MitoLab - General/ Members Folders/Janell/R_Code/2024_08_29_fb_gtex"
setwd(path)

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


Death_causes <- annotations_subject_phenotypes %>%
  select(SUBJID, AGE, WGHT, HGHT, DTHMNNR, DTHCOD, DTHCODD, DTHLUCOD, DTHFUCOD)

# Create combined_string column
Death_causes <- Death_causes %>%
  mutate(combined_string = paste(DTHMNNR, DTHCOD, DTHLUCOD, DTHFUCOD, sep = "_"))

# Define the function to determine COD
determine_cod <- function(combined_string) {
  if (grepl("suicide|Suicide", combined_string)) {
    return("suicide")
  } else if (grepl("homicide|Homicide", combined_string)) {
    return("homicide")
  } else if (grepl("accident|Accident", combined_string) & grepl("drug|OD|toxi|overdose", combined_string)) {
    return("OD")
  } else if (grepl("accident|Accident", combined_string) & grepl("MVA|mva|motor|vehicle", combined_string) &  !grepl("non|NON", combined_string)) {
    return("accident")
  } else if (grepl("accident|Accident", combined_string) & grepl("fall|Fall|fell", combined_string)) {
    return("accident")
  } else if (grepl("Undetermined", combined_string) & grepl("overdose", combined_string)) {
    return("undetermined")
  } else if (grepl("Natural|natural", combined_string) & grepl("disease|Disease", combined_string)) {
    return("natural_disease")
  } else if (grepl("Natural|natural", combined_string) & grepl("unknown cause of death", combined_string)) {
    return("natural_unknown")
    # } else if (grepl("ex8", combined_string) & !grepl("ex9", combined_string)) {
    #   return("ex10")
  } else {
    return("natural") # or any default value you prefer
  }
}


# Create COD column using the function
Death_causes <- Death_causes %>%
  mutate(COD = sapply(combined_string, determine_cod))

COD_counts <- table(Death_causes$COD)
print(COD_counts)

Death_causes_selected <- Death_causes %>%
  select(SUBJID, COD)

# Join the selected columns from Death_causes to annotations_merged by SUBJID
annotations_merged <- annotations_merged %>%
  inner_join(Death_causes_selected, by = "SUBJID")



unique_death_causes_DTHFUCOD <- Death_causes %>%
  distinct(DTHCOD, DTHFUCOD)




annotations_merged <- annotations_merged %>%
  mutate(cardiac = case_when(
    grepl("cardiac|acute mi|myocardial|myocardiac|arrest; cardiac|bradycardia|heart|cardio|chf|aortic dissection|probable mi", 
          DTHCOD, ignore.case = TRUE) |
      grepl("cardiac|acute mi|myocardial|myocardiac|arrest; cardiac|bradycardia|heart|cardio|chf|aortic dissection|probable mi", 
            DTHFUCOD, ignore.case = TRUE) ~ "cardiac",
    TRUE ~ "non cardiac"
  ))









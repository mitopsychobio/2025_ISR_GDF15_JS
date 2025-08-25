


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
  select(SUBJID, AGE, WGHT, HGHT, DTHMNNR, DTHCOD, DTHLUCOD, DTHFUCOD)

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


Death_causes_selected <- Death_causes %>%
  select(SUBJID, COD)

# Join the selected columns from Death_causes to annotations_merged by SUBJID
annotations_merged <- annotations_merged %>%
  inner_join(Death_causes_selected, by = "SUBJID")




# Updated July 15th 2024

#####################################    for DTHCOD    #########################################################
unique_death_causes <- Death_causes %>%
  select(DTHCOD)

unique_death_causes <- unique_death_causes %>%
  distinct(DTHCOD)

unique_death_causes_new <- unique_death_causes %>%
  mutate(New_DTHCOD = case_when(
    grepl("Cardiac|cardiac|myocardial|Myocardial|arrest|cardiovascular|heart failure|acute mi|Massive Heart Attack|probable mi|heart attack", DTHCOD, ignore.case = TRUE) ~ "cardiac failure",
    grepl("stroke|Stroke|cerebrovascular|CVA|cva|Cva|Cerebral Vascular|cerebral vascular", DTHCOD, ignore.case = TRUE) ~ "cerebrovascular accident",
    grepl("COPD|copd|pulmonary edema", DTHCOD, ignore.case = TRUE) ~ "COPD",
    grepl("alcoholism", DTHCOD, ignore.case = TRUE) ~ "alcoholism",
    grepl("anoxia|Anoxia|ANOXIA|anoxic encephalopathy", DTHCOD, ignore.case = TRUE) ~ "anoxia",
    grepl("als|ALS", DTHCOD, ignore.case = TRUE) ~ "ALS",
    grepl("allergic reaction|allergy|allergic|allergies", DTHCOD, ignore.case = TRUE) ~ "allergic reaction",
    grepl("brain cancer", DTHCOD, ignore.case = TRUE) ~ "brain cancer",
    grepl("cancer", DTHCOD, ignore.case = TRUE) ~ "cancer",
    grepl("head trauma|Head trauma|Head Trauma|blunt injury|blunt force trauma|trauma|trauma due to fall|head trama", DTHCOD, ignore.case = TRUE) ~ "head trauma",
    grepl("heart disease|CAD|congestive failure heart", DTHCOD, ignore.case = TRUE) ~ "heart disease",
    grepl("Dementia", DTHCOD, ignore.case = TRUE) ~ "dementia",
    grepl("end stage liver disease|ESLD|ESRD|renal|liver disease", DTHCOD, ignore.case = TRUE) ~ "liver disease",
    grepl("kidney failure| kidney diseases", DTHCOD, ignore.case = TRUE) ~ "kidney disease",
    grepl("sirs|smoke inhalation-respiratory disease|respiratory diseases|lung disease", DTHCOD, ignore.case = TRUE) ~ "respiratory disease",
    grepl("unknown death|unknown", DTHCOD, ignore.case = TRUE) ~ "unknown",
    grepl("toxic effect of unspecified substance|poisoning|poisoning by overdose of substance|overdose", DTHCOD, ignore.case = TRUE) ~ "poison",
    grepl("suicide by hanging|suicide-hanging|strangulation|asphyiation due to hanging", DTHCOD, ignore.case = TRUE) ~ "suicide by hanging",
    grepl("seizure|epilepsy", DTHCOD, ignore.case = TRUE) ~ "seizures",
    TRUE ~ DTHCOD  # Default case to keep the original value if no match
  ))

unique_death_causes_new <- unique_death_causes_new %>%
  select(New_DTHCOD) %>%
  distinct(New_DTHCOD)
#####################################    for DTHCOD    #########################################################


#--------------------------------------
# What to Keep or remove!
#--------------------------------------
all_objects <- ls()

# Specify the objects you want to keep
objects_to_remove <- c("Death_causes_selected", "Death_causes", "annotations_sample_attributes", "annotations_subject_phenotypes", "determine_cod")

# Create a list of objects to remove
objects_to_keep <- setdiff(all_objects, objects_to_remove)

# Remove all objects except the ones specified
rm(list = objects_to_remove)
rm(all_objects)
rm(objects_to_keep, objects_to_remove)

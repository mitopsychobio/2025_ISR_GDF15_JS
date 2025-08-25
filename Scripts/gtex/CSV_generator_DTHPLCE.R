


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

# Set the base path


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

folder_path <- here("Results", "gtex", "CSV_generator", Folder_name)
if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)


# Read necessary input files

script_path <- here("Scripts", "gtex", "Attibutes_Phenos_Merged_plus_COD.R")
source(script_path)

condition_column_names <- c("DTHPLCE")
cond_csv_name <- paste0(condition_column_names, "_csv_generator_for_all_tissues")

# folder_path <- here("Results", "gtex", Folder_name, cond_csv_name)
# if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)



# Path to the folder containing individual tissue folders
#########################################################################################################
individual_tissues_folder <- here("Data", "gtex", "All_Tissues_Indiv_Folders")
#########################################################################################################

# Debug: print the path to check if it is correct
print(individual_tissues_folder)

# Initialize variables

alldata <- data.frame()


# Get a list of all subdirectories within the main folder
subfolders <- list.dirs(individual_tissues_folder, recursive = FALSE, full.names = TRUE)

# Loop through each subfolder
for (folder in subfolders) {
  # Print the folder name
  print(paste("Processing folder:", basename(folder)))
  
  last_component <- basename(folder)
  
  # Set the path to the current folder
  current_folder <- folder
  
  correlations_readin <- here("Data", "Fibroblast_lifespan", "Processed", "20250526_Spearman_rho_Gene_vs_DaysGrown_correlations_for_scaling.csv") 
  correlations <- read.csv(correlations_readin)
  
  
  
  # Print the last component
  print(last_component)
  
  
  # new_folder_path <- paste0(path, "/3_Outputs/gtex/", Folder_name, "/", last_component)
  # 
  # # Check if the folder exists
  # if (!dir.exists(new_folder_path)) {
  #   # If it doesn't exist, create the folder
  #   dir.create(new_folder_path, recursive = TRUE)
  #   cat("Folder created at:", new_folder_path, "\n")
  # } else {
  #   cat("Folder already exists at:", new_folder_path, "\n")
  # }
  # 
  
  # Get a list of all .gz files in the tissue folder
  gz_files <- list.files(folder, pattern = "\\.gz$", full.names = TRUE) # CHANGE BACK TO THIS AFTER DONE WITH BRAIN
  # gz_files <- list.files(individual_tissues_folder, pattern = "\\.gz$", full.names = TRUE)
  
  # Initialize dataframes for each tissue
  df_tissue_list <- list()
  gene_symbols <- NULL
  
  # Loop through each .gz file and read data
  for (gz_file in gz_files) {
    file_name <- basename(gz_file)
    tissue <- sub(".*_v8_(.*?)\\.gct\\.gz", "\\1", file_name)
    
    # Read data from the current .gz file
    raw_gene_counts <- read.delim(gz_file, skip = 2)
    
    # Extract gene symbols if not already done
    if (is.null(gene_symbols)) {
      gene_symbols <- raw_gene_counts[, 2:3]
    }
    
    # Extract sample IDs and gene counts
    gene_counts <- raw_gene_counts %>%
      select(-id, -Description) %>%
      as.data.frame()
    
    # Store the gene counts in the dataframe list with tissue name
    df_tissue_list[[tissue]] <- gene_counts
  }
  
  # Combine dataframes if there are multiple files
  combined_gene_counts <- gene_symbols
  for (tissue in names(df_tissue_list)) {
    combined_gene_counts <- combined_gene_counts %>%
      left_join(df_tissue_list[[tissue]], by = "Name")
  }
  
  # Remove the "Description" column if it exists
  if ("Description" %in% colnames(combined_gene_counts)) {
    combined_gene_counts <- combined_gene_counts %>%
      select(-Description)
  }
  
  combined_gene_counts <- combined_gene_counts %>%
    column_to_rownames("Name")
  
  annotations_sample_attributes_filtered <- annotations_merged %>%
    filter(SAMPID %in% colnames(combined_gene_counts))
  
  # TMM normalization
  dge <- DGEList(combined_gene_counts[, annotations_sample_attributes_filtered$SAMPID], group = factor(annotations_sample_attributes_filtered$COD))
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  exprs <- cpm(dge, log = TRUE) # changed from FALST to TRUE 9/27/2024 
  
  gene_symbols <- gene_symbols %>%
    mutate(Description = ifelse(Description == "KIAA0141", "DELE1", Description),  # changing so that DELE1 is the gene symbol
           Description = ifelse(Description == "WARS", "WARS1", Description),
           Description = ifelse(Description == "NARS", "NARS1", Description))
  
  
  # Join with gene symbols
  exprs <- exprs %>%
    as.data.frame() %>%
    rownames_to_column("Name") %>%
    full_join(gene_symbols, by = "Name") %>%
    na.omit() %>%
    select(-Name) %>%
    rename("Gene" = Description)
  
  # Filter genes of the ISR only into our exprs
  exprs_ISR <- exprs %>%
    filter(Gene %in% Total_ISR_List$Gene) %>%
    column_to_rownames("Gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("SAMPID") %>%
    unique()
  
  
  #------------------ finding missing genes ----------------------
  
  # Step 1: Extract gene names
  exprs_genes <- colnames(exprs_ISR)[-1]
  list_genes <- Total_ISR_List$Gene
  
  # Step 2: Standardize gene names
  exprs_genes_clean <- trimws(toupper(exprs_genes))
  list_genes_clean <- trimws(toupper(list_genes))
  
  # Remove duplicates if necessary
  exprs_genes_clean <- unique(exprs_genes_clean)
  list_genes_clean <- unique(list_genes_clean)
  
  # Step 3: Find missing genes
  genes_in_exprs_not_in_list <- setdiff(exprs_genes_clean, list_genes_clean)
  genes_in_list_not_in_exprs <- setdiff(list_genes_clean, exprs_genes_clean)
  
  # Step 4: Print results
  if (length(genes_in_exprs_not_in_list) > 0) {
    cat("Genes in exprs_ISR but not in Total_ISR_List:\n")
    print(genes_in_exprs_not_in_list)
  } else {
    cat("All genes in exprs_ISR are present in Total_ISR_List.\n")
  }
  
  if (length(genes_in_list_not_in_exprs) > 0) {
    cat("Genes in Total_ISR_List but not in exprs_ISR:\n")
    print(genes_in_list_not_in_exprs)
  } else {
    cat("All genes in Total_ISR_List are present in exprs_ISR.\n")
  }
  
  # Check if there are missing genes
  if (length(genes_in_list_not_in_exprs) > 0 || length(genes_in_exprs_not_in_list) > 0 ) {
    cat(paste("Missing genes detected in folder:", last_component, "- Skipping this folder.\n"))
    # Skip to the next iteration of the loop
    next
  }
  
  #-------------------------------------------------------------------------------------------------------  
  
  # CHANGING THIS mean centering and scaling needs to be performed using fibroblast data
  ###################################################################################################
  pathtoMEANofX <- here("Data", "Fibroblast_lifespan", "Processed", "20250526_Spearman_rho_Gene_vs_DaysGrown_meanx_for_scaling.csv")
  mean_x_fbdata <- read.csv(pathtoMEANofX)
  
  pathtoSDofX <- here("Data", "Fibroblast_lifespan", "Processed", "20250526_Spearman_rho_Gene_vs_DaysGrown_sdx_for_scaling.csv") 
  sd_x_fbdata <- read.csv(pathtoSDofX)
  ###################################################################################################
  
  # Rename columns using dplyr
  sd_x_fbdata <- sd_x_fbdata %>%
    rename(
      Gene = X,
      sd = x
    )
  
  # Rename columns using dplyr
  mean_x_fbdata <- mean_x_fbdata %>%
    rename(
      Gene = X,
      mean = x
    )
  
  sd_x_fbdata <- sd_x_fbdata %>%
    left_join(mean_x_fbdata %>% select(Gene, mean), by = "Gene")
  
  
  # Step 1: Extract gene names from exprs_ISR (excluding SampID)
  gene_columns <- setdiff(names(exprs_ISR), "SAMPID")
  
  # Step 2: Ensure that gene names match between exprs_ISR and sd_x_fbdata
  common_genes <- intersect(gene_columns, sd_x_fbdata$Gene)
  # view(common_genes)
  
  
  #Finding any missing genes
  # Step 1: Extract gene names
  exprs_genes <- setdiff(names(exprs_ISR), "SAMPID")
  sd_genes <- sd_x_fbdata$Gene
  
  # Step 2: Standardize gene names
  exprs_genes_clean <- trimws(toupper(exprs_genes))
  sd_genes_clean <- trimws(toupper(sd_genes))
  
  # Step 3: Find missing genes
  genes_in_exprs_not_in_sd <- setdiff(exprs_genes_clean, sd_genes_clean)
  genes_in_sd_not_in_exprs <- setdiff(sd_genes_clean, exprs_genes_clean)
  
  # Step 4: Print missing genes
  if (length(genes_in_exprs_not_in_sd) > 0) {
    cat("Genes in exprs_ISR but not in sd_x_fbdata:\n")
    print(genes_in_exprs_not_in_sd)
  } else {
    cat("All genes in exprs_ISR are present in sd_x_fbdata.\n")
  }
  
  if (length(genes_in_sd_not_in_exprs) > 0) {
    cat("Genes in sd_x_fbdata but not in exprs_ISR:\n")
    print(genes_in_sd_not_in_exprs)
  } else {
    cat("All genes in sd_x_fbdata are present in exprs_ISR.\n")
  }
  
  
  
  
  #Scaling the data
  # Step 3: Subset exprs_ISR to include only the common genes
  exprs_ISR_subset <- exprs_ISR[, c("SAMPID", common_genes)]
  
  # Step 4: Create named vectors for mean and sd
  mean_values <- sd_x_fbdata$mean
  names(mean_values) <- sd_x_fbdata$Gene
  
  sd_values <- sd_x_fbdata$sd
  names(sd_values) <- sd_x_fbdata$Gene
  
  # Subset mean and sd to include only common genes
  mean_values <- mean_values[common_genes]
  sd_values <- sd_values[common_genes]
  
  # Check for zero standard deviations
  if (any(sd_values == 0)) {
    stop("Standard deviation is zero for some genes. Cannot perform scaling.")
  }
  
  # Step 5: Manually scale the data
  exprs_ISR_scaled <- exprs_ISR_subset
  
  for (gene in common_genes) {
    exprs_ISR_scaled[[gene]] <- (exprs_ISR_subset[[gene]] - mean_values[gene]) / sd_values[gene]
  }
  
  
  exprs_gdf15 <- exprs_ISR_scaled %>%
    select(SAMPID, GDF15)
  
  # # Join with subject phenotypes
  exprs_gdf15 <- exprs_gdf15 %>%
    full_join(annotations_sample_attributes_filtered, by = "SAMPID") %>%
    select(GDF15, SAMPID, AGE, SEX) %>%
    unique() %>%
    na.omit()
  #
  # # Join with subject phenotypes
  exprs_gdf15_tissue <- exprs_gdf15 %>%
    full_join(annotations_sample_attributes_filtered, by = "SAMPID") %>%
    select(GDF15, SMTSD, SAMPID) %>%
    unique() %>%
    na.omit()
  
  ########################################################################################
  
  scaled_exprs_ISR <- exprs_ISR_scaled %>%
    full_join(annotations_sample_attributes_filtered, by = "SAMPID") %>%
    unique()
  # 
  # scaled_exprs_ISR <- scaled_exprs_ISR %>% # not sure how this is different
  #   select(SMTSD, SUBJID, everything()) %>%
  #   unique()
  
  tissue_results <- data.frame(PC = character(), Spearman_Rho = numeric(), P_Value = numeric(), regression_coefficient = numeric(), reg_p_value = numeric(), N = numeric(), stringsAsFactors = FALSE)
  
  
  loadings_readin <- here("Data", "Fibroblast_lifespan", "Processed", "20250526_Spearman_rho_Gene_vs_DaysGrown_FacLoads_12_factors.csv")
  loadings <- read.csv(loadings_readin)
  loadings <- loadings %>%
    select(X, Factor1)
  
  
  selected_loadings <- loadings %>%
    rename(Gene = X, Loading = Factor1) %>%
    as.data.frame()
  
  common_genes <- intersect(selected_loadings$Gene, colnames(scaled_exprs_ISR))
  selected_loadings <- selected_loadings %>%
    filter(Gene %in% common_genes)
  
  # Step 1: Prepare the correlations matrix
  rownames(loadings) <- loadings$X
  loadings$X <- NULL
  loadings_matrix <- as.matrix(loadings)
  
  scaled_data <- scaled_exprs_ISR %>%
    select(SAMPID, all_of(common_genes))
  
  rownames(scaled_data) <- scaled_data$SAMPID
  scaled_data$SAMPID <- NULL
  scaled_data_matrix <- as.matrix(scaled_data)
  
  
  # t(loadings) %*% solve(correlations)
  # Step 1: Prepare the correlations matrix
  rownames(correlations) <- correlations$X
  correlations$X <- NULL
  correlations_matrix <- as.matrix(correlations)
  
  W <- t(loadings_matrix) %*% solve(correlations_matrix)
  Wtransposed <- t(W)
  final <- scaled_data_matrix %*% Wtransposed
  
  final_df <- as.data.frame(final)
  
  final_df <- final_df %>%
    tibble::rownames_to_column(var = "SAMPID")
  
  final_data <- final_df %>%
    full_join(annotations_sample_attributes_filtered, by = "SAMPID") %>%
    unique()
  
  
  
  # Prepare boxplot data
  boxplot_data <- final_data %>%
    select(SAMPID, Factor1, AGE, SEX, DTHPLCE) %>%
    left_join(exprs_ISR_scaled %>% select(SAMPID, GDF15), by = "SAMPID") %>%
    filter(DTHPLCE %in% c("Hospital inpatient", "Emergency room")) %>%
    pivot_longer(cols = c(Factor1, GDF15), names_to = "fa_vs_gdf15", values_to = "value") %>%
    mutate(
      group = paste(DTHPLCE, fa_vs_gdf15, sep = "_"),
      type = ifelse(grepl("Factor1$", group), "fa_Scores", "GDF15_scaled"),
      value = as.numeric(value)
    ) %>%
    mutate(group = factor(group, levels = c(
      "Hospital inpatient_Factor1", "Emergency room_Factor1",
      "Hospital inpatient_GDF15", "Emergency room_GDF15"
    )))
  
  boxplot_data <- boxplot_data %>%
    mutate(Tissue = last_component)
  
  alldata <- rbind(alldata, boxplot_data)
  
}

write.csv(alldata, paste0(folder_path, "/All_Tissue_data_",condition_column_names, ".csv"))

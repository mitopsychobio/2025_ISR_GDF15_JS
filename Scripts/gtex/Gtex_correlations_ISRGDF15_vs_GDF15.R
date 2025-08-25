


#This code currently gets the log2 TMM expression of GDF15, then centers and scales it, and then compared the pc2 score (which also gets centerd and scaled) for each tissue

rm(list = ls())
# 
# BiocManager::install("ComplexHeatmap")
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



Plot_Save = "ON"
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
  sub("\\.[^.]*$", "", basename(script_path))
} else {
  format(Sys.Date(), "%Y_%m_%d")
}


#------------------------------------------------------

correlation_path <- here("Results", "gtex", "correlations")
if (!dir.exists(correlation_path)) dir.create(correlation_path, recursive = TRUE)

folder_path <- here("Results", "gtex", "correlations", Folder_name)
if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)



# Read necessary input files

script_path <- here("Scripts", "gtex", "Attibutes_Phenos_Merged_plus_COD.R")
source(script_path)



# Path to the folder containing individual tissue folders
#########################################################################################################
individual_tissues_folder <- here("Data", "gtex", "All_Tissues_Indiv_Folders")
#########################################################################################################

# Debug: print the path to check if it is correct
print(individual_tissues_folder)

correlation_results <- data.frame(
  Tissue = character(),
  rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Initialize an empty dataframe to store the results
results_df <- data.frame(
  Tissue = character(),
  rho = numeric(),
  p_adjust = numeric(),
  asterisks = character(),
  stringsAsFactors = FALSE
)



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



# new_folder_path <- here("Results", "gtex", Folder_name, last_component)

# # Check if the folder exists
# if (!dir.exists(new_folder_path)) {
#   # If it doesn't exist, create the folder
#   dir.create(new_folder_path, recursive = TRUE)
#   cat("Folder created at:", new_folder_path, "\n")
# } else {
#   cat("Folder already exists at:", new_folder_path, "\n")
# }


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
  
  #### TMM normalization
  ##########       Option 1
  
  dge <- DGEList(combined_gene_counts[, annotations_sample_attributes_filtered$SAMPID], group = factor(annotations_sample_attributes_filtered$COD))
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  exprs <- cpm(dge, log = TRUE)
  
  
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
  
  
  loadings_readin <- here("Data", "Fibroblast_lifespan", "Processed", "20250526_Spearman_rho_Gene_vs_DaysGrown_FacLoads_12_factors.csv") # includes gdf15, is from folder 2024_09_23
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
     

    
     boxplot_data <- final_data %>%
       select(SAMPID, Factor1, AGE, SEX, DTHPLCE)
     
     exprs_GDF15 <- exprs_ISR_scaled %>%
       select(SAMPID, GDF15)
     
     boxplot_data <- boxplot_data %>%
       left_join(exprs_GDF15, by = "SAMPID")
     
     ## 1. Get the tissue name
     Tissue <- last_component  # or use 'last_component' if that's the correct variable
     
     boxplot_data$Tissue <- Tissue
     
     # boxplot_data <- boxplot_data %>%
     #   filter(DTHPLCE %in% c("Hospital inpatient", "Emergency room"))
     
     # Reshape the dataframe
     boxplot_data <- pivot_longer(
       boxplot_data, 
       cols = c(Factor1, GDF15), 
       names_to = "fa_vs_gdf15", 
       values_to = "value"
     )
     
     # Reshape data from long to wide format
     boxplot_data_wide <- boxplot_data %>%
       pivot_wider(
         names_from = fa_vs_gdf15,  # Columns to create from unique values in fa_vs_gdf15
         values_from = value         # Values to fill in the new columns
       )
     
     boxplot_data <- boxplot_data_wide
     
     # Convert Factor1 and GDF15 to numeric
     boxplot_data <- boxplot_data %>%
       mutate(
         Factor1 = as.numeric(as.character(Factor1)),
         GDF15 = as.numeric(as.character(GDF15))
       )

     # boxplot_data$value <- as.numeric(boxplot_data$value)
     # class(boxplot_data$value)

     
     
     
     
  ##### where I pasted
     
     # 
     # # Reshape data from long to wide format
     # boxplot_data_wide <- boxplot_data %>%
     #   pivot_wider(
     #     names_from = fa_vs_gdf15,  # Columns to create from unique values in fa_vs_gdf15
     #     values_from = value         # Values to fill in the new columns
     #   )
     # 
     # boxplot_data <- boxplot_data_wide
     
     # Calculate Spearman correlations per Tissue
     correlation_results <- boxplot_data %>%
       group_by(Tissue) %>%
       summarise(
         rho = cor(Factor1, GDF15, method = "spearman", use = "complete.obs"),
         p_value = cor.test(Factor1, GDF15, method = "spearman", exact = FALSE)$p.value
       ) %>%
       ungroup()
     
     # Adjust p-values for multiple comparisons (Bonferroni correction with 44 tissues)
     correlation_results <- correlation_results %>%
       mutate(
         p_adjust = p_value * 44,          # Bonferroni correction
         p_adjust = ifelse(p_adjust > 1, 1, p_adjust),  # Cap at 1
         asterisks = case_when(
           p_adjust < 0.0001 ~ "****",
           p_adjust < 0.001  ~ "***",
           p_adjust < 0.01   ~ "**",
           p_adjust < 0.05   ~ "*",
           TRUE              ~ ""
         )
       )
     
     # Merge correlation results back to boxplot_data
     boxplot_data <- boxplot_data %>%
       left_join(correlation_results, by = "Tissue")
     
     # Define color mapping for DTHPLCE
     color_mapping <- c(
       "Hospital inpatient" = "orange",
       "Emergency room" = "maroon",
       "Other" = "gray"   # Assign specific color if needed
     )

     # Ensure 'DTHPLCE' is a factor with the correct levels
     boxplot_data$DTHPLCE <- factor(boxplot_data$DTHPLCE, levels = names(color_mapping))
     
     # Replace NA values in DTHPLCE with "Other"
     boxplot_data <- boxplot_data %>%
       mutate(
         DTHPLCE = replace_na(DTHPLCE, "Other")
       )
     sum(is.na(boxplot_data$DTHPLCE))  # Should be 0 or handled in color mapping
     
     sum(is.na(boxplot_data$Factor1))  # Count of NA in Factor1
     sum(is.na(boxplot_data$GDF15))    # Count of NA in GDF15
     boxplot_data$Factor1 <- as.numeric(as.character(boxplot_data$Factor1))
     boxplot_data$GDF15 <- as.numeric(as.character(boxplot_data$GDF15))
     
       
       # Retrieve correlation results
       rho_value <- unique(boxplot_data$rho)
       asterisks <- unique(boxplot_data$asterisks)
       
       # Handle potential multiple entries
       rho_value <- ifelse(length(rho_value) > 1, rho_value[1], rho_value)
       asterisks <- ifelse(length(asterisks) > 1, asterisks[1], asterisks)
       
       # Create the plot
       p <- ggplot(boxplot_data, aes(x = Factor1, y = GDF15, color = DTHPLCE)) +
         geom_point(size = 3, alpha = 0.7) +  # Scatter points
         geom_smooth(method = "lm", se = TRUE, color = "red") + 
         ggtitle(Tissue) +
         theme_minimal() +
         theme(
           plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
           axis.title = element_text(size = 12)
         ) +
         labs(
           x = "Index Score",
           y = "GDF15"
         ) +
         # # Add Spearman rho and asterisks
         # annotate("text", 
         #          x = Inf, y = Inf, 
         #          label = paste0("rho = ", round(rho_value, 3), " ", asterisks),
         #          hjust = 1.1, vjust = 1.5,
         #          size = 5, 
         #          color = "black")+
         # Apply custom color mapping
         scale_color_manual(values = color_mapping)
    
         
       # Display the plot
       print(p)
       
       save_path <- paste0(folder_path, "/")
       # Create the directory if it doesn't exist
       if (!dir.exists(save_path)) {
         dir.create(save_path, recursive = TRUE)
         cat("Directory created at:", save_path, "\n")
       } else {
         cat("Directory already exists at:", save_path, "\n")
       }
       
       # Save the plot
       ggsave(
         filename = paste0("Spearman_Correlation_", gsub(" ", "_", Tissue), ".png"),
         plot = p,
         path = save_path,
         width = 8,
         height = 6,
         dpi = 300
       )
       
       
       
       # Extract rho and p_adjust from boxplot_data
       rho_value <- unique(boxplot_data$rho)
       p_adjust_value <- unique(boxplot_data$p_adjust)
       
       # Assign asterisks based on adjusted p-value
       asterisks <- unique(boxplot_data$asterisks)
       
       # Create a temporary dataframe for the current Tissue
       temp_df <- data.frame(
         Tissue = Tissue,
         rho = rho_value,
         p_adjust = p_adjust_value,
         asterisks = asterisks,
         stringsAsFactors = FALSE
       )
       
       # Append the temporary dataframe to the results dataframe
       results_df <- rbind(results_df, temp_df)
       
       
       
       cat("Plot saved for tissue:", Tissue, "\n")
     }
     



# Optional: Order tissues by rho for better visualization
results_df$Tissue <- factor(results_df$Tissue, levels = results_df$Tissue[order(results_df$rho)])

# Create the correlation plot
correlation_plot <- ggplot(results_df, aes(x = rho, y = Tissue)) +
  geom_point(color = "black", 
             size = 4, 
             alpha = 0.7) +  # Semi-transparent points
  # geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Reference line at rho = 0
  # geom_text(aes(label = asterisks), 
  #           hjust = -0.3,  # Adjust horizontal position
  #           vjust = 0.5,   # Adjust vertical position
  #           size = 6, 
  #           color = "black") +  # Add asterisks next to significant points
  expand_limits(x = 0) +  # Ensure x-axis includes
  labs(
    title = "Spearman Correlation between Factor1 and GDF15 by Tissue",
    x = "Spearman Rho",
    y = "Tissue"
    # color = "Tissue"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  # Hide legend if colors are self-explanatory
  )

# Display the plot
print(correlation_plot)

# Define the directory where the summary plot will be saved
save_path2 <- paste0(folder_path, "/")
     
 # Save the summary plot
ggsave("Spearman_Correlation_Summary_x0.png", plot = correlation_plot, path = save_path2, width = 10, height = 8, dpi = 300)
     
     
     
     

# # Create the correlation plot
# correlation_plot <- ggplot(results_df, aes(x = rho, y = Tissue)) +
#   geom_point(color = "black", 
#              size = 4, 
#              alpha = 0.7) +  # Semi-transparent points
#   geom_text(aes(label = asterisks),
#             hjust = -0.3,  # Adjust horizontal position
#             vjust = 0.5,   # Adjust vertical position
#             size = 6,
#             color = "black") +  # Add asterisks next to significant points
#   labs(
#     title = "Spearman Correlation between Factor1 and GDF15 by Tissue",
#     x = "Spearman Rho",
#     y = "Tissue"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 12),
#     legend.position = "none"  # Hide legend if colors are self-explanatory
#   )
# 
# # Display the plot
# print(correlation_plot)
#      
#      
     
     
     
     
     
    
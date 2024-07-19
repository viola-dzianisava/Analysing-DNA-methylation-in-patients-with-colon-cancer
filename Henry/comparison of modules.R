# Load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("WGCNA", "data.table", "dplyr", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("foreach", "doParallel", "parallel"))

# Load required libraries
library(WGCNA)
library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(foreach)
library(doParallel)
library(parallel)

# Load the results for tumor, normal, and combined datasets
load("C:/Users/Henry/Desktop/Dissertation/WGCNA_results_tumor.RData")  # Replace with actual path
load("C:/Users/Henry/Desktop/Dissertation/WGCNA_results_normal.RData")  # Replace with actual path
load("C:/Users/Henry/Desktop/Dissertation/WGCNA_results_combined.RData")  # Replace with actual path

# Assume moduleColors and MEs are loaded for tumor, normal, and combined datasets
# They should be in variables like moduleColors_tumor, moduleColors_normal, etc.

# Function to compare modules
compareModules <- function(moduleColors1, moduleColors2, data1, data2) {
  common_modules <- intersect(unique(moduleColors1), unique(moduleColors2))
  significant_genes <- list()
  
  for (module in common_modules) {
    genes1 <- names(moduleColors1)[moduleColors1 == module]
    genes2 <- names(moduleColors2)[moduleColors2 == module]
    common_genes <- intersect(genes1, genes2)
    
    if (length(common_genes) > 0) {
      # Perform differential expression analysis between the two datasets for common genes
      data1_common <- data1[, common_genes]
      data2_common <- data2[, common_genes]
      
      # Here, we'll use a simple t-test for illustration; you may use more advanced methods as needed
      p_values <- sapply(common_genes, function(gene) {
        t.test(data1_common[, gene], data2_common[, gene])$p.value
      })
      
      significant_genes[[module]] <- names(p_values)[p.adjust(p_values, method = "BH") < 0.05]
    }
  }
  
  return(significant_genes)
}

# Compare tumor and normal datasets
significant_genes_tumor_vs_normal <- compareModules(moduleColors_tumor, moduleColors_normal, datMethyl_clean_tumor, datMethyl_clean_normal)

# Compare tumor and combined datasets
significant_genes_tumor_vs_combined <- compareModules(moduleColors_tumor, moduleColors_combined, datMethyl_clean_tumor, datMethyl_clean_combined)

# Compare normal and combined datasets
significant_genes_normal_vs_combined <- compareModules(moduleColors_normal, moduleColors_combined, datMethyl_clean_normal, datMethyl_clean_combined)

# Print and save significant genes
printSignificantGenes <- function(significant_genes, comparison_name) {
  cat("Significant genes for", comparison_name, "comparison:\n")
  for (module in names(significant_genes)) {
    cat("\nModule:", module, "\n")
    cat(significant_genes[[module]], "\n")
    write.csv(significant_genes[[module]], paste0("significant_genes_", comparison_name, "_", module, ".csv"), row.names = FALSE)
  }
}

printSignificantGenes(significant_genes_tumor_vs_normal, "tumor_vs_normal")
printSignificantGenes(significant_genes_tumor_vs_combined, "tumor_vs_combined")
printSignificantGenes(significant_genes_normal_vs_combined, "normal_vs_combined")

# Function to identify hub genes
identifyHubGenes <- function(datExpr, moduleColors, moduleEigengenes) {
  hubs <- list()
  for (module in unique(moduleColors)) {
    module_genes <- names(moduleColors[moduleColors == module])
    if (length(module_genes) > 1) {
      kME <- signedKME(datExpr[, module_genes], moduleEigengenes[[module]])
      hub_genes <- names(sort(kME[, 1], decreasing = TRUE)[1:10])  # Top 10 hub genes
      hubs[[module]] <- hub_genes
    }
  }
  return(hubs)
}

# Identify hub genes for tumor, normal, and combined datasets
hub_genes_tumor <- identifyHubGenes(datMethyl_clean_tumor, moduleColors_tumor, MEs_tumor)
hub_genes_normal <- identifyHubGenes(datMethyl_clean_normal, moduleColors_normal, MEs_normal)
hub_genes_combined <- identifyHubGenes(datMethyl_clean_combined, moduleColors_combined, MEs_combined)

# Print and save hub genes
printHubGenes <- function(hub_genes, dataset_name) {
  cat("Hub genes for", dataset_name, "dataset:\n")
  for (module in names(hub_genes)) {
    cat("\nModule:", module, "\n")
    cat(hub_genes[[module]], "\n")
    write.csv(hub_genes[[module]], paste0("hub_genes_", dataset_name, "_", module, ".csv"), row.names = FALSE)
  }
}

printHubGenes(hub_genes_tumor, "tumor")
printHubGenes(hub_genes_normal, "normal")
printHubGenes(hub_genes_combined, "combined")


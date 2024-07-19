# Load required libraries
library(WGCNA)
library(data.table)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(clusterProfiler)
library(org.Hs.eg.db)
library(foreach)
library(doParallel)
library(parallel)
library(pheatmap)

# Load the methylation data for the tumor samples
methylation_data_tumor <- fread("C:/Users/Henry/Desktop/Dissertation/GDCdata/gene_methylation_data_tumor_35_patients.csv")

# Extract the columns with CpG site names, gene names, and the actual methylation data
cpg_sites_tumor <- methylation_data_tumor$Composite
gene_names_tumor <- methylation_data_tumor$Gene
methylation_values_tumor <- methylation_data_tumor[, 5:ncol(methylation_data_tumor), with = FALSE]

# Ensure all data is numeric
methylation_values_clean_tumor <- as.data.frame(lapply(methylation_values_tumor, function(x) as.numeric(as.character(x))))
methylation_values_clean_tumor <- methylation_values_clean_tumor[, colMeans(is.na(methylation_values_clean_tumor)) < 0.5]

# Select top 5000 genes based on variance in the tumor dataset
gene_variances_tumor <- apply(methylation_values_clean_tumor, 1, var)
top_genes_tumor <- order(gene_variances_tumor, decreasing = TRUE)[1:5000]
top_methylation_data_tumor <- methylation_values_clean_tumor[top_genes_tumor, ]

# Save the CpG sites and gene names of the top 5000 genes
top_cpg_sites_tumor <- cpg_sites_tumor[top_genes_tumor]
top_gene_names_tumor <- gene_names_tumor[top_genes_tumor]

# Flatten the semicolon-separated gene names
flattened_gene_names_tumor <- unlist(strsplit(top_gene_names_tumor, ";"))

# Transpose the cleaned data for WGCNA
datMethyl_clean_tumor <- t(as.matrix(top_methylation_data_tumor))

# Ensure datMethyl_clean_tumor is numeric
datMethyl_clean_tumor <- apply(datMethyl_clean_tumor, 2, as.numeric)

# Check for good samples and genes
gsg_tumor <- goodSamplesGenes(datMethyl_clean_tumor, verbose = 3)
if (!gsg_tumor$allOK) {
  datMethyl_clean_tumor <- datMethyl_clean_tumor[gsg_tumor$goodSamples, gsg_tumor$goodGenes]
}

# Choose a set of soft-thresholding powers
powers <- c(1:10, seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft_tumor <- pickSoftThreshold(datMethyl_clean_tumor, powerVector = powers, verbose = 5)

# Choose the power based on the analysis, for example, if power = 6
softPower_tumor <- 6

# Perform blockwise WGCNA with larger block size
blockwiseModulesResult_tumor <- blockwiseModules(
  datExpr = datMethyl_clean_tumor,
  power = softPower_tumor,
  maxBlockSize = 5000,  # Adjust based on available memory and data size
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "TOM-blockwise",
  verbose = 3
)

# Save module colors and labels
moduleLabels_tumor <- blockwiseModulesResult_tumor$colors
moduleColors_tumor <- labels2colors(moduleLabels_tumor)
MEs_tumor <- blockwiseModulesResult_tumor$MEs

# Calculate adjacency matrix
adjacency_tumor <- adjacency(datMethyl_clean_tumor, power = softPower_tumor)

# Calculate Topological Overlap Matrix (TOM)
TOM_tumor <- TOMsimilarity(adjacency_tumor)

# Calculate gene connectivity
connectivity_tumor <- intramodularConnectivity(TOM_tumor, moduleColors_tumor)

# Identify hub genes
hub_genes_tumor <- chooseTopHubInEachModule(datMethyl_clean_tumor, moduleColors_tumor)

# Map hub gene IDs to indices
hub_gene_ids <- match(hub_genes_tumor, flattened_gene_names_tumor)

# Check the range of hub_gene_ids
cat("Range of hub_gene_ids:", range(hub_gene_ids), "\n")

# Verify the mapping of hub gene IDs to gene names
mapped_gene_names <- sapply(hub_gene_ids, function(id) {
  if (!is.na(id) && id <= length(flattened_gene_names_tumor)) {
    return(flattened_gene_names_tumor[id])
  } else {
    return(NA)
  }
})

# Print the top hub genes for each module with mapped names
for (i in seq_along(hub_genes_tumor)) {
  module <- names(hub_genes_tumor)[i]
  gene_id <- hub_genes_tumor[i]
  gene_name <- mapped_gene_names[i]
  cat(module, ":", gene_id, "-", gene_name, "\n")
}

# Save the results for the tumor dataset
save_path_tumor <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/TT3/WGCNA_results_tumor.RData"
save(moduleLabels_tumor, moduleColors_tumor, MEs_tumor, datMethyl_clean_tumor, connectivity_tumor, hub_genes_tumor, file = save_path_tumor)
cat("WGCNA results for the tumor dataset saved to:", save_path_tumor, "\n")







# Print the hub gene IDs
cat("Hub Gene IDs:\n")
print(hub_genes_tumor)

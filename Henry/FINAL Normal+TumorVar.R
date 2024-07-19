

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

# Load the methylation data for the tumor samples to get top 5000 genes
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
top_cpg_sites_tumor <- cpg_sites_tumor[top_genes_tumor]

# Load the methylation data for the normal samples
methylation_data_normal <- fread("C:/Users/Henry/Desktop/Dissertation/GDCdata/gene_methylation_data_normal.csv")

# Extract the columns with CpG site names, gene names, and the actual methylation data
cpg_sites_normal <- methylation_data_normal$Composite
gene_names_normal <- methylation_data_normal$Gene
methylation_values_normal <- methylation_data_normal[, 5:ncol(methylation_data_normal), with = FALSE]

# Ensure all data is numeric
methylation_values_clean_normal <- as.data.frame(lapply(methylation_values_normal, function(x) as.numeric(as.character(x))))
methylation_values_clean_normal <- methylation_values_clean_normal[, colMeans(is.na(methylation_values_clean_normal)) < 0.5]

# Select the same genes (CpG sites) identified in the tumor dataset
selected_genes_normal <- methylation_values_clean_normal[which(cpg_sites_normal %in% top_cpg_sites_tumor), ]

# Transpose the cleaned data for WGCNA
datMethyl_clean_normal <- t(as.matrix(selected_genes_normal))

# Ensure datMethyl_clean_normal is numeric
datMethyl_clean_normal <- apply(datMethyl_clean_normal, 2, as.numeric)

# Check for good samples and genes
gsg_normal <- goodSamplesGenes(datMethyl_clean_normal, verbose = 3)
if (!gsg_normal$allOK) {
  datMethyl_clean_normal <- datMethyl_clean_normal[gsg_normal$goodSamples, gsg_normal$goodGenes]
}

# Choose a set of soft-thresholding powers
powers <- c(1:10, seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft_normal <- pickSoftThreshold(datMethyl_clean_normal, powerVector = powers, verbose = 5)

# Plot the results of the soft-thresholding analysis
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft_normal$fitIndices[, 1], -sign(sft_normal$fitIndices[, 3]) * sft_normal$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", 
     main = "Scale independence")
text(sft_normal$fitIndices[, 1], -sign(sft_normal$fitIndices[, 3]) * sft_normal$fitIndices[, 2], 
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")

plot(sft_normal$fitIndices[, 1], sft_normal$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", 
     main = "Mean connectivity")
text(sft_normal$fitIndices[, 1], sft_normal$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

# Choose the power based on the analysis, for example, if power = 6
softPower_normal <- 6

# Perform blockwise WGCNA with larger block size
blockwiseModulesResult_normal <- blockwiseModules(
  datExpr = datMethyl_clean_normal,
  power = softPower_normal,
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
moduleLabels_normal <- blockwiseModulesResult_normal$colors
moduleColors_normal <- labels2colors(moduleLabels_normal)
MEs_normal <- blockwiseModulesResult_normal$MEs

# Plot the dendrogram and module colors
sizeGrWindow(8, 6)
plotDendroAndColors(blockwiseModulesResult_normal$dendrograms[[1]], moduleColors_normal[blockwiseModulesResult_normal$blockGenes[[1]]], 
                    "Module colors", dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)

# Calculate adjacency matrix
adjacency_normal <- adjacency(datMethyl_clean_normal, power = softPower_normal)

# Calculate Topological Overlap Matrix (TOM)
TOM_normal <- TOMsimilarity(adjacency_normal)

# Calculate gene connectivity
connectivity_normal <- intramodularConnectivity(TOM_normal, moduleColors_normal)

# Identify hub genes
hub_genes_normal <- chooseTopHubInEachModule(datMethyl_clean_normal, moduleColors_normal)

# Save the results for the normal dataset
save_path_normal <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/NT/WGCNA_results_normal_top5000_tumor.RData"
save(moduleLabels_normal, moduleColors_normal, MEs_normal, datMethyl_clean_normal, connectivity_normal, hub_genes_normal, file = save_path_normal)
cat("WGCNA results for the normal dataset saved to:", save_path_normal, "\n")

# Print the top hub genes for each module
print("Top hub genes for each module in the normal dataset:")
print(hub_genes_normal)

# Example code to print the most connected genes in a specific module, e.g., blue module
blue_module_genes_normal <- names(connectivity_normal$kWithin[moduleColors_normal == "blue"])
blue_module_gene_connectivity_normal <- connectivity_normal$kWithin[moduleColors_normal == "blue"]
top_connected_genes_blue_normal <- blue_module_genes_normal[order(blue_module_gene_connectivity_normal, decreasing = TRUE)[1:10]]

print("Top connected genes in the blue module (normal dataset):")
print(top_connected_genes_blue_normal)

# Ensure the output directory exists
output_dir <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/NT/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create and save heatmaps for all modules
unique_modules_normal <- unique(moduleColors_normal)
for (module in unique_modules_normal) {
  selectedGenes_normal <- colnames(datMethyl_clean_normal)[moduleColors_normal == module]
  selectedData_normal <- datMethyl_clean_normal[, selectedGenes_normal, drop = FALSE]
  
  # Ensure there are no NA/NaN/Inf values
  selectedData_normal[!is.finite(selectedData_normal)] <- 0
  
  # Create and save the heatmap for the selected module
  heatmap_file <- paste0(output_dir, "heatmap_normal_", module, ".png")
  tryCatch({
    png(heatmap_file, width = 1200, height = 900)
    pheatmap(selectedData_normal, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, 
             main = paste("Heatmap of", module, "Module"), show_rownames = FALSE, show_colnames = FALSE)
    dev.off()
    cat("Saved heatmap for module:", module, "\n")
  }, error = function(e) {
    cat("Failed to save heatmap for module:", module, "\nError message:", e$message, "\n")
  })
}

cat("Heatmaps for all modules saved.\n")

# Create a WGCNA heatmap
sizeGrWindow(12, 9)
TOMplot(dissim = 1 - TOM_normal, dendro = blockwiseModulesResult_normal$dendrograms[[1]], Colors = moduleColors_normal, main = "Network heatmap plot (Normal Dataset)")

# Create a heatmap of the eigengene networks (clusters)
MEs_col <- orderMEs(MEs_normal)
sizeGrWindow(8, 6)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", marDendro = c(0, 4, 2, 0), marHeatmap = c(3, 4, 2, 1), plotDendrograms = TRUE, xLabelsAngle = 90)

# Save the heatmap plots to files
heatmap_tom_file <- paste0(output_dir, "WGCNA_heatmap_normal.png")
heatmap_clusters_file <- paste0(output_dir, "WGCNA_clusters_heatmap_normal.png")

png(heatmap_tom_file, width = 1200, height = 900)
TOMplot(dissim = 1 - TOM_normal, dendro = blockwiseModulesResult_normal$dendrograms[[1]], Colors = moduleColors_normal, main = "Network heatmap plot (Normal Dataset)")
dev.off()

png(heatmap_clusters_file, width = 800, height = 600)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", marDendro = c(0, 4, 2, 0), marHeatmap = c(3, 4, 2, 1), plotDendrograms = TRUE, xLabelsAngle = 90)
dev.off()

cat("WGCNA heatmap and cluster heatmap saved to:", heatmap_tom_file, "and", heatmap_clusters_file, "\n")

# Save the tumor dataset results into NT directory as well
file.copy("C:/Users/Henry/Desktop/Dissertation/WGCNA_results_tumor.RData", "C:/Users/Henry/Desktop/Dissertation/FinalFindings/NT/WGCNA_results_tumor.RData")
cat("Tumor dataset results copied to NT directory.\n")

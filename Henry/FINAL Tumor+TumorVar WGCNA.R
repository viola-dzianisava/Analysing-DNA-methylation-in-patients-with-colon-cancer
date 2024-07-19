
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

# Save the CpG sites of the top 5000 genes
top_cpg_sites_tumor <- cpg_sites_tumor[top_genes_tumor]

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

# Plot the results of the soft-thresholding analysis
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft_tumor$fitIndices[, 1], -sign(sft_tumor$fitIndices[, 3]) * sft_tumor$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", 
     main = "Scale independence")
text(sft_tumor$fitIndices[, 1], -sign(sft_tumor$fitIndices[, 3]) * sft_tumor$fitIndices[, 2], 
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")

plot(sft_tumor$fitIndices[, 1], sft_tumor$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", 
     main = "Mean connectivity")
text(sft_tumor$fitIndices[, 1], sft_tumor$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

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

# Plot the dendrogram and module colors
sizeGrWindow(8, 6)
plotDendroAndColors(blockwiseModulesResult_tumor$dendrograms[[1]], moduleColors_tumor[blockwiseModulesResult_tumor$blockGenes[[1]]], 
                    "Module colors", dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)

# Calculate adjacency matrix
adjacency_tumor <- adjacency(datMethyl_clean_tumor, power = softPower_tumor)

# Calculate Topological Overlap Matrix (TOM)
TOM_tumor <- TOMsimilarity(adjacency_tumor)

# Calculate gene connectivity
connectivity_tumor <- intramodularConnectivity(TOM_tumor, moduleColors_tumor)

# Identify hub genes
hub_genes_tumor <- chooseTopHubInEachModule(datMethyl_clean_tumor, moduleColors_tumor)

output_dir <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/NN/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the results for the tumor dataset
save_path_tumor <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/NN/WGCNA_results_tumor.RData"
save(moduleLabels_tumor, moduleColors_tumor, MEs_tumor, datMethyl_clean_tumor, connectivity_tumor, hub_genes_tumor, file = save_path_tumor)
cat("WGCNA results for the tumor dataset saved to:", save_path_tumor, "\n")

# Print the top hub genes for each module
print("Top hub genes for each module in the tumor dataset:")
print(hub_genes_tumor)

# Example code to print the most connected genes in a specific module, e.g., blue module
blue_module_genes_tumor <- names(connectivity_tumor$kWithin[moduleColors_tumor == "blue"])
blue_module_gene_connectivity_tumor <- connectivity_tumor$kWithin[moduleColors_tumor == "blue"]
top_connected_genes_blue_tumor <- blue_module_genes_tumor[order(blue_module_gene_connectivity_tumor, decreasing = TRUE)[1:10]]

print("Top connected genes in the blue module (tumor dataset):")
print(top_connected_genes_blue_tumor)

output_dir <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create and save heatmaps for all modules
unique_modules_tumor <- unique(moduleColors_tumor)
for (module in unique_modules_tumor) {
  selectedGenes_tumor <- colnames(datMethyl_clean_tumor)[moduleColors_tumor == module]
  selectedData_tumor <- datMethyl_clean_tumor[, selectedGenes_tumor, drop = FALSE]
  
  # Ensure there are no NA/NaN/Inf values
  selectedData_tumor[!is.finite(selectedData_tumor)] <- 0
  
  # Create and save the heatmap for the selected module
  heatmap_file <- paste0("C:/Users/Henry/Desktop/Dissertation/FinalFindings/heatmap_tumor_", module, ".png")
  png(heatmap_file, width = 1200, height = 900)
  pheatmap(selectedData_tumor, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, 
           main = paste("Heatmap of", module, "Module"), show_rownames = FALSE, show_colnames = FALSE)
  dev.off()
}

cat("Heatmaps for all modules saved.\n")

# Create a WGCNA heatmap
sizeGrWindow(12, 9)
TOMplot(dissTOM = TOM_tumor, dendro = blockwiseModulesResult_tumor$dendrograms[[1]], Colors = moduleColors_tumor, main = "Network heatmap plot (Tumor Dataset)")

# Create a heatmap of the eigengene networks (clusters)
MEs_col <- orderMEs(MEs_tumor)
sizeGrWindow(8, 6)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", marDendro = c(0, 4, 2, 0), marHeatmap = c(3, 4, 2, 1), plotDendrograms = TRUE, xLabelsAngle = 90)

# Save the heatmap plots to files
heatmap_tom_file <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/WGCNA_heatmap_tumor.png"
heatmap_clusters_file <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/WGCNA_clusters_heatmap_tumor.png"

png(heatmap_tom_file, width = 1200, height = 900)
TOMplot(dissim = 1 - TOM_tumor, dendro = blockwiseModulesResult_tumor$dendrograms[[1]], Colors = moduleColors_tumor, main = "Network heatmap plot (Tumor Dataset)")
dev.off()

png(heatmap_clusters_file, width = 800, height = 600)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", marDendro = c(0, 4, 2, 0), marHeatmap = c(3, 4, 2, 1), plotDendrograms = TRUE, xLabelsAngle = 90)
dev.off()

cat("WGCNA heatmap and cluster heatmap saved to:", heatmap_tom_file, "and", heatmap_clusters_file, "\n")

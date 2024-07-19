# Load the tumor data
methylation_data_tumor <- fread("/Users/henry/Desktop/Dissertation/GDCdata/gene_methylation_data_tumor_35_patients.csv")

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
top_cpg_sites <- cpg_sites_tumor[top_genes_tumor]


# Load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("WGCNA", "data.table", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("foreach", "doParallel", "parallel", "pheatmap"))

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

# Specify download directory
download_dir <- "C:/Users/Henry/Desktop/Dissertation/GDCdata"
if (!dir.exists(download_dir)) {
  dir.create(download_dir, recursive = TRUE)
}

# Query and download DNA methylation data for normal samples in TCGA-COAD
query_TCGA_COAD_normal <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  workflow.type = "SeSAMe Methylation Beta Estimation",
  sample.type = "Solid Tissue Normal",
  platform = "Illumina Human Methylation 450"
)

GDCdownload(query = query_TCGA_COAD_normal, directory = download_dir)

# Prepare and clean data
tcga_data_coad_normal <- GDCprepare(query = query_TCGA_COAD_normal, directory = download_dir)
gene_methylation_data_normal <- as.data.frame(assay(tcga_data_coad_normal))
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot_df <- as.data.frame(annot)
merged_data_normal <- merge(gene_methylation_data_normal, annot_df, by.x = "row.names", by.y = "Name", all.x = TRUE)

# Select and rename necessary columns
tcga_columns_normal <- grep("^TCGA", colnames(merged_data_normal), value = TRUE)
final_data_normal <- merged_data_normal[, c("Row.names", "chr", "pos", "UCSC_RefGene_Name", tcga_columns_normal)]
colnames(final_data_normal)[1:4] <- c("Composite", "Chr", "Coordinate", "Gene")

# Save the final data
write.csv(final_data_normal, "C:/Users/Henry/Desktop/Dissertation/GDCdata/gene_methylation_data_normal.csv", row.names = FALSE)

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
selected_genes_normal <- methylation_values_clean_normal[which(cpg_sites_normal %in% top_cpg_sites), ]

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

# List all available modules
unique_modules_normal <- unique(moduleColors_normal)
print(unique_modules_normal)

# Create and save heatmaps for all modules
for (module in unique_modules_normal) {
  selectedGenes_normal <- colnames(datMethyl_clean_normal)[moduleColors_normal == module]
  selectedData_normal <- datMethyl_clean_normal[, selectedGenes_normal, drop = FALSE]
  
  # Ensure there are no NA/NaN/Inf values
  selectedData_normal[!is.finite(selectedData_normal)] <- 0
  
  # Create and save the heatmap for the selected module
  pheatmap(selectedData_normal, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, 
           main = paste("Heatmap of", module, "Module"), show_rownames = FALSE, show_colnames = FALSE)
}

# Print key information from the WGCNA analysis
cat("Module labels:\n")
print(moduleLabels_normal)
cat("\nModule colors:\n")
print(moduleColors_normal)
cat("\nModule Eigengenes (first few rows):\n")
print(head(MEs_normal))

# Save the module colors and labels
write.csv(moduleLabels_normal, "moduleLabels_normal.csv", row.names = FALSE)
write.csv(moduleColors_normal, "moduleColors_normal.csv", row.names = FALSE)
write.csv(MEs_normal, "moduleEigengenes_normal.csv", row.names = FALSE)

# Define the path where you want to save the .RData file
save_path_combined <- "C:/Users/Henry/Desktop/Dissertation/WGCNA_results_combined.RData"

# Save the necessary objects
save(moduleLabels_normal, moduleColors_normal, MEs_normal, datMethyl_clean_normal, file = save_path_combined)

cat("WGCNA results for the combined dataset saved to:", save_path_combined, "\n")

library(TCGAbiolinks)
library(limma)
library(edgeR)
library(SummarizedExperiment)
library(caret)
library(doParallel)
library(randomForest)
library(gplots)
library(RColorBrewer)

# Define download directory
download_dir <- "/Users/linxijun/Desktop/GDCdata"
if (!dir.exists(download_dir)) {
  dir.create(download_dir, recursive = TRUE)
}

# Query and prepare DNA methylation data for TCGA-COAD
query_TCGA_COAD <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  workflow.type = "SeSAMe Methylation Beta Estimation",
  platform = "Illumina Human Methylation 450"
)
tcga_data_coad <- GDCprepare(query = query_TCGA_COAD, directory = download_dir)

# Define differential expression analysis function
limma_pipeline_na_permissive <- function(tcga_data_coad, condition_variable, reference_group=NULL) {
  if (!condition_variable %in% colnames(colData(tcga_data_coad))) {
    stop("The condition variable is not found in colData of the SummarizedExperiment object.")
  }
  
  group <- factor(colData(tcga_data_coad)[[condition_variable]])
  if (!is.null(reference_group) && reference_group %in% levels(group)) {
    group <- relevel(group, ref = reference_group)
  }
  
  design <- model.matrix(~ group)
  counts <- assay(tcga_data_coad)
  valid_rows <- rowSums(is.na(counts)) == 0
  counts_non_na <- counts[valid_rows, ]
  
  dge <- DGEList(counts = counts_non_na)
  dge <- calcNormFactors(dge)
  
  voomObj <- voom(dge, design, plot = FALSE)
  fit <- lmFit(voomObj, design)
  fit <- eBayes(fit)
  
  topGenes <- topTable(fit, coef = ncol(design), number = 500, sort.by = "p", p.value = 0.1)
  
  return(list(voomObj = voomObj, fit = fit, topGenes = topGenes))
}

# Run limma analysis pipeline
limma_result <- limma_pipeline_na_permissive(tcga_data_coad, "definition", "Solid Tissue Normal")

# Check the number of genes in the results
if (nrow(limma_result$topGenes) < 10) {
  stop("Not enough relevant genes to plot the heatmap. Consider adjusting DE analysis settings.")
}

###==============PCA Tumor stage========
plot_PCA <- function(voomObj, data, condition_variable) {
  if (!condition_variable %in% colnames(data)) {
    stop("condition_variable not found in the provided data")
  }
  
  group <- factor(data[[condition_variable]])
  
  pca <- prcomp(t(voomObj$E))
  
  plot(pca$x[, 1:2], col=as.numeric(group), pch=19, xlab="PC1", ylab="PC2", main=paste("PCA Plot based on", condition_variable),)
  legend("topright", legend=levels(group), col=1:length(levels(group)), pch=19)
  
  return(pca)
}

# Perform PCA
if (!is.null(limma_result) && !is.null(limma_result$voomObj)) {
  res_pca <- plot_PCA(limma_result$voomObj, colData(tcga_data_coad), "definition")
  plot_PCA(limma_result$voomObj, colData(tcga_data_coad), "gender")
} else {
  print("Voom object or necessary metadata is not available.")
}

variance_explained <- summary(res_pca)$importance[2, 1:2]
print(variance_explained)
loadings <- res_pca$rotation[, 1:2]
head(loadings)

# Proportion of variance explained by PC1 and PC2
variance_explained <- summary(res_pca)$importance[2, 1:2]
print(paste("Variance explained by PC1:", round(variance_explained[1] * 100, 2), "%"))
print(paste("Variance explained by PC2:", round(variance_explained[2] * 100, 2), "%"))

# Loadings of the top genes on PC1 and PC2
loadings <- res_pca$rotation[, 1:2]
top_genes_PC1 <- head(loadings[order(abs(loadings[, 1]), decreasing = TRUE), 1])
top_genes_PC2 <- head(loadings[order(abs(loadings[, 2]), decreasing = TRUE), 2])

print("Top genes contributing to PC1:")
print(top_genes_PC1)

print("Top genes contributing to PC2:")
print(top_genes_PC2)

gene_info <- data.frame(CpG_id = rownames(rowData(tcga_data_coad)), gene = rowData(tcga_data_coad)$gene)
methylation_data <- data.frame(CpG_id = rownames(assay(tcga_data_coad)), assay(tcga_data_coad))

# Merge the data frames based on CpG_id
annotated_methylation_data <- merge(gene_info, methylation_data, by = "CpG_id")

# Set rownames to gene names for easier access later
rownames(annotated_methylation_data) <- annotated_methylation_data$gene

# Drop the CpG_id and gene columns if they are no longer needed
annotated_methylation_data <- annotated_methylation_data[, -c(1,2)]

# Inspect the first few rows of the annotated data
head(annotated_methylation_data)

merged_data <- cbind(gene_info, gene_methylation_data)

##============PCA race========
library(RColorBrewer)

plot_PCA <- function(voomObj, data, condition_variable, legend_size = 0.6) {
  if (!condition_variable %in% colnames(data)) {
    stop(paste("Condition variable", condition_variable, "not found in the provided data"))
  }
  
  # Gets grouping information for the specified condition variable, explicitly excluding 'not reported'
  group <- factor(data[[condition_variable]], levels = c("american indian or alaska native", "asian", "black or african american", "white"))
  
  valid_indices <- which(!is.na(group))
  group <- group[valid_indices]
  valid_data <- t(voomObj$E)[, valid_indices]
  
  pca <- prcomp(valid_data, scale. = TRUE)

  num_colors <- length(levels(group))
  color_palette <- brewer.pal(num_colors, "Set1")
  
  plot_colors <- color_palette[as.numeric(group)]
  
  plot(pca$x[, 1:2], col=plot_colors, pch=19, xlab="PC1", ylab="PC2", 
       main=paste("PCA Plot based on", condition_variable))
  

  legend("topright", legend=levels(group), col=color_palette, pch=19, cex=legend_size)
  
  return(pca)
}


res_pca_race_filtered <- plot_PCA(limma_result$voomObj, colData(tcga_data_coad), "race")


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
  
  plot(pca$x[, 1:2], col=as.numeric(group), pch=19, xlab="PC1", ylab="PC2", main=paste("PCA Plot based on", condition_variable))
  legend("topright", legend=levels(group), col=1:length(levels(group)), pch=19)
  
  return(pca)
}

# Perform PCA for the "definition" variable
if (!is.null(limma_result) && !is.null(limma_result$voomObj)) {
  res_pca <- plot_PCA(limma_result$voomObj, colData(tcga_data_coad), "definition")
} else {
  print("Voom object or necessary metadata is not available.")
}

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

# Extract the gene names from the rowData of the tcga_data_coad object
gene_info <- data.frame(CpG_id = rownames(rowData(tcga_data_coad)), gene = rowData(tcga_data_coad)$gene)

# Extract the methylation data from the assay of the tcga_data_coad object
methylation_data <- data.frame(CpG_id = rownames(assay(tcga_data_coad)), assay(tcga_data_coad))

# Merge the gene names with the methylation data based on the CpG_id
annotated_methylation_data <- merge(gene_info, methylation_data, by = "CpG_id")

# Inspect the first few rows of the annotated methylation data
head(annotated_methylation_data)

# Create function to get top genes
get_top_genes <- function(top_genes) {
  top_gene_ids <- names(top_genes)  # 提取CpG_id
  top_gene_names <- gene_info[gene_info$CpG_id %in% top_gene_ids, ]
  return(top_gene_names)
}


top_genes_PC1_names <- get_top_genes(top_genes_PC1)

top_genes_PC2_names <- get_top_genes(top_genes_PC2)


print("Top genes contributing to PC1:")
print(top_genes_PC1_names)

print("Top genes contributing to PC2:")
print(top_genes_PC2_names)


# Extract the top 15 genes contributing to PC1 and PC2
top_genes_PC1 <- head(loadings[order(abs(loadings[, 1]), decreasing = TRUE), 1], 15)
top_genes_PC2 <- head(loadings[order(abs(loadings[, 2]), decreasing = TRUE), 2], 15)

# Function to get gene names for top contributing CpG sites
get_top_genes <- function(top_genes) {
  top_gene_ids <- names(top_genes)  # Extract CpG_id
  top_gene_names <- gene_info[gene_info$CpG_id %in% top_gene_ids, ]
  return(top_gene_names)
}

# Get top genes for PC1
top_genes_PC1_names <- get_top_genes(top_genes_PC1)
# Get top genes for PC2
top_genes_PC2_names <- get_top_genes(top_genes_PC2)

# Print the top genes contributing to PC1 and PC2
print("Top genes contributing to PC1:")
print(top_genes_PC1_names)

print("Top genes contributing to PC2:")
print(top_genes_PC2_names)


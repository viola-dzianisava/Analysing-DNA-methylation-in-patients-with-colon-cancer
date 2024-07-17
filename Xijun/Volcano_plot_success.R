library(TCGAbiolinks)
library(limma)
library(edgeR)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(ggrepel)

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
GDCdownload(query = query_TCGA_COAD, directory = download_dir)
tcga_data_coad <- GDCprepare(query = query_TCGA_COAD, directory = download_dir)

# Extract the genetic information and set it to the line name of the DataFrame
gene_info <- rowData(tcga_data_coad)[, "gene", drop = FALSE]
rownames(gene_info) <- rownames(rowData(tcga_data_coad))

# Switch the methylation data from SummarizedExperiment into DataFrame
gene_methylation_data <- as.data.frame(assay(tcga_data_coad))
rownames(gene_methylation_data) <- rownames(rowData(tcga_data_coad))

# merged gene info and methylation data
merged_data <- cbind(gene_info, gene_methylation_data)

# Normal and tumor sample ids were extracted
normal_samples <- colData(tcga_data_coad)$submitter_id[colData(tcga_data_coad)$sample_type == "Solid Tissue Normal"]
tumor_samples <- colData(tcga_data_coad)$submitter_id[colData(tcga_data_coad)$sample_type %in% c("Metastatic", "Primary Tumor", "Recurrent Tumor")]

cat("Number of normal samples:", length(normal_samples), "\n")
cat("Number of tumor samples:", length(tumor_samples), "\n")

# Check that the column name matches the sample ID
sample_columns <- colnames(merged_data)
sample_columns <- sapply(sample_columns, function(x) substr(x, 1, 12))

# Match
normal_expression_data <- merged_data[, sample_columns %in% normal_samples]
tumor_expression_data <- merged_data[, sample_columns %in% tumor_samples]

cat("Columns in normal_expression_data:", ncol(normal_expression_data), "\n")
cat("Columns in tumor_expression_data:", ncol(tumor_expression_data), "\n")

# Combine normal and tumor data
combined_data <- cbind(normal_expression_data, tumor_expression_data)
group <- factor(c(rep("Normal", ncol(normal_expression_data)), rep("Tumor", ncol(tumor_expression_data))))

# Create design matrix
design <- model.matrix(~ group)

# Perform differential expression analysis using limma
fit <- lmFit(combined_data, design)
fit <- eBayes(fit)

# Extract results
results <- topTable(fit, coef = 2, sort.by = "P", adjust.method = "BH", number = Inf)

# Add gene name in results
results$gene <- gene_info$gene[match(rownames(results), rownames(gene_info))]

# Set thresholds and focus on genes
results$Category <- ifelse(results$adj.P.Val > 0.05, "Not significant",
                           ifelse(results$logFC > log2(1.06), "Up regulated",
                                  ifelse(results$logFC < -log2(1.06), "Down regulated", "Not significant")))

# Set thresholds and gene labels
significance_threshold <- -log10(0.05)
logFC_threshold_high <- log2(1.06)
logFC_threshold_low <- -log2(1.06)

# Define special gene
genes_to_label <- c("DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "MBD1", "MBD2", "MBD3", "MBD4", "UHRF1", "UHRF2", "ZBTB4", "ZBTB38")

# Label special gene
results$SpecialHighlight <- ifelse(results$gene %in% genes_to_label, "Special", "Regular")


significant_special_genes <- results %>%
  filter(SpecialHighlight == "Special" & Category != "Not significant")

# Check
table(results$SpecialHighlight)

# Generate the volcano plot
p <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Category)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_point(data = significant_special_genes, aes(color = "Special"), size = 3, shape = 21) +
  geom_text_repel(data = significant_special_genes, aes(label = gene), size = 5, color = "black", box.padding = 0.35, point.padding = 0.5, max.overlaps = 50) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-log2(1.06), log2(1.06)), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Down regulated" = "blue", 
                                "Not significant" = "grey", 
                                "Up regulated" = "red", 
                                "Special" = "green")) +  # Combine all color definitions into one scale
  labs(title = "Enhanced Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_minimal()

# Print the plot
print(p)

# Set thresholds and focus on genes
results$Category <- ifelse(results$adj.P.Val > 0.03, "Not significant",
                           ifelse(results$logFC > log2(1.03), "Up regulated",
                                  ifelse(results$logFC < -log2(1.03), "Down regulated", "Not significant")))

# Set thresholds and gene labels
significance_threshold <- -log10(0.03)
logFC_threshold_high <- log2(1.03)
logFC_threshold_low <- -log2(1.03)
# Define special gene
genes_to_label <- c("DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "MBD1", "MBD2", "MBD3", "MBD4", "UHRF1", "UHRF2", "ZBTB4", "ZBTB38")

results$SpecialHighlight <- ifelse(results$gene %in% genes_to_label, "Special", "Regular")

# Filter the special genes
significant_special_genes <- results %>%
  filter(SpecialHighlight == "Special" & Category != "Not significant")k
# Double checj
table(results$SpecialHighlight)
# Filter special genes and ensure they are significant by adjusting p-value criteria
significant_special_genes <- results %>%
  filter(SpecialHighlight == "Special" & adj.P.Val < 0.03)

# Generate the volcano plot focusing only on significant special genes
p2 <- ggplot(significant_special_genes, aes(x = logFC, y = -log10(adj.P.Val), color = Category)) +
  geom_point(alpha = 0.75, size = 3.5, shape = 21) +
  geom_text_repel(
    aes(label = gene),
    size = 5,
    box.padding = 0.35,
    point.padding = 0.5,
    max.overlaps = 10
  ) +
  geom_hline(yintercept = -log10(0.03), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-log2(1.03), log2(1.03)), linetype = "dashed", color = "black") +
  scale_color_manual(values = c(
    "Down regulated" = "blue", 
    "Not significant" = "grey", 
    "Up regulated" = "red"
  )) +
  labs(
    title = "Enhanced Volcano Plot: Special Genes",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  ) +
  coord_cartesian(ylim = c(0, 10))  # Adjust the y-axis limits to zoom into significant areas more

# Print the plot
print(p2)


# Load required packages
library("sesameData")
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("dplyr")
library("pheatmap")
library("ggplot2")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# Specify download directory
download_dir <- "/Users/linxijun/Desktop/GDCdata"
if (!dir.exists(download_dir)) {
  dir.create(download_dir, recursive = TRUE)
}

# Query DNA methylation data for TCGA-COAD (Solid Tissue Normal)
query_TCGA_COAD_nm <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  workflow.type = "SeSAMe Methylation Beta Estimation",
  sample.type = "Solid Tissue Normal",
  platform = "Illumina Human Methylation 450"
)

# Download data to the specified directory
GDCdownload(query = query_TCGA_COAD_nm, directory = download_dir)

# Prepare data and specify the directory
tcga_data_coad_nm <- GDCprepare(query = query_TCGA_COAD_nm, directory = download_dir)

# Extract assay data and convert to data frame
gene_methylation_data_nm <- as.data.frame(assay(tcga_data_coad_nm))

# Load annotation data from the package
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Retrieve the detailed annotation data
annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Convert the annotation data to a data frame
annot_df <- as.data.frame(annot)

# Merge methylation data with annotation
merged_data <- merge(gene_methylation_data_nm, annot_df, by.x = "row.names", by.y = "Name", all.x = TRUE)

# Select and rename necessary columns using base R
tcga_columns <- grep("^TCGA", colnames(merged_data), value = TRUE)
final_data <- merged_data[, c("Row.names", "chr", "pos", "UCSC_RefGene_Name", tcga_columns)]
colnames(final_data)[1:4] <- c("Composite", "Chr", "Coordinate", "Gene")

# Function to remove duplicate gene names separated by semicolons
clean_gene_names <- function(gene) {
  if (!is.na(gene) && gene != "") {
    unique_genes <- unique(unlist(strsplit(gene, ";")))
    return(paste(unique_genes, collapse = ";"))
  }
  return(gene)
}

# Apply the function to clean up gene names
final_data$Gene <- sapply(final_data$Gene, clean_gene_names)

# Remove rows where all TCGA columns are NA or NULL
final_data <- final_data[rowSums(is.na(final_data[tcga_columns]) | final_data[tcga_columns] == "") < length(tcga_columns), ]

# Display the final data
print(head(final_data))

# Save the final data to a CSV file (optional)
write.csv(final_data, "/Users/linxijun/Desktop/gene_methylation_data_nm.csv", row.names = FALSE)

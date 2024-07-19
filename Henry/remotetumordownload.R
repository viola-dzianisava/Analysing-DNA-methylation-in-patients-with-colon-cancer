# Install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("sesame")
BiocManager::install("sesameData")
BiocManager::install("TCGAbiolinks")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("SummarizedExperiment")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# List of required packages
required_packages <- c(
  "sesameData", "TCGAbiolinks", "limma", "edgeR", "SummarizedExperiment", 
  "IlluminaHumanMethylation450kanno.ilmn12.hg19", "glmnet", "factoextra", 
  "FactoMineR", "caret", "gplots", "survival", "survminer", "RColorBrewer", 
  "gProfileR", "genefilter", "dplyr", "pheatmap", "ggplot2", "WGCNA"
)

# Function to check and install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      if (pkg %in% BiocManager::available() || pkg %in% rownames(installed.packages())) {
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }
}

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install missing packages
install_if_missing(required_packages)

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
library("sesame")
sesameDataCache()

# Specify download directory
download_dir <- "C:/Users/Henry/Desktop/Dissertation/GDCdata"
if (!dir.exists(download_dir)) {
  dir.create(download_dir, recursive = TRUE)
}

# Query and download DNA methylation data for tumor samples in TCGA-COAD
query_TCGA_COAD_tumor <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  workflow.type = "SeSAMe Methylation Beta Estimation",
  sample.type = "Primary Tumor",
  platform = "Illumina Human Methylation 450"
)

GDCdownload(query = query_TCGA_COAD_tumor, directory = download_dir)

# Prepare and clean data
tcga_data_coad_tumor <- GDCprepare(query = query_TCGA_COAD_tumor, directory = download_dir)
gene_methylation_data_tumor <- as.data.frame(assay(tcga_data_coad_tumor))

# Check the dimensions of the data
cat("Dimensions of gene methylation data (tumor):", dim(gene_methylation_data_tumor), "\n")

# Randomly select 35 patients
set.seed(123) # Set seed for reproducibility
selected_patients <- sample(colnames(gene_methylation_data_tumor), 35)
gene_methylation_data_tumor <- gene_methylation_data_tumor[, selected_patients]

# Check the dimensions after selecting patients
cat("Dimensions of gene methylation data after selecting 35 patients:", dim(gene_methylation_data_tumor), "\n")

# Load annotation data
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot_df <- as.data.frame(annot)

# Merge gene methylation data with annotation
merged_data <- merge(gene_methylation_data_tumor, annot_df, by.x = "row.names", by.y = "Name", all.x = TRUE)

# Check the first few rows of the merged data
cat("First few rows of the merged data:\n")
print(head(merged_data))


# Extract relevant columns for gene annotation
annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot_df <- as.data.frame(annot)

# Verify the unique identifiers in the annotation data
cat("Unique gene IDs in annotation data:\n")
print(head(annot_df$UCSC_RefGene_Name))

# Merge gene methylation data with annotation data
merged_data <- merge(gene_methylation_data_tumor, annot_df, by.x = "row.names", by.y = "Name", all.x = TRUE)

# Select and rename necessary columns
tcga_columns <- grep("^TCGA", colnames(merged_data), value = TRUE)
final_data <- merged_data[, c("Row.names", "chr", "pos", "UCSC_RefGene_Name", tcga_columns)]
colnames(final_data)[1:4] <- c("Composite", "Chr", "Coordinate", "Gene")

# Save the final data
write.csv(final_data, "C:/Users/Henry/Desktop/Dissertation/GDCdata/gene_methylation_data_tumor_35_patients.csv", row.names = FALSE)

# Load necessary libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Getting Clinical Data
clinical_coad <- GDCquery_clinic("TCGA-COAD")

# Entire data query
query_coad_all <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  experimental.strategy = "Methylation Array",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)

# Download and prepare data
GDCdownload(query_coad_all)
data_coad_all <- GDCprepare(query_coad_all)

# Extract assay data and convert to data frame
gene_methylation_data <- as.data.frame(assay(data_coad_all))

# Add necessary annotation
annot <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

# Merge methylation data with annotation
merged_data <- merge(gene_methylation_data, annot, by.x = "row.names", by.y = "Name", all.x = TRUE)

# Check column names
colnames(merged_data)

# Select and rename necessary columns using base R
tcga_columns <- grep("^TCGA", colnames(merged_data), value = TRUE)
final_data <- merged_data[, c("Row.names", "chr", "pos", "UCSC_RefGene_Name", tcga_columns)]
colnames(final_data)[1:4] <- c("Composite", "Chr", "Coordinate", "Gene")

# Display the final data
print(head(final_data))

# Save the final data to a CSV file (optional)
write.csv(final_data, "gene_methylation_data.csv", row.names = FALSE)





# ================Query DNA methylation data======
query_TCGA_COAD <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  workflow.type = "SeSAMe Methylation Beta Estimation",
  sample.type = c("Metastatic", "Primary Tumor", "Recurrent Tumor", "Solid Tissue Normal"),
  platform = "Illumina Human Methylation 450"
)

# Specify download directory
download_dir <- "/Users/linxijun/Desktop/GDCdata"
if (!dir.exists(download_dir)) {
  dir.create(download_dir, recursive = TRUE)
}

# Download data to the specified directory
GDCdownload(query = query_TCGA_COAD, directory = download_dir)

# Cache SeSAMe auxiliary data
sesameDataCacheAll()

# Prepare data
tcga_data_coad <- GDCprepare(query = query_TCGA_COAD)

# Convert DFrame to a regular data.frame
data_df_coad <- as.data.frame(colData(tcga_data_coad))

# Convert list columns to character
list_cols <- sapply(data_df_coad, is.list)
data_df_coad[list_cols] <- lapply(data_df_coad[list_cols], function(x) sapply(x, toString))

# Set working directory to Desktop
setwd("/Users/linxijun/Desktop")

# Export the data to a CSV file
write.csv(data_df_coad, file = "tcga_coad_data450.csv", row.names = FALSE)

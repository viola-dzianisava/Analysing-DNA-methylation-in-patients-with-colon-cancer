# Load the clinical data
clinical_data <- read.csv("/Users/henry/Desktop/Dissertation/Combine/clinical.csv")

# Check the number of samples in each group
table(clinical_data$sample_type)  # Replace sample_type with the relevant column indicating groups

# Load the methylation data
methylation_data <- read.csv("/Users/henry/Desktop/Dissertation/gene_methylation_data.csv", row.names = 1)

# Extract sample columns (assume columns starting with "TCGA" are sample columns)
sample_columns <- colnames(methylation_data)[grepl("^TCGA", colnames(methylation_data))]

# Number of samples
num_samples <- length(sample_columns)
cat("Number of samples:", num_samples, "\n")


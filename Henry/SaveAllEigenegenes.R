# Load necessary libraries
library(WGCNA)

# Function to print all eigengenes for each module
printAllEigengenes <- function(MEs, dataset_name) {
  cat("Eigengenes for each module in the", dataset_name, "dataset:\n")
  for (module in colnames(MEs)) {
    cat("\nModule:", module, "\n")
    print(MEs[, module])
  }
}

# Print eigengenes for tumor dataset
print("Tumor Dataset Eigengenes:")
printAllEigengenes(MEs_tumor, "tumor")

# Print eigengenes for normal dataset with top 5000 normal variance genes
print("\nNormal Dataset (Top 5000 Normal Variance Genes) Eigengenes:")
printAllEigengenes(MEs_normal, "normal_top5000_normal")

# Print eigengenes for normal dataset with top 5000 tumor variance genes
print("\nNormal Dataset (Top 5000 Tumor Variance Genes) Eigengenes:")
printAllEigengenes(MEs_normal_top5000_tumor, "normal_top5000_tumor")

# Print eigengenes for tumor dataset with top 5000 normal variance genes
print("\nTumor Dataset (Top 5000 Normal Variance Genes) Eigengenes:")
printAllEigengenes(MEs_tumor_top5000_normal, "tumor_top5000_normal")


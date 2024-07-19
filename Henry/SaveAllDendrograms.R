# Ensure necessary libraries are loaded
library(WGCNA)

# Function to save dendrograms
saveDendrograms <- function(blockwiseModulesResult, dataset_name, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save each dendrogram as a PNG file
  for (i in seq_along(blockwiseModulesResult$dendrograms)) {
    dendrogram <- blockwiseModulesResult$dendrograms[[i]]
    moduleColors <- labels2colors(blockwiseModulesResult$colors)
    
    dendrogram_file <- paste0(output_dir, "/", dataset_name, "_dendrogram_", i, ".png")
    png(dendrogram_file, width = 1200, height = 900)
    plotDendroAndColors(dendrogram, moduleColors[blockwiseModulesResult$blockGenes[[i]]], 
                        main = paste("Dendrogram and Module Colors -", dataset_name, "Block", i), 
                        dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
    dev.off()
  }
  
  cat("Dendrograms for", dataset_name, "saved in", output_dir, "\n")
}

# Save dendrograms for tumor dataset (Tumor Tumor)
output_dir_tt <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/TT/"
saveDendrograms(blockwiseModulesResult_tumor, "Tumor Tumor", output_dir_tt)

# Save dendrograms for normal dataset (Normal Normal)
output_dir_nn <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/NN/"
saveDendrograms(blockwiseModulesResult_normal, "Normal Normal", output_dir_nn)

# Save dendrograms for tumor dataset with top 5000 normal variance genes (Tumor Normal Variance)
output_dir_tn <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/TN/"
saveDendrograms(blockwiseModulesResult_tumor_top5000_normal, "Tumor Normal Variance", output_dir_tn)

# Save dendrograms for normal dataset with top 5000 tumor variance genes (Normal Tumor Variance)
output_dir_nt <- "C:/Users/Henry/Desktop/Dissertation/FinalFindings/NT/"
saveDendrograms(blockwiseModulesResult_normal_top5000_tumor, "Normal Tumor Variance", output_dir_nt)

# Load libraries
library(ggplot2)

# Set working directory
setwd("C:/Users/Administrator/Downloads/CMB_Project/CMB_copy")

# Load the prediction and ground truth CSVs
preds <- read.csv("./Results/testdata_with_predictions/preds.csv", row.names = 1, check.names = FALSE)
ki02 <- read.csv("./Data/KI_02_gene_expression_TPM_Ln_Norm.csv", row.names = 1, check.names = FALSE)

# Clean and standardize Patient IDs by removing "-1" and "-2" suffixes
rownames(preds) <- sub("-1$", "", trimws(rownames(preds)))
rownames(ki02)  <- sub("-2$", "", trimws(rownames(ki02)))

# Add PatientID column
preds$PatientID <- rownames(preds)
ki02$PatientID <- rownames(ki02)

# Merge datasets on PatientID
merged <- merge(preds, ki02, by = "PatientID")

# Debug: show how many patients matched
cat("Number of patients after merge:", nrow(merged), "\n")

# Set output directory for plots
out_dir <- "./Results/Correlation_Plots"
dir.create(out_dir, showWarnings = FALSE)

# Get gene names that are shared (predicted col ends with "_average_predicted")
pred_genes <- sub("_average_predicted$", "", grep("_average_predicted$", colnames(preds), value = TRUE))
shared_genes <- pred_genes[pred_genes %in% colnames(ki02)]

cat("Number of shared genes:", length(shared_genes), "\n")

# Loop over shared genes and generate scatter plots
for (gene in shared_genes) {
  pred_col <- paste0(gene, "_average_predicted")
  gt_col <- gene
  
  df <- merged[, c(pred_col, gt_col)]
  
  cat("\nProcessing gene:", gene, "\n")
  cat(" - Prediction column:", pred_col, "\n")
  cat(" - Ground truth column:", gt_col, "\n")
  cat(" - Data points:", nrow(df), "\n")
  print(head(df))
  
  if (nrow(df) >= 3) {
    cor_test <- cor.test(df[[1]], df[[2]], method = "spearman")
    
    p <- ggplot(df, aes(x = .data[[pred_col]], y = .data[[gt_col]])) +
      geom_point(color = "#2c3e50") +
      geom_smooth(method = "lm", se = FALSE, color = "#e74c3c") +
      labs(title = gene, x = "Prediction", y = "Ground Truth") +
      annotate("text", 
               x = min(df[[1]]), y = max(df[[2]]), 
               label = sprintf("Spearman r = %.3f\np = %.3g", cor_test$estimate, cor_test$p.value),
               hjust = 0, vjust = 1, size = 4.5, fontface = "bold") +
      theme_minimal()
    
    # Save plot
    ggsave(filename = paste0(out_dir, "/", gene, "_correlation_plot.png"),
           plot = p, width = 5, height = 4, dpi = 300)
    
    cat(" --> Plot saved for", gene, "\n")
  } else {
    cat(" --> Skipping", gene, "- not enough data points\n")
  }
}

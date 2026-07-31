###############################################################
# Project : TCGA Prostate Cancer Gene Expression Analysis
# Script  : 09_Heatmap.R
#
# Purpose:
# Generate heatmap of target gene expression
###############################################################

library(pheatmap)


#--------------------------------------------------------------
# Load data
#--------------------------------------------------------------

load("results/Processed_Expression_Data.RData")


#--------------------------------------------------------------
# Prepare matrix
#--------------------------------------------------------------

heatmap_data <- filtered_log2


# Scale genes

heatmap_scaled <- t(scale(t(heatmap_data)))


#--------------------------------------------------------------
# Generate heatmap
#--------------------------------------------------------------

png(
  "figures/Gene_Expression_Heatmap.png",
  width = 2000,
  height = 1500,
  res = 300
)


pheatmap(
  heatmap_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 12,
  border_color = NA,
  main = NULL
)


dev.off()


cat("\nHeatmap generated successfully.\n")

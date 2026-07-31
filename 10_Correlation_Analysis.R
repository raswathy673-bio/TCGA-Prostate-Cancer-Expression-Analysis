###############################################################
# Project : TCGA Prostate Cancer Gene Expression Analysis
# Script  : 10_Correlation_Analysis.R
#
# Purpose:
# Perform correlation analysis among target genes
###############################################################

library(corrplot)


#--------------------------------------------------------------
# Load data
#--------------------------------------------------------------

load("results/Processed_Expression_Data.RData")


#--------------------------------------------------------------
# Calculate Pearson correlation
#--------------------------------------------------------------

cor_matrix <- cor(
  t(filtered_log2),
  method = "pearson"
)


print(cor_matrix)


#--------------------------------------------------------------
# Save correlation values
#--------------------------------------------------------------

write.csv(
  cor_matrix,
  "results/Gene_Correlation_Matrix.csv"
)


#--------------------------------------------------------------
# Plot correlation heatmap
#--------------------------------------------------------------

png(
  "figures/Gene_Correlation_Heatmap.png",
  width = 2000,
  height = 1800,
  res = 300
)


corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  addCoef.col = "black",
  number.cex = 0.8,
  tl.col = "black",
  tl.srt = 45
)


dev.off()


cat("\nCorrelation analysis completed successfully.\n")

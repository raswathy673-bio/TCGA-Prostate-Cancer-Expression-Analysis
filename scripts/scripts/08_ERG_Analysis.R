###############################################################
# Project : TCGA Prostate Cancer Gene Expression Analysis
# Script  : 08_ERG_Analysis.R
#
# Purpose:
# Analyze association between gene expression and
# TMPRSS2-ERG fusion status
###############################################################

library(dplyr)
library(tidyr)
library(tibble)

source("scripts/functions.R")


#--------------------------------------------------------------
# Load data
#--------------------------------------------------------------

load("results/Processed_Expression_Data.RData")


#--------------------------------------------------------------
# Prepare dataframe
#--------------------------------------------------------------

expr_df <- as.data.frame(t(filtered_log2))


expr_df$ERG_status <- clinical_data$paper_ERG_status


#--------------------------------------------------------------
# Keep only fusion and none groups
#--------------------------------------------------------------

expr_df <- expr_df %>%
  filter(
    ERG_status %in% c("fusion", "none")
  )


expr_df$ERG_status <- factor(
  expr_df$ERG_status,
  levels = c("none", "fusion")
)


table(expr_df$ERG_status)


#--------------------------------------------------------------
# Long format
#--------------------------------------------------------------

expr_long <- expr_df %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    cols = c(KLK2, KLK3, TMPRSS2, AKR1C3),
    names_to = "Gene",
    values_to = "Expression"
  )


#--------------------------------------------------------------
# Wilcoxon test
#--------------------------------------------------------------

stats <- run_wilcox(
  expr_long,
  "ERG_status"
)


print(stats)


write.csv(
  stats,
  "results/ERG_Statistics.csv",
  row.names = FALSE
)


#--------------------------------------------------------------
# Plot
#--------------------------------------------------------------

p <- plot_violin(
  expr_long,
  group = "ERG_status",
  method = "wilcox.test",
  xlab = "ERG Fusion Status"
)


print(p)


save_plot(
  p,
  "ERG_Status"
)


cat("\nERG analysis completed successfully.\n")

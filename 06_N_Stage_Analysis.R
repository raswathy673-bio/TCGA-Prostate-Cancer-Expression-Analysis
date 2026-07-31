###############################################################
# Project : TCGA Prostate Cancer Gene Expression Analysis
# Script  : 06_N_Stage_Analysis.R
# Purpose : Association of gene expression with lymph node stage
#           (N0 vs N1)
###############################################################

#--------------------------------------------------------------
# Load libraries
#--------------------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)

source("scripts/functions.R")

#--------------------------------------------------------------
# Load processed data
#--------------------------------------------------------------

load("results/Processed_Expression_Data.RData")

#--------------------------------------------------------------
# Prepare expression dataframe
#--------------------------------------------------------------

expr_df <- as.data.frame(t(filtered_log2))

# Add N stage information
expr_df$N_stage <- clinical_data$ajcc_pathologic_n

#--------------------------------------------------------------
# Remove missing and NX samples
#--------------------------------------------------------------

expr_df <- expr_df %>%
  filter(
    !is.na(N_stage),
    N_stage %in% c("N0", "N1")
  )

#--------------------------------------------------------------
# Convert N stage as factor
#--------------------------------------------------------------

expr_df$N_stage <- factor(
  expr_df$N_stage,
  levels = c("N0", "N1")
)

# Check sample numbers
table(expr_df$N_stage)

#--------------------------------------------------------------
# Convert expression data into long format
#--------------------------------------------------------------

expr_long <- expr_df %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    cols = c(KLK2, KLK3, TMPRSS2, AKR1C3),
    names_to = "Gene",
    values_to = "Expression"
  )

#--------------------------------------------------------------
# Statistical analysis
#--------------------------------------------------------------

stats <- run_wilcox(
  expr_long,
  "N_stage"
)

print(stats)

write.csv(
  stats,
  "results/N_stage_Statistics.csv",
  row.names = FALSE
)

#--------------------------------------------------------------
# Generate violin plot
#--------------------------------------------------------------

p <- plot_violin(
  data = expr_long,
  group = "N_stage",
  method = "wilcox.test",
  xlab = "Lymph Node Stage"
)

print(p)

#--------------------------------------------------------------
# Save figures
#--------------------------------------------------------------

save_plot(
  p,
  "N_stage"
)

cat("\nN stage analysis completed successfully.\n")

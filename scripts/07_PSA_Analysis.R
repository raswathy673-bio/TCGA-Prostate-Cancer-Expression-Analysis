###############################################################
# Project : TCGA Prostate Cancer Gene Expression Analysis
# Script  : 07_PSA_Analysis.R
#
# Purpose:
# Analyze association between gene expression and
# pre-operative PSA levels
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


# Add PSA information

expr_df$PSA <- clinical_data$paper_PSA_preop


#--------------------------------------------------------------
# Convert PSA character values to numeric
#--------------------------------------------------------------

expr_df$PSA <- gsub(",", ".", expr_df$PSA)

expr_df$PSA <- as.numeric(expr_df$PSA)


#--------------------------------------------------------------
# Remove missing PSA values
#--------------------------------------------------------------

expr_df <- expr_df %>%
  filter(!is.na(PSA))


#--------------------------------------------------------------
# Create PSA groups
#--------------------------------------------------------------

expr_df$PSA_group <- case_when(

  expr_df$PSA < 10 ~ "Low (<10)",

  expr_df$PSA >= 10 &
    expr_df$PSA <= 20 ~ "Intermediate (10-20)",

  expr_df$PSA > 20 ~ "High (>20)"

)


expr_df$PSA_group <- factor(
  expr_df$PSA_group,
  levels = c(
    "Low (<10)",
    "Intermediate (10-20)",
    "High (>20)"
  )
)


# Check sample distribution

table(expr_df$PSA_group)


#--------------------------------------------------------------
# Convert to long format
#--------------------------------------------------------------

expr_long <- expr_df %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    cols = c(KLK2, KLK3, TMPRSS2, AKR1C3),
    names_to = "Gene",
    values_to = "Expression"
  )


#--------------------------------------------------------------
# Kruskal-Wallis test
#--------------------------------------------------------------

stats <- run_kruskal(
  expr_long,
  "PSA_group"
)


print(stats)


write.csv(
  stats,
  "results/PSA_Statistics.csv",
  row.names = FALSE
)


#--------------------------------------------------------------
# Generate violin plot
#--------------------------------------------------------------

p <- plot_violin(
  data = expr_long,
  group = "PSA_group",
  method = "kruskal.test",
  xlab = "Pre-operative PSA Level (ng/mL)"
)


print(p)


#--------------------------------------------------------------
# Save figures
#--------------------------------------------------------------

save_plot(
  p,
  "PSA"
)


cat("\nPSA analysis completed successfully.\n")

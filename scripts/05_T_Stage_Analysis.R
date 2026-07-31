###############################################################
# Project : TCGA Prostate Cancer Gene Expression Analysis
# Script  : 05_T_Stage_Analysis.R
# Purpose : Association of gene expression with pathological
#           T stage
###############################################################

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

expr_df$T_stage <- clinical_data$ajcc_pathologic_t

expr_df <- expr_df %>%
  filter(!is.na(T_stage))

#--------------------------------------------------------------
# Merge stages
#--------------------------------------------------------------

expr_df$T_stage <- as.character(expr_df$T_stage)

expr_df$T_stage[expr_df$T_stage %in%
                  c("T2","T2a","T2b","T2c")] <- "T2"

expr_df$T_stage[expr_df$T_stage %in%
                  c("T3","T3a","T3b")] <- "T3"

expr_df$T_stage <- factor(
  expr_df$T_stage,
  levels = c("T2","T3","T4")
)

expr_df <- expr_df %>%
  filter(T_stage %in% c("T2","T3","T4"))

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
# Statistics
#--------------------------------------------------------------

stats <- run_kruskal(
  expr_long,
  "T_stage"
)

print(stats)

write.csv(
  stats,
  "results/T_stage_Statistics.csv",
  row.names = FALSE
)

#--------------------------------------------------------------
# Plot
#--------------------------------------------------------------

p <- plot_violin(
  data = expr_long,
  group = "T_stage",
  method = "kruskal.test",
  xlab = "Pathological T Stage"
)

print(p)

save_plot(
  p,
  "T_stage"
)

cat("\nT stage analysis completed successfully.\n")

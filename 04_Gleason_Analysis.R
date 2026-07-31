###############################################################
# Project : TCGA Prostate Cancer Gene Expression Analysis
# Script  : 04_Gleason_Analysis.R
# Author  : Your Name
# Date    : August 2026
#
# Purpose:
# Compare expression of KLK2, KLK3, TMPRSS2 and AKR1C3
# across Gleason score groups.
###############################################################

#--------------------------------------------------------------
# Load libraries
#--------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tibble)

#--------------------------------------------------------------
# Load processed data
#--------------------------------------------------------------

load("results/Processed_Expression_Data.RData")

#--------------------------------------------------------------
# Prepare expression dataframe
#--------------------------------------------------------------

expr_df <- as.data.frame(t(filtered_log2))

expr_df$Gleason <- clinical_data$gleason_score

#--------------------------------------------------------------
# Remove missing values
#--------------------------------------------------------------

expr_df <- expr_df %>%
  filter(!is.na(Gleason))

#--------------------------------------------------------------
# Merge Gleason groups
#--------------------------------------------------------------

expr_df$Gleason <- as.character(expr_df$Gleason)

expr_df$Gleason[expr_df$Gleason == "10"] <- "9"

expr_df$Gleason <- factor(
  expr_df$Gleason,
  levels = c("6","7","8","9"),
  labels = c("6","7","8","9-10")
)

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

stats <- expr_long %>%
  group_by(Gene) %>%
  summarise(
    p_value = kruskal.test(Expression ~ Gleason)$p.value
  )

print(stats)

write.csv(
  stats,
  "results/Gleason_Statistics.csv",
  row.names = FALSE
)

#--------------------------------------------------------------
# Publication-quality violin plot
#--------------------------------------------------------------

p <- ggplot(
  expr_long,
  aes(Gleason,
      Expression,
      fill = Gleason)
) +

  geom_violin(
    trim = FALSE,
    alpha = 0.8,
    color = "black"
  ) +

  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.shape = NA
  ) +

  stat_compare_means(
    method = "kruskal.test",
    label.y = max(expr_long$Expression) + 0.5
  ) +

  facet_wrap(
    ~Gene,
    ncol = 2,
    scales = "free_y"
  ) +

  labs(
    x = "Gleason Score",
    y = expression(log[2](TPM + 1))
  ) +

  theme_classic(base_size = 15) +

  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

print(p)

#--------------------------------------------------------------
# Save figures
#--------------------------------------------------------------

ggsave(
  "figures/Gleason.png",
  p,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  "figures/Gleason.pdf",
  p,
  width = 10,
  height = 8
)

ggsave(
  "figures/Gleason.tiff",
  p,
  width = 10,
  height = 8,
  dpi = 600,
  compression = "lzw"
)

cat("\nGleason analysis completed successfully.\n")

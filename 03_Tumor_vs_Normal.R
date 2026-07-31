###############################################################
# Project : TCGA Prostate Cancer Gene Expression Analysis
# Script  : 03_Tumor_vs_Normal.R
# Author  : Your Name
# Date    : August 2026
#
# Purpose:
# Compare expression of KLK2, KLK3, TMPRSS2 and AKR1C3
# between normal and tumor tissues.
###############################################################

#--------------------------------------------------------------
# 1. Load required libraries
#--------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tibble)

#--------------------------------------------------------------
# 2. Load processed data
#--------------------------------------------------------------

load("results/Processed_Expression_Data.RData")

#--------------------------------------------------------------
# 3. Extract sample groups
#--------------------------------------------------------------

sample_type <- data_se$shortLetterCode

group <- ifelse(sample_type == "NT",
                "Normal",
                "Tumor")

group <- factor(group,
                levels = c("Normal", "Tumor"))

#--------------------------------------------------------------
# 4. Prepare expression dataframe
#--------------------------------------------------------------

expr_df <- as.data.frame(t(filtered_log2))

expr_df$Group <- group

#--------------------------------------------------------------
# 5. Convert to long format
#--------------------------------------------------------------

expr_long <- expr_df %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    cols = c(KLK2, KLK3, TMPRSS2, AKR1C3),
    names_to = "Gene",
    values_to = "Expression"
  )

#--------------------------------------------------------------
# 6. Wilcoxon test
#--------------------------------------------------------------

stats <- expr_long %>%
  group_by(Gene) %>%
  summarise(
    p_value = wilcox.test(Expression ~ Group)$p.value
  )

print(stats)

write.csv(
  stats,
  "results/Tumor_vs_Normal_Statistics.csv",
  row.names = FALSE
)

#--------------------------------------------------------------
# 7. Publication-quality violin plot
#--------------------------------------------------------------

p <- ggplot(
  expr_long,
  aes(Group,
      Expression,
      fill = Group)
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
    method = "wilcox.test",
    label = "p.format"
  ) +

  facet_wrap(
    ~Gene,
    ncol = 2,
    scales = "free_y"
  ) +

  scale_fill_manual(
    values = c("#4DBBD5", "#E64B35")
  ) +

  labs(
    x = "",
    y = expression(log[2](TPM + 1))
  ) +

  theme_bw(base_size = 15) +

  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  )

print(p)

#--------------------------------------------------------------
# 8. Save figures
#--------------------------------------------------------------

ggsave(
  "figures/Tumor_vs_Normal.png",
  p,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  "figures/Tumor_vs_Normal.pdf",
  p,
  width = 10,
  height = 8
)

ggsave(
  "figures/Tumor_vs_Normal.tiff",
  p,
  width = 10,
  height = 8,
  dpi = 600,
  compression = "lzw"
)

cat("\nTumor vs Normal analysis completed successfully.\n")

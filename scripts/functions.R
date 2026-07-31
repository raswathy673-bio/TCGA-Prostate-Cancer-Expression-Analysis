###############################################################
# Project : TCGA Prostate Cancer Gene Expression Analysis
# File    : functions.R
# Purpose : Reusable functions for statistical analysis and
#           publication-quality visualization
###############################################################

library(ggplot2)
library(ggpubr)
library(dplyr)

#--------------------------------------------------------------
# Kruskal-Wallis Test
#--------------------------------------------------------------

run_kruskal <- function(data, group){

  data %>%
    group_by(Gene) %>%
    summarise(
      p_value = kruskal.test(
        as.formula(
          paste("Expression ~", group)
        )
      )$p.value
    )

}

#--------------------------------------------------------------
# Wilcoxon Test
#--------------------------------------------------------------

run_wilcox <- function(data, group){

  data %>%
    group_by(Gene) %>%
    summarise(
      p_value = wilcox.test(
        as.formula(
          paste("Expression ~", group)
        )
      )$p.value
    )

}

#--------------------------------------------------------------
# Violin Plot
#--------------------------------------------------------------

plot_violin <- function(data,
                        group,
                        method,
                        xlab){

  ggplot(
    data,
    aes_string(
      x = group,
      y = "Expression",
      fill = group
    )
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
      method = method,
      label = "p.format"
    ) +

    facet_wrap(
      ~Gene,
      scales = "free_y",
      ncol = 2
    ) +

    labs(
      x = xlab,
      y = expression(log[2](TPM+1))
    ) +

    theme_classic(base_size = 15) +

    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    )

}

#--------------------------------------------------------------
# Save Plot
#--------------------------------------------------------------

save_plot <- function(plot, filename){

  ggsave(
    paste0("figures/", filename, ".png"),
    plot,
    width = 10,
    height = 8,
    dpi = 300
  )

  ggsave(
    paste0("figures/", filename, ".pdf"),
    plot,
    width = 10,
    height = 8
  )

  ggsave(
    paste0("figures/", filename, ".tiff"),
    plot,
    width = 10,
    height = 8,
    dpi = 600,
    compression = "lzw"
  )

}

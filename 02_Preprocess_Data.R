###############################################################
# Project : TCGA Prostate Cancer Gene Expression Analysis
# Script  : 02_Preprocess_Data.R
# Author  : Your Name
# Date    : August 2026
#
# Purpose:
# Extract TPM expression matrix, map gene symbols,
# prepare expression matrix for downstream analysis.
###############################################################

#--------------------------------------------------------------
# 1. Load required packages
#--------------------------------------------------------------

library(SummarizedExperiment)

#--------------------------------------------------------------
# 2. Load downloaded TCGA object
#--------------------------------------------------------------

load("data/TCGA_PRAD_SummarizedExperiment.RData")

#--------------------------------------------------------------
# 3. Display available assays
#--------------------------------------------------------------

assayNames(data_se)

#--------------------------------------------------------------
# 4. Extract TPM matrix
#--------------------------------------------------------------

tpm_matrix <- assay(data_se, "tpm_unstrand")

#--------------------------------------------------------------
# 5. Extract gene annotation
#--------------------------------------------------------------

gene_info <- as.data.frame(rowData(data_se))

head(gene_info)

#--------------------------------------------------------------
# 6. Replace Ensembl IDs with gene symbols
#--------------------------------------------------------------

rownames(tpm_matrix) <- gene_info$gene_name

# Remove genes without names
tpm_matrix <- tpm_matrix[rownames(tpm_matrix) != "", ]

#--------------------------------------------------------------
# 7. Remove duplicated gene symbols
#--------------------------------------------------------------

tpm_matrix <- tpm_matrix[
  !duplicated(rownames(tpm_matrix)),
]

#--------------------------------------------------------------
# 8. Target genes
#--------------------------------------------------------------

target_genes <- c(
  "KLK2",
  "KLK3",
  "TMPRSS2",
  "AKR1C3"
)

filtered_tpm <- tpm_matrix[
  rownames(tpm_matrix) %in% target_genes,
]

#--------------------------------------------------------------
# 9. Log2 transformation
#--------------------------------------------------------------

filtered_log2 <- log2(filtered_tpm + 1)

#--------------------------------------------------------------
# 10. Check extracted genes
#--------------------------------------------------------------

cat("Genes extracted:\n")
print(rownames(filtered_log2))

cat("\nExpression matrix dimensions:\n")
print(dim(filtered_log2))

#--------------------------------------------------------------
# 11. Save processed data
#--------------------------------------------------------------

#--------------------------------------------------------------
# 11. Save processed data
#--------------------------------------------------------------

clinical_data <- as.data.frame(colData(data_se))

save(
  filtered_log2,
  gene_info,
  clinical_data,
  data_se,
  file = "results/Processed_Expression_Data.RData"
)

cat("\nPreprocessing completed successfully.\n")
cat("\nPreprocessing completed successfully.\n")

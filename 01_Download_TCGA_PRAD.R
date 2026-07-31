###############################################################
# Project : TCGA Prostate Cancer Gene Expression Analysis
# Script  : 01_Download_TCGA_PRAD.R
# Author  : Your Name
# Date    : August 2026
#
# Purpose:
# Download TCGA-PRAD RNA-seq (TPM) data and associated clinical
# information using TCGAbiolinks.
###############################################################

#--------------------------------------------------------------
# 1. Load required packages
#--------------------------------------------------------------

required_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment"
)

for(pkg in required_packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    BiocManager::install(pkg, ask = FALSE)
  }
}

library(TCGAbiolinks)
library(SummarizedExperiment)

#--------------------------------------------------------------
# 2. Create project folders
#--------------------------------------------------------------

dirs <- c("data", "results", "figures")

for(i in dirs){
  if(!dir.exists(i)){
    dir.create(i)
  }
}

#--------------------------------------------------------------
# 3. Query TCGA-PRAD RNA-seq TPM data
#--------------------------------------------------------------

query <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - TPM"
)

#--------------------------------------------------------------
# 4. Download data
#--------------------------------------------------------------

GDCdownload(query)

#--------------------------------------------------------------
# 5. Prepare SummarizedExperiment object
#--------------------------------------------------------------

data_se <- GDCprepare(query)

#--------------------------------------------------------------
# 6. Display dataset summary
#--------------------------------------------------------------

cat("Project:", unique(data_se$project_id), "\n")
cat("Total Samples:", ncol(data_se), "\n")
cat("Total Genes:", nrow(data_se), "\n")

cat("\nSample Types\n")
print(table(data_se$shortLetterCode))

#--------------------------------------------------------------
# 7. Save object
#--------------------------------------------------------------

save(
  data_se,
  file = "data/TCGA_PRAD_SummarizedExperiment.RData"
)

cat("\nDownload completed successfully.\n")

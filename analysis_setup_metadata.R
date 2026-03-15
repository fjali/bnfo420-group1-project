# ============================================================
# BNFO 420 Capstone
# Script 01: Metadata and Data Organization
# Goal: Load expression data, organize metadata, and separate by tissue
# ============================================================

# ----------------------------
# 1. Install packages if needed
# ----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma")
}

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

# ----------------------------
# 2. Load libraries
# ----------------------------
library(limma)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(GEOquery)

# ----------------------------
# 3. Create output folders
# ----------------------------
date_stamp <- format(Sys.Date(), "%Y-%m-%d")
figures_dir <- paste0("outputs/", date_stamp, "/figures")
data_dir <- paste0("outputs/", date_stamp, "/data")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

cat("Outputs will be saved to: outputs/", date_stamp, "\n")

# ----------------------------
# 4. Load expression data
# ----------------------------
# capstonedata.txt.gz is a compressed tab-delimited text file
# row.names = 1 assumes the first column contains gene IDs
expression_data <- read.delim("capstonedata.txt.gz", row.names = 1)

# Check that the data loaded correctly
dim(expression_data)
head(expression_data[, 1:5])

# Save a copy of the raw expression data
write.csv(expression_data, file.path(data_dir, "expression_data_raw.csv"))

# ----------------------------
# 5. Build metadata table
# ----------------------------
metadata <- data.frame(
  gsm_id = c(
    # Glp1r fl/fl - Vehicle - Liver (5 samples)
    "GSM8879750","GSM8879751","GSM8879752","GSM8879753","GSM8879754",
    # Glp1r fl/fl - Vehicle - iWAT (5 samples)
    "GSM8879765","GSM8879766","GSM8879767","GSM8879768","GSM8879769",
    # Glp1r fl/fl - Drug - Liver (5 samples)
    "GSM8879760","GSM8879761","GSM8879762","GSM8879763","GSM8879764",
    # Glp1r fl/fl - Drug - iWAT (5 samples)
    "GSM8879775","GSM8879776","GSM8879777","GSM8879778","GSM8879779",
    # Glp1r Wnt1-/- - Vehicle - Liver (5 samples)
    "GSM8879755","GSM8879756","GSM8879757","GSM8879758","GSM8879759",
    # Glp1r Wnt1-/- - Vehicle - iWAT (5 samples)
    "GSM8879770","GSM8879771","GSM8879772","GSM8879773","GSM8879774",
    # Gipr fl/fl - Vehicle - Liver (5 samples)
    "GSM8879780","GSM8879782","GSM8879786","GSM8879791","GSM8879806",
    # Gipr fl/fl - Vehicle - iWAT (5 samples)
    "GSM8879792","GSM8879794","GSM8879796","GSM8879801","GSM8879804",
    # Gipr fl/fl - Drug - Liver (5 samples)
    "GSM8879784","GSM8879787","GSM8879789","GSM8879805","GSM8879807",
    # Gipr fl/fl - Drug - iWAT (5 samples)
    "GSM8879795","GSM8879797","GSM8879799","GSM8879802","GSM8879808",
    # Gipr Syn-/- - Vehicle - Liver (5 samples)
    "GSM8879781","GSM8879783","GSM8879785","GSM8879788","GSM8879790",
    # Gipr Syn-/- - Vehicle - iWAT (5 samples)
    "GSM8879793","GSM8879798","GSM8879800","GSM8879803","GSM8879809"
  ),
  genotype = c(
    rep("Glp1r_fl_fl", 20),
    rep("Glp1r_Wnt1KO", 10),
    rep("Gipr_fl_fl", 20),
    rep("Gipr_SynKO", 10)
  ),
  treatment = c(
    rep("Vehicle", 5),
    rep("Vehicle", 5),
    rep("Drug", 5),
    rep("Drug", 5),
    rep("Vehicle", 10),
    rep("Vehicle", 5),
    rep("Vehicle", 5),
    rep("Drug", 5),
    rep("Drug", 5),
    rep("Vehicle", 10)
  ),
  tissue = c(
    rep("Liver", 5), rep("iWAT", 5),
    rep("Liver", 5), rep("iWAT", 5),
    rep("Liver", 5), rep("iWAT", 5),
    rep("Liver", 5), rep("iWAT", 5),
    rep("Liver", 5), rep("iWAT", 5),
    rep("Liver", 5), rep("iWAT", 5)
  ),
  replicate = rep(1:5, 12),
  stringsAsFactors = FALSE
)

# Check metadata
dim(metadata)
head(metadata)

# ----------------------------
# 6. Make sure sample names match metadata
# ----------------------------
all(metadata$gsm_id %in% colnames(expression_data))
setdiff(metadata$gsm_id, colnames(expression_data))

# Reorder metadata to match expression data columns
metadata <- metadata %>%
  filter(gsm_id %in% colnames(expression_data)) %>%
  arrange(match(gsm_id, colnames(expression_data)))

# Confirm order matches
identical(metadata$gsm_id, colnames(expression_data))

# ----------------------------
# 7. Separate metadata by tissue
# ----------------------------
metadata_liver <- metadata %>% filter(tissue == "Liver")
metadata_iWAT  <- metadata %>% filter(tissue == "iWAT")

# ----------------------------
# 8. Separate expression data by tissue
# ----------------------------
expression_liver <- expression_data[, metadata_liver$gsm_id]
expression_iWAT  <- expression_data[, metadata_iWAT$gsm_id]

# Check dimensions
dim(expression_liver)
dim(expression_iWAT)

# ----------------------------
# 9. Save metadata and tissue-specific datasets
# ----------------------------
write.csv(metadata, file.path(data_dir, "metadata_all.csv"), row.names = FALSE)
write.csv(metadata_liver, file.path(data_dir, "metadata_liver.csv"), row.names = FALSE)
write.csv(metadata_iWAT, file.path(data_dir, "metadata_iWAT.csv"), row.names = FALSE)

write.csv(expression_liver, file.path(data_dir, "expression_liver.csv"))
write.csv(expression_iWAT, file.path(data_dir, "expression_iWAT.csv"))

cat("Metadata and tissue-specific expression files saved successfully.\n")

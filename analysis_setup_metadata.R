# ============================================================
# BNFO 420 Capstone
# Script 01: Metadata and Data Organization
# Goal: Load expression data, organize metadata, and separate by tissue
# ============================================================

#install packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("limma",       quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("tidyverse",   quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("pheatmap",    quietly = TRUE)) install.packages("pheatmap")

#load libraries
library(limma)
library(tidyverse)
library(ggplot2)
library(pheatmap)

#create output folders
date_stamp  <- format(Sys.Date(), "%Y-%m-%d")
figures_dir <- paste0("outputs/", date_stamp, "/figures")
data_dir    <- paste0("outputs/", date_stamp, "/data")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,    recursive = TRUE, showWarnings = FALSE)

cat("Outputs will be saved to: outputs/", date_stamp, "\n")

#build metadata table
# IMPORTANT: This must come BEFORE loading expression data
# so we can use it to rename columns

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
    rep("Glp1r_fl_fl",  20),
    rep("Glp1r_Wnt1KO", 10),
    rep("Gipr_fl_fl",   20),
    rep("Gipr_SynKO",   10)
  ),
  treatment = c(
    rep("Vehicle", 5),
    rep("Vehicle", 5),
    rep("Drug",    5),
    rep("Drug",    5),
    rep("Vehicle", 10),
    rep("Vehicle", 5),
    rep("Vehicle", 5),
    rep("Drug",    5),
    rep("Drug",    5),
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

cat("Metadata built:", nrow(metadata), "rows\n")
head(metadata)

# load all 4 expression files
# Each file = one tissue + one receptor group = 15 samples
# Total: 4 files x 15 samples = 60 samples

load_expr <- function(filename) {
  df <- read.delim(filename,
                   header      = TRUE,
                   check.names = FALSE,
                   sep         = "\t")
  rownames(df) <- df[, 1]
  df <- df[, -1]
  return(df)
}

expr_liver_glp1r <- load_expr("/Users/kbandukwala/bnfo420/GSE293291_FPKM_Liver_GLP1R.txt")
expr_liver_gipr  <- load_expr("/Users/kbandukwala/bnfo420/GSE293291_FPKM_Liver_GIPR.txt")
expr_iwat_glp1r  <- load_expr("/Users/kbandukwala/bnfo420/GSE293291_FPKM_ING_GLP1R.txt")
expr_iwat_gipr   <- load_expr("/Users/kbandukwala/bnfo420/GSE293291_FPKM_ING_GIPR.txt")

cat("Liver GLP1R:", ncol(expr_liver_glp1r), "samples\n")
cat("Liver GIPR: ", ncol(expr_liver_gipr),  "samples\n")
cat("iWAT GLP1R: ", ncol(expr_iwat_glp1r),  "samples\n")
cat("iWAT GIPR:  ", ncol(expr_iwat_gipr),   "samples\n")

#combine all 4 files
# Order must match metadata order:
# Liver GLP1R + iWAT GLP1R + Liver GLP1R Drug + iWAT GLP1R Drug +
# Liver Wnt1KO + iWAT Wnt1KO + Liver GIPR + iWAT GIPR +
# Liver GIPR Drug + iWAT GIPR Drug + Liver SynKO + iWAT SynKO
expression_data <- cbind(expr_liver_glp1r, expr_liver_gipr,
                         expr_iwat_glp1r,  expr_iwat_gipr)

cat("Total samples after combining:", ncol(expression_data), "\n")

#rename columns to the gsm ids
# The numeric column names need to become GSM IDs
# We assign them in the same order as our metadata table
colnames(expression_data) <- metadata$gsm_id

# Confirm renaming worked
cat("First few column names after renaming:\n")
print(head(colnames(expression_data)))
cat("Should show GSM IDs like GSM8879750, GSM8879751...\n")

#separate by tissue
metadata_liver <- metadata %>% filter(tissue == "Liver")
metadata_iWAT  <- metadata %>% filter(tissue == "iWAT")

expression_liver <- expression_data[, metadata_liver$gsm_id]
expression_iWAT  <- expression_data[, metadata_iWAT$gsm_id]

cat("Liver dataset:", nrow(expression_liver), "genes x", ncol(expression_liver), "samples\n")
cat("iWAT dataset: ", nrow(expression_iWAT),  "genes x", ncol(expression_iWAT),  "samples\n")

# Should say 30 samples each

#save everything
write.csv(metadata,       file.path(data_dir, "metadata_all.csv"),    row.names = FALSE)
write.csv(metadata_liver, file.path(data_dir, "metadata_liver.csv"),  row.names = FALSE)
write.csv(metadata_iWAT,  file.path(data_dir, "metadata_iWAT.csv"),   row.names = FALSE)

write.table(expression_liver, file.path(data_dir, "expression_liver.csv"),
            sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(expression_iWAT,  file.path(data_dir, "expression_iWAT.csv"),
            sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

cat("All files saved successfully.\n")
cat("Script 01 complete. Run Script 02 next.\n")

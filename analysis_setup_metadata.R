# ============================================================
# BNFO 420 Capstone
# Script 02: Normalize and Filter Expression Data
# Goal: Log2 transform, filter low-expression genes,
#       and produce a gene filtering summary figure
# run after the first script - metadata organization

#load libraries
library(limma)
library(tidyverse)
library(ggplot2)

#2. recreate output folders
date_stamp  <- format(Sys.Date(), "%Y-%m-%d")
figures_dir <- paste0("outputs/", date_stamp, "/figures")
data_dir    <- paste0("outputs/", date_stamp, "/data")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,    recursive = TRUE, showWarnings = FALSE)

# 3. load the data saved by script 1
expression_liver <- read.csv(file.path(data_dir, "expression_liver.csv"), 
                             row.names   = 1,
                             check.names = FALSE)
expression_iWAT  <- read.csv(file.path(data_dir, "expression_iWAT.csv"),  
                             row.names   = 1,
                             check.names = FALSE)
metadata_liver   <- read.csv(file.path(data_dir, "metadata_liver.csv"),
                             check.names = FALSE)
metadata_iWAT    <- read.csv(file.path(data_dir, "metadata_iWAT.csv"),
                             check.names = FALSE)

# Confirm they loaded correctly
cat("Liver:", nrow(expression_liver), "genes x", ncol(expression_liver), "samples\n")
cat("iWAT: ", nrow(expression_iWAT),  "genes x", ncol(expression_iWAT),  "samples\n")

# 4. log2 transform
# We use log2(value + 1) to:
# - Bring large values down to a workable scale
# - The +1 prevents log(0) which is undefined
# - Makes the data more normally distributed for statistics

expression_liver_log <- log2(expression_liver + 1)
expression_iWAT_log  <- log2(expression_iWAT  + 1)

cat("Log2 transformation complete.\n")

# 5. filter low expression genes
# We remove genes that are barely expressed across all samples
# because they are just noise and can mess up our statistics
#
# Rule: keep a gene only if it has expression > 1 (log2 scale)
#       in at least 5 samples (= one full biological replicate group)

min_expression <- 1   # log2(FPKM + 1) must be above this
min_samples    <- 5   # in at least this many samples

# Record gene counts before filtering
genes_before_liver <- nrow(expression_liver_log)
genes_before_iWAT  <- nrow(expression_iWAT_log)

# Apply filter to Liver
keep_liver           <- rowSums(expression_liver_log > min_expression) >= min_samples
expression_liver_log <- expression_liver_log[keep_liver, ]

# Apply filter to iWAT
keep_iWAT           <- rowSums(expression_iWAT_log > min_expression) >= min_samples
expression_iWAT_log <- expression_iWAT_log[keep_iWAT, ]

cat("--- Liver filtering ---\n")
cat("Genes before:", genes_before_liver, "\n")
cat("Genes after: ", nrow(expression_liver_log), "\n")
cat("Genes removed:", genes_before_liver - nrow(expression_liver_log), "\n\n")

cat("--- iWAT filtering ---\n")
cat("Genes before:", genes_before_iWAT, "\n")
cat("Genes after: ", nrow(expression_iWAT_log), "\n")
cat("Genes removed:", genes_before_iWAT - nrow(expression_iWAT_log), "\n")

#6. create gene filtering bar chart
filter_summary <- data.frame(
  tissue = c("Liver", "Liver", "iWAT", "iWAT"),
  step   = c("Before Filtering", "After Filtering",
             "Before Filtering", "After Filtering"),
  genes  = c(
    genes_before_liver, nrow(expression_liver_log),
    genes_before_iWAT,  nrow(expression_iWAT_log)
  )
)

ggplot(filter_summary, aes(x = step, y = genes, fill = step)) +
  geom_bar(stat = "identity", width = 0.5, color = "white") +
  geom_text(aes(label = format(genes, big.mark = ",")), 
            vjust = -0.6, size = 4, fontface = "bold", color = "gray30") +
  scale_fill_manual(values = c(
    "Before Filtering" = "#f4a7b9",   # light pink
    "After Filtering"  = "#a8d8ea"    # light blue
  )) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.15))
  ) +
  facet_wrap(~ tissue, labeller = labeller(tissue = c(
    "Liver" = "Liver Tissue",
    "iWAT"  = "Inguinal White Adipose Tissue (iWAT)"
  ))) +
  labs(
    title    = "Gene Expression Filtering Summary",
    subtitle = paste("Genes retained after removing low-expression genes\n",
                     "Threshold: log2(FPKM+1) >", min_expression,
                     "in at least", min_samples, "samples"),
    x        = "Filtering Step",
    y        = "Number of Genes",
    fill     = "Filtering Step"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title.x    = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y    = element_text(face = "bold", margin = margin(r = 10)),
    strip.text      = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = file.path(figures_dir, "figure1_gene_filtering.png"),
  width = 8, height = 6, dpi = 300
)
cat("Figure 1 saved: gene filtering bar chart\n")

# 7. save filtered gene data
write.table(expression_liver_log, file.path(data_dir, "expression_liver_filtered.csv"),
            sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(expression_iWAT_log,  file.path(data_dir, "expression_iWAT_filtered.csv"),
            sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

cat("Filtered expression data saved successfully.\n")
cat("Script 02 complete. Run Script 03 (PCA) next.\n")

#metadata and data organization
#load expression data, organize metadata, and separate by tissue
#install and load packages 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")        # differential expression
install.packages("tidyverse")        # data organization + ggplot2
install.packages("pheatmap")         # heatmaps
install.packages("ggplot2")          # figures
install.packages("GEOquery")         # download GEO data

#load libraries
library(limma)
library(tidyverse)
library(ggplot2)
library(pheatmap)

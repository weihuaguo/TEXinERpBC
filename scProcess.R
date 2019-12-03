# Single-cell RNA sequencing data standard processing procedure
# Weihua Guo, Ph.D.
# Start date: 12/02/2019
# Functions: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

rm(list = ls())

suppressMessages(library(Seurat))

dataDir = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/scRNAseqTex/"
ctsFile = "Tumor_LN_extracted_raw_counts.txt"
cellAnnFile = "Tumor_LN_cell_annotation.txt"
resDir = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_results/scRNAseqTexResults/"

expID = "tumor_ln_wo_scimp_test"
dir.create(file.path(resDir, expID), showWarnings = TRUE)
expDirPre = paste(resDir, expID, "/", expID, "_", sep = "")

cat("Start to read counts...\n")
ctsDf = read.table(paste(dataDir, ctsFile, sep = ""), header = TRUE, row.names = 1, sep = "\t")


# scDataPrep.R
# Weihua Guo, Ph.D.
# Start date: 12/02/2019
# Check date 1: 12/30/2019
# Functions:
# 1) Read raw count matrix with cell annotation
# 2) Extract the cells of interest
# 3) Imputation (optional)
# 4) Batch correction
# 5) Create a Seurat object for further analysis

suppressMessages(library(stringr))

dataDir = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/scRNAseqTex/"
rawCtsFile = "04222019_wta_all_merged_v5_cleaned.txt" # RAW
cellAnnFile = "04252019_flow_ann_tissue_only.csv" # RAW
# resDir = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_results/scRNAseqTexResults/"


extractFlag = TRUE
saveExtract = TRUE
# tisOi = c("Tumor", "LN", "Normal", "PBMC", "Colorectal", "Melanoma")
tisOi = c("Tumor")
nKclust = 5
tisOiStr = paste(c(tisOi,"scImpute", nKclust), collapse = "_")
dir.create(file.path(dataDir, tisOiStr), showWarnings = TRUE)
expDir = paste(dataDir, tisOiStr, "/", sep = "")

scImpFlag = TRUE

cat("Start to read raw count matrix...\n")
st = Sys.time()
rawCtsDf = read.table(paste(dataDir, rawCtsFile, sep = ""), header = TRUE, row.names = 1, sep = "\t")
cat("\tGene:", dim(rawCtsDf)[1], "\t", "Cell:", dim(rawCtsDf)[2], "\n")
cat("\t")
print(Sys.time()-st)

if (extractFlag) {
	cat("Only analyze cells from specific tissue types\t", tisOi, "\n")
	cellAnnDf = read.csv(paste(dataDir, cellAnnFile, sep = ""), header = TRUE, row.names = 1)
#	print(cellAnnDf[1:9,1:4])
	# Add a column to merge the annotation of mutliple columns
	for (itis in 1:length(tisOi)) {
		cat("\t", tisOi[itis], "\n")
		if (itis == 1) {
			tmpCellAnn = as.data.frame(t(cellAnnDf))
			tmpCellAnn$filCol = tmpCellAnn[,tisOi[itis]]
			tmpCellAnn$tissue = NA
		} else {
			tmpCellAnn$filCol = tmpCellAnn[,tisOi[itis]] + tmpCellAnn[,"filCol"]
		}
		tmpCellAnn[tmpCellAnn[,tisOi[itis]] == 1, "tissue"] = tisOi[itis]
	}
	selectCellAnn = tmpCellAnn[tmpCellAnn$filCol > 0,]
	cat("\tAfter extraction,", nrow(selectCellAnn), "cells left\n")
	# Add plate and cell column NOTE: REPLACE DOT is a nightmare!!
	selectCellAnn$plate = str_split_fixed(rownames(selectCellAnn), "\\.", n = 2)[,1]
	proCtsDf = rawCtsDf[,rownames(selectCellAnn)]
	tisOiStr = paste(tisOi, collapse = "_")
	if (saveExtract) {
		write.table(selectCellAnn, 
			    paste(expDir, tisOiStr, "_cell_annotation.txt", sep = ""))
		write.table(proCtsDf,
			    paste(expDir, tisOiStr, "_extracted_raw_counts.txt", sep = ""))
	}
	print(dim(proCtsDf))

} else {
	# Directly use extracted file
	cat("Analyze all the cells...\n")
	selectCellAnn = read.csv(paste(dataDir, cellAnnFile, sep = ""), header = TRUE, row.names = 1)
	proCtsDf = rawCtsDf
}
print(proCtsDf[1:9, 1:6])

if (scImpFlag) {
	st = Sys.time()
	cat("Start to impute with scImpute...\n")
	suppressMessages(library(scImpute))
	proCtsTxt = paste(expDir, tisOiStr, "_extracted_raw_counts.txt", sep = "")
	scimpute(count_path = proCtsTxt, infile = "txt", outfile = "txt", out_dir = expDir,
		 Kcluster = 5, ncores = 12)
	print(Sys.time()-st)
}

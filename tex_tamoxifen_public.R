# GSE82721/2 
# Weihua Guo, Ph.D.
# 07/09/2020
rm(list = ls())

cat("Loading packages...\n")
suppressMessages(library(genefu))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(readxl))
suppressMessages(library(ggplot2))
suppressMessages(library(GEOquery))

data_dir <- "/home/weihua/mnts/group_plee/Weihua/tamoxifen_public_data/"
geoq_read <- FALSE

if (geoq_read) {
	cat("Loading GSE82713...\n")
	geoqst <- Sys.time()
	gse82173 <- getGEO('GSE82173', GSEMatrix = T)
	saveRDS(gse82173, paste(data_dir, "GSE82173_from_GEOQuery.RDS", sep = ""))
	print(Sys.time()-geoqst)
} else {
	cat("Reading RDS file for GSE82173...\n")
	gse82173 <- readRDS(paste(data_dir, "GSE82173_from_GEOQuery.RDS", sep = ""))
}
print(show(gse82173))
file_names <- names(gse82173)
expr_df <- experimentData(gse82173[[1]])
print(dim(expr_df))
q(save = "no")

data_dir <- "/home/weihua/mnts/group_plee/Weihua/tamoxifen_public_data/"
gse1_file <- "GSE82171_series_matrix.txt"
gse2_file <- "GSE82172_series_matrix.txt"
gpl1_file <- "GPL570_55999.xlsx"

cat("Reading series table...\n")
gse1_df <- read.table(paste(data_dir, gse1_file, sep = ""), skip = 60, header = TRUE, 
		      row.names = 1, sep = "\t", comment.char = "!")
print(dim(gse1_df))
print(gse1_df[1:9, 1:6])
gpl1_df <- read_excel(paste(data_dir, gpl1_file, sep = ""))
gpl1_df <- as.data.frame(gpl1_df)
print(dim(gpl1_df))
rownames(gpl1_df) <- gpl1_df$ID
print(gpl1_df[1:9, 1:6])



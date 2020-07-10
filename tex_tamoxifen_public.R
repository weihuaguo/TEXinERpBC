# GSE82721/2 
# Weihua Guo, Ph.D.
# 07/09/2020
rm(list = ls())

cust_sig <- function(gse, gpl, sig, match_col = "Gene Symbol") {
	ssst <- Sys.time()
	cat("Preparing for the input matrices for sig.score...\n")
	x_df <- as.data.frame(matrix(ncol = 3, nrow = nrow(sig)))
	colnames(x_df) <- c("probe", "EntrezGene.ID", "coefficient")
	rownames(x_df) <- rownames(sig)
	x_df$probe <- rownames(sig)
	x_df$coefficient <- sig$logfc
	x_df$EntrezGene.ID <- rownames(x_df)
	clean_gpl <- gpl[!is.na(gpl[,match_col]),]
	clean_gpl$EntrezGene.ID <- clean_gpl$ENTREZ_GENE_ID
	cat("Preparing x...\n")
	for (ir in rownames(x_df)) {
		tmp_mask <- clean_gpl[,match_col] == ir
		if (sum(tmp_mask) > 0) {
			clean_gpl$EntrezGene.ID[tmp_mask] <- ir
		} else {
			detc_mask <- str_detect(clean_gpl[,match_col], ir)
			clean_gpl$EntrezGene.ID[detc_mask] <- ir
	#		print(clean_gpl1[tmp_mask, "ENTREZ_GENE_ID"])
			if (sum(detc_mask) == 0) {
			cat("\t", ir, " is not found!\n")
			} 
		}
	}
	# print(x_df)
	annot_df <- clean_gpl
	annot_df$probe <- rownames(annot_df)
	data_df <- gse[rownames(annot_df),]
	data_df <- t(data_df)
	print(dim(data_df))
	print(dim(annot_df))
	print(dim(x_df))
	sig_score <- sig.score(x = x_df, data = data_df, annot = annot_df, 
				do.mapping = T, signed = T, verbose = T)
	cat("Time cost ", Sys.time()-ssst, "\n")
	return(sig_score)
}

cat("Loading packages...\n")
suppressMessages(library(genefu))
suppressMessages(library(limma))
suppressMessages(library(ggpubr))
suppressMessages(library(edgeR))
suppressMessages(library(readxl))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(GEOquery))

data_dir <- "/home/weihua/mnts/group_plee/Weihua/tamoxifen_public_data/"
geoq_read <- FALSE

if (FALSE) {
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
}

data_dir <- "/home/weihua/mnts/group_plee/Weihua/tamoxifen_public_data/"
gse1_file <- "GSE82171_series_matrix.txt"
gse2_file <- "GSE82172_series_matrix.txt"
gpl1_file <- "GPL570_55999.xlsx"
gpl2_file <- "gpl13158_5065.xlsx"
# sign_file <- "tex_signature_colt_c2.txt"
sign_file <- "interferon_gamma_signature.txt"
spl_file <- "manual_annotation_GSE82171_2_v1.xlsx"

signame <- "IFNg"
png_res <- 300

cat("Reading signatures...\n")
sign_df <- read.table(paste(data_dir, sign_file, sep = ""), header = TRUE, row.names = 1, 
		      sep = "\t")

cat("Reading the sample information...\n")
spl_df <- as.data.frame(read_excel(paste(data_dir, spl_file, sep = ""), sheet = "sum"))
rownames(spl_df) <- spl_df$GSM
spl_df$Res <- "Responsed"
ns_mask <- str_detect(spl_df$Response, "D")
spl_df$Res[ns_mask] <- "Non-responsed"
print(head(spl_df))

cat("Reading series table -- GSE82171\n")
gse1_df <- read.table(paste(data_dir, gse1_file, sep = ""), skip = 60, header = TRUE, 
		      row.names = 1, sep = "\t", comment.char = "!")
# print(dim(gse1_df))
# print(gse1_df[1:9, 1:6])
gpl1_df <- read_excel(paste(data_dir, gpl1_file, sep = ""))
gpl1_df <- as.data.frame(gpl1_df)
# print(dim(gpl1_df))
rownames(gpl1_df) <- gpl1_df$ID
# print(gpl1_df[1:9,])
# print(colnames(gpl1_df))
sig_res1 <- cust_sig(gse = gse1_df, gpl = gpl1_df, sig = sign_df)

cat("Reading series table -- GSE82172\n")
gse2_df <- read.table(paste(data_dir, gse2_file, sep = ""), skip = 60, header = TRUE, 
		      row.names = 1, sep = "\t", comment.char = "!")
# print(dim(gse2_df))
# print(gse2_df[1:9, 1:6])
gpl2_df <- read_excel(paste(data_dir, gpl2_file, sep = ""))
gpl2_df <- as.data.frame(gpl2_df)
# print(dim(gpl2_df))
rownames(gpl2_df) <- gpl2_df$ID
# print(gpl2_df[1:9, 1:6])
sig_res2 <- cust_sig(gse = gse2_df, gpl = gpl2_df, sig = sign_df)

sig_score1 <- as.data.frame(sig_res1$score)
colnames(sig_score1) <- c(signame)
sig_score2 <- as.data.frame(sig_res2$score)
colnames(sig_score2) <- c(signame)

sig_scores <- rbind(sig_score1, sig_score2)
spl_df <- spl_df[rownames(sig_scores),]
plot_scores <- cbind(sig_scores, spl_df)
print(head(plot_scores))

boxGG <- ggplot(plot_scores, aes_string(x = "Response", y = signame)) +
	geom_boxplot(aes_string(color = "Response"), fill = NA) +
	geom_jitter(position = position_jitter(0.2), aes_string(color = "Response", shape = "GSE"), 
		    size = 2.4, alpha = 0.6) +
	scale_color_manual(values = c("firebrick1","dodgerblue1",  "firebrick4", "dodgerblue4")) +
	stat_compare_means(data = plot_scores, mapping = aes_string(group = "Response"), 
			   label = "p.signif", label.x = 1.35, vjust = 1, method = "anova") +
	labs(y = paste(signame, "signature"), color = "Tamoxifen\nresponses", shape = "Dataset") +
#	facet_wrap(~ gene, ncol = 7, scales = "free") +
	theme_bw() +
	theme(axis.title.x = element_blank(),
	      legend.position = "right")
ggsave(paste(data_dir, signame, "_boxplot_tex_specific_response_all.png", sep = ""), boxGG,
       dpi = png_res, width = 4, height = 4.5)

boxGG <- ggplot(plot_scores, aes_string(x = "Res", y = signame)) +
	geom_boxplot(aes_string(color = "Res"), fill = NA) +
	geom_jitter(position = position_jitter(0.2), aes_string(color = "Res", shape = "GSE"), 
		    size = 2.4, alpha = 0.6) +
	scale_color_manual(values = c("firebrick2","dodgerblue2")) +
	stat_compare_means(data = plot_scores, mapping = aes_string(group = "Res"), 
			   label = "p.signif", label.x = 1.35, vjust = 1, method = "kruskal.test") +
	labs(y = paste(signame, "signature"), color = "Tamoxifen\nresponses", shape = "Dataset") +
#	facet_wrap(~ gene, ncol = 7, scales = "free") +
	theme_bw() +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1),
	      legend.position = "right")
ggsave(paste(data_dir, signame,  "_boxplot_tex_simple_response_all.png", sep = ""), boxGG,
       dpi = png_res, width = 3, height = 4.5)

boxGG <- ggplot(plot_scores, aes_string(x = "Res", y = signame)) +
	geom_boxplot(aes_string(color = "Res"), fill = NA) +
	geom_jitter(position = position_jitter(0.2), aes_string(color = "Res", shape = "GSE"), 
		    size = 2.4, alpha = 0.6) +
	scale_color_manual(values = c("firebrick2","dodgerblue2")) +
	stat_compare_means(data = plot_scores, mapping = aes_string(group = "Res"), 
			   label = "p.signif", label.x = 1.35, vjust = 1, method = "kruskal.test") +
	labs(y = paste(signame, "signature"), color = "Tamoxifen\nresponses", shape = "Dataset") +
	facet_wrap(~ GSE, ncol = 2, scales = "free") +
	theme_bw() +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1),
	      legend.position = "right")
ggsave(paste(data_dir, signame,  "_boxplot_tex_simple_response_split.png", sep = ""), boxGG,
       dpi = png_res, width = 5, height = 4.5)


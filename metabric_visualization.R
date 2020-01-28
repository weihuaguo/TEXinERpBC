# This script is for visualizing the results from METABRIC database
# Heatmap of transcriptomics
# CIBERSORT results with statistical test
# Weihua Guo, Ph.D.
# 01/27/2020

rm(list = ls())
suppressMessages(library(RColorBrewer))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(circlize))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))


dataDir = "/home/weihua/mnts/group_plee/Weihua/PD1_CD39_Public_Data/final_used_data/"
resDir = "/home/weihua/mnts/group_plee/Weihua/PD1_CD39_Public_Data/final_outputs/"

deCsv = "_sym25_tex_ER_cbpt_de_limma.csv"
exprCsv = "data_expression_median.csv"
texCsv = "sub_scres_group_ER_cbpt.csv"
hallGMT = "h.all.v7.0.symbols.gmt"
cbstTxt = "cbpt_brca_cibersort_annot_merged.txt"
geneWNames = c("GZMB", "IL7R", "HLA-DRB1", "CCL5", "HLA-DRA", "IDO1", "CD79A", "CD79B", "HLA-DQA1", "CXCL13", "STAT1", "CXCL10")

cat("Reading the METABRIC data...")
st = Sys.time()
deResDf = read.table(paste(dataDir, deCsv, sep = ""), header = TRUE, row.names = 1, sep = ",")
# exprDf = read.table(paste(dataDir, exprCsv, sep = ""), header = TRUE, row.names = 1, sep = ",")
# exprDf = readRDS(paste(dataDir, "data_expression_median.RDS", sep = ""))
texGroup = read.table(paste(dataDir, texCsv, sep = ""), header = TRUE, row.names = 1, sep = ",")
hallDf = as.data.frame(read.table(paste(dataDir, hallGMT, sep = ""), sep = "\t", row.names = 1, fill = TRUE, stringsAsFactors = FALSE))
cbstDf = read.table(paste(dataDir, cbstTxt, sep = ""), sep = "\t", row.names = 2, header = TRUE)

print(Sys.time()-st)

# print(head(deResDf))
# print(exprDf[1:9,1:6])
# print(head(texGroup))
# print(head(cbstDf))

pidOi = rownames(texGroup)[texGroup$group != "Medium"]
clinOiAnn = texGroup[texGroup$group != "Medium",]
# print(pidOi)

cbstDf[,1] = NULL
cbstERDf = cbstDf[pidOi,]
cbstERDf = cbind(cbstERDf, clinOiAnn)
cbstCellTypes = colnames(cbstERDf)[1:(match("P.value", colnames(cbstERDf))-1)]
cbstSprERDf = gather(cbstERDf, "cell_type", "relative_abundance", cbstCellTypes)
print(head(cbstSprERDf))

cbstGG = ggboxplot(cbstSprERDf, x = "cell_type", y = "relative_abundance", color = "group", 
		   add = "jitter", add.params = list(alpha = 0.2, size = 0.9),
		   palette = c("#D25565", "#0B9BC6"), outlier.shape = NA,
		   ylab = "Relative abundance", xlab = "Immune cell type", legend = "right") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif.., group = group), hide.ns = TRUE, label.y = 0.72) +
	rotate_x_text(45)
cbstHGG = ggpar(cbstGG, legend.title = "Tex group")
ggsave(cbstHGG, filename = paste(resDir, "cbst_er_pos_metabric_boxplot_horz.tiff", sep = ""), dpi = 300, width = 9, height = 6)
cbstGG = ggboxplot(cbstSprERDf, x = "cell_type", y = "relative_abundance", color = "group", 
		   add = "jitter", add.params = list(alpha = 0.2, size = 0.9),
		   palette = c("#D25565", "#0B9BC6"), outlier.shape = NA,
		   ylab = "Relative abundance", xlab = "Immune cell type", legend = "right") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif.., group = group), hide.ns = TRUE, label.y = 0.72) +
	coord_flip()
cbstVGG = ggpar(cbstGG, legend.title = "Tex group")
ggsave(cbstVGG, filename = paste(resDir, "cbst_er_pos_metabric_boxplot_vert.tiff", sep = ""), dpi = 300, width = 6, height = 9)

sigUpGenes = rownames(deResDf)[deResDf$logFC >= 1.0 & deResDf$adj.P.Val <= 0.10]
sigDwGenes = rownames(deResDf)[deResDf$logFC <= -1.0 & deResDf$adj.P.Val <= 0.10]
# print(sigUpGenes)
# print(sigDwGenes)

# Prepare hallmark GMT
hallDf = hallDf[,-1]
hallDf = t(hallDf)

# Overlap pathway-related genes and significantly differentially expressed genes 
hallOlDf = as.data.frame(matrix(ncol = 7, nrow = ncol(hallDf)))
rownames(hallOlDf) = colnames(hallDf)
colnames(hallOlDf) = c("hall.gene", "in.up", "in.down", "up.per", "down.per", "per.up", "per.down")
for (icol in colnames(hallDf)) {
	tmpGenes = hallDf[hallDf[,icol] != "", icol]
	hallOlDf[icol, "hall.gene"] = length(tmpGenes)
	hallOlDf[icol, "in.up"] = length(intersect(tmpGenes, sigUpGenes))
	hallOlDf[icol, "in.down"] = length(intersect(tmpGenes, sigDwGenes))
	hallOlDf[icol, "up.per"] = length(intersect(sigUpGenes, tmpGenes))/length(tmpGenes)
	hallOlDf[icol, "down.per"] = length(intersect(sigDwGenes, tmpGenes))/length(tmpGenes)
	hallOlDf[icol, "per.up"] = length(intersect(sigUpGenes, tmpGenes))/length(sigUpGenes)
	hallOlDf[icol, "per.down"] = length(intersect(sigDwGenes, tmpGenes))/length(sigDwGenes)

}

# Generate gene annotation for upregulated genes
hallUpOlDf = hallOlDf[order(hallOlDf$per.up),]
hallUpOlDf = hallUpOlDf[hallUpOlDf$per.up != 0,]
# print(head(hallUpOlDf))
upGeneAnn = as.data.frame(matrix(ncol = 2, nrow = length(sigUpGenes)))
colnames(upGeneAnn) = c("Pathway1", "FC")
rownames(upGeneAnn) = sigUpGenes
for (irow in rownames(hallUpOlDf)) {
#	print(hallOlDf[irow,])
	tmpOl = intersect(hallDf[, irow], sigUpGenes)
	upGeneAnn[rownames(upGeneAnn) %in% tmpOl, "Pathway1"] = irow
#	print(tmpOl)
}

# Generate gene annotation for upregulated genes
hallDwOlDf = hallOlDf[order(hallOlDf$per.down),]
hallDwOlDf = hallDwOlDf[hallDwOlDf$per.down != 0,]
# print(head(hallDwOlDf))
dwGeneAnn = as.data.frame(matrix(ncol = 2, nrow = length(sigDwGenes)))
colnames(dwGeneAnn) = c("Pathway1", "FC")
rownames(dwGeneAnn) = sigDwGenes
for (irow in rownames(hallDwOlDf)) {
#	print(hallOlDf[irow,])
	tmpOl = intersect(hallDf[, irow], sigDwGenes)
	dwGeneAnn[rownames(dwGeneAnn) %in% tmpOl, "Pathway1"] = irow
#	print(tmpOl)
}

geneAnn = rbind(upGeneAnn, dwGeneAnn)
geneAnn = geneAnn[!is.na(geneAnn$Pathway1),]
geneAnn$Pathway1 = str_replace(geneAnn$Pathway1, "HALLMARK_", "")
subDE = deResDf[rownames(geneAnn),]
geneAnn$FC = subDE$logFC

## Remove the pathways with <3 genes
keepPaths = names(which(table(geneAnn$Pathway1) > 3))
print(keepPaths)
geneAnn = geneAnn[geneAnn$Pathway1 %in% keepPaths,]
print(dim(geneAnn))
sigGenes = rownames(geneAnn)

write.csv(hallOlDf, paste(resDir, "hallmark_row_annotation_overlap_result.csv", sep = ""))


# print(length(intersect(c(sigUpGenes, sigDwGenes), exprDf$Hugo_Symbol)))
exprOiData = exprDf[exprDf$Hugo_Symbol %in% sigGenes,]
exprOiData = exprOiData[,c("Hugo_Symbol", pidOi)]
# print(length(c(sigUpGenes, sigDwGenes)))
# print(dim(exprOiData))
rownames(exprOiData) = exprOiData$Hugo_Symbol
exprOiData$Hugo_Symbol = NULL


backClinAnn = clinOiAnn
clinAnn = backClinAnn[order(backClinAnn$group),]
print(head(clinAnn))

exprOiData = exprOiData[,rownames(backClinAnn)[order(backClinAnn$group)]]

scaleExprOiData = scale(t(exprOiData), center = TRUE, scale = TRUE)
hmScaleExprData = t(scaleExprOiData)
hmScaleExprData = hmScaleExprData[, rownames(backClinAnn)[order(backClinAnn$group)]]
hmScaleExprData = hmScaleExprData[rownames(geneAnn)[order(geneAnn$Pathway1)],]
# print(exprOiData[1:9,1:6])

geneHa = rowAnnotation(ig = anno_mark(at = match(geneWNames, rownames(hmScaleExprData)), labels = geneWNames))

geneAnn$Pathway1 = factor(geneAnn$Pathway1)
pathCols = colorRampPalette(brewer.pal(6,"Spectral"))(length(levels(geneAnn$Pathway1)))
names(pathCols) = levels(geneAnn$Pathway1)
pathHa = rowAnnotation(Pathway = geneAnn$Pathway1, 
		       col = list(Pathway = pathCols),
		       annotation_legend_param = list(Pathway = list(title = "Pathways (HALLMARK V7.0)", 
								     legend_direction = "horizontal", 
								     nrow = 3)))
names(pathHa) = c("")

texHa = HeatmapAnnotation(texg = clinAnn$group,
			  col = list(texg = c("Low" = "#0B9BC6", "High" = "#D25565")),
			  annotation_legend_param = list(texg = list(title = "Tex group", legend_direction = "horizontal", nrow = 2)))
names(texHa) = c("Tex group")
scaleCol = colorRamp2(seq(-3, 3, length = 3), c("purple", "black", "yellow"))
## Unscaled version
exprHm = Heatmap(as.matrix(exprOiData), name = "Expression",
#		     col = scaleCol,
                     show_column_names = FALSE,
                     cluster_columns = FALSE,
                     cluster_rows = TRUE,
#		     column_order = order(clinAnn$TexFACSGroup),
#		     row_order = order(colnames(rawScore)),
                     top_annotation = texHa,
                     heatmap_legend_param = list(direction = "horizontal"))#, col_fun = scaleCol))

tiff(paste(resDir, "normalized_expression_heatmap_unscale.tiff", sep = ""), res = 300, height = 12, width = 9, units = "in")
draw(exprHm, heatmap_legend_side = "bottom", merge_legend = TRUE)
gar = dev.off()

## Scaled version
exprHm = Heatmap(as.matrix(hmScaleExprData), name = "Expression",
		     col = scaleCol,
                     show_column_names = FALSE,
		     show_row_names = FALSE,
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
		     row_order = order(geneAnn$Pathway1),
                     top_annotation = texHa,
		     left_annotation = pathHa,
		     right_annotation = geneHa,
                     heatmap_legend_param = list(direction = "horizontal", col_fun = scaleCol))

tiff(paste(resDir, "normalized_expression_heatmap_scaled_by_gene.tiff", sep = ""), res = 300, height = 6, width = 9, units = "in")
draw(exprHm, heatmap_legend_side = "top", merge_legend = TRUE)
gar = dev.off()

## Scaled version with gene names
exprHm = Heatmap(as.matrix(hmScaleExprData), name = "Expression",
		     col = scaleCol,
                     show_column_names = FALSE,
		     show_row_names = TRUE,
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
		     row_order = order(geneAnn$Pathway1),
                     top_annotation = texHa,
		     left_annotation = pathHa,
                     heatmap_legend_param = list(direction = "horizontal", col_fun = scaleCol))

tiff(paste(resDir, "normalized_expression_heatmap_scaled_by_gene_w_genename.tiff", sep = ""), 
     res = 300, height = 18, width = 9, units = "in")
draw(exprHm, heatmap_legend_side = "top", merge_legend = TRUE)
gar = dev.off()


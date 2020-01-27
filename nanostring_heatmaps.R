# This script is for visualizing NanoString results for PD1+CD39+ T exhausted cells
# Weihua Guo, Ph.D.
# 11/13/2019

suppressMessages(library(ComplexHeatmap))
suppressMessages(library(dplyr))
suppressMessages(library(readxl))
suppressMessages(library(tibble))
suppressMessages(library(circlize))

# Directory in which you save the SI excels
dataDir = "/home/weihua/mnts/group_plee/Weihua/PD1_CD39_nanostring/final_used_data/"
resDir = "/home/weihua/mnts/group_plee/Weihua/PD1_CD39_nanostring/final_results/"

clinAnnXlsx = "facs_ihc_clinical_annot.xlsx"
rawScoreXlsx = "ns_cell_score_raw_log2_no_TIL.xlsx"
rltScoreXlsx = "ns_cell_score_relative_log2_vs_TIL.xlsx"

clinAnn = read_excel(paste(dataDir, clinAnnXlsx, sep = ""))
rawScore = read_excel(paste(dataDir,rawScoreXlsx, sep = ""))
rltScore = read_excel(paste(dataDir, rltScoreXlsx, sep = ""))

clinAnn = clinAnn %>% column_to_rownames("pid")
rawScore = rawScore %>% column_to_rownames("pid")
rltScore = rltScore %>% column_to_rownames("pid")

# print(rownames(clinAnn))
# print(head(clinAnn))
# print(rownames(rltScore)) 


scaleRawScore = scale(rawScore, center = TRUE, scale = TRUE)
scaleRLTScore = scale(rltScore, center = TRUE, scale = TRUE)

write.csv(scaleRawScore, file = paste(resDir, "raw_score_no_TIL_log2_scaled_by_cell_type.csv", sep = ""))
write.csv(scaleRLTScore, file = paste(resDir, "relative_score_no_TIL_log2_scaled_by_cell_type.csv", sep = ""))

backClinAnn = clinAnn
clinAnn = backClinAnn[order(backClinAnn$TexFACS),]

hmRawScore = t(rawScore)
hmRawScore = as.data.frame(hmRawScore[,rownames(backClinAnn)[order(backClinAnn$TexFACS)]])

hmScaleRawScore = t(scaleRawScore)
hmScaleRawScore = as.data.frame(hmScaleRawScore[,rownames(backClinAnn)[order(backClinAnn$TexFACS)]])

hmRltScore = t(rltScore)
hmRltScore = as.data.frame(hmRltScore[,rownames(backClinAnn)[order(backClinAnn$TexFACS)]])

hmScaleRltScore = t(scaleRLTScore)
hmScaleRltScore = as.data.frame(hmScaleRltScore[,rownames(backClinAnn)[order(backClinAnn$TexFACS)]])

clinHa = HeatmapAnnotation(FACS = anno_barplot(clinAnn[,"TexFACS"], gp = gpar(fill = "#D25565", col = "#D25565")),
                           CD8 = anno_barplot(clinAnn[,"CD8"], gp = gpar(fill = "cyan", col = "cyan")),
                           CD20 = anno_barplot(clinAnn[,"CD20"], gp = gpar(fill = "#7030A0", col = "#7030A0")),
                           PDL1 = anno_barplot(clinAnn[,"PDL1"], gp = gpar(fill = "#996633", col = "#996633")),
                           gap = unit(2.4, "mm")
                           )
names(clinHa) = c("Tex (FACS%)", "CD8 (IHC%)", "CD20 (IHC%)", "PD-L1(IHC%)")
scaleCol = colorRamp2(seq(-3, 3, length = 3), c("blue", "white", "red"))
# print(scaleCol(c(-4,-3.5,3.5, 4)))

rawScoreHm = Heatmap(hmRawScore, name = "Cell score (log2)",
                     show_column_names = TRUE,
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     column_order = order(clinAnn$TexFACS),
                     row_order = order(colnames(rawScore)),
                     top_annotation = clinHa,
		     cell_fun = function(j, i, x, y, width, height, fill) {
			     grid.text(sprintf("%.2f", hmRawScore[i, j]), x, y, gp = gpar(fontsize = 9))},
                     heatmap_legend_param = list(direction = "horizontal"))

tiff(paste(resDir, "raw_score_log2_no_TIL_heatmap_test.tiff", sep = ""), res = 300, height = 6, width = 15, units = "in")
draw(rawScoreHm, heatmap_legend_side = "bottom")
gar = dev.off()

rltScoreHm = Heatmap(hmRltScore, name = "Relative abundance",
                     show_column_names = TRUE,
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     column_order = order(clinAnn$TexFACS),
                     row_order = order(colnames(rltScore)),
                     top_annotation = clinHa,
		     cell_fun = function(j, i, x, y, width, height, fill) {
			     grid.text(sprintf("%.2f", hmRltScore[i, j]), x, y, gp = gpar(fontsize = 9))},
                     heatmap_legend_param = list(direction = "horizontal"))

tiff(paste(resDir, "relative_score_log2_heatmap_test.tiff", sep = ""), res = 300, height = 6, width = 15, units = "in")
draw(rltScoreHm, heatmap_legend_side = "bottom")
gar = dev.off()

rawSScoreHm = Heatmap(hmScaleRawScore, name = "Scaled cell score (log2)",
		     col = scaleCol,
                     show_column_names = TRUE,
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     column_order = order(clinAnn$TexFACS),
                     row_order = order(colnames(rawScore)),
                     top_annotation = clinHa,
		     cell_fun = function(j, i, x, y, width, height, fill) {
			     grid.text(sprintf("%.2f", hmScaleRawScore[i, j]), x, y, gp = gpar(fontsize = 9))},
                     heatmap_legend_param = list(direction = "horizontal", col_fun = scaleCol))

tiff(paste(resDir, "raw_score_log2_no_TIL_heatmap_scaled_by_cell_type_test.tiff", sep = ""), res = 300, height = 6, width = 15, units = "in")
draw(rawSScoreHm, heatmap_legend_side = "bottom")
gar = dev.off()


rltSScoreHm = Heatmap(hmScaleRltScore, name = "Scaled relative abundance",
		     col = scaleCol,
                     show_column_names = TRUE,
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     top_annotation = clinHa,
                     row_order = order(colnames(rltScore)),
		     cell_fun = function(j, i, x, y, width, height, fill) {
			     grid.text(sprintf("%.2f", hmScaleRltScore[i, j]), x, y, gp = gpar(fontsize = 9))},
                     heatmap_legend_param = list(direction = "horizontal"))

tiff(paste(resDir, "relative_score_log2_heatmap_scaled_by_cell_type_test.tiff", sep = ""), res = 300, height = 6, width = 15, units = "in")
draw(rltSScoreHm, heatmap_legend_side = "bottom")
gar = dev.off()

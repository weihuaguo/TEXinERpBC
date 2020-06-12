# Compare tumor cells between patients with and without Tex
# Weihua Guo, Ph.D.
# 06/11/2020

rm(list = ls())

scQuickPro <- function(srsc, plotPf, res = 300, ndim = 15, rdsSave = FALSE, jsFlag = FALSE, plotFeat = NULL) {
	cat("Start to QC and pre-processing Seurat object...\n")

	mito.genes <- grep(pattern = "^MT-", x = rownames(srsc@assays[["RNA"]]), value = TRUE)
	srsc[["percent.mt"]] <- PercentageFeatureSet(srsc, pattern = "^MT-")
#	print(head(srsc@meta.data))

	## NOTE: need some basic filter to remove totally empty droplet!!!
	qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), 
	       dpi = res, width = 9, height = 6)

	qcsc1 <- FeatureScatter(srsc, feature1 = "nCount_RNA", feature2 = "percent.mt")
	qcsc2 <- FeatureScatter(srsc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	cqcgg <- CombinePlots(plots = list(qcsc1, qcsc2))
	ggsave(plot = cqcgg, filename = paste(plotPf, "QCPlot2.png", sep = ""), 
	       dpi = res, width = 9, height = 6)
	srsc <- subset(srsc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

	###################
	srsc <- NormalizeData(srsc, assay = "RNA", normalization.method = "LogNormalize", 
			     scale.factor = 10000)
	################
	srsc <- FindVariableFeatures(srsc, selection.method = "vst", nfeatures = 2000)
	top10Features <- head(VariableFeatures(srsc), 10)
	varGG <- VariableFeaturePlot(srsc)
	labVargg <- LabelPoints(plot = varGG, points = top10Features, repel = TRUE)
	ggsave(plot = labVargg, filename = paste(plotPf, "variable_features.png", sep = ""), 
	       dpi = res, width = 9, height = 6)

	##############
	allGenes <- rownames(srsc[["RNA"]])
	srsc <- ScaleData(srsc, features = allGenes, assay = "RNA")

	#############
	srsc <- RunPCA(srsc, features = VariableFeatures(object = srsc))

	pcaGG1 <- VizDimLoadings(srsc, dim = 1:6, reduction = "pca")
	ggsave(plot = pcaGG1, filename = paste(plotPf, "PCA1.png", sep = ""), dpi = res, width = 9, height = 9)
	
	pcaGG2 <- DimPlot(srsc, reduction = "pca")
	ggsave(plot = pcaGG2, filename = paste(plotPf, "PCA2.png", sep = ""), dpi = res, width = 9, height = 6)

	png(paste(plotPf, "PCA_HEATMAP.png", sep = ""), res = res, width = 9, height = 16, units = "in")
	DimHeatmap(srsc, dims = 1:15, cells = 500, balanced = TRUE)
	gar <- dev.off()

	# TODO: JUST FOR SAVING TIME
	if (jsFlag) {
		srsc <- JackStraw(srsc, num.replicate = 100)
		srsc <- ScoreJackStraw(srsc, dims = 1:20)
		jsGG <- JackStrawPlot(srsc, dims = 1:20)
		ggsave(plot = jsGG, filename = paste(plotPf, "JackStraw.png", sep = ""), dpi = res, width = 9, heigh = 6)
	}

	elGG <- ElbowPlot(srsc, ndims = 50)
	ggsave(plot = elGG, filename = paste(plotPf, "Elbow.png", sep = ""), dpi = res, width = 9, heigh = 6)

	srsc <- FindNeighbors(srsc, dims = 1:ndim)
	srsc <- FindClusters(srsc, resolution = 0.5)

	srsc <- RunUMAP(srsc, dims = 1:ndim)
	umapGG <- DimPlot(srsc, reduction = "umap")
	ggsave(plot = umapGG, filename = paste(plotPf, "UMAP.png", sep = ""), dpi = res, width = 9, heigh = 6)

	outData <- FetchData(srsc, c("ident", "nGene", "UMAP1", "UMAP2"))
	write.csv(outData, file = paste(plotPf, "useful_results.csv", sep = ""))

	if (!is.null(plotFeat)) {
		# TODO: Change width and heigh automatically
		cat("Start to plot feature plots and violin plots\n")
		vlnGG <- VlnPlot(srsc, features = plotFeat, slot = "counts", log = TRUE, 
				ncol = 3, pt.size = 0.5)
		ggsave(plot = vlnGG, 
		       filename = paste(plotPf, "violin_plot_for_cell_annotation.png", sep = ""), 
		       dpi = res, width = 16, heigh = 12)

		featGG <- FeaturePlot(srsc, features = plotFeat, ncol = 3, pt.size = 0.2, sort.cell = TRUE)
		ggsave(plot = featGG, 
		       filename = paste(plotPf, "umap_feature_plot_for_cell_annotation.png", sep = ""),
		       dpi = res, width = 16, heigh = 10)

		dotGG <- DotPlot(srsc, features = plotFeat) + RotatedAxis()
		ggsave(plot = dotGG,
		       filename = paste(plotPf, "dotplot_for_cell_annotation.png", sep = ""),
		       dpi = res, width = 6, heigh = 9)
	}

	if (length(levels(outData$ident)) >= 2) {
		srsc.markers <- FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.0)
		topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 18, wt = avg_logFC)

		print(head(topMarkers))

		write.csv(srsc.markers, file <- paste(plotPf, "ALL_POS_MARKERS.csv", sep = ""))
		write.csv(topMarkers, file <- paste(plotPf, "top_pos_markers.csv", sep = ""))

		png(paste(plotPf, "marker_heatmap.png", sep = ""), res = res, width = 16, height = 12, units = "in")
		print(DoHeatmap(object = srsc, features = topMarkers$gene) + NoLegend())
	gar <- dev.off()
	}

	if (rdsSave) {
		saveRDS(srsc, file = paste(plotPf, "Seurat_Objects_Clustered.RDS", sep = ""))
	}
	print(srsc)
	return(srsc)
}

cat("Loading packages...\n")
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(readxl))
suppressMessages(library(MAST))

dataDir <- "/home/weihua/mnts/data_plee/Group/scmodb/ALL10X/"
resDir <- "/home/weihua/mnts/group_plee/Weihua/scrnaseq_results/tex_tumor/"

expID <- "tex_tumor_comparison"
inDataFormat <- "10X"
inDataType <- "filtered"
saSelect <- c("Tumor") #c("ContralateralNormal", "IpsilateralNormal", "LymphNode")
rnaSelect <- NULL # c("five")
inteFLAG <- FALSE # Run integration with Seurat
inteFilter <- 100
saveListFlag <- TRUE # Save the input list for Seurat integration
postFlag <- FALSE # Run experiments for integrated objects!
visFlag <- FALSE # Visualize results
subClusterFlag <- FALSE # Re-cluster for each general cell types

dir.create(file.path(resDir, expID), showWarnings = FALSE)
expDir <- paste(resDir, expID, "/", sep = "")

git_scfolder <- "/home/weihua/git_repo/TEXinERpBC"
rscript_file <- list.files(git_scfolder, "tumor_comp.R")
file.copy(rscript_file, expDir)

if (inDataFormat == "10X") {
	if (inDataType == "raw") {
		inFiles <- list.files(dataDir, pattern = "_raw_feature_bc_matrices")
		saName <- str_split_fixed(inFiles, "_raw", n = 2)[,1]
	}
	if (inDataType == "filtered") {
		inFiles <- list.files(dataDir, pattern = "_filtered_gene_bc_matrices")
		saName <- str_split_fixed(inFiles, "_fil", n = 2)[,1]
	}

	if (!is.null(saSelect)) {
		splCter <- 1
		for (ispl in saSelect) {
			tmpFiles <- inFiles[grepl(ispl, inFiles)]
			tmpName <- saName[grepl(ispl, saName)]
			if (splCter == 1) {
				filterFiles <- tmpFiles
				filterNames <- tmpName
			} else {
				filterFiles <- c(filterFiles, tmpFiles)
				filterNames <- c(filterNames, tmpName)
			}
			splCter <- splCter + 1
		}
		inFiles <- filterFiles
		saName <- filterNames
	}

	if (!is.null(rnaSelect)) {
		print(inFiles)
		inFiles <- inFiles[grepl(rnaSelect, inFiles)]
		saName <- saName[grepl(rnaSelect, saName)]
	}


	cat("All the samples:", saName, "\n")
	## Initital a dataframe to record the basic process results
	recColName <- c("tissue", "patient", "before_cell", "after_cell", "before_nFeature", "after_nFeature", 
			"median_mtper", "mean_mtper", "cell_kept")
	inteRecDf <- as.data.frame(matrix(ncol = length(recColName), nrow = length(saName)))
	rownames(inteRecDf) <- saName
	colnames(inteRecDf) <- recColName
	## Initial list for integration
	srData <- vector(mode = "list", length = length(saName))
	srRawObjs <- vector(mode = "list", length = length(saName))
	srCleanObjs <- vector(mode = "list", length = length(saName))
	names(srData) <- saName
	names(srRawObjs) <- saName
	names(srCleanObjs) <- saName
	for (isa in 1:length(inFiles)) {
		st <- Sys.time()
		cat(saName[isa], ":" , inFiles[isa], "\n")
		tmpDir <- paste(dataDir, inFiles[isa], sep = "")
		tmpData <- try(Read10X(tmpDir))
		if (class(tmpData) == "try-error") {
			stop("STILL NOT WORKING!!!")
		}

		if (is.null(names(tmpData))) {
			## NOTE: need some basic filter to remove totally empty droplet!!!
			tmpObj <- CreateSeuratObject(counts = tmpData, project = saName[isa], 
						     min.cells = 3, min.features = 5)
		} else {
			tmpObj <- CreateSeuratObject(counts = tmpData$`Gene Expression`, 
						     project = saName[isa], min.cells = 3, min.features = 5)
		}
		print(tmpObj)

		tisName <- str_split_fixed(saName[isa], "_", n = 4)[,3]
		rnaDirc <- str_split_fixed(saName[isa], "_", n = 4)[,2]
		patName <- str_split_fixed(saName[isa], "_", n = 4)[,1]
		cat(tisName, patName, rnaDirc, "\n")
		tmpObj@meta.data$tissue <- tisName
		tmpObj@meta.data$patient <- patName
		tmpObj@meta.data$sampleID <- saName[isa]
		tmpObj@meta.data$rna <- rnaDirc

		srData[saName[isa]] <- tmpData
		srRawObjs[saName[isa]] <- tmpObj

		plotPrefix <- paste(expDir, expID, "_", saName[isa], "_", sep = "")
		tmpProObj <- scQuickPro(tmpObj, plotPf = plotPrefix, 
					plotFeat = c("EPCAM", "PTPRC", "CD3D", "CD8A", "PDCD1", "CXCL13"))
		srCleanObjs[saName[isa]] <- tmpProObj

		inteRecDf[saName[isa], "tissue"] <- tisName
		inteRecDf[saName[isa], "patient"] <- patName
		inteRecDf[saName[isa], "before_cell"] <- length(colnames(tmpObj))
		inteRecDf[saName[isa], "after_cell"] <- length(colnames(tmpProObj))
		inteRecDf[saName[isa], "before_nFeature"] <- median(tmpObj@meta.data$nFeature_RNA)
		inteRecDf[saName[isa], "after_nFeature"] <- median(tmpProObj@meta.data$nFeature_RNA)
		inteRecDf[saName[isa], "median_mtper"] <- median(tmpProObj@meta.data$percent.mt)
		inteRecDf[saName[isa], "mean_mtper"] <- mean(tmpProObj@meta.data$percent.mt)

		print(Sys.time() - st)
		cat("+++++++++++++++++++++++++++++++++++++++++++++\n\n")
	}
	inteRecDf$cell_kept <- inteRecDf$after_cell/inteRecDf$before_cell
	write.csv(inteRecDf, paste(expDir, expID, "_integration_records.csv", sep = ""))

	if (saveListFlag) {
		st <- Sys.time()
		cat("Saving Seurat object list to a RDS file...\n")
		if (inDataType == "raw") {
			dataRDS <- paste(expDir, expID, "_all_raw_10X_data_list.RDS", sep = "")
			rawObjRDS <- paste(expDir, expID, "_all_raw_raw_10X_seurat_objects_list.RDS", sep = "")
			cleanObjRDS <- paste(expDir, expID, "_all_clean_raw_10X_seurat_objects_list.RDS", sep = "")
		} else {
			dataRDS <- paste(expDir, expID, "_all_filtered_10X_data_list.RDS", sep = "")
			rawObjRDS <- paste(expDir, expID, "_all_raw_filtered_10X_seurat_objects_list.RDS", sep = "")
			cleanObjRDS <- paste(expDir, expID, "_all_clean_filtered_10X_seurat_objects_list.RDS", sep = "")
		}
		saveRDS(srData, file = dataRDS)
		saveRDS(srRawObjs, file = rawObjRDS)
		saveRDS(srCleanObjs, file = cleanObjRDS)
		cat("Save to RDS cost")
		print(Sys.time()-st)
	}

}

# This script is designed to analyze scRNA-seq data ONLY in batch with Seurat
# Weihua Guo, Ph.D.
# 03/27/2020

rm(list = ls())
options(warn = -1)

srQCPrePro = function(srsc, plotPf, res = 300, ndim = 15, rdsSave = FALSE, basicQCOnly = TRUE, 
		      jsFlag = FALSE, hashtagFlag = TRUE, plotFeat = NULL, atbFlag = FALSE) {
	cat("Start to QC and pre-processing Seurat object...\n")
	suppressMessages(library(Seurat))
	suppressMessages(library(dplyr))
	suppressMessages(library(stringr))
	mito.genes = grep(pattern = "^MT-", x = rownames(srsc@assays[["RNA"]]), value = TRUE)
	srsc[["percent.mt"]] = PercentageFeatureSet(srsc, pattern = "^MT-")
#	print(head(srsc@meta.data))

	if (atbFlag) {
		cat("Search for Hashtag...\n")
		allProteins = srsc[["protein"]]@data
		print(dim(allProteins))
		protFeas = rownames(allProteins)
		ratioThr = 0.05
		if (any(grepl("hashtag", tolower(protFeas)))) {
			srsc@meta.data$hashtagCheck = "duplicate"
			cat("Start to QC based on hashtag...\n")
			protMat = t(as.matrix(allProteins))
			hashtagMat = protMat[,grepl("hashtag", tolower(colnames(protMat)))]
			strict2 = rownames(protMat)[hashtagMat[,1] == 0 & hashtagMat[,2] != 0]
			strict1 = rownames(protMat)[hashtagMat[,1] != 0 & hashtagMat[,2] == 0]
			loose1 = rownames(protMat)[hashtagMat[,1]*ratioThr > hashtagMat[,2]]
			loose2 = rownames(protMat)[hashtagMat[,2]*ratioThr > hashtagMat[,1]]
			srsc@meta.data[loose2, "hashtagCheck"] = "S2"
			srsc@meta.data[loose1, "hashtagCheck"] = "S1"
			rmID = rownames(srsc@meta.data)[srsc@meta.data$hashtagCheck == "duplicate"]
			cat(length(rmID), "cells are removed based on Hastag\n")
		} else {
			srsc@meta.data$hashtagCheck = "NA"
		}
		srsc = subset(srsc, subset = hashtagCheck == "duplicate", invert=TRUE)
		print(srsc)
	}
	## NOTE: need some basic filter to remove totally empty droplet!!!
	qcgg = VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), 
	       dpi = res, width = 9, height = 6)

	qcsc1 = FeatureScatter(srsc, feature1 = "nCount_RNA", feature2 = "percent.mt")
	qcsc2 = FeatureScatter(srsc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	cqcgg = CombinePlots(plots = list(qcsc1, qcsc2))
	ggsave(plot = cqcgg, filename = paste(plotPf, "QCPlot2.png", sep = ""), 
	       dpi = res, width = 9, height = 6)
	srsc = subset(srsc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

	###################
	srsc = NormalizeData(srsc, assay = "RNA", normalization.method = "LogNormalize", 
			     scale.factor = 10000)
	if (atbFlag) {
		srsc <- NormalizeData(srsc, assay = "protein", normalization.method = "CLR")
	}

	################
	srsc = FindVariableFeatures(srsc, selection.method = "vst", nfeatures = 2000)
	top10Features = head(VariableFeatures(srsc), 10)
	varGG = VariableFeaturePlot(srsc)
	labVargg = LabelPoints(plot = varGG, points = top10Features, repel = TRUE)
	ggsave(plot = labVargg, filename = paste(plotPf, "variable_features.png", sep = ""), 
	       dpi = res, width = 9, height = 6)

	if (basicQCOnly) {
		cat("Basic QC is finished and ready for integration...\n")
		if (rdsSave) {
			saveRDS(srsc, file = paste(plotPf, "basicQC_Seurat_Objects_Clustered.RDS", sep = ""))
		}
		return(srsc)
	}
	##############
	allGenes = rownames(srsc[["RNA"]])
	srsc = ScaleData(srsc, features = allGenes, assay = "RNA")

	#############
	srsc = RunPCA(srsc, features = VariableFeatures(object = srsc))

	pcaGG1 = VizDimLoadings(srsc, dim = 1:6, reduction = "pca")
	ggsave(plot = pcaGG1, filename = paste(plotPf, "PCA1.png", sep = ""), dpi = res, width = 9, height = 9)
	
	pcaGG2 = DimPlot(srsc, reduction = "pca")
	ggsave(plot = pcaGG2, filename = paste(plotPf, "PCA2.png", sep = ""), dpi = res, width = 9, height = 6)

	png(paste(plotPf, "PCA_HEATMAP.png", sep = ""), res = res, width = 9, height = 16, units = "in")
	DimHeatmap(srsc, dims = 1:15, cells = 500, balanced = TRUE)
	gar = dev.off()

	# TODO: JUST FOR SAVING TIME
	if (jsFlag) {
		srsc = JackStraw(srsc, num.replicate = 100)
		srsc = ScoreJackStraw(srsc, dims = 1:20)
		jsGG = JackStrawPlot(srsc, dims = 1:20)
		ggsave(plot = jsGG, filename = paste(plotPf, "JackStraw.png", sep = ""), dpi = res, width = 9, heigh = 6)
	}

	elGG = ElbowPlot(srsc, ndims = 50)
	ggsave(plot = elGG, filename = paste(plotPf, "Elbow.png", sep = ""), dpi = res, width = 9, heigh = 6)

	srsc = FindNeighbors(srsc, dims = 1:ndim)
	srsc = FindClusters(srsc, resolution = 0.5)

	srsc = RunUMAP(srsc, dims = 1:ndim)
	umapGG = DimPlot(srsc, reduction = "umap")
	ggsave(plot = umapGG, filename = paste(plotPf, "UMAP.png", sep = ""), dpi = res, width = 9, heigh = 6)

	outData = FetchData(srsc, c("ident", "nGene", "UMAP1", "UMAP2"))
	write.csv(outData, file = paste(plotPf, "useful_results.csv", sep = ""))
	
	if (length(levels(outData$ident)) >= 2) {
		srsc.markers = FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.0)
		topMarkers = srsc.markers %>% group_by(cluster) %>% top_n(n = 18, wt = avg_logFC)

		print(head(topMarkers))

		write.csv(srsc.markers, file = paste(plotPf, "ALL_POS_MARKERS.csv", sep = ""))
		write.csv(topMarkers, file = paste(plotPf, "top_pos_markers.csv", sep = ""))

		png(paste(plotPf, "marker_heatmap.png", sep = ""), res = res, width = 16, height = 12, units = "in")
		print(DoHeatmap(object = srsc, features = topMarkers$gene) + NoLegend())
	gar = dev.off()
	}

	if (!is.null(plotFeat)) {
		# TODO: Change width and heigh automatically
		cat("Start to plot feature plots and violin plots\n")
		vlnGG = VlnPlot(srsc, features = plotFeat, slot = "counts", log = TRUE, 
				ncol = 3, pt.size = 0.5)
		ggsave(plot = vlnGG, 
		       filename = paste(plotPf, "violin_plot_for_cell_annotation.png", sep = ""), 
		       dpi = res, width = 16, heigh = 12)

		featGG = FeaturePlot(srsc, features = plotFeat, ncol = 3, pt.size = 0.2, sort.cell = TRUE)
		ggsave(plot = featGG, 
		       filename = paste(plotPf, "umap_feature_plot_for_cell_annotation.png", sep = ""),
		       dpi = res, width = 16, heigh = 10)
	}

	if (rdsSave) {
		saveRDS(srsc, file = paste(plotPf, "Seurat_Objects_Clustered.RDS", sep = ""))
	}
	print(srsc)
	return(srsc)
}

inteSrPro <- function(srsc, plotPf, batchName = NULL, res = 300, ndim = 18, 
		      rdsSave = FALSE, plotFeat = NULL, jsFlag = FALSE, 
		      clusterRes = 0.5, atbFlag = FALSE, abtFeats = NULL) {
	suppressMessages(library(Seurat))
	suppressMessages(library(dplyr))
	suppressMessages(library(stringr))
	##############
	qcgg = VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
		       ncol = 3, group.by = batchName)
	ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), 
	       dpi = res, width = 9, height = 6)

	srsc = ScaleData(srsc)
	if (atbFlag) {
		srsc <- ScaleData(srsc, assay = "protein")
	}
	## TODO: Need to scale the protein data as well

	#############
	srsc = RunPCA(srsc)
	pcaGG1 = VizDimLoadings(srsc, dim = 1:6, reduction = "pca")
	ggsave(plot = pcaGG1, filename = paste(plotPf, "PCA1.png", sep = ""), 
	       dpi = res, width = 9, height = 9)
	
	pcaGG2 = DimPlot(srsc, reduction = "pca")
	ggsave(plot = pcaGG2, filename = paste(plotPf, "PCA2.png", sep = ""), 
	       dpi = res, width = 9, height = 6)

	png(paste(plotPf, "PCA_HEATMAP.png", sep = ""), res = res, width = 9, height = 16, units = "in")
	DimHeatmap(srsc, dims = 1:15, cells = 500, balanced = TRUE)
	gar = dev.off()

	if (jsFlag) {
		jst = Sys.time()
		cat("Running Jack Straw method to determine the PC dimension to use...\n")
		srsc = JackStraw(srsc, num.replicate = 100)
		srsc = ScoreJackStraw(srsc, dims = 1:20)
		jsGG = JackStrawPlot(srsc, dims = 1:20)
		ggsave(plot = jsGG, filename = paste(plotPf, "JackStraw.png", sep = ""), 
		       dpi = res, width = 9, heigh = 6)
		cat("Jack Straw costs")
		print(Sys.time()-jst)
	}

	elGG = ElbowPlot(srsc, ndims = 50)
	ggsave(plot = elGG, filename = paste(plotPf, "Elbow.png", sep = ""), dpi = res, width = 9, heigh = 6)

	cat("Start to cluster...\n")
	srsc = FindNeighbors(srsc, dims = 1:ndim)
	srsc = FindClusters(srsc, resolution = clusterRes)

	cat("Start to reduce dimension...\n")
	srsc = RunUMAP(srsc, dims = 1:ndim)
	srsc = RunTSNE(srsc, dims = 1:ndim)
	
	umapGG = DimPlot(srsc, reduction = "umap", label = TRUE)
	ggsave(plot = umapGG, filename = paste(plotPf, "UMAP.png", sep = ""), 
	       dpi = res, width = 9, heigh = 6)

	tsneGG = DimPlot(srsc, reduction = "tsne", label = TRUE)
	ggsave(plot = tsneGG, filename = paste(plotPf, "TSNE.png", sep = ""), 
	       dpi = res, width = 9, heigh = 6)

	if (atbFlag) {
		cat("Start to analyze protein data...\n")
		pst <- Sys.time()
		if (FALSE) {
		proteinMarkers <- FindAllMarkers(srsc, assay = "protein", only.pos = TRUE, min.pct = 0.0)
		write.csv(proteinMarkers, file = paste(plotPf, "ALL_POS_PROTEIN_MARKERS.csv", sep = ""))
		}
		allProteins <- srsc[["protein"]]@data
		print(dim(allProteins))
		protFeas <- rownames(allProteins)
		protFeats = paste0("protein_", rownames(srsc[["protein"]]))
#		print(protFeas)
		protFeatureGG = FeaturePlot(srsc, features = protFeats, reduction = "umap", 
					    min.cutoff = "q05", max.cutoff = "q95", 
					    ncol = 3, cols = c("orange1", "grey", "purple"), 
					    sort.cell = TRUE)
		ggsave(plot = protFeatureGG, filename = paste(plotPf, "TotalSeq_UMAP_overlay.png", sep = ""),
		       dpi = res, width = 24, heigh = 30)

		protDotGG = DotPlot(srsc, features = protFeas, 
				    cols = c("honeydew", "purple"), assay = "protein") + coord_flip()
		ggsave(plot = protDotGG, filename = paste(plotPf, "protein_dotplot.png", sep = ""), 
		       dpi = res, width = 9, heigh = 6)

		png(paste(plotPf, "protein_marker_heatmap.png", sep = ""), 
		    res = res, width = 12, height = 9, units = "in")
		print(DoHeatmap(object = srsc, features = protFeas, assay = "protein"))
		gar = dev.off()

		ridgeGG = RidgePlot(srsc, features = protFeas, assay = "protein", ncol = 4)
		ggsave(plot = ridgeGG, filename = paste(plotPf, "protein_marker_RidgePlot.png", sep = ""), 
		       dpi = res, width = 24, height = 30)
	}


	if (!is.null(plotFeat)) {
		cat("Start to plot feature plots and violin plots\n")
		vlnGG = VlnPlot(srsc, features = plotFeat, slot = "counts", 
				log = TRUE, ncol = 3, pt.size = 0.5)
		ggsave(plot = vlnGG, 
		       filename = paste(plotPf, "violin_plot_for_cell_annotation.png", sep = ""), 
		       dpi = res, width = 16, height = 30)

		featGG = FeaturePlot(srsc, features = plotFeat, ncol = 3, pt.size = 0.2, sort.cell = TRUE)
		ggsave(plot = featGG, 
		       filename = paste(plotPf, "umap_feature_plot_for_cell_annotation.png", sep = ""), 
		       dpi = res, width = 16, height = 27)

		dotGG = DotPlot(srsc, features = plotFeat) + RotatedAxis()
		ggsave(plot = dotGG, filename = paste(plotPf, "DotPlot_for_cell_annotation.png", sep = ""),
		       dpi = res, width = 9, height = 6)
	}

	if (!is.null(batchName)) {
		batchCheckGG = DimPlot(srsc, reduction = "umap", group.by = batchName)
		ggsave(plot = batchCheckGG, 
		       filename = paste(plotPf, "_batch_effect_check_dimPlot.png", sep = ""),
		       dpi = 180, width = 18, height = 6)
	}
	
	srsc.markers = FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
	topMarkers = srsc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

	write.csv(srsc.markers, file = paste(plotPf, "ALL_POS_MARKERS.csv", sep = ""))
	write.csv(topMarkers, file = paste(plotPf, "top_pos_markers.csv", sep = ""))

	png(paste(plotPf, "marker_heatmap.png", sep = ""), res = res, width = 16, height = 18, units = "in")
	print(DoHeatmap(object = srsc, features = topMarkers$gene) + NoLegend())
	gar = dev.off()

	topDotMarkers = srsc.markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_logFC)
	geneMarkers = unique(as.character(topDotMarkers$gene))
	geneMarkers = unique(c(as.character(topDotMarkers$gene), "CD14", "CD68", "HLA-DRA", "CD163", "STAB1"))

	dotMarkerGG = DotPlot(srsc, features = geneMarkers) + RotatedAxis() + coord_flip()
	ggsave(plot = dotMarkerGG, filename = paste(plotPf, "TopMarker_DotPlot.png", sep = ""), 
	       dpi = res, width = 9, height = 16)

	
	if (rdsSave) {
		saveRDS(srsc, file = paste(plotPf, "Seurat_Objects_Clustered.RDS", sep = ""))
	}
	print(srsc)
	return(srsc)
}

cellCompVis <- function(srsc, ctrlItem, itemOi, plotPf, res = 300) {
	cat("Visualize the composition of each cluster...\n")
	## TODO: normalize the cell counts to the total cell counts of its own category (tissue/patient)!
	ctsTable <- as.data.frame(table(srsc@meta.data[,c(ctrlItem, itemOi)]))
#	print(head(ctsTable))
#	print(dim(ctsTable))
	allComb <- combn(c(ctrlItem, itemOi), 2)
	print(allComb)
	for (i in 1:ncol(allComb)) {
		cat("\tPlotting for", allComb[1,i], allComb[2,i], "\n")
		tmpAbsGG <- ggplot(ctsTable, aes_string(x = allComb[1,i], y = "Freq", fill = allComb[2,i])) +
			labs(y = "Cell numbers") +
			geom_bar(stat = "identity", color = "white") +
			theme_classic()
#		ggsave(plot = tmpAbsGG, filename = paste(plotPf, allComb[2,i], "_abs_cell_count_", allComb[1,i], ".png", sep = ""),
#		       dpi = res, width = 9, height = 6)

		tmpRelGG <- ggplot(ctsTable, aes_string(x = allComb[1,i], y = "Freq", fill = allComb[2,i])) +
			labs(y = "Relative abundance") +
			geom_bar(stat = "identity", position = "fill", color = "white") +
			theme_classic()
#		ggsave(plot = tmpRelGG, filename = paste(plotPf, allComb[2,i], "_rel_cell_count_", allComb[1,i], ".png", sep = ""),
#		       dpi = res, width = 9, height = 6)

		combFig <- ggarrange(tmpAbsGG, tmpRelGG, nrow = 2, ncol = 1, 
				     common.legend = TRUE, legend = "right")
		ggexport(combFig, 
			 filename = paste(plotPf, allComb[2,i], "_comb_cell_count_", 
					  allComb[1,i], ".png", sep = ""),
			 res = res, width = 1200, height = 1200, verbose = FALSE)
	}
}

get_earliest_principal_node <- function(cds, clusterOi) {
	cell_ids <- which(colData(cds)[, "seurat_clusters"] == clusterOi)
	closet_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
	closet_vertex <- as.matrix(closet_vertex[colnames(cds),])
	keyV <- as.numeric(names(which.max(table(closet_vertex[cell_ids,]))))
	root_pr_nodes <- paste("Y", keyV, sep = "_")
	return(root_pr_nodes)
}


subMonocle3 <- function(cts, metadata, plotPf, alignBatch="None", pngres=300, rdsSave=FALSE, 
			pseudoTimeFlag=FALSE, rootCluster="None", vdjFlag = FALSE) {
	cat("Start to analyze with Monocle 3...\n")
	mncSt <- Sys.time()
	mncCDS <- new_cell_data_set(cts, cell_metadata = metadata)
	mncCDS <- preprocess_cds(mncCDS, num_dim = 100)
#	mncCDS <- preprocess_cds(mncCDS, num_dim = 100, method = "LSI")
	if (alignBatch != "None") {
		mncCDS <- align_cds(mncCDS, alignment_group = alignBatch)
	}
	mncCDS <- reduce_dimension(mncCDS)
	mncCDS <- cluster_cells(mncCDS)
	mncCDS <- learn_graph(mncCDS)
	
	patientGG <- plot_cells(mncCDS, color_cells_by = "patient", label_groups_by_cluster=FALSE, 
				label_leaves=FALSE, label_branch_points=FALSE)
	ggsave(plot = patientGG, filename = paste(plotPf, "trajectory_umap_with_patient.png", sep = ""), 
	       dpi = pngres, width = 6, height = 4)
	tissueGG <- plot_cells(mncCDS, color_cells_by = "tissue", label_groups_by_cluster=FALSE, 
			       label_leaves=FALSE, label_branch_points=FALSE)
	ggsave(plot = tissueGG, filename = paste(plotPf, "trajectory_umap_with_tissue.png", sep = ""), 
	       dpi = pngres, width = 6, height = 4)
	seuratGG <- plot_cells(mncCDS, color_cells_by = "seurat_clusters", 
			       label_groups_by_cluster = FALSE, label_leaves = TRUE, 
			       label_branch_points = TRUE, reduction_method = "UMAP") + 
		theme(legend.position="right")
	ggsave(plot = seuratGG, 
	       filename = paste(plotPf, "trajectory_umap_with_seurat_clusters.png", sep = ""), 
	       dpi = pngres, width = 6, height = 4)

	if (pseudoTimeFlag) {
		cat("Start to calculate pseudotime...\n")
		rootNodes <- get_earliest_principal_node(mncCDS, clusterOi = rootCluster)
		#clusterOimncMetaData <- colData(mncCDS)
		#colDatarootCells <- rownames(mncMetaData)[mncMetaData$seurat_clusters == 3]
		#mncMetaDataprint(head(mncMetaData))
		#headprint(length(rootCells))
		mncCDS <- order_cells(mncCDS, root_pr_nodes = rootNodes, verbose = TRUE)
		pseudoTGG <- plot_cells(mncCDS, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
					label_leaves = FALSE, label_branch_points = FALSE, 
					label_roots = TRUE, graph_label_size = 1.5)
		ggsave(plot = pseudoTGG, 
		       filename = paste(plotPf, "trajectory_umap_pseudotime.png", sep = ""), 
		       dpi = pngres, width = 9, height = 6)
	}
	if (vdjFlag) {
		cloneSizeGG <- plot_cells(mncCDS, color_cells_by = "clone_size", alpha = 0.45,
			       label_groups_by_cluster = FALSE, label_leaves = TRUE, 
			       label_cell_groups = FALSE, label_branch_points = TRUE, 
			       reduction_method = "UMAP") + 
		theme(legend.position="right")
		ggsave(plot = cloneSizeGG, 
		       filename = paste(plotPf, "trajectory_umap_with_clone_size.png", sep = ""), 
		       dpi = pngres, width = 6, height = 4)

		cloneSpaceGG <- plot_cells(mncCDS, color_cells_by = "clone_space", alpha = 0.45,
					   label_groups_by_cluster = FALSE, label_leaves = TRUE, 
					   label_cell_groups = FALSE, label_branch_points = TRUE, 
					   reduction_method = "UMAP") + 
		theme(legend.position="right")
		ggsave(plot = cloneSpaceGG, 
		       filename = paste(plotPf, "trajectory_umap_with_clone_space.png", sep = ""), 
		       dpi = pngres, width = 6, height = 4)

	}
	## TODO: Protein? Plot_gene_in_peudotime?
	if (rdsSave) {saveRDS(mncCDS, paste(plotPf, "monocle3_object_with_trajectory.RDS", sep = ""))}
	cat("Time cost\t")
	print(Sys.time()-mncSt)
	return(mncCDS)
}

tst = Sys.time()
suppressMessages(library(Seurat))
suppressMessages(library(monocle3))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(readxl))
suppressMessages(library(MAST))

dataDir <- "/home/weihua/mnts/data_plee/Group/scmodb/ALL10X/"
# rdsDir = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_results/10x_10312019/10X_tumor_11072019/TAM_results/10X_tumor_11072019_subsetted_seurat_object_c4_21.RDS"
# dataDir = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/triomics_attempt_BC394_BC401_042020/"
resDir = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_results/tex_tumor/"

# resDir <- "/home/weihua/mnts/group_plee/Weihua/scrnaseq_results/ALL10X_04032020/" ## For Colt all the samples
# resDir <- "/home/weihua/mnts/group_plee/Weihua/scrnaseq_results/AllThreePrimeRNAResults/" # For Lei, TAM

# sigFile = "/home/weihua/mnts/data_plee/Group/scmodb/signatures/lee_lab_signature_vis_list_v0.1.xlsx"
sigMarker <- c("PTPRC", "EPCAM", "COL1A1", "CD3D", "CD8A", "CD4", "PDCD1", "TCF7", "CD274", "MS4A1", "CD27", "CD19", 
	      "CD68", "FCER1A", "CD14", "NCAM1", "TPSAB1", "GNLY", "FOXP3", "HLA-DRA", "FCGR3A")
sigMarker <- c("PTPRC", "EPCAM", "COL1A1", "CD3D", "CD8A", "CD4", "PDCD1", "TCF7", "CD274", "MS4A1", "CD27", "CD19", 
	      "CD68", "FCER1A", "KRT19", "NCAM1", "TPSAB1", "GNLY", "FOXP3", "HLA-DRA", "FCGR3A")

geneOi <- c("CD3D", "CD8A", "CD4", "PDCD1", "CD28", "CTLA4", "ICOS", "BTLA", "LAG3", "HAVCR2", "TNFRSF9", "TNFRSF4", "CD274", "PDCD1LG2", "CD80", "SIGLEC15", "IL10", "CD44", "VDR", "RXRA", "RXRB", "ENTPD1", "CD69", "ITGAE")

expID <- "tex_tumor_integrate_comparison"
# expID <- "clean3prime_tumor_02262020" # For Lei
inDataFormat <- "10Xx"
inDataType <- "filtered"
saSelect <- c("Tumor") #c("ContralateralNormal", "IpsilateralNormal", "LymphNode")
rnaSelect <- NULL # c("five")
atbFlag <- TRUE # Antibody capture flag (TRUE, includes the totalseq data)
vdjFlag <- FALSE
# hashtagFlag <- TRUE
inteFLAG <- FALSE # Run integration with Seurat
inteFilter <- 100
saveListFlag <- TRUE # Save the input list for Seurat integration
postFlag <- FALSE # Run experiments for integrated objects!
visFlag <- FALSE # Visualize results
subClusterFlag <- FALSE # Re-cluster for each general cell types
tcrAnaFlag <- FALSE
monocle3TCR <- FALSE
specDEFlag <- FALSE
geneDEFlag <- TRUE
pngRes <- 300

dir.create(file.path(resDir, expID), showWarnings = FALSE)
expDir <- paste(resDir, expID, "/", sep = "")

git_scfolder <- "/home/weihua/git_repo/scrna_seq_pipeline"
rscript_file <- list.files(git_scfolder, "scrna_seurat_batch_v1_41.R")
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
#		print(names(tmpData))

		if (is.null(names(tmpData))) {
			## NOTE: need some basic filter to remove totally empty droplet!!!
			tmpObj <- CreateSeuratObject(counts = tmpData, project = saName[isa], 
						     min.cells = 3, min.features = 5)
			finalAtbFlag <- FALSE
		} else {
			tmpObj <- CreateSeuratObject(counts = tmpData$`Gene Expression`, 
						     project = saName[isa], min.cells = 3, min.features = 5)
			if (atbFlag) {
				cat("Including TotalSeq data\n")
				## TODO: Standarize the input antibody name!!!
				print(dim(tmpData$`Antibody Capture`))
				tmpRawRow <- rownames(tmpData$`Antibody Capture`)
#				print(tmpRawRow)
				tmpNewRow <- str_replace_all(tmpRawRow, " ", "_")
				tmpNewRow <- str_split_fixed(tmpNewRow, "man_", n = 2)[,2]
				tmpNewRow <- str_split_fixed(tmpNewRow, "_A", n = 2)[,1]
#				print(tmpNewRow)
				tmpNewRow <- str_split_fixed(tmpNewRow, "-", n = 2)[,1]
#				print(str_detect(tmpNewRow, "\\("))
				if (any(str_detect(tmpNewRow, "\\("))) {
					tmpNewRow <- str_split_fixed(tmpNewRow, "_\\(", n = 2)[,1]
				}
#				print(tmpNewRow)
				rownames(tmpData$`Antibody Capture`) <- tmpNewRow
				tmpObj[["protein"]] = CreateAssayObject(counts = tmpData$`Antibody Capture`)
			}
			finalAtbFlag <- atbFlag
		}
		print(tmpObj)
		if (vdjFlag) {
			## NOTE: contig annotation is REQUIRED to merge before adding to seurat data
			cat("TCR data is available!\n")
			tcrCsv = paste(dataDir, saName[isa],"_filtered_new_contig_annotations.csv", sep = "")
			rawTcrTable = read.csv(tcrCsv, row.names = 1, header = TRUE)
			tcrTable = rawTcrTable
			tcrTable = tcrTable[tcrTable$new_clonotype_id != "None", ]
			rownames(tcrTable) = tcrTable$barcode
			cellInRNA = rownames(tmpObj@meta.data)
			tmpCellID = str_split_fixed(cellInRNA, "-", n = 2)
			if (nrow(tmpCellID[tmpCellID[,2] != 1,]) > 0) {error("Mismatched TCRs!")}

			useTCRTable = tcrTable[, c("new_clonotype_id", "cdr3s_aa")]
			rownames(useTCRTable) = paste0(rownames(useTCRTable), "-1", sep = "")
			metaTCRTable = useTCRTable[cellInRNA,]
			rownames(metaTCRTable) = cellInRNA
			tmpObj@meta.data$clonotype_id = metaTCRTable$new_clonotype_id
			tmpObj@meta.data$cdr3s_aa = metaTCRTable$cdr3s_aa
#			print(head(metaTCRTable))
#			print(metaTCRTable[!is.na(metaTCRTable$cdr3s_aa),])
#			print(head(tmpObj@meta.data))
		}
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
		tmpProObj <- srQCPrePro(tmpObj, plotPf = plotPrefix, atbFlag = finalAtbFlag)
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
#	q(save = "no")

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

intePf <- paste(expDir, expID, "_Integrated_Results_", sep = "")
if (inteFLAG) {
	srCleanObjs <- readRDS(paste(expDir, expID, 
				     "_all_clean_filtered_10X_seurat_objects_list.RDS", sep = ""))
	#####################
	## Note: Batch correction seems to automated here!!!
	st <- Sys.time()
	cat("Start to merge the samples...\n")

	integrateList <- srCleanObjs
	for (inm in names(srCleanObjs)) {
		if(length(colnames(srCleanObjs[[inm]])) <= inteFilter) {integrateList[[inm]] <- NULL}
	}
	integrateAnchors <- FindIntegrationAnchors(object.list = integrateList, 
						   dims = 1:50, k.filter = inteFilter)
	inteObjs <- IntegrateData(anchorset = integrateAnchors, dims = 1:50)
	DefaultAssay(inteObjs) <- "integrated"
	print(inteObjs)
	saveRDS(inteObjs, file = paste(intePf, "_unprocessed_integrated_seurat_object.RDS", sep = ""))
	print(Sys.time() - st)

#	rnaMarkers <- paste0("rna_", sigMarker)
#	inteProObj <- inteSrPro(inteObjs, batchName = c("patient", "tissue"), plotPf = intePf, rdsSave = TRUE, plotFeat = rnaMarkers)
}

if (postFlag) {
	cat("Analyzing integrated Seurat object...\n")
	pst <- Sys.time()
	inteObjs <- readRDS(paste(intePf, "_unprocessed_integrated_seurat_object.RDS", sep = ""))
	rnaMarkers <- paste0("rna_", sigMarker)
	inteProObj <- inteSrPro(inteObjs, batchName = c("patient", "tissue"), 
				plotPf = intePf, rdsSave = TRUE, ndim = 20,
				plotFeat = rnaMarkers, jsFlag = TRUE, atbFlag = atbFlag)
	cat("Total")
	print(Sys.time()-pst)
}

if (visFlag) {# Highly customized for macrophage!
	cat("Visualize the processed Seurat object...\n")
	proObj <- readRDS(paste(intePf, "Seurat_Objects_Clustered.RDS", sep = ""))
	print(proObj)
	print(table(Idents(proObj)))
	print(table(Idents(proObj), proObj@meta.data$patient))
	ctsIdentPat <- as.data.frame.matrix(table(Idents(proObj), proObj@meta.data$patient))
	write.csv(ctsIdentPat, paste(intePf, "counts_idents_patient.csv", sep = ""))
}

if (subClusterFlag) {
	cat("Re-cluster for each general cell types...\n")
#	cellTypeOi <- c("Macrophage", "T cell", "Fibroblast")
	clsRes <- 0.8
#	cellTypeOi <- c("Macrophage_Res100")
	cellTypeOi <- c("CD8_T_cell", "CD4_T_cell")
	clusterOiList <- vector(mode = "list", length = length(cellTypeOi))
	names(clusterOiList) <- cellTypeOi

	## Define the cluster for each general cell type
	clusterOiList[["CD8_T_cell"]] <- c(0, 5, 15)
	clusterOiList[["CD4_T_cell"]] <- c(7)
#	clusterOiList[["B cell"]] <- c(14)
#	clusterOiList[["Fibroblast"]] <- c(1, 11)

	dimList <- vector(mode = "list", length = length(cellTypeOi))
	names(dimList) <- cellTypeOi
	dimList[["CD8_T_cell"]] <- 16
	dimList[["CD4_T_cell"]] <- 16
#	dimList[["B cell"]] <- 45
#	dimList[["Fibroblast"]] <- 12

	subInteFilter <- 30

	## Read the Seurat objects
	proObj <- readRDS(paste(intePf, "Seurat_Objects_Clustered.RDS", sep = ""))

	cellCompVis(srsc = proObj, ctrlItem = "seurat_clusters", 
		    itemOi = c("tissue", "sampleID", "patient"), plotPf = intePf)
	
	# Initialize the record matrix for each cell types
	colOi <- c("final_cell_num", "cell_cluster_num", "ndim", 
		   paste0(unique(proObj@meta.data$patient), "_cell_num"))
	cellTypeRecDf <- as.data.frame(matrix(ncol = length(colOi), nrow = length(cellTypeOi)))
	colnames(cellTypeRecDf) <- colOi
	rownames(cellTypeRecDf) <- cellTypeOi

	for (ict in cellTypeOi) {
		ist <- Sys.time()
		cat("Re-cluster", ict, ":", clusterOiList[[ict]], "\n")
		tmpSubSrsc <- subset(proObj, idents = clusterOiList[[ict]])

		tmpRSObj <- CreateSeuratObject(counts = tmpSubSrsc[["RNA"]]@counts, 
					       meta.data = tmpSubSrsc@meta.data)
		if (atbFlag) {
			tmpRSObj[["protein"]] <- CreateAssayObject(counts = tmpSubSrsc[["protein"]]@data)
		}
		tmpSeuratList <- SplitObject(tmpRSObj, split.by = "sampleID")
		print(tmpSeuratList)
		
		subFolder <- paste(ict, "_subtype_analysis", sep = "")
		dir.create(file.path(expDir, subFolder), showWarnings = FALSE)
		subDir <- paste(expDir, subFolder, "/", sep = "")
		subPrefix <- paste(subDir, expID, "_", ict, "_", sep = "")
		
		for (i in names(tmpSeuratList)) {
			cat("Patient", i, ":", length(colnames(tmpSeuratList[[i]])), "\n")
			###################
			tmpSeuratList[[i]] <- NormalizeData(tmpSeuratList[[i]], assay = "RNA", 
							    normalization.method = "LogNormalize", 
							    scale.factor = 10000)
			if (atbFlag) {
				tmpSeuratList[[i]] <- NormalizeData(tmpSeuratList[[i]], 
								    assay = "protein", 
								    normalization.method = "CLR")
			}
#			tmpSeuratList[[i]] <- NormalizeData(tmpSeuratList[[i]], verbose = FALSE)
			tmpSeuratList[[i]] <- FindVariableFeatures(tmpSeuratList[[i]], 
								   selection.method = "vst", 
								   nfeatures = 2000, verbose = FALSE)
			cellTypeRecDf[ict, paste(i, "_cell_num", sep = "")] <- length(colnames(tmpSeuratList[[i]]))
			if (length(colnames(tmpSeuratList[[i]])) <= subInteFilter) {
				## NOTE: remove the patients with few cells (<100 cells) for integration!
				tmpSeuratList[[i]] <- NULL
				cat("\tPatient", i, "is removed due to small cell number!!!\n")
			}
		}
		integrateAnchors <- FindIntegrationAnchors(object.list = tmpSeuratList, 
							   k.filter = subInteFilter)
		inteObjs <- IntegrateData(anchorset = integrateAnchors)
		DefaultAssay(inteObjs) <- "integrated"
		saveRDS(inteObjs, paste(subPrefix, "_Integrated_results.RDS", sep = ""))
		cellTypeRecDf[ict, "final_cell_num"] <- length(colnames(inteObjs))
		cat("Integration time cost:")
		print(Sys.time()-ist)

		pst <- Sys.time()
		# TODO: automatic select ndim from elbow?
		plotGenes <- paste0("rna_", geneOi)
		inteProObj <- inteSrPro(inteObjs, plotPf = subPrefix, rdsSave = TRUE, jsFlag = FALSE, 
					ndim = dimList[[ict]], batchName = c("tissue", "patient"), 
					plotFeat = plotGenes, clusterRes = clsRes, atbFlag = atbFlag)
		ctsIdentPat <- as.data.frame.matrix(table(Idents(inteProObj), inteProObj@meta.data$patient))
		write.csv(ctsIdentPat, paste(subPrefix, "counts_idents_patient.csv", sep = ""))
		cellTypeRecDf[ict, "cell_cluster_num"] <- length(unique(inteProObj@meta.data$seurat_clusters))
		write.csv(cellTypeRecDf, paste(subPrefix, "cell_type_analysis_record.csv", sep = ""))
		cat("Analysis time cost:")
		print(Sys.time()-pst)
		cat("++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
	}
}

if (geneDEFlag) {
	cat("Run differential expression between specific cohorts...\n")
	rawObj <- readRDS(paste(intePf, "Seurat_Objects_Clustered.RDS", sep = ""))
	proObj <- subset(rawObj, idents = c(1, 6, 8, 10, 16))
	intePf <- paste(intePf, "ER_only_", sep = "")
	print(proObj)
	## List for annotate PD-L1 status on TAM
	texTypes <- c("TEXhigh", "TEXlow")
	texAnnList <- vector(mode = "list", length = length(texTypes))
	names(texAnnList) <- texTypes
	texAnnList[["TEXhigh"]] <- c("BC258", "BC394", "BC397")
	texAnnList[["TEXlow"]] <- c("BC389", "BC392", "BC393", "BC401")

	proObj@meta.data$TEX <- "high"
	for (itt in names(texAnnList)) {
		tmpMask <- proObj@meta.data$patient %in% texAnnList[[itt]]
		proObj@meta.data[tmpMask, "TEX"] <- itt
	}

	proObj <- subset(proObj, subset = patient == "BC394", invert = TRUE)

	ifnrs <- c("IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2")
	ifnrGG <- FeaturePlot(proObj, features = ifnrs, split.by = "TEX", order = TRUE)
	ggsave(paste(intePf, "ifnrs_umap_split.png", sep = ""), ifnrGG,
	       dpi = pngRes, width = 9, height = 12)

	ifnrGG <- FeaturePlot(proObj, features = ifnrs, order = TRUE)
	ggsave(paste(intePf, "ifnrs_umap.png", sep = ""), ifnrGG,
	       dpi = pngRes, width = 9, height = 6)

	q(save = "no")

	exprOi <- FetchData(proObj, c("MX1", "IFI27"), slot = "data")
	gathExpr <- gather(exprOi, "gene", "counts", colnames(exprOi))
	histGG <- ggplot(gathExpr, aes_string(x = "counts", color = "gene")) +
		geom_histogram(fill = "white", alpha = 0.5, position = "identity", bins = 100) +
		theme_classic()
	ggsave(paste(intePf, "histogram_data.png", sep = ""), histGG,
	       dpi = pngRes, width = 9, height = 6)

	exprOi <- FetchData(proObj, c("MX1", "IFI27"), slot = "counts")
	gathExpr <- gather(exprOi, "gene", "counts", colnames(exprOi))
	histGG <- ggplot(gathExpr, aes_string(x = "counts", color = "gene")) +
		geom_histogram(fill = "white", alpha = 0.5, position = "identity", bins = 100) +
		theme_classic()
	ggsave(paste(intePf, "histogram_counts.png", sep = ""), histGG,
	       dpi = pngRes, width = 9, height = 6)


	proObj@meta.data$mx1 <- "N"
	proObj@meta.data$mx1[exprOi[,1] > 0] <- "MX1"
	proObj@meta.data$ifi27 <- "N"
	proObj@meta.data$ifi27[exprOi[,2] > 0] <- "IFI27"

	print(table(proObj@meta.data$TEX, proObj@meta.data$ifi27))
	print(table(proObj@meta.data$TEX, proObj@meta.data$mx1))

	cat("MX1\n")
	de_res <- FindMarkers(proObj, ident.1 = "MX1", ident.2 = "N", group.by = "mx1", test.use = "MAST", verbose = TRUE)
	write.csv(de_res, paste(intePf, "mx1_count_de_all.csv", sep = ""))
	sig_res <- de_res %>% filter(p_val_adj <= 0.10)

	cat("\tVisualize the top differential expression genes...\n")
	subObj <- proObj
	Idents(subObj) <- "mx1"
	cat("\t\tDot plot...\n")
	top_features <- c(rownames(sig_res %>% filter(avg_logFC > 0) %>% top_n(20, avg_logFC)), 
			  rownames(sig_res %>% filter(avg_logFC < 0)%>% top_n(-20, avg_logFC)))
	dot_gg <- DotPlot(subObj, features = top_features, cols = c("goldenrod", "purple"), 
			  dot.scale = 8, split.by = "mx1") + RotatedAxis()
	ggsave(paste(intePf, "mx1_count_dotplot_de_top20.png", sep = ""), dot_gg, dpi = pngRes, width = 12, height = 6)
	cat("\t\tScatter plot...\n")
	top_less_features <- c(rownames(sig_res %>% filter(avg_logFC > 0) %>% top_n(5, avg_logFC)), 
			       rownames(sig_res %>% filter(avg_logFC < 0)%>% top_n(-5, avg_logFC)))
	avg_expr <- log1p(AverageExpression(subObj, verbose = TRUE)$RNA)
	avg_expr$logFC <- NA
	for (ir in rownames(de_res)) {
		avg_expr[ir, "logFC"] <- de_res[ir, "avg_logFC"]
	}
	avg_expr <- avg_expr %>% arrange(desc(is.na(logFC)))
	sct_gg <- ggplot(avg_expr, aes(x = MX1, y = N)) + 
		geom_point(aes(color = logFC)) + 
		theme_classic() + 
		scale_color_viridis_c() +
		ggtitle("MX1")
	sct_gg <- LabelPoints(plot = sct_gg, points = top_less_features, repel = TRUE)
	ggsave(paste(intePf, "mx1_count_scatter_de_top5.png", sep = ""), sct_gg, dpi = pngRes, width = 4.5, height = 4.5)

	cat("IFI27\n")
	de_res <- FindMarkers(proObj, ident.1 = "IFI27", ident.2 = "N", group.by = "ifi27", test.use = "MAST", verbose = TRUE)
	write.csv(de_res, paste(intePf, "ifi27_count_de_all.csv", sep = ""))
	sig_res <- de_res %>% filter(p_val_adj <= 0.10)

	cat("\tVisualize the top differential expression genes...\n")
	subObj <- proObj
	Idents(subObj) <- "ifi27"
	cat("\t\tDot plot...\n")
	top_features <- c(rownames(sig_res %>% filter(avg_logFC > 0) %>% top_n(20, avg_logFC)), 
			  rownames(sig_res %>% filter(avg_logFC < 0)%>% top_n(-20, avg_logFC)))
	dot_gg <- DotPlot(subObj, features = top_features, cols = c("goldenrod", "purple"), 
			  dot.scale = 8, split.by = "ifi27") + RotatedAxis()
	ggsave(paste(intePf, "ifi27_count_dotplot_de_top20.png", sep = ""), dot_gg, dpi = pngRes, width = 12, height = 6)
	cat("\t\tScatter plot...\n")
	top_less_features <- c(rownames(sig_res %>% filter(avg_logFC > 0) %>% top_n(5, avg_logFC)), 
			       rownames(sig_res %>% filter(avg_logFC < 0)%>% top_n(-5, avg_logFC)))
	avg_expr <- log1p(AverageExpression(subObj, verbose = TRUE)$RNA)
	avg_expr$logFC <- NA
	for (ir in rownames(de_res)) {
		avg_expr[ir, "logFC"] <- de_res[ir, "avg_logFC"]
	}
	avg_expr <- avg_expr %>% arrange(desc(is.na(logFC)))
	sct_gg <- ggplot(avg_expr, aes(x = IFI27, y = N)) + 
		geom_point(aes(color = logFC)) + 
		theme_classic() + 
		scale_color_viridis_c() +
		ggtitle("IFI27")
	sct_gg <- LabelPoints(plot = sct_gg, points = top_less_features, repel = TRUE)
	ggsave(paste(intePf, "ifi27_count_scatter_de_top5.png", sep = ""), sct_gg, dpi = pngRes, width = 4.5, height = 4.5)

}

if (specDEFlag) {
	## Run DE analysis between PD-L1+/PD-L1- TAMs (Highly customized!)
	cat("Run differential expression between specific cohorts...\n")
	geneOi <- c("IFIT3", "IFI44L", "IFI44", "ISG15", "MX1", "IFIT1", "IFI6", 
		    "IFIT2", "IFI27", "STAT1", "CXCL10", "CXCL9", "OAS1", "OASL", "MX1")
	## Read the Seurat objects
	proObj <- readRDS(paste(intePf, "Seurat_Objects_Clustered.RDS", sep = ""))

	## List for annotate PD-L1 status on TAM
	texTypes <- c("TEXhigh", "TEXlow")
	texAnnList <- vector(mode = "list", length = length(texTypes))
	names(texAnnList) <- texTypes
	texAnnList[["TEXhigh"]] <- c("BC258", "BC394", "BC397")
	texAnnList[["TEXlow"]] <- c("BC389", "BC392", "BC393", "BC401")

	proObj@meta.data$TEX <- "high"
	for (itt in names(texAnnList)) {
		tmpMask <- proObj@meta.data$patient %in% texAnnList[[itt]]
		proObj@meta.data[tmpMask, "TEX"] <- itt
	}
	proObj@meta.data$tumor <- "Non-tumor cells"
	tumor_mask <- proObj@meta.data$seurat_clusters %in% c(1, 6, 8, 10, 16)
	proObj@meta.data[tumor_mask, "tumor"] <- "Tumor cells"
	tumorFGG <- FeaturePlot(proObj, features = geneOi[1:7], sort.cell = T, split.by = "tumor")
	ggsave(paste(intePf, "umap_ifn_gene_tumor_only_part1.png", sep = ""), tumorFGG,
	       dpi = pngRes, width = 9, height = 21)
	tumorFGG <- FeaturePlot(proObj, features = geneOi[8:14], sort.cell = T, split.by = "tumor")
	ggsave(paste(intePf, "umap_ifn_gene_tumor_only_part2.png", sep = ""), tumorFGG,
	       dpi = pngRes, width = 9, height = 21)



	## Subset 
	rawObj <- proObj
	proObj <- subset(rawObj, idents = c(1, 6, 8, 10, 16))
	proObj <- subset(proObj, subset = patient == "BC394", invert = TRUE)
	intePf <- paste(intePf, "ER_only_", sep = "")
	texCol <- c("dodgerblue", "firebrick")

	cat("Vln plot visualization for gene of interest...\n")
	vlnGG <- VlnPlot(proObj, features = geneOi, group.by = "TEX", ncol = 5,
			 pt.size = 0.01, cols = c("firebrick", "dodgerblue"))
	ggsave(paste(intePf, "vlnplot_tex_high_vs_low_tumor_cells.png", sep = ""), vlnGG,
	       dpi = pngRes, width = 9, height = 7.5)
	print(proObj)

	tumorFGG <- FeaturePlot(proObj, features = geneOi[1:7], sort.cell = T, split.by = "TEX")
	ggsave(paste(intePf, "umap_ifn_gene_tumor_only_tex_part1.png", sep = ""), tumorFGG,
	       dpi = pngRes, width = 9, height = 21)
	tumorFGG <- FeaturePlot(proObj, features = geneOi[8:14], sort.cell = T, split.by = "TEX")
	ggsave(paste(intePf, "umap_ifn_gene_tumor_only_tex_part2.png", sep = ""), tumorFGG,
	       dpi = pngRes, width = 9, height = 21)

	exprDf <- as.matrix(proObj[["RNA"]]@data[geneOi,])
	exprDf <- cbind(t(exprDf), proObj@meta.data)
	gathDf <- gather(exprDf, "gene", "expr", geneOi)

	sumGathDf <- gathDf %>% 
		group_by(gene, TEX) %>%
		summarise(
			  sd = sd(expr, na.rm=T),
			  expr = mean(expr),
			  n = n()
			  )
	sumGathDf$sem <- sumGathDf$sd/(sumGathDf$n)^0.5
#	print(head(gathDf))
	print(head(sumGathDf))

	barPtGG <- ggplot(gathDf, aes_string(x = "TEX", y = "expr")) +
		geom_bar(stat = "identity", data = sumGathDf, fill = NA, aes_string(color = "TEX")) +
		geom_jitter(position = position_jitter(0.2), aes_string(color = "TEX"), size = 0.2, alpha = 0.6) +
		geom_errorbar(aes(ymin = expr-sem, ymax = expr+sem), data = sumGathDf, width = 0.2) +
		scale_color_manual(values = c("firebrick", "dodgerblue")) +
		stat_compare_means(data = gathDf, mapping = aes_string(group = "TEX"), 
				   label = "p.signif", label.x = 1.35, vjust = 1, method = "wilcox.test") +
		labs(y = "Expression levels", color = "Group") +
		facet_wrap(~ gene, ncol = 7, scales = "free") +
		theme_bw() +
		theme(axis.title.x = element_blank(),
		      legend.position = "top")
	ggsave(paste(intePf, "manual_barplot_tex_high_vs_low_tumor_cells.png", sep = ""), barPtGG,
	       dpi = pngRes, width = 9, height = 6)


	boxGG <- ggplot(gathDf, aes_string(x = "TEX", y = "expr")) +
		geom_boxplot(fill = NA, aes_string(color = "TEX"), outlier.size = 0.5, outlier.shape = 16) +
		geom_jitter(position = position_jitter(0.2), aes_string(color = "TEX"), size = 0.2, alpha = 0.6) +
		scale_color_manual(values = c("firebrick", "dodgerblue")) +
		stat_compare_means(data = gathDf, mapping = aes_string(group = "TEX"), 
				   label = "p.signif", label.x = 1.35, vjust = 1, method = "wilcox.test") +
		labs(y = "Expression levels", color = "Group") +
		facet_wrap(~ gene, ncol = 7, scales = "free") +
		theme_bw() +
		theme(axis.title.x = element_blank(),
		      legend.position = "top")
	ggsave(paste(intePf, "manual_boxplot_tex_high_vs_low_tumor_cells.png", sep = ""), boxGG,
	       dpi = pngRes, width = 9, height = 6)

	vlnGG <- ggplot(gathDf, aes_string(x = "TEX", y = "expr")) +
		geom_violin(fill = NA, aes_string(color = "TEX"), trim = FALSE) +
		geom_jitter(position = position_jitter(0.2), aes_string(color = "TEX"), size = 0.2, alpha = 0.6) +
		scale_color_manual(values = c("firebrick", "dodgerblue")) +
		stat_compare_means(data = gathDf, mapping = aes_string(group = "TEX"), 
				   label = "p.signif", label.x = 1.35, vjust = 1, method = "wilcox.test") +
		labs(y = "Expression levels", color = "Group") +
		facet_wrap(~ gene, ncol = 7, scales = "free") +
		theme_bw() +
		theme(axis.title.x = element_blank(),
		      legend.position = "top")
	ggsave(paste(intePf, "manual_vlnplot_tex_high_vs_low_tumor_cells.png", sep = ""), vlnGG,
	       dpi = pngRes, width = 9, height = 6)

	q(save = "no")

	cat("Start to run DE analysis -- Wilcox\n")
	deRes <- FindMarkers(proObj, ident.1 = "TEXhigh", ident.2 = "TEXlow", group.by = "TEX", test.use = "wilcox")
	write.csv(deRes, paste(intePf, "post_wilcox_high_vs_low_TEX_res.csv", sep = ""))
	deRes <- deRes[order(deRes$avg_logFC),]
	sigMask <- deRes$p_val_adj <= 0.10 &(deRes$avg_logFC >= 0.5 | deRes$avg_logFC <= -0.5)
	sigFeat <- rownames(deRes)[sigMask]
	png(paste(intePf, "post_wilcox_TEX_high_low_marker_heatmap.png", sep = ""), res = 300, width = 9, height = 6, units = "in")
	print(DoHeatmap(object = proObj, features = sigFeat, group.by = "TEX", label = FALSE, 
			group.colors = texCol))
	gar <- dev.off()

	cat("Start to run DE analysis -- MAST\n")
	deRes <- FindMarkers(proObj, ident.1 = "TEXhigh", ident.2 = "TEXlow", group.by = "TEX", test.use = "MAST")
	write.csv(deRes, paste(intePf, "post_mast_high_vs_low_TEX_res.csv", sep = ""))
	deRes <- deRes[order(deRes$avg_logFC),]
	sigMask <- deRes$p_val_adj <= 0.10 &(deRes$avg_logFC >= 0.5 | deRes$avg_logFC <= -0.5)
	sigFeat <- rownames(deRes)[sigMask]
	png(paste(intePf, "post_mast_TEX_high_low_marker_heatmap.png", sep = ""), res = 300, width = 9, height = 6, units = "in")
	print(DoHeatmap(object = proObj, features = sigFeat, group.by = "TEX", label = FALSE, 
			group.colors = texCol))
	gar <- dev.off()

	proObj <- subset(rawObj, idents = c(1, 6, 8, 16))
	proObj <- subset(proObj, subset = patient == "BC394", invert = TRUE)
	texCol <- c("dodgerblue", "firebrick")

	cat("Start to run DE analysis -- Wilcox\n")
	deRes <- FindMarkers(proObj, ident.1 = "TEXhigh", ident.2 = "TEXlow", group.by = "TEX", test.use = "wilcox")
	write.csv(deRes, paste(intePf, "wo10_post_wilcox_high_vs_low_TEX_res.csv", sep = ""))
	deRes <- deRes[order(deRes$avg_logFC),]
	sigMask <- deRes$p_val_adj <= 0.10 &(deRes$avg_logFC >= 0.5 | deRes$avg_logFC <= -0.5)
	sigFeat <- rownames(deRes)[sigMask]
	png(paste(intePf, "wo10_post_wilcox_TEX_high_low_marker_heatmap.png", sep = ""), res = 300, width = 9, height = 6, units = "in")
	print(DoHeatmap(object = proObj, features = sigFeat, group.by = "TEX", label = FALSE, 
			group.colors = texCol))
	gar <- dev.off()

	cat("Start to run DE analysis -- MAST\n")
	deRes <- FindMarkers(proObj, ident.1 = "TEXhigh", ident.2 = "TEXlow", group.by = "TEX", test.use = "MAST")
	write.csv(deRes, paste(intePf, "wo10_post_mast_high_vs_low_TEX_res.csv", sep = ""))
	deRes <- deRes[order(deRes$avg_logFC),]
	sigMask <- deRes$p_val_adj <= 0.10 &(deRes$avg_logFC >= 0.5 | deRes$avg_logFC <= -0.5)
	sigFeat <- rownames(deRes)[sigMask]
	png(paste(intePf, "wo10_post_mast_TEX_high_low_marker_heatmap.png", sep = ""), res = 300, width = 9, height = 6, units = "in")
	print(DoHeatmap(object = proObj, features = sigFeat, group.by = "TEX", label = FALSE, 
			group.colors = texCol))
	gar <- dev.off()

}

if (tcrAnaFlag) {
	cat("Start to analyze the TCR data for the Seurat object...\n")
	cellTypeOi <- c("CD8_T_cell", "CD4_T_cell")
	for (ict in cellTypeOi) {
		subFolder <- paste(ict, "_subtype_analysis", sep = "")
		subDir <- paste(expDir, subFolder, "/", sep = "")
		subPrefix <- paste(subDir, expID, "_", ict, "_", sep = "")
		cat("Reading Seurat object...")
		st = Sys.time()
		proObj = readRDS(paste(subPrefix, "Seurat_Objects_Clustered.RDS", sep = ""))
		print(proObj)
		proObj@meta.data$tcrFlag = "No TCR"
		proObj@meta.data$tcrFlag[!is.na(proObj@meta.data$clonotype_id)] = "TCR"
		proObj@meta.data$pid_clonotype = paste0(proObj@meta.data$patient, "_", 
							proObj@meta.data$clonotype_id, sep = "")
		proObj@meta.data$clone_size = 0
		proObj@meta.data$total_top_prop = "None"
		proObj@meta.data$clone_space = "None"

		print(head(proObj@meta.data))
		tcrFolder = paste(ict, "tcr_results", sep = "_")
		dir.create(file.path(subDir, tcrFolder), showWarnings = FALSE)
		tcrDir = paste(subDir, tcrFolder, "/", sep = "")
		plotPrefix = paste(tcrDir, expID, "_", ict, "_", sep = "")

		for (ipat in unique(proObj@meta.data$patient)) {
			cat("\t", ipat, "\n")
			print(Sys.time()-st)
			# Put all the plots and results for TCR in one folder
			tmpMetaData <- proObj@meta.data[proObj@meta.data$patient == ipat,]
			## For UMAP with TCR data
			avaCtTable = as.data.frame.matrix(table(tmpMetaData$pid_clonotype, 
								tmpMetaData$tissue)) 
			# Available clone type table
			avaCtTable$Total = rowSums(avaCtTable)
			avaCtTable$Prop2Sum = avaCtTable$Total/sum(avaCtTable$Total)
			## Top proportions:
			powerBase = 5
			cateThr = 5^(0:ceiling(nrow(avaCtTable)^(1/powerBase)))
		#	print(cateThr)
			avaCtTable = avaCtTable[order(avaCtTable$Total, decreasing = TRUE),]
			for (icate in 2:length(cateThr)) {
				if (icate == 2) {
					bg = cateThr[icate-1]
					ed = cateThr[icate]
				} else {
					bg = cateThr[icate-1]+1
					ed = cateThr[icate]
				}
				if (icate == length(cateThr)) {
					bg = cateThr[icate-1]+1
					ed = nrow(avaCtTable)
				}
		#		cat(head(bg:ed), "\n")
				avaCtTable[bg:ed, "totalTop"] = paste("Top group ", icate-1, 
								      ": Top", bg, "~", "Top", ed, sep = "")
			}
			## Clonal space homeostasis:
			spcThr = 10^(floor(log(min(avaCtTable$Prop2Sum),10)):0)
			for (ispc in 2:length(spcThr)) {
				bg = spcThr[ispc-1]
				ed = spcThr[ispc]
				tmpMask = avaCtTable$Prop2Sum > bg & avaCtTable$Prop2Sum <= ed
				avaCtTable[tmpMask, "cloneSpace"] = paste("Space ", ispc-1, ": ", 
									  bg, "~", ed, sep = "")
			}
#			print(head(avaCtTable))
			## Map the TCR parameters back to each cell
			for (irow in rownames(avaCtTable)) {
				tmpMask = proObj@meta.data$pid_clonotype == irow & !is.na(proObj@meta.data$clonotype_id)
				proObj@meta.data[tmpMask, "clone_size"] = avaCtTable[irow, "Total"]
				proObj@meta.data[tmpMask, "total_top_prop"] = avaCtTable[irow, "totalTop"]
				proObj@meta.data[tmpMask, "clone_space"] = avaCtTable[irow, "cloneSpace"]
			}
		}
		## Per cell type analysis
		saveRDS(proObj, paste(subPrefix, "Seurat_clustered_TCR_analyzed.RDS", sep = ""))
		
		cat("Visualize TCR data on 2D dimension...\n")
		print(head(proObj@meta.data))
		proObj@meta.data$clone_size = as.numeric(as.character(proObj@meta.data$clone_size))

		pngdpi <- 300
		pngW <- 12
		pngH <- 4
		gnlTCRGG = DimPlot(proObj, reduction = "umap", group.by = "tcrFlag", split.by = "patient",
				   cols = c("grey", "#5BC4FF"))
		ggsave(plot = gnlTCRGG, filename = paste(plotPrefix, "_tcrFlag_dimPlot.png", sep = ""),
		       dpi = pngdpi, width = pngW, height = pngH)

		totalPropGG = DimPlot(proObj, reduction = "umap", group.by = c("total_top_prop"),
				      split.by = "patient")
		ggsave(plot = totalPropGG, 
		       filename = paste(plotPrefix, "_tcr_total_prop_dimPlot.png", sep = ""),
		       dpi = pngdpi, width = pngW, height = pngH)

		cloneSpcGG = DimPlot(proObj, reduction = "umap", group.by = c("clone_space"), 
				     split.by = "patient", cols = "YlOrRd")
		ggsave(plot = cloneSpcGG, 
		       filename = paste(plotPrefix, "_tcr_clone_space_dimPlot.png", sep = ""),
		       dpi = pngdpi, width = pngW, height = pngH)

		cloneSizeGG = FeaturePlot(proObj, reduction = "umap", features = "clone_size", 
					  split.by="patient", cols = c("#E0E0E0", "#4C0099"), 
					  sort.cell = TRUE)
		ggsave(plot = cloneSizeGG, 
		       filename = paste(plotPrefix, "_clone_size_dimPlot.png", sep = ""),
		       dpi = pngdpi, width = pngW, height = pngH)
	}
}

if (monocle3TCR) {
	cat("Start to analyze the TCR data for the Seurat object...\n")
	cellTypeOi <- c("CD8_T_cell", "CD4_T_cell")
	rootList <- vector(mode = "list", length = length(cellTypeOi))
	names(rootList) <- cellTypeOi
	rootList[["CD8_T_cell"]] <- 0
	rootList[["CD4_T_cell"]] <- 0
	ctmnList <- vector(mode = "list", length = length(cellTypeOi))
	names(ctmnList) <- cellTypeOi
	ctmnList[["CD8_T_cell"]] <- c(11, 13)	
	ctmnList[["CD4_T_cell"]] <- c(2, 7)

	for (ict in cellTypeOi) {
		ist <- Sys.time()
		subFolder <- paste(ict, "_subtype_analysis", sep = "")
		subDir <- paste(expDir, subFolder, "/", sep = "")
		subPrefix <- paste(subDir, expID, "_", ict, "_", sep = "")
		cat("Reading Seurat object...")
		st = Sys.time()
		proObj = readRDS(paste(subPrefix, "Seurat_clustered_TCR_analyzed.RDS", sep = ""))
		print(proObj)
		print(table(proObj@meta.data$seurat_clusters, proObj@meta.data$tissue))
		print(table(proObj@meta.data$tissue, proObj@meta.data$patient))
		q(save = "no")
		cleanProObj <- subset(proObj, idents = ctmnList[[ict]], invert=TRUE)
		trajFolder <- paste(ict, "trajectory_analysis", sep = "_")
		dir.create(file.path(subDir, trajFolder), showWarnings = FALSE)
		trajDir <- paste(subDir, trajFolder, "/", sep = "")
		trajPrefix <- paste(trajDir, expID, "_", ict, "_", sep = "")
		## NOTE: Reset the object!!!
		tmpCtsMat <- as.matrix(cleanProObj[["RNA"]]@counts)
		tmpMNCObj <- subMonocle3(cts = tmpCtsMat, metadata = cleanProObj@meta.data, 
					 plotPf = trajPrefix, alignBatch = "patient", rdsSave = TRUE, 
					 pseudoTimeFlag=TRUE, rootCluster = rootList[[ict]], 
					 vdjFlag = vdjFlag)
		cat(ict, "cost")
		print(Sys.time()-ist)
		cat("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n")
	}
}




cat("Total time cost:")
print(Sys.time()-tst)

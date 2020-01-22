# Single-cell RNA sequencing data standard processing procedure
# Weihua Guo, Ph.D.
# Start date: 12/02/2019
# NOTE: 01/21/2018 Add cell filter function, count-meta.data cell alignment, 
# Functions: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

rm(list = ls())

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

dataDir = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/scRNAseqTex/Tumor_LN_Normal_PBMC_Colorectal_Melanoma_scImpute_5/"
# ctsFile = "Tumor_LN_extracted_raw_counts.txt"
ctsFile = "scimpute_count.txt"
cellAnnFile = "Tumor_LN_Normal_PBMC_Colorectal_Melanoma_cell_annotation.txt"
facsTagFile = "FACS_TAG_ALL.txt"
resDir = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_results/scRNAseqTexResults"

# expID = "tumor_crc_melanoma_scimp_k5_rmkeep_proenc_custMT_ndim_10"
expID = "new_scimpute_breast_tissue_only"
tisFilter = c("Tumor", "Normal", "LN", "PBMC")
mtThr = 20
custMTGeneFlag = TRUE # TRUE: use manually refined mitochondrial genes to calculate percent_mt
rmFlag = TRUE # TRUE: remove ribosome genes and mitochondrial genes
keepFlag = TRUE # TRUE: only use characterized protein-encoding genes
tifres = 180
ndim = 12
nFeat = 1000
jsFlag = TRUE
dir.create(file.path(resDir, expID), showWarnings = TRUE)
expDirPre = paste(resDir, "/", expID, "/", expID, "_", sep = "")

##########################################################################
# Save this Rscript to the result folder
git_scfolder = "/home/weihua/git_repo/TEXinERpBC"
rscript_file = list.files(git_scfolder, "scProcess.R")
file.copy(rscript_file, paste(resDir, "/", expID, "/", sep = ""))

cat("Start to read counts...\n")
ctsDf = read.table(paste(dataDir, ctsFile, sep = ""), header = TRUE, row.names = 1, sep = " ")
cellAnnDf = read.table(paste(dataDir, cellAnnFile, sep = ""), header = TRUE, row.names = 1, sep = " ",
		       stringsAsFactors = FALSE)
facsDf = read.table(paste(dataDir, facsTagFile, sep = ""), header = TRUE, row.names = 1, sep = "\t")

## TODO: move this out side finally!!!
if (is.null(tisFilter)) {
	cat("Use all the tissue:", unique(cellAnnDf$tissue),"\n")
	cellAnnDf = cellAnnDf[complete.cases(cellAnnDf$tissue),]
	ctsDf = ctsDf[,colnames(ctsDf) %in% rownames(cellAnnDf)]
	ctsDf = ctsDf[,rownames(cellAnnDf)]
} else {
	cat("Only use the following tissue:", tisFilter, "\n")
	cellAnnDf = cellAnnDf[complete.cases(cellAnnDf$tissue),]
	cellAnnDf = cellAnnDf[cellAnnDf$tissue %in% tisFilter,]
	ctsDf = ctsDf[,colnames(ctsDf) %in% rownames(cellAnnDf)]
	ctsDf = ctsDf[,rownames(cellAnnDf)]
}

## Add flow-cytometry tag
facsDf = facsDf[rownames(facsDf) %in% colnames(ctsDf),]
singleFACS = colnames(facsDf)
# PD1-CD39 double positive
facsDf$PD1CD39 = facsDf$PD1 + facsDf$CD39
facsDf[facsDf$PD1CD39 == 1,"PD1CD39"] = 0
facsDf[facsDf$PD1CD39 == 2,"PD1CD39"] = 1

metaDf = cbind(facsDf, cellAnnDf[,c("tissue", "plate")])
metaDf[] = lapply(metaDf, factor)

srsc = CreateSeuratObject(counts = ctsDf, project = expID, min.cells = 3, min.features = 200, meta.data = metaDf)
# NOTE: before adding meta.data, align the rows of meta.data to ctsDf is NECCESSARY
# srsc@meta.data$tissue = as.factor(cellAnnDf[rownames(srsc@meta.data), "tissue"]) 
# srsc@meta.data$plate = as.factor(cellAnnDf[rownames(srsc@meta.data), "plate"])
cat("Raw Seurat Object:")
print(srsc)
# stop("TT")

if (custMTGeneFlag) {
	cat("Use manually refined mitochondrial gene list to calculate percentage of mitochondrial gene expression\n")
	mtGeneFile = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/scRNAseqTex/mt_gene_list.txt"
	mtGenes = scan(mtGeneFile, what="", sep="\n", quiet = TRUE)
	avaMTGenes = intersect(rownames(srsc[["RNA"]]@data), mtGenes)
	srsc[["percent.mt"]] = PercentageFeatureSet(srsc, features = avaMTGenes)
} else {
	cat("Use genes starting with MT- to calculate percentage of mitochondrial gene expression\n")
	srsc[["percent.mt"]] = PercentageFeatureSet(srsc, pattern = "^MT-")
}

qcVln = VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plot = qcVln, filename = paste(expDirPre, "QC_Vln_1.tiff", sep = ""), width = 9, height = 6, dpi = tifres)

plot1 = FeatureScatter(srsc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(srsc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
qcComb = CombinePlots(plots = list(plot1, plot2))
ggsave(plot = qcComb, filename = paste(expDirPre, "QC_Correlation_2.tiff", sep = ""), 
       width = 12, height = 6, dpi = tifres)

ctsData = srsc[["RNA"]]@data
metaData = srsc@meta.data[,c("plate", "tissue", "percent.mt")]

if (rmFlag) {
	cat("Remove the non-informative genes (ribosomal or mitochondrial genes...)\n")
	# NOTE: Key QC step to remove dead or dying cells
	srsc = subset(srsc, subset = percent.mt < 20)
	##########################################################################
	# Pre_filter
	rmGeneFile = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/scRNAseqTex/remove_list.txt"
	rmGenes = scan(rmGeneFile, what="", sep="\n", quiet = TRUE)
	
	olRmGenes = intersect(rownames(ctsData), rmGenes)
	keepGenes = setdiff(rownames(ctsData), olRmGenes)
	cleanCtsData = ctsData[keepGenes,]
	
	srsc = CreateSeuratObject(counts = cleanCtsData, project = expID, min.cells = 3, min.features = 200, 
				  meta.data = metaData)
	srsc = subset(srsc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
} else {
	# NOTE: Key QC step to remove low quality cells
	srsc = subset(srsc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
}

if (keepFlag) {
	cat("Only keep the characterized protein-encoding genes for further analysis...\n")
	# NOTE: Key QC step to remove dead or dying cells
	if (!rmFlag) {
		srsc = subset(srsc, subset = percent.mt < 20)
	} else {
		ctsData = srsc[["RNA"]]@data
		metaData = srsc@meta.data[,c("plate", "tissue", "percent.mt")]
	}
	##########################################################################
	# Pre_filter
	kpGeneFile = "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/scRNAseqTex/keep_protein_coding_gene.txt"
	kpGenes = scan(kpGeneFile, what="", sep="\n", quiet = TRUE)
	
	olKPGenes = intersect(rownames(ctsData), kpGenes)
	cleanCtsData = ctsData[olKPGenes,]
	
	srsc = CreateSeuratObject(counts = cleanCtsData, project = expID, min.cells = 3, min.features = 200, 
				  meta.data = metaData)
	srsc = subset(srsc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
} else {
	# NOTE: Key QC step to remove low quality cells
	srsc = subset(srsc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
}

cat("After data clean:")
print(head(colnames(srsc)))
srsc[["FACS"]] = CreateAssayObject(counts = t(facsDf[colnames(srsc),]))
print(srsc)
print(head(srsc@meta.data))

# NOTE: Normalization
srsc = NormalizeData(srsc, normalization.method = "LogNormalize", scale.factor = 10000)
# NOTE: Identify highly variable features
srsc = FindVariableFeatures(srsc, selection.method = "vst", nfeatures = 1000)
topVarGene = head(VariableFeatures(srsc), 18)
plot1 = VariableFeaturePlot(srsc)
plot2 = LabelPoints(plot = plot1, points = topVarGene, repel = TRUE)
varFeatPlot = CombinePlots(plots = list(plot1, plot2))
ggsave(plot = varFeatPlot, filename = paste(expDirPre, "top_variated_genes.tiff", sep = ""), 
       width = 12, height = 6, dpi = tifres)

# NOTE: Scale data
all.genes = rownames(srsc)
srsc = ScaleData(srsc, features = all.genes)
# NOTE: RunPCA
srsc = RunPCA(srsc, features = VariableFeatures(object = srsc))
print(srsc[["pca"]], dims = 1:5, nfeatures = 5)

visDimPlot = VizDimLoadings(srsc, dims = 1:9, reduction = "pca", ncol = 3)
ggsave(plot = visDimPlot, filename = paste(expDirPre, "PCA_VisDimLoadings.tiff", sep = ""), 
       width = 15, height = 15, dpi = tifres)

pcaDimPlot = DimPlot(srsc, reduction = "pca")
ggsave(plot = pcaDimPlot, filename = paste(expDirPre, "PCA_DimPlot.tiff", sep = ""), 
       width = 7.5, height = 6, dpi = tifres)

tiff(paste(expDirPre, "PCA_DimHeatmap.tiff", sep = ""), width = 9, height = 12, res = tifres, units = 'in')
DimHeatmap(srsc, dims = 1:15, cells = 500, balanced = TRUE)
gar = dev.off()

if (jsFlag) {
	# NOTE: Determine what components to use for further dimension reduction
	srsc = JackStraw(srsc, num.replicate = 100)
	srsc = ScoreJackStraw(srsc, dims = 1:20)

	jsPlot = JackStrawPlot(srsc, dims = 1:15)
	ggsave(plot = jsPlot, filename = paste(expDirPre, "JackStrawPlot.tiff", sep = ""), 
	       width = 7.5, height = 6, dpi = tifres)
}
elbowPlot = ElbowPlot(srsc)
ggsave(plot = elbowPlot, filename = paste(expDirPre, "ElbowPlot.tiff", sep = ""), 
       width = 7.5, height = 6, dpi = tifres)

# NOTE: Find clusters...
srsc = FindNeighbors(srsc, dims = 1:ndim)
srsc = FindClusters(srsc, resolution = 0.5) # NOTE: higher resolution, more clusters

# NOTE: further dimension reduction
srsc = RunUMAP(srsc, dims = 1:ndim)
srsc = RunTSNE(srsc, dims = 1:ndim, tsne.method = "Rtsne", check_duplicates = FALSE)

umapDim = DimPlot(srsc, reduction = "umap")
ggsave(plot = umapDim, filename = paste(expDirPre, "UMAP.tiff", sep = ""), 
       width = 7.5, height = 6, dpi = tifres)

tsneDim = DimPlot(srsc, reduction = "tsne")
ggsave(plot = tsneDim, filename = paste(expDirPre, "TSNE.tiff", sep = ""), 
       width = 7.5, height = 6, dpi = tifres)

# NOTE: Find cluster markers
clusterMarkers = FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
posCsv = paste(expDirPre, "positive_cluster_markers_only.csv", sep = "")
write.csv(clusterMarkers, file = posCsv)

topMarkers = clusterMarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
tiff(paste(expDirPre, "cluster_marker_heatmap.tiff", sep = ""), width = 9, height = 12, res = tifres, units = 'in')
DoHeatmap(srsc, features = topMarkers$gene)
gar = dev.off()

platePlot = DimPlot(srsc, reduction = "umap", group.by = "plate")
ggsave(plot = platePlot, filename = paste(expDirPre, "plate_umap_check_plot.tiff", sep = ""), 
       width = 7.5, height = 6, dpi = tifres)

tisPlot = DimPlot(srsc, reduction = "umap", group.by = "tissue")
ggsave(plot = tisPlot, filename = paste(expDirPre, "tissue_umap_check_plot.tiff", sep = ""), 
       width = 7.5, height = 6, dpi = tifres)

print(head(srsc@meta.data))
facsPlot = FeaturePlot(srsc, reduction = "umap", features = singleFACS)
ggsave(plot = facsPlot, filename = paste(expDirPre, "facs_single_marker_umap_check_plot.tiff", sep = ""), 
       width = 15, height = 12, dpi = tifres)

pd1cd39Plot = FeaturePlot(srsc, reduction = "umap", features = "PD1CD39")
ggsave(plot = pd1cd39Plot, filename = paste(expDirPre, "exhausted_marker_umap_check_plot.tiff", sep = ""), 
       width = 7.5, height = 6, dpi = tifres)

mtPlot = FeaturePlot(srsc, reduction = "umap", features = "percent.mt")
ggsave(plot = mtPlot, filename = paste(expDirPre, "mt_percent_umap_check_plot.tiff", sep = ""), 
       width = 7.5, height = 6, dpi = tifres)

saveRDS(srsc, file = paste(expDirPre, "seurat_object.RDS", sep = ""))
print(srsc)


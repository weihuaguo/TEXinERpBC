# This R script is designed based on Seurat 3.0.0
# Weihua Guo, Ph.D.
# 04/17/2019

# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear workspace
rm(list=ls())

library(dplyr)
library(Seurat)
library(tcltk)
library(ggplot2)

exp_id <- "tex_brtissue_scim_k5_res45_pptest" # Experimental ID, a new folder with this name will be created
sample_id <- exp_id
# input_dir <- "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/trm_tex/04222019_wta_all_merged_v5_cleaned.txt" # Before imputation
input_dir <- "/home/weihua/mnts/group_plee/Weihua/scrnaseq_results/trm_tex_bd/trm_tex_scim_k5/imputation_res/scimpute_count.csv"

select_file <- "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/trm_tex/04252019_flow_ann_tissue_only.csv"
meta_file <- "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/trm_tex/01222020_flow_ann_all_merged_v2.txt"
main_res_dir <- "/home/weihua/mnts/group_plee/Weihua/scrnaseq_results/scRNAseqTexResults/"

dir.create(file.path(main_res_dir, exp_id), showWarnings = TRUE)
res_dir <- paste(c(main_res_dir, exp_id, "/"),collapse="")

mt_gene_file <- "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/trm_tex/mt_gene_list.txt"
rm_gene_file <- "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/trm_tex/remove_list.txt"
kp_gene_file <- "/home/weihua/mnts/group_plee/Weihua/scrnaseq_leelab/trm_tex/keep_protein_coding_gene.txt"

mt_genes <- scan(mt_gene_file, what="", sep="\n")
rm_genes <- scan(rm_gene_file, what="", sep="\n")
kp_genes <- scan(kp_gene_file, what="", sep="\n")

##########################################################################
# Save this Rscript to the result folder
git_scfolder <- "/home/weihua/git_repo/TEXinERpBC/"
rscript_file <- list.files(git_scfolder, "seurat_upplat_v1.r")
file.copy(rscript_file, res_dir)

# PARAMETER SECTION!!!
selfil <- TRUE # TRUE: will filter by some values from select_file
fil_items <- c("Tumor","Normal","LN","PBMC") # Breast only
facs_markers <- c("CCR7", "CD39", "CD45RA", "CD69", "CD103", "CD137", "PD1")
meta_cols <- c("tissue", "patient", "plate", "PD1CD39_Status")
num_dim <- 12 # dim in PCA, tSNE
mt_thr <- 20 # threshold for percentage.mito
mt_cust <- TRUE # Use customized mitochondrial genes (mt_gene_file)
rm_flag <- TRUE # Remove the genes in rm_gene_file
kp_flag <- TRUE # Only keep the genes in kp_gene_file
impu_app <- "--" # sci: scImpute, --: no imputation
sci_kc <- 5 # kcluster in scImpute
default_assay <- "RNA" # Default assay name
mat_sep <- "," # Separator for input matrix file
low_rna_thr <- 200 # If the cell has less than this number of genes, it will be removed
up_rna_thr <- 2500 # If the cell has more than this number of genes, it will be removed
fv_thr <- 1000
mt_pat <- "^MT."
clRes = 0.45

############################################################################
## Generate a list to store the output plot names
plot_names = c("qc1", "qc2", "var_gene", "pca_2d", "pca_hm", "pca_viz", "JSP", "elbow", 
	       "umap", "tsne", "batch")
plot_save <- vector(mode = "list", length=length(plot_names))
names(plot_save) <- plot_names
plot_save["qc1"] <- paste(c(res_dir, sample_id, "_", "qcplot1.tiff"),collapse="")
plot_save["qc2"] <- paste(c(res_dir, sample_id, "_", "qcplot2.tif"),collapse="")
plot_save["var_gene"] <- paste(c(res_dir, sample_id, "_", "variable_genes.tif"),collapse="")
plot_save["pca_2d"] <- paste(c(res_dir, sample_id, "_", "pca_2dcluster.tif"),collapse="")
plot_save["pca_hm"] <- paste(c(res_dir, sample_id, "_", "pca_heatmap.tif"),collapse="")
plot_save["pca_viz"] <- paste(c(res_dir, sample_id, "_", "pca_viz.tif"),collapse="")
plot_save["JSP"] <- paste(c(res_dir, sample_id, "_", "JackStrawPlot.tif"),collapse="")
plot_save["elbow"] <- paste(c(res_dir, sample_id, "_", "ElbowPlot.tif"),collapse="")
plot_save["umap"] <- paste(c(res_dir, sample_id, "_", "umap.tif"),collapse="")
plot_save["tsne"] <- paste(c(res_dir, sample_id, "_", "tsne.tif"),collapse="")
plot_save["batch"] <- paste(c(res_dir, sample_id, "_", "batch_effect_check.tiff"),collapse="")


## Read the data from cell-expression matrix
cat("Start to reading the matrix...\n")
sc.data <- read.table(file = input_dir, header = TRUE, sep = mat_sep, row.names=NULL)
sc.data <- sc.data[!duplicated(sc.data$gene),]
rownames(sc.data) <- sc.data$gene

############################################################################
## Filter by markers/tissues
meta.data = read.table(file = meta_file, header = TRUE, sep = "\t", row.names=1)
if (selfil) {
	cat("Filter the cells based on the tissue type: \n")
	cat("\t", paste(fil_items, sep=" "), "\n")
	cell_oi = rownames(meta.data)[meta.data$tissue %in% fil_items]
	print(length(cell_oi))

#	sc.filter <- read.table(file = select_file, header=TRUE, sep=",",row.names=1)
#	temp.fildf <- as.data.frame(t(sc.filter))
#	temp.mask <- temp.fildf[,fil_items] == 1
#	temp.sum <- apply(temp.mask,1,sum)
#	temp.smask <- temp.sum == 1
	cat("Selected cell numbers: ", length(cell_oi),"\n")
#	print(head(temp.sum))
#	temp.seldf <- temp.fildf[temp.smask,]
#	dim(temp.seldf)
#	sel_cells <- row.names(temp.seldf)
#	print(setdiff(sel_cells, cellOi))
	sc.data <- sc.data[,cell_oi]
	sc.data[duplicated(sc.data$rownames),]
	print(dim(sc.data))
}

############################################################################
## Imputation section and generate Seurat objects
if (impu_app == "sci") { # Can we do basic filter first and then do the imputation?
	library(scImpute)
	temp_file <- paste(c(res_dir,"temp_sci.csv"),collapse="")
	dir.create(file.path(res_dir, "imputation_res"), showWarnings = TRUE)
	imp_dir <- paste(c(res_dir,"imputation_res","/"),collapse="")
	write.csv(sc.data, file=temp_file,row.names=TRUE)
	scimpute(count_path=temp_file, infile="csv",outfile="csv",out_dir=imp_dir,
		 Kcluster=sci_kc,ncore=12)
	sc.data <- read.table(file=paste(c(imp_dir,"scimpute_count.csv"),collapse=""),
			      header=TRUE, sep=",",row.names=1)
	srsc <- CreateSeuratObject(counts=sc.data, project=sample_id, min.cells=5, min.features=200)
	print("Printing basic info of Seurat object...",quote=FALSE)
#	GetAssayData(object = srsc, slot = 'data')[1:3, 1:9]
}

if (impu_app == "--") {
	srsc <- CreateSeuratObject(counts=sc.data, project=sample_id, min.cells=5, min.features=200)
#	GetAssayData(object = srsc, slot = 'data')[1:3,1:9]
}
print(srsc)
cat(paste(c(rep("+++",33),"\n"),collapse=""),quote=FALSE)

############################################################################
## QC 
# Calculate the mitochondrial gene percentages
if (mt_cust){
	cur_genes <- rownames(srsc[[default_assay]]@data)
	used_mt <- intersect(cur_genes, mt_genes)
	srsc[["percent.mt"]] <- PercentageFeatureSet(object = srsc, features = used_mt, assay = default_assay)
} else {
	srsc[["percent.mt"]] <- PercentageFeatureSet(object = srsc, pattern = mt_pat)
}
tiff(plot_save[["qc1"]], res = 300, width = 16, heigh = 10, units = 'in')
VlnPlot(object = srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
garbage <- dev.off()

# Plot the QC correlations 
plot1 <- FeatureScatter(object = srsc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = srsc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
tiff(plot_save[["qc2"]], res = 300, width = 16, heigh = 9, units = 'in')
CombinePlots(plots = list(plot1, plot2))
garbage <- dev.off()

# QC control remove the cells with 
# 1) too many RNA [possible cell doublets or mutliplets],
# 2) too less RNA [low-quality or empty cells],
# 3) dead cells [more than 5% expression are from mitochondrial genes]
srsc <- subset(x = srsc, subset = nFeature_RNA > low_rna_thr & nFeature_RNA < up_rna_thr & percent.mt < mt_thr)
cat("After QC, basic info of Seurat object:")
print(srsc)

## Manual QC to keep or remove some genes
orig_srsc <- srsc
cur_genes <- rownames(srsc[[default_assay]]@data)
clean.data <- srsc[[default_assay]]@data
# Remove genes
if (rm_flag){
	used_rm <- setdiff(cur_genes, rm_genes)
	clean.data <- clean.data[used_rm,]
	cat("After removing pseudogenes: ", ncol(clean.data), "\n")
	srsc <- CreateSeuratObject(counts=clean.data, project=sample_id, min.cells=5, min.features=200)
}
# Extract genes
if (kp_flag){
	cur_genes <- rownames(clean.data)
	used_kp <- intersect(cur_genes, kp_genes)
	clean.data <- clean.data[used_kp,]
	cat("After keeping protein-encoding genes: ", ncol(clean.data), "\n")
	srsc <- CreateSeuratObject(counts=clean.data, project=sample_id, min.cells=5, min.features=200)
}

## Add meta.data and FACS data
ava_cell = colnames(srsc)
if (length(intersect(ava_cell, rownames(meta.data))) != length(ava_cell)) {
	warning("Cell lost in meta.data!!!")
}
ava_metadata = meta.data[intersect(ava_cell, rownames(meta.data)),]
facs_data = ava_metadata[, facs_markers]
srsc[["FACS"]] = CreateAssayObject(counts = t(facs_data[colnames(srsc),]))
srsc@meta.data = ava_metadata[colnames(srsc), meta_cols]
cat("After manual QC, basic info of Seurat object:\n")
print(srsc)

srsc <- NormalizeData(object = srsc, normalization.method = "LogNormalize", scale.factor = 10000)

## Find variable genes
srsc <- FindVariableFeatures(object = srsc, selection.method = "vst", nfeatures = fv_thr)
# Identify the 20 most highly variable genes
top_vg <- head(x = VariableFeatures(object = srsc), 20)

plot1 <- VariableFeaturePlot(object = srsc)
plot2 <- LabelPoints(plot = plot1, points = top_vg, repel = TRUE)
tiff(plot_save[["var_gene"]], res = 300, width = 16, heigh = 9, units = 'in')
CombinePlots(plots = list(plot1, plot2))
garbage <- dev.off()

## Scaling the data
all.genes <- rownames(x = srsc)
srsc <- ScaleData(object = srsc, features = all.genes)

## Perform linear dimensional reduction
srsc <- RunPCA(object = srsc, features = VariableFeatures(object = srsc))
# print(x = srsc[["pca"]], dims = 1:10, nfeatures = 5)
tiff(plot_save[["pca_2d"]], res = 300, width = 4, heigh = 3, units = 'in')
DimPlot(object = srsc, reduction = "pca")
garbage <- dev.off()

tiff(plot_save[["pca_hm"]],res = 300, width = 9, heigh = 16, units = 'in')
DimHeatmap(object = srsc, dims = 1:15, cells = 500, balanced = TRUE)
garbage <- dev.off()

## Determine how many pricinple components we should use for tSNE
# JackStrawPlot method
srsc <- JackStraw(object = srsc, num.replicate = 100)
srsc <- ScoreJackStraw(object = srsc, dims = 1:20)

st = Sys.time()
tiff(plot_save[["JSP"]], res = 300, width = 7.5, heigh = 4.5, units = 'in')
JackStrawPlot(object = srsc, dims = 1:20)
garbage <- dev.off()
print(Sys.time()-st)

tiff(plot_save[["elbow"]], res = 300, width = 7.5, heigh = 4.5, units = 'in')
ElbowPlot(object = srsc)
garbage <- dev.off()
# IMPORTANT NOTE: how many components [dims = 1:n] used for clustering need to be determined
# by JackStrawPlot and ElbowPlot!!!

## Cluster the cells
srsc <- FindNeighbors(object = srsc, dims = 1:num_dim)
srsc <- FindClusters(object = srsc, resolution = clRes)

### Non-linear dimensional reduction!!!

## UMAP
srsc <- RunUMAP(object = srsc, dims = 1:num_dim)
tiff(plot_save[["umap"]], res = 360, width = 7.5, heigh = 4.5, units = 'in')
DimPlot(object = srsc, reduction = "umap")
garbage <- dev.off()

## tSNE
srsc <- RunTSNE(object = srsc, dims = 1:num_dim, check_duplicates = FALSE)
tiff(plot_save[["tsne"]], res = 360, width = 7.5, heigh = 4.5, units = 'in')
DimPlot(object = srsc, reduction = "tsne")
garbage <- dev.off()

# NOTE: Find cluster markers
clusterMarkers <- FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
clusterMarkers$sigFlag <- 0
clusterMarkers[clusterMarkers$avg_logFC > 1 & clusterMarkers$p_val_adj <0.10, "sigFlag"] <- 2
posCsv <- paste(res_dir, sample_id, "_positive_cluster_markers_only.csv", sep = "")
write.csv(clusterMarkers, file = posCsv)

topMarkers <- clusterMarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
tiff(paste(res_dir, sample_id, "_cluster_marker_heatmap.tiff", sep = ""), width = 9, height = 12, res = 300, units = 'in')
DoHeatmap(srsc, features = topMarkers$gene)
gar = dev.off()

## Plots
# FACS
bePlot <- DimPlot(srsc, reduction = "umap", group.by = c("tissue", "plate", "patient"), ncol = 3)
ggsave(plot = bePlot, filename = plot_save[["batch"]], width = 7.5, height = 6, dpi = 300)

print(head(srsc@meta.data))
facsPlot <- FeaturePlot(srsc, reduction = "umap", features = facs_markers)
ggsave(plot = facsPlot, filename = paste(expDirPre, "facs_single_marker_umap_check_plot.tiff", sep = ""), 
       width = 15, height = 12, dpi = tifres)

pd1cd39Plot <- FeaturePlot(srsc, reduction = "umap", features = "PD1CD39")
ggsave(plot = pd1cd39Plot, filename = paste(expDirPre, "exhausted_marker_umap_check_plot.tiff", sep = ""), 
       width = 7.5, height = 6, dpi = tifres)


## Save to RDS
rds_file = paste(c(res_dir, sample_id, "_knn_umap_tsne.rds"),collapse="")
saveRDS(srsc, file = rds_file)

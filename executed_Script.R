library(Seurat)
dataset1 <- Read10X(data.dir = paste0("counts/",  "GSM6204109"))
dataset1 <- CreateSeuratObject(counts = dataset1)
dataset1@meta.data$disease = "PDAC liver met"
dataset1@meta.data$orig.ident = "P1"
dataset1@meta.data$stage = "Metastatic"
dataset1@meta.data$treatment = "Untreated"
dataset2 <- Read10X(data.dir = paste0("counts/",  "GSM6204110"))
dataset2 <- CreateSeuratObject(counts = dataset2)
dataset2@meta.data$disease = "PDAC liver met"
dataset2@meta.data$orig.ident = "P2"
dataset2@meta.data$stage = "Metastatic"
dataset2@meta.data$treatment = "Untreated"
dataset3 <- Read10X(data.dir = paste0("counts/",  "GSM6204111"))
dataset3 <- CreateSeuratObject(counts = dataset3)
dataset3@meta.data$disease = "Primary PDAC"
dataset3@meta.data$orig.ident = "P3"
dataset3@meta.data$stage = "Resectable"
dataset3@meta.data$treatment = "Treated "
dataset4 <- Read10X(data.dir = paste0("counts/",  "GSM6204112"))
dataset4 <- CreateSeuratObject(counts = dataset4)
dataset4@meta.data$disease = "Primary PDAC"
dataset4@meta.data$orig.ident = "P4"
dataset4@meta.data$stage = "Resectable"
dataset4@meta.data$treatment = "Untreated"
dataset5 <- Read10X(data.dir = paste0("counts/",  "GSM6204113"))
dataset5 <- CreateSeuratObject(counts = dataset5)
dataset5@meta.data$disease = "Primary PDAC"
dataset5@meta.data$orig.ident = "P5"
dataset5@meta.data$stage = "Resectable"
dataset5@meta.data$treatment = "Untreated"
dataset6 <- Read10X(data.dir = paste0("counts/",  "GSM6204114"))
dataset6 <- CreateSeuratObject(counts = dataset6)
dataset6@meta.data$disease = "Primary PDAC"
dataset6@meta.data$orig.ident = "P6"
dataset6@meta.data$stage = "Locally Advanced"
dataset6@meta.data$treatment = "Treated"
dataset7 <- Read10X(data.dir = paste0("counts/",  "GSM6204115"))
dataset7 <- CreateSeuratObject(counts = dataset7)
dataset7@meta.data$disease = "Primary PDAC"
dataset7@meta.data$orig.ident = "P7"
dataset7@meta.data$stage = "Resectable"
dataset7@meta.data$treatment = "Untreated"
dataset8 <- Read10X(data.dir = paste0("counts/",  "GSM6204116"))
dataset8 <- CreateSeuratObject(counts = dataset8)
dataset8@meta.data$disease = "Primary PDAC"
dataset8@meta.data$orig.ident = "P8"
dataset8@meta.data$stage = "Borderline"
dataset8@meta.data$treatment = "Treated"
dataset9 <- Read10X(data.dir = paste0("counts/",  "GSM6204117"))
dataset9 <- CreateSeuratObject(counts = dataset9)
dataset9@meta.data$disease = "Primary PDAC"
dataset9@meta.data$orig.ident = "P9"
dataset9@meta.data$stage = "Metastatic"
dataset9@meta.data$treatment = "Untreated"
dataset10 <- Read10X(data.dir = paste0("counts/",  "GSM6204118"))
dataset10 <- CreateSeuratObject(counts = dataset10)
dataset10@meta.data$disease = "Primary PDAC"
dataset10@meta.data$orig.ident = "P10"
dataset10@meta.data$stage = "Borderline"
dataset10@meta.data$treatment = "Treated"
dataset11 <- Read10X(data.dir = paste0("counts/",  "GSM6204119"))
dataset11 <- CreateSeuratObject(counts = dataset11)
dataset11@meta.data$disease = "PDAC liver met"
dataset11@meta.data$orig.ident = "P11"
dataset11@meta.data$stage = "Metastatic"
dataset11@meta.data$treatment = "Untreated"
dataset12 <- Read10X(data.dir = paste0("counts/",  "GSM6204120"))
dataset12 <- CreateSeuratObject(counts = dataset12)
dataset12@meta.data$disease = "Primary PDAC"
dataset12@meta.data$orig.ident = "P12"
dataset12@meta.data$stage = "Metastatic"
dataset12@meta.data$treatment = "Treated"
dataset13 <- Read10X(data.dir = paste0("counts/",  "GSM6204121"))
dataset13 <- CreateSeuratObject(counts = dataset13)
dataset13@meta.data$disease = "Primary PDAC"
dataset13@meta.data$orig.ident = "P13"
dataset13@meta.data$stage = "Metastatic"
dataset13@meta.data$treatment = "Untreated"
dataset14 <- Read10X(data.dir = paste0("counts/",  "GSM6204122"))
dataset14 <- CreateSeuratObject(counts = dataset14)
dataset14@meta.data$disease = "Primary PDAC"
dataset14@meta.data$orig.ident = "P14"
dataset14@meta.data$stage = "Metastatic"
dataset14@meta.data$treatment = "Treated"
dataset15 <- Read10X(data.dir = paste0("counts/",  "GSM6204123"))
dataset15 <- CreateSeuratObject(counts = dataset15)
dataset15@meta.data$disease = "Primary PDAC"
dataset15@meta.data$orig.ident = "P15"
dataset15@meta.data$stage = "Metastatic"
dataset15@meta.data$treatment = "Untreated"
dataset16 <- Read10X(data.dir = paste0("counts/",  "GSM6204124"))
dataset16 <- CreateSeuratObject(counts = dataset16)
dataset16@meta.data$disease = "PDAC liver met"
dataset16@meta.data$orig.ident = "P16"
dataset16@meta.data$stage = "Metastatic"
dataset16@meta.data$treatment = "Untreated"
dataset17 <- Read10X(data.dir = paste0("counts/",  "GSM6204125"))
dataset17 <- CreateSeuratObject(counts = dataset17)
dataset17@meta.data$disease = "PDAC liver met"
dataset17@meta.data$orig.ident = "P17"
dataset17@meta.data$stage = "Metastatic"
dataset17@meta.data$treatment = "Treated"
dataset18 <- Read10X(data.dir = paste0("counts/",  "GSM6204126"))
dataset18 <- CreateSeuratObject(counts = dataset18)
dataset18@meta.data$disease = "PDAC liver met"
dataset18@meta.data$orig.ident = "P18"
dataset18@meta.data$stage = "Metastatic"
dataset18@meta.data$treatment = "Untreated"
dataset19 <- Read10X(data.dir = paste0("counts/",  "GSM6204127"))
dataset19 <- CreateSeuratObject(counts = dataset19)
dataset19@meta.data$disease = "Primary PDAC"
dataset19@meta.data$orig.ident = "P19"
dataset19@meta.data$stage = "Resectable"
dataset19@meta.data$treatment = "Untreated"
dataset20 <- Read10X(data.dir = paste0("counts/",  "GSM6204128"))
dataset20 <- CreateSeuratObject(counts = dataset20)
dataset20@meta.data$disease = "Primary PDAC"
dataset20@meta.data$orig.ident = "P20"
dataset20@meta.data$stage = "Resectable"
dataset20@meta.data$treatment = "Untreated"
dataset21 <- Read10X(data.dir = paste0("counts/",  "GSM6204129"))
dataset21 <- CreateSeuratObject(counts = dataset21)
dataset21@meta.data$disease = "PDAC liver met"
dataset21@meta.data$orig.ident = "P21"
dataset21@meta.data$stage = "Metastatic"
dataset21@meta.data$treatment = "Untreated"
dataset22 <- Read10X(data.dir = paste0("counts/",  "GSM6204130"))
dataset22 <- CreateSeuratObject(counts = dataset22)
dataset22@meta.data$disease = "Primary PDAC"
dataset22@meta.data$orig.ident = "P22"
dataset22@meta.data$stage = "Locally Advanced"
dataset22@meta.data$treatment = "Untreated"
dataset23 <- Read10X(data.dir = paste0("counts/",  "GSM6204131"))
dataset23 <- CreateSeuratObject(counts = dataset23)
dataset23@meta.data$disease = "Primary PDAC"
dataset23@meta.data$orig.ident = "P23"
dataset23@meta.data$stage = "Resectable"
dataset23@meta.data$treatment = "Untreated"
dataset24 <- Read10X(data.dir = paste0("counts/",  "GSM6204132"))
dataset24 <- CreateSeuratObject(counts = dataset24)
dataset24@meta.data$disease = "PDAC liver met"
dataset24@meta.data$orig.ident = "P24"
dataset24@meta.data$stage = "Metastatic"
dataset24@meta.data$treatment = "Untreated"
dataset25 <- Read10X(data.dir = paste0("counts/",  "GSM6204133"))
dataset25 <- CreateSeuratObject(counts = dataset25)
dataset25@meta.data$disease = "PDAC liver met"
dataset25@meta.data$orig.ident = "P25"
dataset25@meta.data$stage = "Metastatic"
dataset25@meta.data$treatment = "Untreated"
dataset26 <- Read10X(data.dir = paste0("counts/",  "GSM6204134"))
dataset26 <- CreateSeuratObject(counts = dataset26)
dataset26@meta.data$disease = "Primary PDAC"
dataset26@meta.data$orig.ident = "P26"
dataset26@meta.data$stage = "Metastatic"
dataset26@meta.data$treatment = "Untreated"
dataset27 <- Read10X(data.dir = paste0("counts/",  "GSM6204135"))
dataset27 <- CreateSeuratObject(counts = dataset27)
dataset27@meta.data$disease = "PDAC liver met"
dataset27@meta.data$orig.ident = "P27"
dataset27@meta.data$stage = "Metastatic"
dataset27@meta.data$treatment = "Untreated"
# Merge dataset1 and dataset2
merged1 <- merge(dataset1, dataset2)
# Remove dataset1 and dataset2
rm(dataset1, dataset2)
# Merge merged1 with dataset3
merged2 <- merge(merged1, dataset3)
# Remove merged1 and dataset3
rm(merged1, dataset3)
# Merge merged2 with dataset4
merged3 <- merge(merged2, dataset4)
# Remove merged2 and dataset4
rm(merged2, dataset4)
# Merge merged3 with dataset5
merged4 <- merge(merged3, dataset5)
# Remove merged3 and dataset5
rm(merged3, dataset5)
# Merge merged4 with dataset6
merged5 <- merge(merged4, dataset6)
# Remove merged4 and dataset6
rm(merged4, dataset6)
# Merge merged5 with dataset7
merged6 <- merge(merged5, dataset7)
# Remove merged5 and dataset7
rm(merged5, dataset7)
# Merge merged6 with dataset8
merged7 <- merge(merged6, dataset8)
# Remove merged6 and dataset8
rm(merged6, dataset8)
# Merge merged7 with dataset9
merged8 <- merge(merged7, dataset9)
# Remove merged7 and dataset9
rm(merged7, dataset9)
# Merge merged8 with dataset10
merged9 <- merge(merged8, dataset10)
# Remove merged8 and dataset10
rm(merged8, dataset10)
# Merge merged9 with dataset11
merged10 <- merge(merged9, dataset11)
# Remove merged9 and dataset11
rm(merged9, dataset11)
# Merge merged10 with dataset12
merged11 <- merge(merged10, dataset12)
# Remove merged10 and dataset12
rm(merged10, dataset12)
# Merge merged11 with dataset13
merged12 <- merge(merged11, dataset13)
# Remove merged11 and dataset13
rm(merged11, dataset13)
# Merge merged12 with dataset14
merged13 <- merge(merged12, dataset14)
# Remove merged12 and dataset14
rm(merged12, dataset14)
# Merge merged13 with dataset15
merged14 <- merge(merged13, dataset15)
# Remove merged13 and dataset15
rm(merged13, dataset15)
# Merge merged14 with dataset16
merged15 <- merge(merged14, dataset16)
# Remove merged14 and dataset16
rm(merged14, dataset16)
# Merge merged15 with dataset17
merged16 <- merge(merged15, dataset17)
# Remove merged15 and dataset17
rm(merged15, dataset17)
# Merge merged16 with dataset18
merged17 <- merge(merged16, dataset18)
# Remove merged16 and dataset18
rm(merged16, dataset18)
# Merge merged17 with dataset19
merged18 <- merge(merged17, dataset19)
# Remove merged17 and dataset19
rm(merged17, dataset19)
# Merge merged18 with dataset20
merged19 <- merge(merged18, dataset20)
# Remove merged18 and dataset20
rm(merged18, dataset20)
# Merge merged19 with dataset21
merged20 <- merge(merged19, dataset21)
# Remove merged19 and dataset21
rm(merged19, dataset21)
# Merge merged20 with dataset22
merged21 <- merge(merged20, dataset22)
# Remove merged20 and dataset22
rm(merged20, dataset22)
# Merge merged21 with dataset23
merged22 <- merge(merged21, dataset23)
# Remove merged21 and dataset23
rm(merged21, dataset23)
# Merge merged22 with dataset24
merged23 <- merge(merged22, dataset24)
# Remove merged22 and dataset24
rm(merged22, dataset24)
# Merge merged23 with dataset25
merged24 <- merge(merged23, dataset25)
# Remove merged23 and dataset25
rm(merged23, dataset25)
# Merge merged24 with dataset26
merged25 <- merge(merged24, dataset26)
# Remove merged24 and dataset26
rm(merged24, dataset26)
# Merge merged25 with dataset27
final_merged <- merge(merged25, dataset27)
# Remove merged25 and dataset27
rm(merged25, dataset27)
# Verify dimensions and metadata of the final merged dataset
dim(final_merged)
head(final_merged@meta.data)
#loading in the libraries that are required for scRNA-seq analysis
#using Seurat for (QC, normalization, scaling, features selection, cluster)
#using scType for cell annotation
library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(DropletUtils)
library(patchwork)
library(harmony)
library(SoupX)
library(glmnet)
library(ggplot2)
library(ComplexHeatmap)
##scRNA-seq analysis pipeline
pc_combined_data <- final_merged
dim(pc_combined_data)
#Quality Control (identification of dead cells and subsequent removal)
vln_plot1 <- VlnPlot(pc_combined_data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
vln_plot1
ggsave("feature_count_RNA.png", vln_plot1, dpi = 600)
pc_combined_data[["percent.mt"]] <- PercentageFeatureSet(pc_combined_data, pattern = "^MT-")
vln_plot2 <- VlnPlot(pc_combined_data, features = c("percent.mt"), ncol = 1)
ggsave("percent_mt.png", vln_plot2, dpi = 600)
#filtering dead cells out, highly expressive cells and lowly expressive cells are being filtered out
##nFeature_RNA is the number of genes detected in each cell.
##nCount_RNA is the total number of molecules detected within a cell.
##Low nFeature_RNA for a cell indicates that it may be dead/dying or an empty droplet.
##High nCount_RNA and/or nFeature_RNA indicates that the "cell" may in fact be a doublet (or multiplet).
##In combination with %mitochondrial reads, removing outliers from these groups removes most doublets/dead cells/empty droplets, hence why filtering is a common pre-processing step.
pc_combined_data = subset(pc_combined_data, subset = nFeature_RNA > 200 &
nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.mt < 10)
pc_combined_data[["percent.mt"]] <- PercentageFeatureSet(pc_combined_data, pattern = "^MT-")
vln_plot3 <- VlnPlot(pc_combined_data, features = c("percent.mt"), ncol = 1)
vln_plot3
ggsave("after-qc-percent_mt.png", vln_plot3, dpi = 600)
#Normalization of the data
##The NormalizeData step is basically just ensuring expression values across cells are on a comparable scale.
##By default, it will divide counts for each gene by the total counts in the cell, multiply that value for each gene by the scale.factor (10,000 by default), and then natural log-transform them.
pc_combined_data = NormalizeData(object = pc_combined_data, normalization.method = "LogNormalize", scale.factor = 1e4)
##Finding variable features (2000 highly variable genes that have very high expression in some cells and very low in some)
pc_combined_data = FindVariableFeatures(pc_combined_data, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pc_combined_data), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pc_combined_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
ggsave("highly-variable-features.png", plot2, dpi = 600, bg = "white")
#Scaling the data
pc_combined_data = ScaleData(pc_combined_data, verbose = F)
#High dimension reduction (PCA, Heatmap, JackStarwPlot)
pc_combined_data = RunPCA(pc_combined_data, npcs = 50, features=VariableFeatures(object = pc_combined_data), verbose = F)
print(pc_combined_data[["pca"]], dims = 1:5, nfeatures = 10)
pca_plot <- DimPlot(pc_combined_data, reduction = "pca", dims = (c(1,2)))
pca_plot
ggsave("pca_plot.png", pca_plot, dpi = 600, bg = "white")
##plot not saved
DimHeatmap(pc_combined_data, dims = 1:2, cells = 500, balanced = TRUE)
ggsave("basic_heatmap.png", last_plot(), dpi = 600, bg = "white")
#Cell clustering based on gene expression data based on nearest-neighbor method implemented in Seurat
pc_combined_data = FindNeighbors(pc_combined_data, dims = 1:20, verbose = F)
pc_combined_data = FindClusters(pc_combined_data, resolution = 0.6, verbose = F)
pc_combined_data = RunUMAP(pc_combined_data, dims = 1:20, verbose = F)
DimPlot(pc_combined_data, label=T)
ggsave("umap_unannotated.png", last_plot(), dpi = 600, bg = "white")
DimPlot(pc_combined_data, group.by = "orig.ident", label = F)
ggsave("umap_unannotated_normal_tumor.png", last_plot(), dpi = 600, bg = "white")
pc_combined_data = RunTSNE(pc_combined_data, dims = 1:20, verbose = F)
DimPlot(object = pc_combined_data, reduction = "tsne", label = T)
ggsave("tsne_unannotated.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "orig.ident")
ggsave("tsne_unannotated_normal_tumor-original-identity.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "disease")
ggsave("tsne_unannotated_normal_tumor-disease.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "stage")
ggsave("tsne_unannotated_normal_tumor-stage.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "treatment")
ggsave("tsne_unannotated_normal_tumor-treatment.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "orig.ident", split.by = "disease", shape.by = "stage")
ggsave("tsne_unannotated_normal_tumor-treatment-diease-stage.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "orig.ident", split.by = "stage", shape.by = "disease")
ggsave("tsne_1.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "disease", split.by = "stage")
ggsave("tsne_2.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "stage", split.by = "disease")
ggsave("tsne_3.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "orig.ident", split.by = "stage")
ggsave("tsne_4.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "orig.ident", split.by = "disease")
ggsave("tsne_5.png", last_plot(), dpi = 600, bg = "white")
#optional
VizDimLoadings(pc_combined_data, dims = 1:2, reduction = "pca")
ggsave("pca_genes.png", last_plot(), dpi = 600, bg = "white")
#Cell annotation using a supervised method that uses experimentally validated marker database, scType
library(HGNChelper)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = c("Immune system", "Liver", "Pancreas") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = pc_combined_data[["RNA"]]$scale.data, scaled = TRUE,
gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix.
# In case Seurat is used, it is either pc[["RNA"]]@scale.data (default), pc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(pc_combined_data@meta.data$seurat_clusters), function(cl){
es.max.cl = sort(rowSums(es.max[ ,rownames(pc_combined_data@meta.data[pc_combined_data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pc_combined_data@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
as.data.frame(cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores))[,2]
# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
scores <- sctype_scores[,1:3]
write.csv(scores, "cell_annotations.csv")
pc_combined_data@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
cl_type = sctype_scores[sctype_scores$cluster==j,];
pc_combined_data@meta.data$customclassif[pc_combined_data@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(pc_combined_data, reduction = "umap", label = TRUE, group.by = 'customclassif', repel = T)
ggsave("umap_annotated_normal_tumor_cells_removed.png", last_plot(), dpi = 600, bg = "white")
DimPlot(pc_combined_data, reduction = "tsne", label = TRUE, group.by = 'customclassif', repel = T)
ggsave("tsne_annotated_normal_tumor_cells_removed.png", last_plot(), dpi = 600, bg = "white")
DimPlot(pc_combined_data, reduction = "tsne", label = TRUE, group.by = 'customclassif', split.by = "disease", repel = T)
ggsave("tsne_annotated_1.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "customclassif", split.by = "stage", repel = T)
ggsave("tsne_annotated_2.png", last_plot(), dpi = 600, bg = "white")
DimPlot(object = pc_combined_data, reduction = "tsne", group.by = "customclassif", split.by = "stage", label = T)
ggsave("tsne_annotated_3.png", last_plot(), dpi = 600, bg = "white")
#Differential gene expression analysis and subsequent visualizations
levels(pc_combined_data)
pc_combined_data <- RenameIdents(object = pc_combined_data,
"0" = "Naive CD4+ T cells",
"1" = "Non-classical monocytes",
"2" = "Effector CD8+ T cells",
"3" = "Liver progenitor cell",
"4" = "Natural killer cells",
"5" = "Hepatocytes",
"6" = "Hepatocytes",
"7" = "Liver progenitor cell",
"8" = "Non-classical monocytes",
"9" = "Hepatocytes",
"10" = "Hematopoietic cell",
"11" = "Hepatocytes",
"12" = "Cancer stem cells",
"13" = "Endothelial cells",
"14" = "Memory B cells",
"15" = "Liver progenitor cell",
"16" = "Effector CD8+ T cells",
"17" = "Plasmacytoid Dendritic cells",
"18" = "Hepatocytes",
"19" = "Hepatocytes",
"20" = "Endothelial cells",
"21" = "Effector CD8+ T cells",
"22" = "Memory B cells",
"23" = "Liver progenitor cell",
"24" = "Unknown",
"25" = "HSC/MPP cells",
"26" = "Effector CD8+ T cells",
"27" = "Platelets",
"28" = "Basophils",
"29" = "Non-classical monocytes",
"30" = "Non-classical monocytes",
"31" = "Plasmacytoid Dendritic cells",
"32" = "Non-classical monocytes",
"33" = "Memory B cells",
"34" = "Beta cells",
"35" = "Non-classical monocytes",
"36" = "Hepatocytes")
Idents(pc_combined_data)
delete(final_merged)
del(final_merged)
rm(final_merged)
# Find differentially expressed features between
effoctor.markers <- FindMarkers(pc_combined_data, ident.1 = "Effector CD8+ T cells", ident.2 = NULL)
pc_combined_data<-JoinLayers(pc_combined_data)
gc()
save.image(file = "pancreatic-us.RData")
?FindAllMarkers()
##finding all markers, computing altogether
markers_genes_up <- FindAllMarkers(pc_combined_data, log2FC.threshold = 0.2, test.use = "wilcox",
min.pct = 0.1, min.diff.pct = 0.2,
assay = "RNA")
write.csv("marker_genes.csv", markers_genes_up)
write.csv(markers_genes_up, "marker.genes.csv")
# Filter top 10 genes for each cell type
top10 <- markers_genes_up %>%
group_by(cluster) %>%
top_n(5, avg_log2FC)
down_top5 <- markers_genes %>%
group_by(cluster) %>%
top_n(-5, avg_log2FC)
down_top5 <- markers_genes_up %>%
group_by(cluster) %>%
top_n(-5, avg_log2FC)
write.csv(down_top5, "down.top.csv")
write.csv(up_top5, "up.top.csv")
write.csv(top10, "up.top.csv")
DoHeatmap(object = pc_combined_data, features = top10$gene)
ggsave("heatmap_top5_genes.png", last_plot(), dpi = 900)

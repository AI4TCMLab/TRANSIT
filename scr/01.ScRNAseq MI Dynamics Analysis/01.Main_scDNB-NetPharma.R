# ==============================================================================
# Script: Pipeline_01_Global_Integration_and_Comprehensive_Exploration.R
# Purpose: Raw 10X data processing, QC, Harmony integration, global annotation, 
#          and exhaustive exploratory marker visualizations.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment Initialization & Dependencies
# ------------------------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE, future.globals.maxSize = 1e9)
setwd("/home/tsh/E-MTAB-7895/")

library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(tidyverse)
library(stringr)
library(harmony)
library(slingshot)
library(batchelor)
library(ggplot2)
library(patchwork)
library(clustree)
library(RColorBrewer)
library(wesanderson)
library(circlize)
library(pheatmap)
library(gplots)
library(corrplot)
library(gridExtra)
library(igraph)
library(dynamicTreeCut)

# ------------------------------------------------------------------------------
# 2. Raw Data Ingestion & Merging
# ------------------------------------------------------------------------------
rawdata_path <- "/home/tsh/E-MTAB-7895/single_cell/"
filenames <- list.files(rawdata_path)
full_paths <- paste0(rawdata_path, filenames)

sceList <- lapply(full_paths, function(x){
  project_name <- str_split(x, '/')[[1]][6]
  CreateSeuratObject(counts = Read10X(x), project = project_name)
})
names(sceList) <- filenames

sce <- merge(sceList[[1]], sceList[-1], add.cell.ids = names(sceList), project = 'Bobo')

# ------------------------------------------------------------------------------
# 3. Quality Control (QC) & Preprocessing
# ------------------------------------------------------------------------------
sce <- PercentageFeatureSet(sce, pattern = "^mt-", col.name = "percent_mt")
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3, pt.size = 0, group.by = 'orig.ident')

sce <- subset(sce, subset = nCount_RNA < 10000 & nFeature_RNA > 500 & nFeature_RNA < 4000 & percent_mt < 5)

sce <- NormalizeData(sce) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 50, verbose = FALSE)

# ------------------------------------------------------------------------------
# 4. Dimensionality Reduction & Harmony Integration
# ------------------------------------------------------------------------------
sce <- FindNeighbors(sce, dims = 1:30, reduction = "pca")
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 2)){
  sce <- FindClusters(sce, resolution = res)
}
clustree(sce, prefix = "RNA_snn_res.")

sce <- IntegrateLayers(
  object = sce, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE
)

sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:30)
sce <- FindClusters(sce, resolution = 2, cluster.name = "harmony_clusters")
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
sce <- RunTSNE(sce, reduction = "harmony", dims = 1:30, reduction.name = "tsne.harmony")
sce <- RunUMAP(sce, dims = 1:30, reduction.name = "umap")
sce <- JoinLayers(sce)

custom_order <- c("0day", "1day", "3day", "5day", "7day", "14day", "28day")
sce@meta.data$orig.ident <- factor(sce@meta.data$orig.ident, levels = custom_order)

# ------------------------------------------------------------------------------
# 5. Global Cell Type Annotation & Compositional Plotting
# ------------------------------------------------------------------------------
knowmarkers <- factor(c(
  "ZsGreen","Wt1", "Col1a1","Gsn","Dcn", "Cthrc1","Postn","Acta2", "Chil3","Plac8","Ly6c2","Cd14", 
  "C1qa","Cd68","Ms4a7", "Wif1","Dkk3", "Arg1","Itgam","Saa3", "H2-Ab1","Cd74", "Fabp4","Pecam1", 
  "Mt1","Mt2", "Timp1","Prg4", "Cd79a","H2-DMb2", "Clu","Dmkn", "S100a9","S100a8", "Nkg7","Cd3d",
  "Ms4a4b","Gzma", "Rgs5","Kcnj8","P2ry14", "Lyve1","Cldn5", "Plp1","Kcna1"
))

anno_human <- c(
  "Fibroblast I", "Fibroblast I", "Fibroblast I", "Myofibroblasts", "Monocytes I", "Fibroblast I", 
  "Macrophages", "Fibroblast II", "Fibroblast I", "Macrophages", "Macrophages", "Fibroblast I", 
  "Fibroblast I", "Monocytes I", "Macrophages", "Monocytes II", "Macrophages", "DCs", "Endothelium", 
  "Fibroblast III", "Monocytes I", "Macrophages", "Macrophages", "B cells", "Myofibroblasts", 
  "Endothelium", "Myofibroblasts", "Endothelium", "Epicardium", "Epicardium", "Granulocytes", 
  "Endothelium", "Granulocytes", "NK cell", "Endothelium", "NK cell", "SMCs_Pericytes", "Myofibroblasts", 
  "Endothelium", "Macrophages", "Lymphatic endothelium", "Endothelium", "B cells", "DCs", "Fibroblast I", 
  "NK cell", "Schwann cells", "B cells", "Fibroblast I"
)

Idents(sce) <- "seurat_clusters"
names(anno_human) <- levels(sce)
sce <- RenameIdents(sce, anno_human)
sce$anno_human <- Idents(sce)

cell_order <- c("Fibroblast I", "Macrophages","Endothelium", "Monocytes I", "Myofibroblasts",
                "B cells","Fibroblast II","Fibroblast III", "DCs","Monocytes II","NK cell",
                "Granulocytes", "Epicardium","SMCs_Pericytes","Lymphatic endothelium",'Schwann cells')
sce@meta.data$anno_human <- factor(sce@meta.data$anno_human, levels = cell_order)
Idents(sce) <- "anno_human"

# Visualizations
DotPlot(object = sce, features = as_factor(knowmarkers), dot.scale = 6) +
  RotatedAxis() + scale_x_discrete("") + scale_y_discrete("") + 
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) +
  scale_color_continuous(low="#424da7", high = "#dd2b19")

plot2 <- DimPlot(sce, reduction = "tsne.harmony", group.by = "orig.ident") + 
  scale_color_manual(values = c("0day"="#DEB887","1day"="#B22222","3day"="#7FFF00",
                                "5day"="#7FFFD4","7day"="#4169E1","14day"="#191970","28day"="#0000CD")) +
  labs(x = "t-SNE1", y = "t-SNE2", title = "DAY")
print(plot2)

# Composition Stacked Barplot
cell_data <- as.data.frame(prop.table(table(Idents(sce), sce$orig.ident), margin = 2))
colnames(cell_data) <- c("Cell_Type", "Timepoint", "Abundance")
cell_data <- cell_data %>% arrange(Timepoint) %>% mutate(Cell_Type = factor(Cell_Type, levels = cell_order))

color_palette <- c("Fibroblast I" = "#BEFF2F","Macrophages" = "#A0522D", "Myofibroblasts" = "#406400",
                   "Fibroblast III" = "#0F4F4F", "Endothelium"="#808080","Fibroblast II"="#00FA9A",
                   "Monocytes I"="#D2691E","DCs"="#4169E1","Monocytes II"="#DAA520","Epicardium"="#90EE90",
                   "NK cell"="#AFEEEE","SMCs_Pericytes"="#800000","Granulocytes"="#6A5ACD",
                   "Lymphatic endothelium"="#8B4513","B cells"="#191970",'Schwann cells'="#708090")

p_comp <- ggplot(cell_data, aes(x = Timepoint, y = Abundance, fill = fct_rev(Cell_Type))) +
  geom_col() + scale_fill_manual(values = color_palette) + 
  labs(x = "Timepoint", y = "Relative Abundance", fill = "Cell Type") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_text(size = 16, color = "black")) +
  guides(fill = guide_legend(reverse = TRUE))
print(p_comp)

# ==============================================================================
# 6. Primary Subsetting (First Pass Sub-clustering)
# ==============================================================================
# Historical subsetting selections preserved for record
wanted_clusters_history <- list(
  fibro_ext = c("Fibroblast I","Fibroblast II", "Fibroblast III","Myofibroblasts","Epicardium"),
  fibro_core = c("Fibroblast I","Fibroblast III","Myofibroblasts"),
  myeloid = c("Macrophages","Monocytes I","Monocytes II"),
  granulo = c("Granulocytes")
)

# Active Subsetting Execution
subsce <- subset(sce, subset = anno_human %in% wanted_clusters_history$granulo) 
Idents(subsce) <- "orig.ident"

subsce <- RunHarmony(subsce, c("orig.ident")) %>%
  FindVariableFeatures() %>% ScaleData(features = rownames(subsce)) %>% RunPCA() %>%
  FindNeighbors(dims = 1:30, features = VariableFeatures(object = subsce)) %>%
  FindClusters(resolution = 0.36) %>% # 2.5 used for Fibroblast/Macrophages historically
  RunTSNE(dims = 1:30, reduction = "harmony") %>% RunUMAP(dims = 1:30, reduction = "harmony") %>%
  RunTSNE(dims = 1:30, reduction = "pca") %>% RunUMAP(dims = 1:30, reduction = "pca")

subsce <- JoinLayers(subsce)
sce.markers <- FindAllMarkers(subsce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(sce.markers, file = "Granulocytesmarkers.csv")

# ==============================================================================
# 7. Exhaustive Marker Explorations (The "Digital Lab Notebook" Section)
# ==============================================================================

# --- 7.1 Global Features & Pre-defined Vectors ---
Pro_inflammatory <- c("Il1b", "Il6","Cxcl1", "Cxcl2", "Cxcl5", "Cxcl8",'Mt2', 'Mt1','Angptl4',"Gsn","Pdpn")
DotPlot(object = subsce, features = Pro_inflammatory) + RotatedAxis() + scale_x_discrete("") + scale_y_discrete("")

DimPlot(subsce, reduction = "umap", group.by = c("anno_human", "seurat_clusters"), label = TRUE)

FeaturePlot(subsce, features = c("Ly6a","Mt1","Mt2","Angptl4","Ccl2","Prg4","Csf1","Tgfb1","Cxcl1","Cxcl5"), reduction = "umap" )
FeaturePlot(subsce, features = c("Stmn1","Tubb3","Acta2","Tpm2","Cthrc1","H2afz","Ckap4"), reduction = "umap" )
FeaturePlot(subsce, features = c("Pdgfra","Hif1a","Vegfd","Mmp14"), reduction = "umap" )
FeaturePlot(subsce, features = c("Comp","Sfrp2","Col1a1","Dkk3","Postn","Ly6a","Acta2","Ccnb2"), reduction = "umap" )
FeaturePlot(subsce, features = c("Pdgfra","Ly6a","Postn","Ccl2","Angptl4","Cilp","Acta2","Ccnb2","Cdk1","Col1a1","Comp","Sfrp2","Sfrp2","Tgfb1","Scx","Thbs4","Ccn4"), reduction = "umap" )
VlnPlot(subsce, features = c("Pdgfra","Ly6a","Postn","Ccl2","Angptl4","Cilp","Acta2","Ccnb2","Cdk1","Col1a1","Comp","Sfrp2","Tgfb1","Scx","Thbs4","Ccn4"), group.by = "seurat_clusters", pt.size = 0)

# --- 7.2 Fibroblast Sub-cluster Verifications ---
# Clusters 7, 24
VlnPlot(subsce, features = c("Ly6a","Mt1","Mt2","Angptl4","Ccl2","Prg4","Csf1","Tgfb1","Cxcl1","Cxcl5"), group.by = "seurat_clusters", pt.size = 0)
# Clusters 18, 20
VlnPlot(subsce, features = c("Postn","Cilp","Acta2","Ccnb2","Cdk1","Mki67"), group.by = "seurat_clusters", pt.size = 0)
# Clusters 6, 9, 10, 11, 22, 28
VlnPlot(subsce, features = c("Postn","Cilp","Acta2","Cthrc1","Col1a1","Lyz2","Notch2"), group.by = "seurat_clusters", pt.size = 0)
# Cluster 11
VlnPlot(subsce, features = c("Postn","Comp","Sfrp2","Tgfb1","Scx","Thbs4","Ccn4"), group.by = "seurat_clusters", pt.size = 0)
# Cluster 24
VlnPlot(subsce, features = c("Postn","Pdgfra","Hif1a","Vegfd","Mmp14","Angptl4","Angpt2"), group.by = "seurat_clusters", pt.size = 0)

# --- 7.3 Macrophage / Monocyte Polarization Validations ---
# Core Monocytes (Clusters 5, 7, 12, 13, 18)
VlnPlot(subsce, features = c("Ly6c2","Plac8","Ccr2"), group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(subsce, features = c("Ly6c2","Plac8","Ccr2"), reduction = "umap" )

# M1 Subtype (iNOS, Il1b, Il6, Tnfa)
VlnPlot(subsce, features = c("Nos2","Il1b","Il6","Il23","Mmp8","Mmp9","Tnf"), group.by = "seurat_clusters", pt.size = 0)

# M2 Subtype (Arg1, Tgfb, Il10, Vegfa/b/c/d)
VlnPlot(subsce, features = c("Arg1","Tgfb1","Tgfb2","Tgfb3","Il10","Il12","Vegfa"), group.by = "seurat_clusters", pt.size = 0)

# Leukocyte chemotaxis & IFN-stimulated (Clusters 7, 22)
VlnPlot(subsce, features = c("S100a9","S100a8","Ccl7","Rsad2","Isg15","Il18","Irf7","Cxcl10"), group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(subsce, features = c("Ifit1","Ifit3","Isg15", "Rsad2"), reduction = "umap" )

# Anti-inflammatory potential (Clusters 1, 2, 3, 4, 20)
VlnPlot(subsce, features = c("Arg1","Trem2","Spp1","Gpnmb","Ctsd","Cd63","Tgfb1","Tgfb2","Il10","Il12","Hmox1","Fabp4","Fabp5"), group.by = "seurat_clusters", pt.size = 0)
VlnPlot(subsce, features = c("Arg1","Trem2","Gdf15","Psap","Pld3","Nceh1","Tgfb1","Tgfb2","Tgfb3","Il10","Il12","Gpnmb","Hmox1","Fabp4","Fabp5"), group.by = "seurat_clusters", pt.size = 0)
VlnPlot(subsce, features = c("Arg1","Trem2","Prdx1","Tgfb1","Tgfb2","Tgfb3","Il10","Il12","Gpnmb","Hmox1","Fabp4","Fabp5"), group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(subsce, features = c("Trem2","Spp1","Gpnmb","Gdf15","Fabp5","Apoe","Ms4a7"), reduction = "umap" )

# Monocyte-derived antigen-presenting (Clusters 0, 7, 10)
VlnPlot(subsce, features = c("Ccr2","H2-Eb1","H2-DMa","Il1b","Tnfsf9","Tnip3"), group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(subsce, features = c("Ccr2","H2-Aa","H2-Eb1","H2-DMa"), reduction = "umap" )

# Pro-angiogenic (Clusters 6, 22)
VlnPlot(subsce, features = c("Pf4", "Fn1","Slc7a2","Sdc3","Pdgfra","Hif1a","Vegfd","Vegfa","Vegfc","Kdr"), group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(subsce, features = c( "Pf4", "Fn1","Saa3","Fcna"), reduction = "umap" )

# Proliferative Phenotypes (Clusters 14, 11, 21)
VlnPlot(subsce, features = c( "Stmn1","Birc5","Top2a","Ube2c" ), group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(subsce, features = c("Stmn1","Cenpa","Pclaf"), reduction = "umap" )

# Tissue-resident Macrophages (Cluster 8)
VlnPlot(subsce, features = c("Cd14","Mmp12","Cxc3r1","Il1b","Ccr2","H2-Eb1"), group.by = "seurat_clusters", pt.size = 0)
VlnPlot(subsce, features = c("Timd4","Lyve1","Folr2","Mrc1","Cd163","Ccr2","Igf1"), group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(subsce, features = c("Mgl2","Folr2","Adgre1"), reduction = "umap" )

# Polarization transitional states & Collagen subtypes
FeaturePlot(subsce, features = c("Arg1","Il1b","Ccl6","F10"), reduction = "umap" )
FeaturePlot(subsce, features = c("Col5a2","Col1a1","Col3a1","Dcn"), reduction = "umap" )

# ==============================================================================
# 8. Secondary Extraction (Fine-tuning Final Sub-populations)
# ==============================================================================
# Historical refinement selections
refinement_history <- list(
  fibro_final = c("7","24","18","20","9","10","11"),
  myeloid_final = c("5","12","13","18","15","9","1","2","3","4","20","0","10","14","11","21","17")
)

# Execute secondary extraction (Using Myeloid as active example based on the list length)
temp_sce <- subset(subsce, subset = seurat_clusters %in% refinement_history$myeloid_final)

temp_sce <- RunHarmony(temp_sce, c("orig.ident")) %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData(features = rownames(temp_sce)) %>% RunPCA() %>%
  FindNeighbors(dims = 1:30, features = VariableFeatures(object = temp_sce)) %>%
  FindClusters(resolution = 1) %>%
  RunTSNE(dims = 1:30, reduction = "harmony") %>% RunUMAP(dims = 1:30, reduction = "harmony") %>%
  RunTSNE(dims = 1:30, reduction = "pca") %>% RunUMAP(dims = 1:30, reduction = "pca")

DimPlot(temp_sce, reduction = "tsne", group.by = c("anno_human", "seurat_clusters"), label = TRUE)
DimPlot(temp_sce, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), label = TRUE)
DimPlot(temp_sce, reduction = "umap", split.by = c("orig.ident"), label = TRUE)

# -- Dynamic Proportion Area Plot --
cell_data_temp <- as.data.frame(table(Idents(temp_sce), temp_sce$orig.ident))
colnames(cell_data_temp) <- c("Cell_Type", "Timepoint", "Abundance")
ggplot(cell_data_temp, aes(Timepoint, Abundance, fill = Cell_Type)) +
  geom_area() + theme_minimal() + theme(legend.position = c(0.85, 0.85), legend.title = element_blank())

# --- Additional Feature Validations on Final `temp_sce` ---
FeaturePlot(temp_sce, features = c("Ly6a","Mt1","Mt2","Angptl4","Ccl2","Prg4","Csf1","Tgfb1","Cxcl1","Cxcl5"), reduction = "umap")
FeaturePlot(temp_sce, features = c("Postn","Cilp","Acta2","Ccnb2","Cdk1","Mki67"), reduction = "umap")
VlnPlot(temp_sce, features = c("Postn","Cilp","Acta2","Ccnb2","Cdk1","Mki67"), group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(temp_sce, features = c("Postn","Pdgfra","Hif1a","Vegfd","Vegfc","Mmp14","Angptl4","Dpep1","Dcn","Mmp2"), reduction = "umap")
VlnPlot(temp_sce, features = c("Postn","Pdgfra","Hif1a","Vegfd","Vegfc","Mmp14","Angptl4","Dpep1","Dcn","Mmp2"), group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(temp_sce, features = c("Postn","Cilp","Tnc","Fn1","Acta2","Cthrc1","Col1a1"), reduction = "umap")
FeaturePlot(temp_sce, features = c("Postn","Comp","Sfrp2","Tgfb1","Scx","Thbs4","Ccn4","Ccn5","Col8a2"), reduction = "umap")
FeaturePlot(temp_sce, features = c("Runx1","Runx2"), reduction = "umap", split.by = "orig.ident")

temp_sce2.markers <- FindAllMarkers(temp_sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(temp_sce2.markers, file = "newfibroblastmarkers.csv")

# ==============================================================================
# 9. Trajectory Initialization & Final Summary Visualizations
# ==============================================================================
# Minimum Spanning Tree exploration (Requires previously computed JS object)
# FIXME: 'JS' object dependency required for this block to execute.
# sce_exp <- as.SingleCellExperiment(JS)
# by.cluster <- aggregateAcrossCells(sce_exp, ids = colData(sce_exp)$seurat_clusters)
# centroids.cluster <- reducedDim(by.cluster, "UMAP")
# dmat.cluster <- as.matrix(dist(centroids.cluster))
# set.seed(1000)
# trajectory.cluster <- graph.adjacency(dmat.cluster, mode = "undirected", weighted = TRUE)
# plot(minimum.spanning.tree(trajectory.cluster))

# Core Trajectory Drivers Overlay
FeaturePlot(temp_sce, features = c("Stmn1"), split.by = "orig.ident", reduction = "umap")
FeaturePlot(temp_sce, features = c("Ly6c2","Plac8","Ccr2"), reduction = "umap")
FeaturePlot(temp_sce, features = c("Arg1","Il1b","Ccl6","F10"), reduction = "umap")
FeaturePlot(temp_sce, features = c("Stmn1","Cenpa","Pclaf"), reduction = "umap")
FeaturePlot(temp_sce, features = c("Trem2","Spp1","Gpnmb","Fabp5"), reduction = "umap")
FeaturePlot(temp_sce, features = c("Ccr2","H2-Aa","H2-Eb1","H2-DMa"), reduction = "umap")
FeaturePlot(temp_sce, features = c("Col5a2","Col1a1","Col3a1","Dcn"), reduction = "umap")

# Custom Stacked Violin Plot for Trajectory Genes
features_custom <- c("Ly6c2","Ccr2","Il1b","Arg1","Trem2","Gpnmb","Col5a2","Col1a1")
colors_custom <- c("#0072B2","#009E73","#D55E00","#CC79A7","#F0E442",
                   "#56B4E9","#E69F00","#00ADA9","#D0E429","#ED008C","#68217A")

p_custom_vln <- VlnPlot(temp_sce, features = features_custom, stack = TRUE, sort = FALSE, flip = TRUE, cols = colors_custom) +
  theme(legend.position = "none")
print(p_custom_vln)

# saveRDS(temp_sce, file = "temp_sce.rds")

# Genes of interest
features<-c("Ly6c2","Ccr2","Il1b","Arg1","Trem2","Gpnmb",
              "Col5a2","Col1a1")

# Define color scheme for visualization
colors<-c("#0072B2","#009E73","#D55E00","#CC79A7","#F0E442",
            "#56B4E9","#E69F00","#00ADA9","#D0E429","#ED008C","#68217A")
p2<-VlnPlot(temp_sce,features,stack=TRUE,
sort=F,flip=TRUE,cols=colors)+
theme(legend.position="none")
## Save final results
## saveRDS(temp_sce, file = "temp_sce.rds")
## saveRDS(temp_sce, file = "./temp_fibroblast_sce.rds")
## saveRDS(temp_sce, file = "成纤维_sce.rds") # Fibroblasts

temp_sce <- readRDS("/home/tsh/E-MTAB-7895/temp_sce.rds") # Fibroblasts
temp_sce <- readRDS("../E-MTAB-7895/成纤维_sce.rds") # Fibroblasts
temp_sce <- readRDS("../E-MTAB-7895/成纤维大亚群未筛选.rds") # Unfiltered fibroblast macro-population

getwd()
##saveRDS(temp_sce, file = "inflammation_temp_sce.rds") # Monocytes/Macrophages
##saveRDS(temp_sce, file = "filtered_inflammation_temp_sce.rds")
##saveRDS(temp_sce, file = "巨噬_sce.rds")
temp_sce <- readRDS("../E-MTAB-7895/inflammation_temp_sce.rds")
temp_sce <- readRDS("../E-MTAB-7895/filtered_inflammation_temp_sce.rds")
# temp_sce$biotip <- temp_sce$seurat_clusters
##saveRDS(subsce, file = "Granulocytes_temp_sce.rds") # Neutrophils / Granulocytes
##saveRDS(temp_sce, file = "中性粒_sce.rds") # Neutrophils
##saveRDS(temp_sce, file = "Granulocytes_temp_sce2.0.rds")
temp_sce <- readRDS("../E-MTAB-7895/Granulocytes_temp_sce.rds")


######### Immune Subpopulations - Extract Subsets of Interest #################################
# Exclude proliferative and antigen-presenting subpopulations, along with their respective 
# progenitor states, based on prior Monocle 3 pseudotime inference.
Idents(temp_sce) <- "seurat_clusters"
Idents(temp_sce) <- "inflammation_human"
wanted_clusters <- c("6","12","13","3","5","2","0","14","10","4","9")
temp_sce <- subset(temp_sce,subset = seurat_clusters %in% wanted_clusters)
#temp_sce <- subset(temp_sce,subset = seurat_clusters != c("Proliferation"))
temp_sce$"label" <- Idents(temp_sce)

## Filter ribosomal genes if necessary to mitigate technical artifacts
rb.genes <- rownames(temp_sce)[grep("^Rp[sl]",rownames(temp_sce))]
Ct <- GetAssayData(object = temp_sce, layer = "counts")
percent.ribo <- Matrix::colSums(Ct[rb.genes,])/Matrix::colSums(Ct)*100
temp_sce <- AddMetaData(temp_sce, percent.ribo, col.name = "percent.ribo")
VlnPlot(temp_sce, features = "percent.ribo", pt.size = 0.1 ) + NoLegend()
temp_sce <- subset(temp_sce,subset = percent.ribo<25)

temp_sce <- RunHarmony(temp_sce,c( "orig.ident" ))
temp_sce <- FindVariableFeatures(temp_sce)
temp_sce <- ScaleData(temp_sce,features = rownames(temp_sce))
temp_sce <- RunPCA(temp_sce)
temp_sce <- FindNeighbors(temp_sce, dims = 1:30, features = VariableFeatures(object = temp_sce))
temp_sce <- FindClusters(temp_sce, resolution = 1.2)
temp_sce <- FindClusters(temp_sce, resolution = 1)
temp_sce <- FindClusters(temp_sce, resolution = 0.8)
temp_sce <- FindClusters(temp_sce, resolution = 0.36)
temp_sce <- RunTSNE(temp_sce, dims = 1:30, reduction = "harmony")
temp_sce <- RunUMAP(temp_sce, dims = 1:30, reduction = "harmony")
temp_sce <- RunTSNE(temp_sce, dims = 1:30, reduction = "pca")
temp_sce <- RunUMAP(temp_sce, dims = 1:30, reduction = "pca")
DimPlot(temp_sce, reduction = "tsne", group.by = c("subanno_human", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("anno_human", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("inflammation_human", "seurat_clusters","label"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("anno_human","inflammation_human", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", split.by = c("orig.ident"),label = T)
Idents(temp_sce)
head(temp_sce@meta.data)
head(temp_sce@assays$RNA$scale.data)[1:5,1:5]
head(temp_sce[["RNA"]]$scale.data)[1:5,1:5]
getwd()
subset(temp_sce,seurat_clusters == "4")


# Save critical states (Dynamic Network Biomarkers tipping points)
saveRDS(subset(temp_sce,seurat_clusters == "2"),file = "./中性粒早期临界态.rds")
saveRDS(subset(temp_sce,seurat_clusters == "4"),file = "./巨噬后期临界态.rds")
subset_sce <- subset(temp_sce, seurat_clusters == "2")

write.csv( temp_sce %>%
             subset(seurat_clusters == "7") %>%
             `[[`("RNA") %>%
             `$`("data") %>% .[c(DNB_before,DNB_after),],"成纤维cluster7.csv")
which( c(DNB_before,DNB_after),rownames( temp_sce %>% subset(seurat_clusters == "2") %>%`[[`("RNA") %>% `$`("data")))

temp_sce.markers <- FindAllMarkers(temp_sce, only.pos = TRUE, min.pct =0.25, logfc.threshold = 0.25)

write.csv(temp_sce.markers, "Granulocytes2.0.csv", sep = ",", quote = FALSE, row.names = T, col.names = T)
write.csv(temp_sce.markers, "../临界态图/Table2-3成纤维marker.csv", sep = ",", quote = FALSE, row.names = T, col.names = T)
write.csv(temp_sce.markers, "./临界态图/Table2-2中性粒marker.csv", sep = ",", quote = FALSE, row.names = T, col.names = T)
FeaturePlot(temp_sce, features = c("Ly6c2","Arg1","Cxcl3","Ccr1","Trem2","Col3a1"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Angptl4","Ccnb2","Acta2","Sfrp2"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Cxcl3","Fn1","H2-DMa","Ms4a4c"),reduction = "umap" )

FeaturePlot(temp_sce, features =  c("Mt1", "Mt2", "Angptl4", "Ccl2", "Acta2", "Ccnb2", "Cdk1", "Mki67",
  "Fn1", "Acta2", "Chtrc1", "Col1a1", "Ccn5", "Comp", "Sfrp2", "Col8a2"),reduction = "umap" )
table(Idents(temp_sce),temp_sce$orig.ident)
table(temp_sce$inflammation_human,temp_sce$seurat_clusters)

bulk_data <- AverageExpression(temp_sce, assays = "RNA", return.seurat = TRUE, group.by = "seurat_clusters")
write.csv(bulk_data[["RNA"]]$data,"./JS_ALL_THREE_PATH.csv")
# Idents(temp_sce) <- "seurat_clusters"
# wanted_clusters <- c("5")
# temp_sce2 <- subset(temp_sce,subset = seurat_clusters %in% wanted_clusters)
# Idents(temp_sce2) <- "orig.ident"
# wanted_clusters <- c("14day","28day")
# temp_sce2 <- subset(temp_sce2,subset = orig.ident %in% wanted_clusters)
# VariableFeatures(temp_sce2)
# index <- which(temp_sce@meta.data["seurat_clusters"] == "5")
# temp_sce@meta.data[index,]["orig.ident"] <- "1day"

Idents(temp_sce) <- "seurat_clusters"
cell_data <-as.data.frame(prop.table(table(Idents(temp_sce),temp_sce$orig.ident),margin = 2))
cell_data <-as.data.frame(table(Idents(temp_sce),temp_sce$orig.ident),margin = 2)
colnames(cell_data) <- c("细胞类型", "时间", "细胞丰度")
plot <- ggplot(cell_data, aes(x = 时间, y = 细胞丰度, fill = 细胞类型)) +
  geom_col() +
  labs(x = "时间", y = "细胞数量", fill = "细胞类型") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "grey")

  ) +
  facet_wrap(~ 细胞类型, scales = "free_x") +
  geom_area()

plot


Idents(subsce) <- "seurat_clusters"
cell_data <-as.data.frame(prop.table(table(Idents(subsce),subsce$orig.ident),margin = 2))
cell_data <-as.data.frame(table(Idents(subsce),subsce$orig.ident),margin = 2)
colnames(cell_data) <- c("细胞类型", "时间", "细胞丰度")
plot <- ggplot(cell_data, aes(x = 时间, y = 细胞丰度, fill = 细胞类型)) +
  geom_col() +
  labs(x = "时间", y = "细胞数量", fill = "细胞类型") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "grey")

  ) +
  facet_wrap(~ 细胞类型, scales = "free_x") +
  geom_area()

plot


# Export the raw count matrix for all cells
write.table(temp_sce2[["RNA"]]$counts, "967counts.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Export the normalized data matrix for all cells
write.table(temp_sce2[["RNA"]]$data, "967data.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Export the scaled data matrix for all cells
write.table(temp_sce2[["RNA"]]$scale.data, "967scale_data.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Export the metadata matrix for all cells
write.table(temp_sce2@meta.data, "967meta.data.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)



### Further extraction of specific sub-clusters (Optional) ####
### ### Evaluate early-stage Fibroblasts: Stratified into 'before', 'after' critical states, and combined analysis
Idents(temp_sce) <- "seurat_clusters"
wanted_clusters <- c("1","4","9") # DNB 4 FIBROBLAST
log2fc_markers <- FindMarkers(temp_sce,ident.1 = "4",ident.2 = "1",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
DNB_before <- c("Adamts4", "Arc", "Atp1b3", "Bhlhe40", "Birc3", "Btg2", "Cebpb", "Cxcl1", "Egr2", "Egr3",
           "F2rl1", "Fetub", "Fosb", "Fosl2", "Gclc", "Gem", "Gm17268", "Gm26532", "Gm41609", "Has1",
           "Hspa1a", "Ier3", "Kcnip1", "Kdm6b", "Khdrbs3", "Klf4", "Lmo2", "Lypla2", "Maff", "Mpp7",
           "Nfkb1", "Nol12", "Nr4a1", "Nr4a3", "Opa3", "Pde4b", "Phlda1", "Pla1a", "Plekho2", "Ptgs2",
           "Sbno2", "Slc10a6", "Smad3", "Snai1", "Tgm2", "Tiparp", "Tnfaip2", "Tnfaip6", "Trib1",
           "Ubash3b", "Ugdh", "Zfp36")
sub_df <- log2fc_markers[rownames(log2fc_markers) %in% DNB_before, ]
write.table(sub_df, "Fibroblast_log2fc_4_1.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)

log2fc_markers <- FindMarkers(temp_sce,ident.1 = "9",ident.2 = "4",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
DNB_after <- c("Abce1", "Adam12", "Adamts9", "Atf3", "Cxcl2", "Depdc5", "Dot1l", "Ell2", "Fbxo33", "Hsph1",
           "Ier5", "Il6", "Irf3", "Nes", "Nfkbia", "Nr4a2", "Pim1", "Rgs2", "Rnf217", "Rras2", "Serpine1",
           "Sertad1", "Sgk1", "Sms")
sub_df <- log2fc_markers[rownames(log2fc_markers) %in% DNB_after, ]
write.table(sub_df, "Fibroblast_log2fc_9_4.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)

DNB <- unique(c(DNB_before, DNB_after))
log2fc_4_1 <- FindMarkers(temp_sce, ident.1 = "4", ident.2 = "1", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_9_4 <- FindMarkers(temp_sce, ident.1 = "9", ident.2 = "4", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_9_1 <- FindMarkers(temp_sce, ident.1 = "9", ident.2 = "1", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_9_1_DIF <- subset(log2fc_9_1, abs(avg_log2FC) > 2 & p_val_adj < 0.05)
write.table(log2fc_9_1_DIF, "Fibroblast_log2fc_early_DIFF.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)
sub_4_1 <- log2fc_4_1[rownames(log2fc_4_1) %in% DNB, "avg_log2FC", drop=FALSE]
sub_9_4 <- log2fc_9_4[rownames(log2fc_9_4) %in% DNB, "avg_log2FC", drop=FALSE]
sub_9_1 <- log2fc_9_1[rownames(log2fc_9_1) %in% DNB, "avg_log2FC", drop=FALSE]
all_log2fc <- merge(sub_4_1, sub_9_4, by="row.names", all=TRUE, suffixes = c("_4_1", "_9_4"))
all_log2fc <- merge(all_log2fc, sub_9_1, by.x="Row.names", by.y="row.names", all=TRUE)
colnames(all_log2fc) <- c("Gene", "log2fc_4_1", "log2fc_9_4", "log2fc_9_1")
write.table(all_log2fc, "Fibroblast_log2fc_early_all.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)

### ### Evaluate late-stage Fibroblasts: Stratified into 'before', 'after' critical states, and combined analysis
wanted_clusters <- c("9","6","7") # DNB 6 FIBROBLAST
log2fc_markers <- FindMarkers(temp_sce,ident.1 = "6",ident.2 = "9",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
DNB_before <- c("Eif1", "H19", "C1qtnf3", "Slit2", "Ncam1", "Pamr1", "Fbn2", "Olfml2b", "Grb10", "Nxn",
           "Tln2", "Dleu2", "Hs6st2", "Jag1", "Susd6", "Hmgcr", "Nrp2", "A930037H05Rik", "Kcnq3",
           "Fmn1", "Erbin", "Kdm5b", "Tmem263", "Osbpl8", "Ltbp1", "Ednra", "Flrt2", "Slc31a1")
sub_df <- log2fc_markers[rownames(log2fc_markers) %in% DNB_before, ]
write.table(sub_df, "Fibroblast_log2fc_6_9.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)

log2fc_markers <- FindMarkers(temp_sce,ident.1 = "7",ident.2 = "6",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
DNB_after <- c("Itm2a", "Fndc1", "Mdk", "Wfdc18", "Thbs2", "Mfap2", "Rin2", "Piezo2", "Col11a1", "Phldb2",
           "Prickle1", "Tsc22d1", "Scx", "Picalm", "Kcnma1", "4931406P16Rik", "Antxr1", "Nhs", "Plxna4",
           "Sulf2", "Neo1", "Cdh11", "4930523C07Rik", "Col16a1", "Rai14", "Ifitm10", "Serf1", "Strbp",
           "Trabd2b", "Etl4", "Aqp2", "Nuak1", "Mpped2", "Slc7a2", "Lurap1l", "Kcnj15", "Plaat3", "Itga11",
           "Unc119", "Rgl1", "Gm5617", "St3gal5", "Adamts17", "Nlk", "Gm4117", "Cd59a", "Ssc5d", "Acan",
           "Myo1e", "Pik3r3", "Tcea3", "Ramp1", "Jph2")
sub_df <- log2fc_markers[rownames(log2fc_markers) %in% DNB_after, ]
write.table(log2fc_markers, "Fibroblast_log2fc_7_6.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)


DNB <- unique(c(DNB_before, DNB_after))
log2fc_6_9 <- FindMarkers(temp_sce, ident.1 = "6", ident.2 = "9", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_7_6 <- FindMarkers(temp_sce, ident.1 = "7", ident.2 = "6", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_7_9 <- FindMarkers(temp_sce, ident.1 = "7", ident.2 = "9", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_7_9_DIF <- subset(log2fc_7_9, abs(avg_log2FC) > 2 & p_val_adj < 0.05)
write.table(log2fc_7_9_DIF, "Fibroblast_log2fc_late_DIFF.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)
sub_6_9 <- log2fc_6_9[rownames(log2fc_6_9) %in% DNB, "avg_log2FC", drop=FALSE]
sub_7_6 <- log2fc_7_6[rownames(log2fc_7_6) %in% DNB, "avg_log2FC", drop=FALSE]
sub_7_9 <- log2fc_7_9[rownames(log2fc_7_9) %in% DNB, "avg_log2FC", drop=FALSE]
all_log2fc <- merge(sub_6_9, sub_7_6, by="row.names", all=TRUE, suffixes = c("_6_9", "_7_6"))
all_log2fc <- merge(all_log2fc, sub_7_9, by.x="Row.names", by.y="row.names", all=TRUE)
colnames(all_log2fc) <- c("Gene", "log2fc_6_9", "log2fc_7_6", "log2fc_7_9")
write.table(all_log2fc, "Fibroblast_log2fc_late_all.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)


### ### Evaluate early-stage Neutrophils: Stratified into 'before', 'after' critical states, and combined analysis
wanted_clusters <- c("0","2","1") # DNB 2 Granulocyte
wanted_clusters <- c("0","5","7") # DNB 2 Granulocyte 2.0
log2fc_markers <- FindMarkers(temp_sce,ident.1 = "2",ident.2 = "0",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
log2fc_markers <- FindMarkers(temp_sce,ident.1 = "5",ident.2 = "0",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
DNB_before <- c("3-Mar", "Arg1", "Ccl24", "Ccl7", "Cd38", "Cox4i2", "Crem", "Cyfip1", "Errfi1", "Fcgr2b",
           "Inhba", "Mmp19", "Mt2", "Ormdl2", "Pde4d", "Pdpn", "Pgm1", "Ppbp", "Pxk", "Spp1", "Tfrc",
           "Uck2", "Vps8", "Yipf5")
DNB_before <- c ("2010310C07Rik", "4921511C10Rik", "9830107B12Rik", "Acss2", "Ankrd13c", "Arap3", "B430010I23Rik", "Btbd10", "Ccdc138", "Ccdc180", "Ccl21a", "Cfh", "Cox6a2", "Cpd",
  "Dars2", "Eif2ak3", "Esam", "F5", "Frmd5", "Gm36738", "Hdac4", "Hspa4l", "Il1bos", "Itgbl1", "Mfap4", "Mtfr1", "Mxra7", "Nbeal2", "Rflnb", "Rfx2", "Satb1", "Serping1", "Sh2d3c", "Sp2", "Strbp", "Sun2", "Tgds", "Zhx2")

sub_df <- log2fc_markers[rownames(log2fc_markers) %in% DNB_before, ]
write.table(sub_df, "Granulocyte_log2fc_2_0.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(sub_df, "Granulocyte_log2fc_5_0.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)

log2fc_markers <- FindMarkers(temp_sce,ident.1 = "1",ident.2 = "2",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
log2fc_markers <- FindMarkers(temp_sce,ident.1 = "7",ident.2 = "5",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
DNB_after <- c("AU022252", "B4galt1", "Ccl12", "Ccl2", "Ccl9", "Cd86", "Chd9", "Clec4a3", "Dab2", "Dapk1",
               "Ear2", "Ecm1", "Eif2s2", "Emp1", "Exosc3", "F10", "F13a1", "Fam172a", "Fcgr1", "Fn1", "Glul",
               "Gm47283", "Gpr171", "Hes1", "Hmox1", "Imp4", "Klf4", "Lmna", "Lpl", "Mafb", "Man1a", "Mrc1",
               "Ms4a4c", "Ms4a6b", "Ms4a6d", "Mt1", "Myo1e", "Nr4a2", "Nrg1", "Papss2", "Pdia6", "Pf4", "Plcb1",
               "Prdx2", "Rap1gap2", "Slc16a10", "Smpdl3b", "Snx5", "Tbxas1", "Tubb6", "Vcan")
DNB_after <- c("Abhd15", "Adamts2", "Agpat2", "Bloc1s1", "Ccpg1os", "Cdc27", "Cdh5", "Col6a2", "Col8a1", "Csad", "Dpt", "Gcnt7", "Gm38843", "Gmcl1", "Hba-a1", "Hps4", "Khnyn",
               "Lypd6b", "Oplah", "Ralgps1", "Rhof", "Rnf141", "Scai", "Slc25a33", "Stfa2", "Svip", "Tnxb", "Tspan2", "Ttc7")

sub_df <- log2fc_markers[rownames(log2fc_markers) %in% DNB_after, ]
write.table(sub_df, "Granulocyte_log2fc_1_2.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(sub_df, "Granulocyte_log2fc_7_5.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)


DNB <- unique(c(DNB_before, DNB_after))
log2fc_2_0 <- FindMarkers(temp_sce, ident.1 = "2", ident.2 = "0", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_1_2 <- FindMarkers(temp_sce, ident.1 = "1", ident.2 = "2", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_1_0 <- FindMarkers(temp_sce, ident.1 = "1", ident.2 = "0", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_1_0_DIF <- subset(log2fc_1_0, abs(avg_log2FC) > 2 & p_val_adj < 0.05)
write.table(log2fc_1_0_DIF, "Granulocyte_log2fc_early_DIFF.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)

sub_2_0 <- log2fc_2_0[rownames(log2fc_2_0) %in% DNB, "avg_log2FC", drop=FALSE]
sub_1_2 <- log2fc_1_2[rownames(log2fc_1_2) %in% DNB, "avg_log2FC", drop=FALSE]
sub_1_0 <- log2fc_1_0[rownames(log2fc_1_0) %in% DNB, "avg_log2FC", drop=FALSE]
all_log2fc <- merge(sub_2_0, sub_1_2, by="row.names", all=TRUE, suffixes = c("_2_0", "_1_2"))
all_log2fc <- merge(all_log2fc, sub_1_0, by.x="Row.names", by.y="row.names", all=TRUE)
colnames(all_log2fc) <- c("Gene", "log2fc_2_0", "log2fc_1_2", "log2fc_1_0")
write.table(all_log2fc, "Granulocyte_log2fc_early_all.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)


### Evaluate early-stage Monocytes/Macrophages: Stratified into 'before', 'after' critical states, and combined analysis
wanted_clusters <- c("12","5","8") # DNB 5 mono/macrophages
log2fc_markers <- FindMarkers(temp_sce,ident.1 = "5",ident.2 = "12",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
DNB_before <- c("2610037D02Rik", "A530064D06Rik", "Atf4", "Atp9b", "Bst1", "Ccr1", "Cst3", "Ctss", "Dnajb9",
           "Emilin2", "Esyt2", "F10", "F5", "Fam102b", "Fbxo28", "Fcer1g", "Fndc3a", "Gadd45a", "Gm15283",
           "Gpat3", "Gpr141", "Grina", "Herpud1", "Hk1", "Il13ra1", "Inhba", "Lpar1", "Lrrc4", "Lyz2",
           "Maml2", "Manf", "Map2k1", "Nedd9", "Nucb2", "P4hb", "Pcgf5", "Pfkp", "Piwil2", "Prdx6", "Resf1",
           "Rfx3", "Rgcc", "Rnf11", "Rsbn1", "Sc5d", "Sec23b", "Sec61a1", "Serpinb2", "Sgms2", "Slbp",
           "Slc16a3", "Slc23a2", "Slc2a1", "Slc35b1", "Slfn4", "Slpi", "Smox", "Snx18", "Spidr", "Ssr3",
           "Tes", "Tlr4", "Tmem248", "Trem1", "Tuba1a", "Tyrobp", "Upp1", "Uso1", "Zfp292")
sub_df <- log2fc_markers[rownames(log2fc_markers) %in% DNB_before, ]
write.table(sub_df, "Macropahges_log2fc_5_12.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)
log2fc_markers <- FindMarkers(temp_sce,ident.1 = "8",ident.2 = "5",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
DNB_after <- c("Bnip3", "Bnip3l", "Ccdc47", "Ccnd2", "Cxcl3", "Dok2", "Egln3", "Eif4ebp1", "Ero1l",
           "Fam162a", "Fbxl5", "Gbe1", "Gm11290", "Gm47507", "Higd1a", "Impa2", "Klf7", "Lyar",
           "Me2", "Mt2", "Mthfd1l", "Ndrg1", "P4ha1", "Pdia4", "Pgm1", "Pilra", "Rbpms",
           "Selenbp1", "Serpine1", "Tdrd7")
sub_df <- log2fc_markers[rownames(log2fc_markers) %in% DNB_after, ]
write.table(log2fc_markers, "Macropahges_log2fc_8_5.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)

DNB <- unique(c(DNB_before, DNB_after))
log2fc_5_12 <- FindMarkers(temp_sce, ident.1 = "5", ident.2 = "12", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_8_5 <- FindMarkers(temp_sce, ident.1 = "8", ident.2 = "5", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_8_12 <- FindMarkers(temp_sce, ident.1 = "8", ident.2 = "12", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_8_12_DIF <- subset(log2fc_8_12, abs(avg_log2FC) > 2 & p_val_adj < 0.05)
write.table(log2fc_8_12_DIF, "Macropahges_log2fc_early_DIF.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)
sub_5_12 <- log2fc_5_12[rownames(log2fc_5_12) %in% DNB, "avg_log2FC", drop=FALSE]
sub_8_5 <- log2fc_8_5[rownames(log2fc_8_5) %in% DNB, "avg_log2FC", drop=FALSE]
sub_8_12 <- log2fc_8_12[rownames(log2fc_8_12) %in% DNB, "avg_log2FC", drop=FALSE]

all_log2fc <- merge(sub_5_12, sub_8_5, by="row.names", all=TRUE, suffixes = c("_5_12", "_8_5"))
all_log2fc <- merge(all_log2fc, sub_8_12, by.x="Row.names", by.y="row.names", all=TRUE)
colnames(all_log2fc) <- c("Gene", "log2fc_5_12", "log2fc_8_5", "log2fc_8_12")
write.table(all_log2fc, "Macropahges_log2fc_early_all.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)

### Evaluate late-stage Monocytes/Macrophages: Stratified into 'before', 'after' critical states, and combined analysis
wanted_clusters <- c("3","4","1") # DNB 4 mono/macrophages
log2fc_markers <- FindMarkers(temp_sce,ident.1 = "4",ident.2 = "3",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
DNB_before <- c( "Fabp5", "Malat1", "Rgs1", "Lpl", "Pdpn", "Lmna", "Fabp4", "Cd36", "Adam8", "Airn",
           "Il7r", "Atp6v1a", "Tnfsf13", "Gusb", "Uhrf2", "Plk2", "Soat1", "Dennd4c", "Cxcl14", "Psd3",
           "Ankrd33b", "Atp6v1b2", "Ftl1-ps1", "Atp6v1c1", "Bcl2l1", "Fam20c", "Slc48a1", "Uap1l1",
           "Plgrkt", "Samd8", "Ap2m1", "Slc9a3r1", "F7", "Amdhd2", "Mpp1", "Lhfpl2", "Tpm1", "Pcna",
           "Opa1", "Ndufs6")
sub_df <- log2fc_markers[rownames(log2fc_markers) %in% DNB_before, ]
write.table(sub_df, "Macropahges_log2fc_4_3.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)

log2fc_markers <- FindMarkers(temp_sce,ident.1 = "1",ident.2 = "4",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
DNB_after <- c("Gpnmb", "Ccl8", "Cd9", "Syngr1", "Igf1", "Cd72", "Serpinb6a", "Mmp12", "Ckb", "Pld3",
           "Hist1h2bc", "Acp5", "Gdf15", "Fcgr4", "Camk1", "Aopep", "Wwp1", "Stab1", "Slamf9", "Hist1h1c",
           "Lgals3bp", "Myo1e", "Gas6", "Lipa", "Ptms", "Slamf7", "Pmp22", "Abhd12", "Stmn1", "Folr2",
           "Cd84", "Mfge8", "Blnk", "Cyb5a", "Cfb", "Pea15a", "S100a1", "Ifi27", "Gpr137b", "Tmem140",
           "Rhoc", "Selenoh", "Tspan4", "Bhlhe41", "Atp6v0d2", "Cadm1", "Man1c1", "Sptssa", "Selenom",
           "Tle5", "Prdx4", "Fam219a", "Cela1", "Ttyh2", "Mtss1", "B4galt6", "Slc25a4", "Dnmt3a", "Naa38",
           "Pkib", "Hebp1", "Acss1", "Pdgfc", "Npy", "Fbxw4", "Hmgn5", "Glmp", "Nenf", "Il3ra", "Tlk1",
           "Itgb1", "Rabggtb", "Map4k3", "Osbpl11")
sub_df <- log2fc_markers[rownames(log2fc_markers) %in% DNB_after, ]
write.table(sub_df, "Macropahges_log2fc_1_4.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)

DNB <- unique(c(DNB_before, DNB_after))
log2fc_4_3 <- FindMarkers(temp_sce, ident.1 = "4", ident.2 = "3", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_1_4 <- FindMarkers(temp_sce, ident.1 = "1", ident.2 = "4", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_1_3 <- FindMarkers(temp_sce, ident.1 = "1", ident.2 = "3", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
log2fc_1_3_DIF <- subset(log2fc_1_3, abs(avg_log2FC) > 2 & p_val_adj < 0.05)
write.table(log2fc_1_3_DIF, "Macropahges_log2fc_late_DIF.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)
sub_4_3 <- log2fc_4_3[rownames(log2fc_4_3) %in% DNB, "avg_log2FC", drop=FALSE]
sub_1_4 <- log2fc_1_4[rownames(log2fc_1_4) %in% DNB, "avg_log2FC", drop=FALSE]
sub_1_3 <- log2fc_1_3[rownames(log2fc_1_3) %in% DNB, "avg_log2FC", drop=FALSE]
all_log2fc <- merge(sub_4_3, sub_1_4, by="row.names", all=TRUE, suffixes = c("_4_3", "_1_4"))
all_log2fc <- merge(all_log2fc, sub_1_3, by.x="Row.names", by.y="row.names", all=TRUE)
colnames(all_log2fc) <- c("Gene", "log2fc_4_3", "log2fc_1_4", "log2fc_1_3")
write.table(all_log2fc, "Macropahges_log2fc_late_all.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)


### Evaluate early-stage Monocytes/Macrophages: Stratified into 'before', 'after', and combined analysis.
### Discarded: This was ultimately not utilized for early critical state evaluation.
wanted_clusters <- c("5","8","3") # DNB 8 mono/macrophages
log2fc_markers <- FindMarkers(temp_sce,ident.1 = "8",ident.2 = "3",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
write.table(log2fc_markers, "Granulocyte_log2fc_8_3.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)
log2fc_markers <- FindMarkers(temp_sce,ident.1 = "1",ident.2 = "2",only.pos = FALSE, min.pct =0, logfc.threshold = 0)
write.table(log2fc_markers, "Granulocyte_log2fc_3_8.txt", sep = "\t", quote = FALSE, row.names = T, col.names = T)


rownames(log2fc_9_1_DIF)  # Fibro_early
rownames(log2fc_7_9_DIF)  # Fibro_late
rownames(log2fc_1_0_DIF)  # Gran_early
rownames(log2fc_8_12_DIF) # Macro_early
rownames(log2fc_1_3_DIF)  # Macro_late

# Convert vectors into a unified list structure
venn_data <- list(DNB = DNB, Fibro_early = rownames(log2fc_9_1_DIF))

# Render and customize Venn diagram for gene intersection
ggVennDiagram(venn_data, set_size = 6) +
  scale_fill_gradient(low = "lightblue", high = "dodgerblue") +  # Configure gradient fill mapping
  theme(legend.position = "none")  +                               # Hide primary legend
  theme(
    legend.position = "none",
    text = element_text(size = 2),  # Adjust global text size for compact presentation
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Increase margins to prevent clipping
  ) +
  scale_x_continuous(expand = expansion(mult = .5))
##### Secondary Sub-clustering
temp_sce2 <- subset(temp_sce,subset = seurat_clusters %in% wanted_clusters)
# temp_sce2 <- subset(subsce,subset = seurat_clusters %in% wanted_clusters)
temp_sce2$"label" <- temp_sce2$seurat_clusters
Idents(temp_sce2) <- "label"
temp_sce2 <- RunHarmony(temp_sce2,c( "orig.ident" ))
temp_sce2 <- FindVariableFeatures(temp_sce2)
temp_sce2 <- ScaleData(temp_sce2,features = rownames(temp_sce2))
temp_sce2 <- RunPCA(temp_sce2)
temp_sce2 <- FindNeighbors(temp_sce2, dims = 1:30, features = VariableFeatures(object = subsce))
temp_sce2 <- FindClusters(temp_sce2, resolution = 0.36)
temp_sce2 <- RunTSNE(temp_sce2, dims = 1:30, reduction = "pca")
temp_sce2 <- RunUMAP(temp_sce2, dims = 1:30, reduction = "pca")
DimPlot(temp_sce2, reduction = "tsne", group.by = c("label", "seurat_clusters"),label = T)
DimPlot(temp_sce2, reduction = "umap", group.by = c("label", "seurat_clusters","orig.ident"),label = T)
table(Idents(temp_sce2),temp_sce2$orig.ident)
cell_data <- as.data.frame(prop.table(table(Idents(temp_sce2),temp_sce2$orig.ident),margin = 2))
DotPlot(object = temp_sce2, features =Proliferation)+RotatedAxis()+scale_x_discrete("")+scale_y_discrete("")
FeaturePlot(temp_sce2, features = c("Mt1","Mt2","Angptl4","Ccl2","Il6","Mif"),reduction = "umap" )
FeaturePlot(temp_sce2, features = c("Stmn1","Tubb3","Acta2","Tpm2","Cthrc1","H2afz","Ckap4"),reduction = "umap" )
FeaturePlot(temp_sce2, features = c("Pdgfra","Hif1a","Vegfd","Mmp14"),reduction = "umap" )
FeaturePlot(temp_sce2, features = c("Comp","Sfrp2","Col1a1","Dkk3"),reduction = "umap" )
temp_sce2.markers <- FindMarkers(temp_sce, ident.1 = "6", ident.2 = c("9","7"),only.pos = TRUE, min.pct =0.25, logfc.threshold = 0.25)
temp_sce2.markers <- FindAllMarkers(temp_sce2, only.pos = TRUE, min.pct =0.25, logfc.threshold = 0.25)
write.csv(temp_sce2.markers,file = "成纤维967差异基因.csv")
write.csv(temp_sce2.markers,file = "temp_sce2markers.csv")

temp_sce2.markers <- FindAllMarkers(temp_sce, only.pos = FALSE, min.pct =0, logfc.threshold = 0)
write.csv(temp_sce2.markers,file = "all_unfiltered_markers.csv")


############# Slingshot Pseudotime Trajectory Inference #############################
# Convert the Seurat object to a SingleCellExperiment format required by Slingshot
filter_sce_trans <- as.SingleCellExperiment(temp_sce)
Slingshot <- slingshot(filter_sce_trans,
                       clusterLabels = 'seurat_clusters',  # Specify the colData column containing cluster annotations
                       reducedDim = 'UMAP',
                       reweight = FALSE, # Ensure trajectories do not overlap
                       start.clus= "10",  # Define the developmental root (starting node)
                       end.clus = NULL     # Let the algorithm infer the terminal nodes dynamically
)
colnames(colData(Slingshot))
plot.new()
dev.off()
plot(reducedDims(Slingshot)$UMAP, pch=16, asp = 1)
lines(SlingshotDataSet(Slingshot), lwd=1, col=brewer.pal(9,"Set1"))
legend("right",
       legend = paste0("lineage",1:1),
       col = unique(brewer.pal(3,"Set1")),
       inset=0.8,
       pch = 16)
# Construct continuous color gradients along pseudotime
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) # Interpolate 100 color gradients to simulate a continuous transition
plotcol <- colors[cut(Slingshot$slingPseudotime_1, breaks=100)] # Discretize Lineage 1 pseudotime into 100 intervals, mapping cells to corresponding gradient colors
plotcol <- colors[cut(Slingshot$slingPseudotime_2, breaks=100)]


plotcol[is.na(plotcol)] <- "lightgrey" # Assign light grey to points outside the active trajectory (NA values)
Slingshot$plotcol <- plotcol
plot.new()
plot(reducedDims(Slingshot)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(Slingshot), lwd=1, col=brewer.pal(9,"Set1"))
legend("right",
       legend = paste0("lineage",1:1),
       col = unique(brewer.pal(3,"Set1")),
       inset=0.8,
       pch = 16)

#################### Slingshot Sub-trajectory Extraction & Analysis ############################################
# The colData of the Slingshot object is not natively in data.frame format; reconstruction is required.
# Note: Manual extraction for downstream compatibility.
coldata <- colData(Slingshot)
coldata <- data.frame(celltype = coldata@listData$RNA_snn_res.0.36,
                      sampleId = coldata@listData$orig.ident,
                      plotcol = coldata@listData$plotcol)
rownames(coldata) = Slingshot@colData@rownames

# Isolate cells corresponding to Lineage 1
# Filter cell barcodes
filter_cell <- dplyr::filter(coldata, plotcol != "lightgrey")
filter_cell <- rownames(filter_cell)
head(filter_cell)
# Extract subset from count matrix
counts <- Slingshot@assays@data@listData$counts
filter_counts <- counts[,filter_cell]
filter_counts[1:5,1:5]
# Randomly downsample to 2,000 cells (adjustable based on experimental design)
# Mitigates computational overhead and processing time
#set.seed(111)
#scell <- sample(colnames(filter_counts), size = 2000)

# Final subsampled count matrix
filter_counts = filter_counts[, scell]
dim(filter_counts)

# Reconstruct the SingleCellExperiment object with the filtered matrix
filter_sim <- SingleCellExperiment(assays = List(counts = filter_counts))

# Map colData to the filtered object
filter_coldata = colData(Slingshot)[colnames(filter_counts), 1:3]
filter_sim@colData = filter_coldata

# Map dimensionality reduction embeddings
rd = reducedDim(Slingshot)
filter_rd <- rd[colnames(filter_counts),]
reducedDims(filter_sim) <- SimpleList(UMAP = filter_rd)

# K-Means clustering for sub-trajectory anchor points
set.seed(111)
cl <- kmeans(filter_rd, centers = 6)$cluster # Specify the target number of clusters
head(cl)
colData(filter_sim)$kmeans <- cl

## Visualize the K-Means clustering results
library(RColorBrewer)
mycolors = brewer.pal(6,"Set1") # colors

plt_k = data.frame(filter_rd,
                   kmeans_clusters = factor(cl, levels = sort(unique(cl))))
ggplot(plt_k, aes(PC_1, PC_2))+
  geom_point(aes(color = kmeans_clusters), size = .5)+
  scale_color_manual(values = mycolors)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab('PC_1')+
  ylab('PC_2')+
  guides(color = guide_legend(override.aes = list(size = 5)) # Adjust legend point size for clarity
  )



# Perform trajectory inference utilizing the defined K-Means subpopulations
filter_sim <- slingshot(filter_sim,
                        clusterLabels = 'kmeans',
                        reducedDim = 'UMAP',
                        start.clus= "2", # Define biological root/start cluster
                        end.clus = '5' # Define terminal state/end cluster
)
head(colnames(colData(filter_sim)))

plot(reducedDims(filter_sim)$UMAP, pch=16, asp = 1)
lines(SlingshotDataSet(filter_sim), lwd=2, col=mycolors)


#################### Pseudotime-associated Gene Analysis (tradeSeq) ##########################
BiocManager::install("tradeSeq")
library(tradeSeq)
# Fit negative binomial model
counts <- filter_sim@assays@data$counts
crv <- SlingshotDataSet(filter_sim)
# Prior to NB-GAM model fitting, determine the optimal number of basis functions (knots).
# Knots play a critical role in the smoothing spline estimation.
# Evaluate optimal K using evaluateK. Note: Computationally intensive.
# ~16 mins for 2k cells (duration scales with the number of discrete trajectories).
# Automatically generates optimization plots upon completion.
set.seed(111)
icMat <- evaluateK(counts = counts,
                   sds = crv,
                   k = 3:10,    # no more than 12
                   nGenes = 200, # Number of genes randomly sampled for knot evaluation (default 500, reduced here for efficiency)
                   verbose = T)




# we pick nknots = 6.
set.seed(111)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
# fit negative binomial GAM
# 2k cells ~13 min
# Utilize system.time() to profile computational overhead
system.time({
  sce <- fitGAM(counts = counts,
                pseudotime = pseudotime,
                cellWeights = cellWeights,
                nknots = 6,
                verbose = FALSE)
})
table(rowData(sce)$tradeSeq$converged)
assoRes <- associationTest(filter_sce_trans)
head(assoRes)




############# Monocle 3 Trajectory Inference #############################
library(devtools)
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(Seurat)
library(leidenbase)
DefaultAssay(temp_sce) <- "RNA"
DefaultAssay(subsce) <- "RNA"
# Construct data structures required for Monocle 3 CDS
# expression_matrix: Feature (gene) by Cell count matrix
expression_matrix = temp_sce[["RNA"]]$counts
expression_matrix = subsce[["RNA"]]$counts
# cell_metadata: Extract meta.data directly from the Seurat object
cell_metadata = data.frame(temp_sce@meta.data)
cell_metadata = data.frame(subsce@meta.data)
# gene_metadata: Construct feature annotations (Gene symbols)
gene_annotation = data.frame(expression_matrix[,1])
gene_annotation[,1] = row.names(gene_annotation)
colnames(gene_annotation)=c("gene_short_name")
## Instantiate the Monocle 3 cell_data_set (CDS) object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
# Data preprocessing and normalization
cds <- preprocess_cds(cds, num_dim = 50,norm_method = c("log"))

## Dimensionality reduction (Defaulting to UMAP coordinates)
cds <- reduce_dimension(cds,reduction_method="UMAP",cores=5)
plot_cells(cds,color_cells_by = "seurat_clusters")
## Cluster cells (Resolution downscaled intentionally for consolidated broad subpopulation visualization)
#cds <- cluster_cells(cds,resolution = 0.0000001)
cds <- cluster_cells(cds,cluster_method = "louvain", reduction_method ="UMAP")
## Infer trajectory principal graph
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster=FALSE, label_leaves=FALSE,label_branch_points=FALSE)



## Designate specific nodes as the developmental root
## e.g., assuming 'HEpiD' represents the progenitor/starting state
myselect <- function(cds,select.classify,my_select){
  cell_ids <- which(colData(cds)[,select.classify] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes}


cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "6"))
cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "11"))
cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "10"))
## Integrate native Seurat UMAP/Harmony coordinates to ensure spatial consistency across figures
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(temp_sce , reduction = "umap")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(temp_sce , reduction = "harmony")
int.embed <- Embeddings(subsce , reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

#root_group = colnames(cds)[clusters(cds) == 1]
#cds = order_cells(cds, root_cells = root_group)

## Map calculated pseudotime values across distinct cell types
## Note: Higher pseudotime theoretically indicates higher differentiation maturity (demonstrative purpose)
library(ggsci)
plot_cells(cds,color_cells_by = "pseudotime",
            #group_cells_by = "new.labels",
            #color_cells_by = "orig.ident",
            show_trajectory_graph=F) + plot_cells(cds,
                                                  #color_cells_by = "seurat_clusters",
                                                  color_cells_by = "seurat_clusters",
                                                  #color_cells_by = "label",
                                                  #color_cells_by = "orig.ident",
                                                  label_cell_groups=FALSE,
                                                  label_leaves=FALSE,
                                                  label_branch_points=FALSE,
                                                  graph_label_size=1)+ scale_color_manual(values = pal_jco("default", alpha = 0.6)(19))
# Differential gene expression along pseudotime
Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=6)
# Extract top driving genes based on Moran's I autocorrelation metric
Track_genes_sig <- Track_genes %>%top_n(n=30, morans_I) %>% pull(gene_short_name) %>% as.character()
pdf('./临界态图/Figure2-3E.pdf',height = 18,width = 12,onefile = F)
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by="seurat_clusters",min_expr=0.5, ncol= 2,cell_size=1.5) + scale_color_manual(values = pal_jco("default", alpha = 0.6)(11))
dev.off()
# Visualize specific user-defined genes along the trajectory axis
plot_cells(cds, genes=c("Mt2","Mt1","Cxcl3"),show_trajectory_graph=FALSE)#IR

# Scatter plot of gene expression kinetics across pseudotime
plot_genes_in_pseudotime(cds[c("Vegfa"),],
                          color_cells_by="seurat_clusters",
                          min_expr=0.5, ncol= 2,cell_size=1.5)+ scale_color_manual(values = pal_jco("default", alpha = 0.6)(7))


############################ Monocle 2 Trajectory Inference #####################################
# Convert the Seurat count matrix to a sparse matrix format to initiate Monocle 2
library(monocle)
temp_sce2 <- temp_sce
temp_sce2 <- subset(temp_sce,seurat_clusters %in% c("12","5","8","3","4","1"))
temp_sce2 <- subset(temp_sce,seurat_clusters %in% c("1","4","9","6","7","0"))
Idents(temp_sce2) <- "inflammation_human"
sparse_data <-  as(as.matrix(temp_sce2[["RNA"]]$counts),"sparseMatrix")
# Construct the phenotype data matrix
mdata <- new("AnnotatedDataFrame",data = temp_sce2@meta.data)
# Construct the feature data matrix containing gene_short_name
fData <- data.frame(gene_short_name=row.names(sparse_data),row.names =row.names(sparse_data))
fd <- new("AnnotatedDataFrame",data = fData)
# Instantiate the Monocle CellDataSet (CDS)
monocle_cds <- newCellDataSet(cellData = sparse_data,
                              phenoData = mdata,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
# Estimate size factors and empirical dispersions
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
# Filter lowly expressed cells based on defined threshold
monocle_cds <- detectGenes(monocle_cds,min_expr = 3)

# Identify highly variable genes governing the trajectory
# Utilize the robust variable features previously selected by Seurat
expressed_genes <- VariableFeatures(temp_sce2)
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01)) ## Apply strict q-value thresholding (qval < 0.01) instead of the default 0.1.
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)

# Dimensionality Reduction via DDRTree
library(igraph)
install.packages("igraph")
install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz",repos = NULL)
monocle_cds <- reduceDimension(monocle_cds,max_components = 2,reduction_method = "DDRTree")
# Order cells along the inferred pseudotime manifold
monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds,root_state = "2")
monocle_cds <- orderCells(monocle_cds,root_state = "1")
# Visualization
# Note: Color aesthetics can also be mapped directly to 'Pseudotime'
color <- c(
  "1"  = "#4DBBD5",  # Cold blue (Reparative)
  "3"  = "#00A087",  # Teal/Blue-green (Reparative)
  "4"  = "#3C5488",  # Dark blue (Reparative)
  "5"  = "#F39C12",  # Orange (Pro-inflammatory)
  "8"  = "#E64B35",  # Red-orange (Pro-inflammatory)
  "12" = "#D95F02"   # Deep warm orange (Pro-inflammatory)
)
color <- c(
  # Inflammatory (red) - High contrast mapping: Deep red vs Light red
  "1" = "#B2182B",  # Deep red
  "4" = "#EF8A62",  # Light coral
  
  # Proliferative (green)
  "9" = "#1B9E77",  # Standard teal/green (often used in top-tier journals)
  
  # Myofibroblast (blue) - Three-tier gradient (Dark/Mid/Light)
  "6" = "#08306B",  # Very dark blue
  "7" = "#2171B5",  # Mid blue
  "0" = "#6BAED6"   # Light blue
)
p <- plot_cell_trajectory(monocle_cds,color_by = "Pseudotime",cell_size = 0.6,show_backbone = TRUE) + facet_wrap("~orig.ident",nrow=1)
p0 <- plot_cell_trajectory(monocle_cds,color_by = "seurat_clusters",size = 1,show_backbone = TRUE) + facet_wrap("~orig.ident",nrow=1)
p1 <- plot_cell_trajectory(monocle_cds,color_by = "seurat_clusters",cell_size = 0.6,show_backbone = TRUE)+
  theme(legend.position = 'none',panel.border = element_blank())+
  scale_color_manual(values= color) # Apply custom aesthetic mappings
plot_cell_trajectory(monocle_cds,color_by = "Pseudotime",cell_size = 0.6,show_backbone = TRUE)+
  theme(legend.position = 'none',panel.border = element_blank())

# Append branching dendrogram visualization
p2 <- plot_complex_cell_trajectory(monocle_cds,x=1,y=2,color_by = "label")+
      theme(legend.title = element_blank())
      #+scale_color_manual(values= color) # Custom aesthetic mappings
plot_cell_trajectory(monocle_cds,color_by = "State",size = 1,show_backbone = TRUE)
p0|p|p1|p2

plot_cell_trajectory(monocle_cds,color_by = "inflammation_human",size = 1,show_backbone = TRUE) + facet_wrap("~orig.ident",nrow=1)
plot_complex_cell_trajectory(monocle_cds,x=1,y=2,color_by = "inflammation_human")+
  theme(legend.title = element_blank())

#+scale_color_manual(values= color) # Custom aesthetic mappings
## Identify pseudotime-dependent gene expression cascades
Time_diff <- differentialGeneTest(monocle_cds[ordering_genes,],cores =1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
# Filter for genes expressed in a minimum threshold of cells (n >= 10)
Time_diff <- Time_diff %>%
  filter(num_cells_expressed > 10)
# Isolate the top 100 most significantly varying genes along the trajectory
Time_genes <- top_n(Time_diff,n=100,desc(qval)) %>%  pull(gene_short_name) %>% as.character()
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p4 <- plot_pseudotime_heatmap(monocle_cds[Time_genes,],num_clusters = 3,show_rownames = T,return_heatmap = T)
my_color <- colorRampPalette(
  c("#3DAEB7", "#F5EDF0", "#DB643E")
)(100)

plot_pseudotime_heatmap(
  monocle_cds[Time_genes, ],
  num_clusters = 4,
  cores = 1,
  show_rownames = FALSE,
  hmcols = my_color
)

## Extract trajectory-associated genes stratified by hierarchical clustering
clusters <- cutree(p4$tree_row,k=4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_clusters"
table(clustering)

##################### DESTINY (Diffusion Maps) Trajectory #############################
BiocManager::install("destiny")
library(destiny)
library(Biobase)
data(guo_norm)



### Automatically construct Biobase-compatible data structures from the Seurat object
library(Biobase)
ct <-GetAssayData(object = temp_sce)
ct<-ct[VariableFeatures(temp_sce),]
ct <- as.ExpressionSet(as.data.frame(t(ct)))
# Append metadata annotations
#. Annotations can be accessed directly via ct$column and ct[['column']].
ct$celltype <- DPT@meta.data[,c("integrated_merge_cluster")]
dm <- DiffusionMap(ct,k = 10)
palette(cube_helix(4)) # Configure continuous color palette
plot(dm, pch = 20, # pch adjusted for aesthetic rendering
     col_by = "celltype")

## Manually construct the ExpressionSet object
ct <- ExpressionSet(assayData=as.matrix(temp_sce[["RNA"]]$counts),
                     phenoData= mdata,featureData = fd)
dm <- DiffusionMap(ct,k = 3,n_pcs = 50)
plot(dm)
palette(cube_helix(6)) # Utilize cube_helix for a continuous spectral scale
#palette(hue_pal()(6)) # Alternative: utilize default ggplot2 discrete hues
plot(dm, pch = 20, # pch adjusted for aesthetic rendering
     col_by = 'num_cells', # Map variable to color vector
     legend_main = 'Cell stage')
# 2D Projection
plot(dm, 1:2, pch = 20, col_by = 'num_cells',
     legend_main = 'Cell stage')
# 3D Projection
library(rgl)
plot3d(eigenvectors(dm)[, 1:3],
       col = log2(guo_norm$num_cells),
       type = 's', radius = .01)
view3d(theta = 10, phi = 30, zoom = .8)
# Interactive rotation of the 3D plot is enabled via mouse input in the rgl device
rgl.close()

##################### Dyno Ensemble Trajectory Inference #############################
devtools::install_github("dynverse/dyno")
devtools::install_github("dynverse/dynmethods")
library(dyno)
library(tidyverse)
library(Matrix)
library(Seurat)

# Integrate raw counts and normalized expression matrices
# Note: Seurat matrices require transposition (cells as rows, genes as columns) for Dyno compatibility
dataset <- wrap_expression(
  counts = t(temp_sce[["RNA"]]$counts),
  expression = t(temp_sce[["RNA"]]$data)
)

dataset <- wrap_expression(
  counts = t(subsce[["RNA"]]$counts),
  expression = t(subsce[["RNA"]]$data)
)



# Define prior biological information (anchor starting nodes)
# Note: Starting cell IDs can be iteratively refined based on specific trajectory outputs
dataset <- add_prior_information(
  dataset,
  start_id = "0day_ACTTTCAAGGGAACGG-1"
)  # Fibroblast subpopulation anchor
dataset <- add_prior_information(
  dataset,
  start_id = "1day_TACCTATTCCAAACAC-1"
) # Monocyte/Macrophage subpopulation anchor
dataset <- add_prior_information(
  dataset,
  start_id = "0day_ACGATACAGGGTCGAT-1"
) # Neutrophil/Granulocyte anchor
# Map clustering metadata to the dataset (utilizing 'seurat_clusters' annotations)
dataset <- add_grouping(
  dataset,
  temp_sce$seurat_clusters
)
dataset <- add_grouping(
  dataset,
  temp_sce$RNA_snn_res.1
)
dataset <- add_grouping(
  dataset,
  subsce$seurat_clusters
)
## dynwrap::test_docker_installation(detailed = TRUE) # Verify Docker accessibility for backend method execution
## Evaluate and select optimal trajectory algorithms tailored to dataset characteristics
guidelines <- guidelines_shiny(dataset)
methods_selected <- guidelines$methods_selected
## Implementation of the PAGA tree methodology
model_paga <- infer_trajectory(dataset, methods_selected[1])
model_paga <- model_paga %>% add_dimred(dyndimred::dimred_mds, expression_source = dataset$expression)
plot_dendro(model,
            expression_source = dataset$expression, grouping = dataset$grouping)
set.seed(1)
dataset <- dyntoy::generate_dataset(model = "bifurcating", num_cells = 200)
## Preliminary visualization mapping (highly tunable parameters available)
plot_dimred(
  model_paga,
  expression_source = dataset$expression,
  grouping = dataset$grouping
)
system()
dev.off()
## Define rooting genes to anchor the biological directionality of differentiation
model <- model_paga %>%
  add_root_using_expression(c("Mt1","Mt2","Angptl4","Ccl2"), dataset$expression)

model <- model_paga %>%
  add_root_using_expression(c("Ly6c2","Il1b","Stmn1","Trem2","Col5a2"), dataset$expression)
model <- model %>% 
  add_root_using_expression(c("Ccr2"), dataset$expression)
### Map functional marker genes to defined biological milestones
model <- label_milestones_markers(
  model,
  markers = list(
    Pro_inflammatory = c("Angptl4"),
    Prolifmyo = c("Ccnb2"),
    Myofib = c("Acta2"),
    MFC =c("Sfrp2")

  ),
  dataset$expression
)

model <- label_milestones_markers(
  model,
  markers = list(
    Monocytes =c("Ly6c2"),
    Pro_inflammation = c("Il1b"),
    #Proliferation = c("Stmn1"),
    Reparative = c("Trem2"),
    Collagens =c("Col5a2")

  ),
  dataset$expression
)
model_rooted <- model %>% add_root(root_milestone_id = "6")
## Advanced milestone visualization
#model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = dataset$expression)
plot_dimred(
  model,
  color_cells = "pseudotime",
  #color_density = "grouping",
  #color_cells = "grouping",
  expression_source = dataset$expression,
  grouping = dataset$grouping,
  label_milestones = TRUE
)


### Specific Gene Expression Projections

plot_dimred(
  model,
  expression_source = dataset$expression,
  color_cells = "feature",
  feature_oi = "Ly6c2",
  #color_density = "grouping",
  grouping = dataset$grouping,
  label_milestones = TRUE
)

plot_dimred(
  model,
  expression_source = dataset$expression,
  color_cells = "grouping",
  feature_oi = "Ly6c2",
  color_density = "grouping",
  grouping = dataset$grouping,
  label_milestones = TRUE
)

plot_dimred(
  model,
  expression_source = dataset$expression,
  color_cells = "feature",
  feature_oi = "Col3a1",
  color_density = "grouping",
  grouping = dataset$grouping,
  label_milestones = FALSE
)

# Global overview of highly predictive trajectory genes
plot_heatmap(
  model,
  expression_source = dataset$expression,
  grouping = dataset$grouping,
  features_oi = 50
)


## Isolate key features driving trajectory bifurcations

branching_milestone <- model$milestone_network %>% group_by(from) %>% filter(n() > 1) %>% pull(from) %>% dyndimred::first()

branch_feature_importance <- calculate_branching_point_feature_importance(model, expression_source=dataset$expression, milestones_oi = branching_milestone)

branching_point_features <- branch_feature_importance %>% top_n(20, importance) %>% pull(feature_id)

plot_heatmap(
  model,
  expression_source = dataset$expression,
  features_oi = branching_point_features
)

space <- dyndimred::dimred_mds(dataset$expression)
map(branching_point_features[1:20], function(feature_oi) {
  plot_dimred(model, dimred = space, expression_source = dataset$expression, feature_oi = feature_oi, label_milestones = FALSE) +
    theme(legend.position = "none") +
    ggtitle(feature_oi)
}) %>% patchwork::wrap_plots()




################### BioTip Tipping Point Analysis ################################################
# BioTIP Main Functions Walk Through
# Partition matrices by biological samples (samplesL)
subsce_trans <- as.SingleCellExperiment(temp_sce)
subsce_trans <- as.SingleCellExperiment(subsce)
samplesL <- split(rownames(colData(subsce_trans)),f = colData(subsce_trans)$seurat_clusters)
# Transcript Pre-selection Modeling
dec.pois <- modelGeneVarByPoisson(subsce_trans)
hvg <- getTopHVGs(dec.pois, n=4000)
hvg <- getTopHVGs(dec.pois, n=2000)
hvg <- intersect(hvg, rownames(subsce_trans))
dat <- subsce_trans[hvg,]
logmat <- as.matrix(logcounts(dat))
global <- as.matrix(logcounts(subsce_trans))
## Generate scatter plot for global transcript metrics
gene_mean <- apply(global, 1, mean)
gene_sd <- apply(global, 1, sd)
blue_gene_indices <- which(names(gene_mean) %in% hvg)
red_gene_indices <- which(names(gene_mean) %in% cluster_specificgenes)
red_gene_indices <- which(gene_sd >1)
points(gene_mean[red_gene_indices], gene_sd[red_gene_indices], col = "red", pch = 10)
points(gene_mean[blue_gene_indices], gene_sd[blue_gene_indices], col = "lightblue", pch = 10)
# Append figure legend
legend("topright", legend = c("Red Genes"), col = "pink", pch = 16)
plot(gene_mean, gene_sd, xlab = "Mean Expression", ylab = "Standard Deviation")

ggplot(data, aes(x = gene_mean, y = gene_sd)) +
  geom_point(size = data$size)

# Generate smoothed density scatter plot (smoothScatter)
smoothScatter(x = gene_mean, y = gene_sd, nbin = 128,
              colramp = colorRampPalette(c("white", blues9)),
              nrpoints = 100, ret.selection = FALSE,
              transformation = function(x) x^0.25,
              postPlotHook = box,
              xlab = "Mean Expression", ylab = "Standard Deviation",
              xaxs = par("xaxs"), yaxs = par("yaxs"))


################### BioTIP Tipping Point Analysis ################################################
# BioTIP Main Functions Walk Through
# samplesL
subsce_trans <- as.SingleCellExperiment(temp_sce)
subsce_trans <- as.SingleCellExperiment(subsce)
samplesL <- split(rownames(colData(subsce_trans)),f = colData(subsce_trans)$seurat_clusters)

# Transcript Pre-selection
# Select highly variable genes (HVGs)
# Incremental threshold adjustments to balance computational load and HVG retention
cut.preselect = 0.10
cut.preselect = 0.12
cut.preselect = 0.15
cut.preselect = 0.20
# cut.preselect = 0.01 Increased threshold due to computationally prohibitive runtimes
set.seed(2020)
pb <- txtProgressBar(min = 0, max = 100, style = 3)
testres <- optimize.sd_selection(logmat[,unlist(samplesL)], samplesL, B=100, cutoff=cut.preselect, times=.75, percent=0.8)
#save(testres, file=paste0("Fibroblasts_biotip.RData"), compress=TRUE)
#save(testres, file=paste0("macrophages_biotip.RData"), compress=TRUE)
#save(testres, file=paste0("Granulocytes_biotip.RData"), compress=TRUE)
#save(testres, file=paste0("Granulocytes_biotip2.0.RData"), compress=TRUE)
#save(testres, file=paste0("macrophages_biotip2.0.RData"), compress=TRUE)
testres <- load
testres_list<- list()
#### Iterative testing block
nested_list<- list()
#### Iterative testing block
for (i in 1:15) {
  inner_list <- list()
  for (j in 1:1) {
    inner_list[[j]] <- optimize.sd_selection(logmat[,unlist(samplesL)], samplesL, B=100, cutoff=cut.preselect, times=.75, percent=0.8)
  }
  nested_list[[i]] <- inner_list
}
testres <- optimize.sd_selection(logmat[,unlist(samplesL)], samplesL, B=100, cutoff=cut.preselect, times=.75, percent=1)

class(testres$`8`)[,1]
cluster_specificgenes <- rownames(as.data.frame(testres$`0`))
cluster_specificgenes <- append(cluster_specificgenes,rownames(as.data.frame(testres$`14`)))
cluster_specificgenes <- unique(cluster_specificgenes)

# Network Partition
cut.fdr = 0.05
cut.fdr = 0.10
cut.fdr = 0.12
cut.fdr = 0.15
cut.fdr = 0.20
igraphL <- getNetwork(testres, fdr = cut.fdr)
cluster <- getCluster_methods(igraphL)

## Network partition using random walk
dev.off()
par(mfrow=c(3,4))
cluster = getCluster_methods(igraphL)
#i =5
for (i in 1:length(igraphL)){
  tmp = igraphL[[i]]
  E(tmp)$width <- E(tmp)$weight*3
  V(tmp)$community= cluster[[i]]$membership
  mark.groups = table(cluster[[i]]$membership)
  colrs = rainbow(length(mark.groups), alpha = 0.3)
  V(tmp)$label <- NA
  plot(tmp,vertex.color=colrs[V(tmp)$community],vertex.size = 5,
       mark.groups=cluster[[i]])
}
table(V(tmp)$community)
which(V(tmp)$community  == 114)
plot(tmp,vertex.label = V(tmp)$name,vertex.color=colrs[V(tmp)$community],vertex.size =5,mark.groups=cluster[[i]])
nodes_in_community <- V(tmp)$name[V(tmp)$community ==8]
# Extract induced subgraph containing exclusively nodes from the target community
subgraph <- induced_subgraph(tmp, nodes_in_community)
plot(subgraph,vertex.label = V(subgraph)$name,vertex.color=colrs[V(subgraph)$community],vertex.size = 6)
# Calculate node degree centrality
node_degrees <- degree(subgraph)
# Determine minimum and maximum degrees for normalization
min_degree <- min(node_degrees)
max_degree <- max(node_degrees)
# Map node degrees to a continuous cool-to-warm color palette
color_palette <- colorRampPalette(c("#8FD2E6", "white","#ED884C" ))
color_palette <- colorRampPalette(c( "white","#ED884C" ))
node_colors <- color_palette(100)[as.integer(100 * (node_degrees - min_degree) / (max_degree - min_degree)) + 1]

# Map gene expression levels to the continuous color palette
local_sce <- subset(temp_sce ,subset = seurat_clusters == "6")
#local_sce <- subset(plasma_b_cell,subset = anno_clusters == "5")
local_sce <- as.SingleCellExperiment(local_sce)
local_sce <- as.matrix(logcounts(local_sce))
gene_sd <- apply(local_sce, 1, sd)
gene_sd <- gene_sd[V(subgraph)$name]
gene_exp <- apply(local_sce, 1, mean)
gene_exp <- gene_exp[V(subgraph)$name]
# Determine minimum and maximum mean/SD values for dynamic scaling
min_exp <- min(gene_exp)
max_exp <- max(gene_exp)
node_size <- gene_sd*10
node_colors <- color_palette(100)[as.integer(100 * (gene_exp - min_exp) / (max_exp - min_exp)) ]



# Apply Fruchterman-Reingold (FR) force-directed layout
layout <- layout_with_fr(subgraph)
# Render network graph, scaling node color and size by degree centrality and standard deviation
plot(subgraph,
     vertex.label = V(subgraph)$name,
     vertex.color = node_colors,
     vertex.size = node_size*3,
     #vertex.label.dist = 1.5,
     vertex.label.cex = 0.8,
     layout = layout,
     vertex.label.color = "black")
# Export edge list topologies
edges <- as_edgelist(subgraph)
# Fibroblast Lineage Networks
write.table(edges, "7edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "7gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
write.table(edges, "4edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "4gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
# Fibroblast Lineage - Updated Clusters 4 and 6
write.table(edges, "fibroblast_4edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "fibroblast_4gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
write.table(edges, "fibroblast_6edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "fibroblast_6gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)

# Granulocyte/Neutrophil Networks
write.table(edges, "2edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "2gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)

write.table(edges, "Granulocyte_5edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "Granulocyte_5gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)

# Monocyte/Macrophage Networks
write.table(edges, "macrophages_4_2.0edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "macrophages_4gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
# Monocyte/Macrophage Networks
write.table(edges, "macrophages_8edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "macrophages_8gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)

write.table(edges, "macrophages_5edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "macrophages_5gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)






# Map Pearson Correlation Coefficients (PCC) to edge line widths

# Putative Critical Transition Signals (CTSs) Identification by MCI score
membersL <- getMCI(cluster, testres, adjust.size = FALSE, fun='BioTIP')



# Visualize Dynamic Network Biomarker (DNB) scores for modules across all cell clusters
dev.off()
cut.minsize = 60
pdf('./临界态图/Figure2-6A.pdf',height = 12,width = 22,onefile = F)
par(oma = c(0, 0, 0, 0))  # Configure outer margins
par(mar = c(1, 1, 1, 1))  # Configure inner margins
plotBar_MCI(membersL, ylim=c(0,20), minsize = cut.minsize)
dev.off()
# Isolate the network module exhibiting the maximum Module Criticality Index (MCI)
subg <- induced_subgraph(subgraph, vids = V(subgraph)[name %in% V(subgraph)])
edgelist_sub <- as_edgelist(subgraph, names = TRUE)
edge_attrs_sub <- edge_attr(subgraph)
result_internal <- data.frame(edgelist_sub, edge_attrs_sub, check.names = FALSE)
head(result_internal) # Inspect the top rows of the extracted network attributes
write.csv(result_internal,"./大论文图和数据/Macro_DNB_early_network.csv")
write.csv(result_internal,"./大论文图和数据/Macro_DNB_late_network.csv")
write.csv(result_internal,"./大论文图和数据/Fibro_DNB_early_network.csv")
membersL$sd$`5`
topMCI = getTopMCI(membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]], min=cut.minsize, n=10)
maxMCIms <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], min=cut.minsize, n=10)
maxMCI = getMaxStats(membersL[['MCI']], maxMCIms[['idx']])
CTS.Lib = getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])
CTS.Lib.Symbol <- CTS.Lib
maxMCI <- sort(maxMCI, decreasing = TRUE)
maxMCI <- head(maxMCI, 10)
# Export identified Dynamic Network Biomarkers (DNBs)
# Determine maximum vector length for data frame padding
max_length <- max(lengths(CTS.Lib))
df <- data.frame(matrix(NA, nrow = max_length, ncol = length(CTS.Lib)))
for (i in 1:length(CTS.Lib)) {
  col_name <- names(CTS.Lib)[i]
  col_data <- CTS.Lib[[i]]
  df[, i] <- c(col_data, rep(NA, max_length - length(col_data)))
}
colnames(df) <- names(CTS.Lib)
write.csv(df, file = "fibrobalst_DNB.csv", row.names = TRUE)

write.csv(df, file = "Granulocytes_DNB.csv", row.names = TRUE)
write.csv(df, file = "Granulocytes_DNB2.0.csv", row.names = TRUE)

write.csv(df, file = "macropahges_DNB.csv", row.names = TRUE)
write.csv(df, file = "macropahges_DNB2.0.csv", row.names = TRUE)

names(CTS.Lib.Symbol)
# Calculate shrinkage estimation of covariance matrix (IC.shrink.score)
M <- cor.shrink(logmat[,unlist(samplesL)], Y = NULL, MARGIN = 1, shrink = TRUE)
# Statistically validate the significance of DNB scores via permutation testing
C = 200
simuMCI = list()
set.seed(2020)
for (i in 1:length(CTS.Lib)){
  n <- length(CTS.Lib[[i]])
  simuMCI[[i]] <- simulationMCI(n, samplesL, logmat,  B=C, fun="BioTIP", M=M)
}
dev.off()
par(mfrow=c(5,1))
for (i in 1:length(CTS.Lib)){
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las=2,
                      main=paste0("Cluster ", names(maxMCI)[i], "; ",
                                  length(CTS.Lib[[i]]), " genes", "\n","vs. ",
                                  "100 times of gene-permutation"),
                      which2point=names(maxMCI)[i])
}

dev.off()
pdf('./临界态图/Figure2-6B.pdf',height = 10,width = 8,onefile = F)
par(mfrow=c(5,1))
for (i in 1:5){
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las=2,
                      main=paste0("Cluster ", names(maxMCI)[i], "; ",
                                  length(CTS.Lib[[i]]), " genes", "\n","vs. ",
                                  "100 times of gene-permutation"),
                      which2point=names(maxMCI)[i])
}
dev.off()



# Confirm candidate tipping points using Information Criterion (IC) and Delta-IC scores
C= 200
BioTIP_scores <- SimResults_g <- list()
set.seed(101010)
for(i in 1:length(CTS.Lib)){
  CTS <- CTS.Lib.Symbol[[i]]
  n <- length(CTS)
  BioTIP_scores[[i]] <- getIc(logmat[,unlist(samplesL)], samplesL, CTS, fun="BioTIP", shrink=TRUE, PCC_sample.target = 'none' )
  SimResults_g[[i]]  <- simulation_Ic(n, samplesL, logmat, B=C, fun="BioTIP", shrink=TRUE, PCC_sample.target = 'none')
}
names(BioTIP_scores) <- names(SimResults_g) <- names(CTS.Lib)


dev.off()
ylim = 1
par(mfrow=c(2,2))
for(i in 1:length(BioTIP_scores)){
  n = length(CTS.Lib[[i]])
  interesting = which(names(samplesL) == names(BioTIP_scores[i]))
  p = length(which(SimResults_g[[i]][interesting,] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
  p = p/ncol(SimResults_g[[i]])
  p2 = length(which(SimResults_g[[i]] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
  p2 = p2/ncol(SimResults_g[[i]])
  p2 = p2/nrow(SimResults_g[[i]])
}
ylim = 1
par(mfrow=c(2,2))
for(i in 1:length(BioTIP_scores)){
  n <- length(CTS.Lib[[i]])
  interesting = which(names(samplesL) == names(BioTIP_scores[i]))
  plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], las = 2, ylab="Ic.shrink", ylim=c(0,ylim),
                     main=paste("Cluster ",names(CTS.Lib)[i],"_",n, "genes", "\n","vs. ",
                                "500 gene-permutations"),
                     fun="matplot",  #fun="boxplot",
                     which2point= interesting)
  plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]],
                     main = paste("Delta Ic*",n,"genes"), ylab=NULL,
                     xlim=range(c(BioTIP_scores[[i]][names(BioTIP_scores)[i]],
                                  SimResults_g[[i]])))
}


dev.off()
pdf('./临界态图/Figure2-6C.pdf',height = 10,width = 8,onefile = F)
ylim = 0.4
par(mfrow=c(5,2))
for(i in 1:5){
  n = length(CTS.Lib[[i]])
  interesting = which(names(samplesL) == names(BioTIP_scores[i]))
  p = length(which(SimResults_g[[i]][interesting,] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
  p = p/ncol(SimResults_g[[i]])
  p2 = length(which(SimResults_g[[i]] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
  p2 = p2/ncol(SimResults_g[[i]])
  p2 = p2/nrow(SimResults_g[[i]])
}
ylim = 0.4
par(mfrow=c(5,2))
for(i in 1:5){
  n <- length(CTS.Lib[[i]])
  interesting = which(names(samplesL) == names(BioTIP_scores[i]))
  plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], las = 2, ylab="Ic.shrink", ylim=c(0,ylim),
                     main=paste("Cluster ",names(CTS.Lib)[i],"_",n, "genes", "\n","vs. ",
                                "500 gene-permutations"),
                     fun="matplot",  #fun="boxplot",
                     which2point= interesting)
  plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]],
                     main = paste("Delta Ic*",n,"genes"), ylab=NULL,
                     xlim=range(c(BioTIP_scores[[i]][names(BioTIP_scores)[i]],
                                  SimResults_g[[i]])))
}

############################### DNB Expression Trend and Clustering Analysis ################################
# Define custom function for DNB mean expression matrix extraction
getwd()
i=2
h=1
temp_sce2 <- temp_sce
Idents(temp_sce2) <- "seurat_clusters"
DNBs_expression <- function(expression_matrix,target_cluster,DNB_names){
  failed_genes <- c()
  # Remove duplicated transcription factor names
  DNB_names <- unique(na.omit(DNB_names))
  #tryCatch({
  average_expression <- c()
  df <- data.frame(matrix(nrow = length(target_cluster), ncol = length(DNB_names)))
  rownames(df) <- target_cluster
  colnames(df) <- DNB_names
  # Extract cell indices for the target clusters
  for (i in seq_along(target_cluster)){
    cell_indices <- which(temp_sce2$seurat_clusters == target_cluster[i])
    for (h in 1:length(DNB_names)){
      # Verify if the target gene exists in the scRNA-seq expression matrix
      if (DNB_names[h] %in% rownames(expression_matrix)) {
        average_expression <- mean(expression_matrix[DNB_names[h], cell_indices])

      }else {
        failed_genes <<- c(failed_genes, DNB_names[h])
        average_expression <- NA # Computation failed, return NA
      }
      df[i, DNB_names[h]] <- average_expression
  }
  if (length(failed_genes) > 0) {

    print(paste("List of DNBs failing computation:", paste(unique(failed_genes), collapse = ", ")))
    cat("\n")}
  }
  return(df)
  #})
}
DNB_list <- df
DNB_list <- as.list(read.csv("./fibrobalst_DNB.csv"))

DNB_list <- as.list(read.csv("../E-MTAB-7895/Granulocytes_DNB.csv"))
DNB_list <- as.list(read.csv("../E-MTAB-7895/Granulocytes_DNB2.0.csv"))

DNB_list <- as.list(read.csv("./macropahges_DNB.csv"))
DNB_list <- as.list(read.csv("../E-MTAB-7895/macropahges_DNB2.0.csv"))


Idents(temp_sce2) # Verify scRNA-seq object to ensure 'label' is set as the active identity
#expression_matrix = subsce[["RNA"]]$scale.data #normalized_scale
expression_matrix = temp_sce2[["RNA"]]$scale.data #normalized_scale
expression_matrix = temp_sce[["RNA"]]$data #normalized_scale
#expression_matrix = t(scale(t(subsce[["RNA"]]$counts))) #scale
expression_matrix = t(scale(t(temp_sce2[["RNA"]]$counts))) #scale
expression_matrix = t(scale(t(temp_sce2[["RNA"]]$data))) #scale
expression_matrix = filter_macrophage_cells[["RNA"]]$data #normalized_scale
expression_matrix = plasma_b_cell[["RNA"]]$data #normalized_scale
expression_matrix = filter_sce_mesen[["RNA"]]$data #normalized_scale

### Construct Dynamic Network Biomarker (DNB) Protein-Protein Interaction (PPI) Networks
### Step 0: Load required libraries
library(igraph)
library(Seurat)
library(tidyverse)
library(ggraph)
library(viridis)
library(tidyr)

# Map gene symbols to STRING protein IDs
protein_info <- read.delim("10090.protein.info.v12.0.txt", header = TRUE, sep = "\t", comment.char = "#")
colnames(protein_info)[1:2] <- c("protein_id", "gene_symbol")
symbol_to_protein <- setNames(protein_info$protein_id, protein_info$gene_symbol)
id_to_symbol <- setNames(protein_info$gene_symbol, protein_info$protein_id)

# Load empirical PPI interaction links
ppi_links <- read.delim("10090.protein.links.v12.0.txt", header = TRUE, sep = " ")
ppi_links <- ppi_links %>% dplyr::filter(combined_score >= 700)

# Define temporal trajectory sequence
trajectory <- c("12", "5", "8", "3", "4", "1", "10")
trajectory <- c("6", "1", "3")
trajectory <- c("0", "3", "2","1")
trajectory <- c("11", "5", "6","9")
trajectory <- c("5", "6", "12")
# Extract active DNB gene modules
gene_list <- na.omit(DNB_list[["5"]])
gene_list <- na.omit(DNB_list$X5)
gene_list <- na.omit(df[["1"]])
gene_list <- na.omit(df[["3"]])
gene_list <- na.omit(df[["5"]])
gene_list <- na.omit(df[["6"]])
# Load supplementary graphic libraries
library(Seurat)
library(tidyverse)
library(pheatmap)
library(igraph)
library(ggraph)

# Extract first-order PPI neighbor gene symbols
get_ppi_neighbors <- function(gene, 
                              symbol_to_protein, 
                              id_to_symbol, 
                              ppi_links, 
                              filter_combined_score, 
                              allowed_genes = NULL) {
  protein_id <- symbol_to_protein[gene]
  if (is.na(protein_id)) return(character(0))

  neighbors <- ppi_links %>%
    dplyr::filter(combined_score >= filter_combined_score) %>%
    dplyr::filter(protein1 == protein_id | protein2 == protein_id) %>%
    dplyr::mutate(partner = ifelse(protein1 == protein_id, protein2, protein1)) %>%
    dplyr::pull(partner)

  neighbor_symbols <- id_to_symbol[neighbors]
  neighbor_symbols <- neighbor_symbols[!is.na(neighbor_symbols)]
  
  # Constrain neighbors to an externally provided allowed_genes list (e.g., DEGs)
  if (!is.null(allowed_genes)) {
    neighbor_symbols <- intersect(neighbor_symbols, allowed_genes)
  }
  return(unique(neighbor_symbols))
}

# Identify top N co-expressed genes within specific clusters
get_top_correlated <- function(
    gene,
    seurat_obj,
    cluster_id = "5",
    top_n = 20,
    pval_cutoff = 0.05,
    allowed_genes = NULL,
    balance_sign = T,          # Boolean: Enforce equal distribution of positive/negative correlations
    method = c("pearson", "spearman"),
    adjust_p = FALSE,              # Boolean: Apply multiple testing correction
    p_adjust_method = "BH"         # p.adjust methodology
) {
  method <- match.arg(method)
  
  # 1) Extract cell indices and expression matrix for the target cluster
  cells <- colnames(seurat_obj)[seurat_obj$seurat_clusters == cluster_id]
  if (length(cells) == 0L) {
    warning("No cells found for cluster_id = ", cluster_id)
    return(list(genes = gene, correlations = NA, signs = "Center"))
  }
  data_mat <- as.matrix(seurat_obj[["RNA"]]$data[, cells, drop = FALSE])
  
  # 2) Validate target gene presence
  if (!gene %in% rownames(data_mat)) {
    warning("gene not found in RNA assay: ", gene)
    return(list(genes = gene, correlations = NA, signs = "Center"))
  }
  
  # 3) Compute expression correlations and statistical significance
  gene_expr <- data_mat[gene, ]
  # Iteratively calculate correlation metrics across the matrix
  cors  <- apply(data_mat, 1, function(x) suppressWarnings(cor(x, gene_expr, method = method)))
  pvals <- apply(data_mat, 1, function(x) suppressWarnings(cor.test(x, gene_expr, method = method)$p.value))
  
  # 4) Apply multiple testing correction (optional)
  p_use <- if (adjust_p) p.adjust(pvals, method = p_adjust_method) else pvals
  
  # 5) Preliminary filtering: Significance threshold + exclude self + intersect with allowed_genes
  valid_genes <- names(p_use)[p_use < pval_cutoff]
  valid_genes <- setdiff(valid_genes, gene)
  if (!is.null(allowed_genes)) valid_genes <- intersect(valid_genes, allowed_genes)
  
  # Exclude genes yielding NA correlations
  valid_genes <- valid_genes[!is.na(cors[valid_genes])]
  if (length(valid_genes) == 0L) {
    warning("No valid genes passed the filters.")
    return(list(genes = gene, correlations = NA, signs = "Center"))
  }
  
  # 6) Subset top_n candidates
  take_n <- min(top_n, length(valid_genes))
  
  if (!balance_sign) {
    # Rank uniformly by absolute correlation |cor| in descending order
    sorted <- valid_genes[order(abs(cors[valid_genes]), decreasing = TRUE)]
    top_genes <- head(sorted, take_n)
  } else {
    # Balance correlation signs: ceil(n/2) for positive, floor(n/2) for negative correlations. Pad via |cor| if insufficient.
    k_pos <- take_n %/% 2 + (take_n %% 2)  # ceil
    k_neg <- take_n - k_pos                # floor
    
    pos <- valid_genes[cors[valid_genes] > 0]
    neg <- valid_genes[cors[valid_genes] < 0]
    zer <- valid_genes[cors[valid_genes] == 0]
    
    # Sort positive values descending; sort negative values by magnitude |cor| descending
    pos_sorted <- pos[order(cors[pos], decreasing = TRUE)]
    neg_sorted <- neg[order(abs(cors[neg]), decreasing = TRUE)]
    
    pick_pos <- head(pos_sorted, k_pos)
    pick_neg <- head(neg_sorted, k_neg)
    picked   <- c(pick_pos, pick_neg)
    top_genes <- picked
  }
  # 6) Assemble output structure (Central gene / Correlation values / Directionality signs)
  correlation_signs <- ifelse(cors[top_genes] > 0, "Positive",
                              ifelse(cors[top_genes] < 0, "Negative", "Zero"))
  signs <- c("Center", correlation_signs)
  
  return(list(
    genes = c(gene, top_genes),
    correlations = c(NA, cors[top_genes]),
    signs = signs
  ))

  }


# Generate trajectory-ordered expression matrix
get_cluster_expr_matrix <- function(genes, seurat_obj, trajectory) {
  expr_mat <- matrix(NA, nrow = length(genes), ncol = length(trajectory))
  rownames(expr_mat) <- genes
  colnames(expr_mat) <- as.character(trajectory)

  for (i in seq_along(trajectory)) {
    cells <- colnames(seurat_obj)[seurat_obj$seurat_clusters == trajectory[i]]
    avg_expr <- rowMeans(seurat_obj[["RNA"]]$data[genes, cells, drop = FALSE])
    expr_mat[, i] <- avg_expr
  }
  return(expr_mat)
}
#rowMeans(seurat_obj[["RNA"]]$data["Zdhhc15", colnames(seurat_obj)[seurat_obj$seurat_clusters == trajectory[2]], drop = FALSE])
#seurat_obj[["RNA"]]$data["Serpine1", colnames(seurat_obj)[seurat_obj$seurat_clusters == trajectory[2]], drop = FALSE]
# Plot expression heatmap
center_gene <- gene
plot_expression_heatmap <- plot_expression_heatmap <- function(expr_mat, center_gene, correlations, signs, filename)  {
  # Create correlation sign annotations
  correlation_annotations <- data.frame(Sign = signs)
  rownames(correlation_annotations) <- rownames(expr_mat)
  
  # Configure aesthetic palette for annotations
  sign_colors <- list(Sign = c("Positive" = "red", "Negative" = "blue", "Center" = "black"))
  
  # Render pheatmap with row annotations representing correlation dynamics
  pheatmap::pheatmap(expr_mat,
                     scale = "row",
                     cluster_cols = FALSE,
                     annotation_row = correlation_annotations,  # Append row annotations
                     annotation_colors = sign_colors,  # Map specific colors to center genes
                     filename = filename,
                     width = 8, height = 12)

}
dev.off()
# Plot expression trajectory line graph
plot_expression_lineplot <- function(expr_mat, center_gene, filename) {
  expr_long <- expr_mat %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cluster", values_to = "expression")

  # Enforce structural cluster ordering to match the trajectory sequence
  cluster_order <- colnames(expr_mat)
  expr_long$cluster <- factor(expr_long$cluster, levels = cluster_order)
  
  p <- ggplot(expr_long, aes(x = cluster, y = expression, color = gene, group = gene)) +
    geom_line() + geom_point() +
    theme_minimal() +
    labs(title = paste("Expression Trend:", center_gene))

  ggsave(filename, plot = p, width = 7, height = 5)
}
# Plot permutation test null distribution density
plot_permutation_density_to_file <- function(permuted_means,
                                             observed_mean,
                                             p_val,
                                             gene,
                                             direction = c("Positive", "Negative"),
                                             outdir = "cor_analysis") {
  direction <- match.arg(direction)
  # Establish output directory infrastructure
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Formulate output file string
  outfile <- file.path(outdir, paste0(gene, "_", direction, "_permutation.png"))
  
  png(outfile, width = 800, height = 600)
  plot(density(permuted_means),
       main = paste0("Permutation Null Distribution (", direction, ")"),
       xlab = "Mean log2FC of Random Modules",
       col = "darkgreen", lwd = 2)
  
  abline(v = observed_mean, col = "red", lty = 2, lwd = 2)
  
  text(x = observed_mean, 
       y = max(density(permuted_means)$y) * 0.9,
       labels = paste0("Observed = ", round(observed_mean, 3), 
                       "\nP = ", signif(p_val, 3)),
       pos = ifelse(observed_mean >= mean(permuted_means), 4, 2),
       col = "red", cex = 0.9)
  dev.off()
}


# Plot temporal dynamic PPI network topologies
plot_temporal_PPI_network_pdf <- function(center_gene, node_genes, expr_mat, network_edges, trajectory, outpath) {
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(dplyr)

  # Infer linear expression trend topology
  get_trend <- function(expr_vector) {
    x <- seq_along(expr_vector)
    model <- lm(expr_vector ~ x)
    slope <- coef(model)[2]
    if (slope > 0.05) return("Up")
    if (slope < -0.05) return("Down")
    return("Flat")
  }

  # Calculate coordinate layouts enforcing center node centrality
  trends <- sapply(rownames(expr_mat), function(g) get_trend(expr_mat[g, ]))
  node_genes_wo_center <- setdiff(node_genes, center_gene)
  trend_subset <- trends[node_genes_wo_center]

  assign_circular_position_by_trend <- function(genes, trends) {
    df <- data.frame(name = genes, trend = trends, stringsAsFactors = FALSE)
    df$trend <- factor(df$trend, levels = c("Up", "Flat", "Down"))
    df <- df[order(df$trend), ]
    n <- nrow(df)
    angle_step <- 360 / n
    df$angle <- seq(0, 360 - angle_step, by = angle_step)
    theta_rad <- pi * df$angle / 180
    df$x <- cos(theta_rad)
    df$y <- sin(theta_rad)
    return(df)
  }

  layout_df <- assign_circular_position_by_trend(node_genes_wo_center, trend_subset)
  layout_df <- rbind(layout_df, data.frame(name = center_gene, trend = "Center", x = 0, y = 0, angle = NA))

  plots <- list()

  for (i in seq_along(trajectory)) {
    clus <- as.character(trajectory[i])
    expr_vector <- expr_mat[, clus]
    names(expr_vector) <- rownames(expr_mat)

    g <- igraph::graph_from_data_frame(network_edges, vertices = layout_df[, c("name")], directed = FALSE)
    V(g)$expression <- expr_vector[V(g)$name]
    V(g)$trend <- trends[V(g)$name]

    p <- ggraph(g, layout = "manual", x = layout_df$x, y = layout_df$y) +
      geom_edge_link(edge_colour = "grey80", width = 0.8, alpha = 0.4) +
      geom_node_point(aes(size = expression, fill = expression, shape = trend)) +
      geom_node_text(aes(label = name), repel = TRUE, size = 3) +
      scale_fill_gradientn(
        colours = c("navy", "deepskyblue", "yellow", "red"),
        limits = c(0, 5),
        oob = scales::squish
      ) +
      scale_shape_manual(values = c("Up" = 24, "Down" = 25, "Flat" = 21, "Center" = 22)) +
      theme_void() +
      ggtitle(paste("Cluster", clus)) +
      guides(
        fill = guide_colorbar(title = "Expression"),
        shape = guide_legend(title = "Trend")
      )

    plots[[i]] <- p
  }

  ggsave(
    filename = file.path(outpath, paste0(center_gene, "_temporal_PPI_network.pdf")),
    plot = wrap_plots(plots, ncol = length(trajectory)) +
      plot_annotation(title = paste("Temporal PPI Network for", center_gene)),
    width = 5 * length(trajectory), height = 6, device = "pdf"
  )
}

# Plot temporal co-expression correlation network
plot_temporal_correlation_network_pdf <- function(gene, all_genes, expr_mat, edge_list, trajectory, outpath) {
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(viridis)
  library(purrr)
  library(dplyr)
  #global_max_cor <- max(unlist(lapply(edge_list, function(df) abs(df$weight))), na.rm = TRUE)
  get_trend <- function(expr_vector) {
    x <- seq_along(expr_vector)
    model <- lm(expr_vector ~ x)
    slope <- coef(model)[2]
    if (slope > 0.05) return("Up")
    if (slope < -0.05) return("Down")
    return("Flat")
  }

  assign_circular_position_by_trend <- function(genes, trends) {
    if (length(genes) != length(trends)) stop("Length of genes and trends must be equal.")
    df <- data.frame(name = genes, trend = trends, stringsAsFactors = FALSE)
    df$trend <- factor(df$trend, levels = c("Up", "Flat", "Down"))
    df <- df[order(df$trend), ]
    n <- nrow(df)
    angle_step <- 360 / n
    df$angle <- seq(0, 360 - angle_step, by = angle_step)
    theta_rad <- pi * df$angle / 180
    df$x <- cos(theta_rad)
    df$y <- sin(theta_rad)
    return(df)
  }

  trends <- sapply(rownames(expr_mat), function(g) get_trend(expr_mat[g, ]))
  node_genes <- setdiff(all_genes, gene)
  trend_subset <- trends[node_genes]
  layout_df <- assign_circular_position_by_trend(node_genes, trend_subset)
  layout_df <- rbind(layout_df, data.frame(name = gene, trend = "Center", x = 0, y = 0, angle = NA))
  i = 2
  plots <- list()
  #pdf(file = file.path(outpath, paste0(gene, "_temporal_correlation_network.pdf")), width = 6, height = 5)
  for (i in seq_along(trajectory)) {
    cluster_id <- trajectory[i]
    cluster_expr <- expr_mat[, i]
    edge_df <- edge_list[[i]]

    g <- igraph::graph_from_data_frame(edge_df, vertices = layout_df[, c("name")], directed = FALSE)
    V(g)$expression <- cluster_expr[V(g)$name]
    V(g)$trend <- trends[V(g)$name]
    E(g)$sign <- ifelse(E(g)$weight > 0, "Positive", "Negative")
    E(g)$correlation <- E(g)$weight
    #E(g)$strength <- scales::rescale(abs(E(g)$correlation), to = c(0.5, 3))
    E(g)$strength <- abs(E(g)$correlation)
    #E(g)$strength <- abs(E(g)$correlation) / global_max_cor * 3  # Set maximum line width to 3

    p <- ggraph(g, layout = "manual", x = layout_df$x, y = layout_df$y) +
      geom_edge_link(aes(edge_color = sign, edge_width = strength), alpha = 0.4) +
      geom_node_point(aes(size = expression, fill = expression, shape = trend), color = "grey70", stroke = 0.1) +
      geom_node_text(aes(label = name), repel = TRUE, size = 3) +
      scale_edge_color_manual(values = c("Positive" = "#FF6B6B", "Negative" = "#1E90FF")) +
      scale_edge_width(range = c(0.5, 3)) +
      scale_fill_gradientn(
        colours = c("navy", "deepskyblue", "yellow", "red"),
        limits = c(0, 5),     # Absolute expression range boundary
        oob = scales::squish
      ) +
      scale_shape_manual(values = c("Up" = 24, "Down" = 25, "Flat" = 21, "Center" = 22)) +
      theme_void() +
      ggtitle(paste("Cluster", cluster_id)) +
      guides(
        edge_width = guide_legend(title = "Correlation Strength"),
        edge_color = guide_legend(title = "Correlation Sign"),
        fill = guide_colorbar(title = "Expression"),
        shape = guide_legend(title = "Trend")
      )

    plots[[i]] <- p
  }

  # Compile figures and export to PDF
  combined_plot <- wrap_plots(plots, ncol = length(trajectory)) +
    plot_annotation(title = paste("Temporal Correlation Network for", gene))

  ggsave(
    filename = file.path(outpath, paste0(gene, "_temporal_correlation_network.pdf")),
    plot = combined_plot, width = 5 * length(trajectory), height = 6, device = "pdf"
  )
}

# Export integrated PPI and co-expressed gene matrices
save_ppi_coexpressed_genes <- function(gene, neighbors) {
  ppi_coexp_file <- paste0("ppi_analysis/", gene, "_ppi_coexpressed_neighbors.txt")

  write.table(neighbors, ppi_coexp_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
}
# Initialize environments to store downstream PPI effectors and co-expressed targets per DNB
ppi_list <- list()
coexp_list <- list()
# Main wrapper function: Integrated PPI dynamics
analyze_gene_by_ppi <- function(gene, 
                                seurat_obj, 
                                symbol_to_protein, 
                                id_to_symbol, 
                                ppi_links, 
                                trajectory, 
                                filter_combined_score = 700,  
                                use_deg = TRUE, 
                                deg_genes = NULL) {
  
  # Condition downstream targets to statistically significant DEGs if use_deg = TRUE
  if (use_deg && !is.null(deg_genes)) {
    allowed_genes <- deg_genes  # Enforce externally provided DEG constraint
  } else {
    allowed_genes <- NULL  # Null constraint bypass
  }
  
  neighbors <- get_ppi_neighbors(gene, 
                                 symbol_to_protein, 
                                 id_to_symbol, 
                                 ppi_links, 
                                 filter_combined_score, 
                                 allowed_genes)
  
  all_genes <- c(gene, neighbors)
  all_genes <- intersect(all_genes, rownames(seurat_obj))
  if (length(all_genes) < 2) return(NULL)

  expr_mat <- get_cluster_expr_matrix(all_genes, seurat_obj, trajectory)
  # Filter invariant genes or zero-expression vectors
  expr_mat <- expr_mat[apply(expr_mat, 1, function(x) sd(x, na.rm = TRUE) > 0), , drop = FALSE]
  all_genes <- rownames(expr_mat)
  if (length(all_genes) < 2) return(NULL)

  # Update valid neighbors (excluding central DNB target)
  neighbors <- setdiff(all_genes, gene)
  # Quantify downstream PPI functional effectors
  ppi_genes_count <- length(neighbors)
  print(paste("Number of PPI genes for", gene, ":", ppi_genes_count,":",neighbors))
  print(ppi_list)  # Diagnostic tracking of list structures
  # Store extracted downstream PPI targets structurally
  ppi_list[[gene]] <<- neighbors  # Cache validated downstream PPI genes
  # Export integrated object sets
  #save_ppi_coexpressed_genes(gene, neighbors)
  dir.create("ppi_analysis", showWarnings = FALSE)
  plot_expression_heatmap(expr_mat, gene, paste0("ppi_analysis/", gene, "_ppi_heatmap.png"))
  plot_expression_lineplot(expr_mat, gene, paste0("ppi_analysis/", gene, "_ppi_lineplot.png"))

  # # Conceptual static PPI mapping via absolute average expression
  # edge_df <- data.frame(from = gene, to = neighbors)
  # expr_vector <- rowMeans(seurat_obj[["RNA"]]$data[all_genes, , drop = FALSE])
  # #plot_expression_network(gene, all_genes, expr_vector, edge_df, paste0("ppi_analysis/", gene, "_ppi_network.png"))
  #
  # plot_temporal_PPI_network_pdf(
  #   center_gene = gene,
  #   node_genes = all_genes,
  #   expr_mat = expr_mat,
  #   network_edges = edge_df,
  #   trajectory = trajectory,
  #   outpath = paste0("ppi_analysis"))

  }

# Main wrapper function: Integrated co-expression correlation inference
dnb_gene <- gene
plot_quadrant_summary <- function(dnb_gene, expr_mat, signs, direction = "Positive", 
                                  outpath = "quadrant_summary_plots", csv_path = "quadrant_summary_data.csv") {
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
  
  M1_clusters <- c( "5", "8")
  M2_clusters <- c("3", "4", "1")
  
  # Calculate absolute log2FC for DNB targets
  expr_dnb_m1 <- rowMeans(expr_mat[dnb_gene, M1_clusters, drop = FALSE])
  expr_dnb_m2 <- rowMeans(expr_mat[dnb_gene, M2_clusters, drop = FALSE])
  dnb_log2fc <- log2((expr_dnb_m2 + 1e-6) / (expr_dnb_m1 + 1e-6))
  
  # Isolate direction-congruent co-expressed gene modules
  all_genes <- rownames(expr_mat)
  other_genes <- setdiff(all_genes, dnb_gene)
  coexpr_genes <- other_genes[signs[-1] == direction]
  
  if (length(coexpr_genes) == 0) {
    message(paste("No", direction, "co-expressed genes for", dnb_gene))
    return(NULL)
  }
  
  # Calculate aggregate mean log2FC across co-expression module
  log2fc_list <- sapply(coexpr_genes, function(g) {
    expr_m1 <- mean(expr_mat[g, M1_clusters])
    expr_m2 <- mean(expr_mat[g, M2_clusters])
    log2((expr_m2 + 1e-6) / (expr_m1 + 1e-6))
  })
  avg_coexpr_log2fc <- mean(log2fc_list)
  
  # Categorize coordinates into logical quadrants
  quadrant <- case_when(
    dnb_log2fc > 0 & avg_coexpr_log2fc > 0 ~ "I",
    dnb_log2fc < 0 & avg_coexpr_log2fc > 0 ~ "II",
    dnb_log2fc < 0 & avg_coexpr_log2fc < 0 ~ "III",
    dnb_log2fc > 0 & avg_coexpr_log2fc < 0 ~ "IV",
    TRUE ~ "Border"
  )
  
  # Persist categorical data
  df_save <- data.frame(
    gene = dnb_gene,
    direction = direction,
    dnb_log2fc = round(dnb_log2fc, 3),
    avg_coexpr_log2fc = round(avg_coexpr_log2fc, 3),
    quadrant = quadrant
  )
  if (!file.exists(csv_path)) {
    write.csv(df_save, csv_path, row.names = FALSE)
  } else {
    write.table(df_save, csv_path, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  }
  
  # Structure object frame for plotting aesthetics
  df_plot <- data.frame(
    x = dnb_log2fc,
    y = avg_coexpr_log2fc,
    label = dnb_gene
  )
  
  # Define background quadrants
  bg_df <- expand.grid(
    xmin = c(-Inf, 0), xmax = c(0, Inf),
    ymin = c(-Inf, 0), ymax = c(0, Inf),
    fill = c("lightblue", "mistyrose", "lightgray", "lightyellow")
  )
  # Automatically generate symmetrical axis boundaries + margin padding (0.5)
  x_min <- floor(dnb_log2fc - 0.5)
  x_max <- ceiling(dnb_log2fc + 0.5)
  
  max_abs_val <- max(abs(dnb_log2fc), abs(avg_coexpr_log2fc))
  axis_limit <- ceiling(max_abs_val + 0.5)  # Enforce padding constraints + ceiling function mapping
  p <- ggplot() +
    geom_rect(data = bg_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(data = df_plot, aes(x = x, y = y), color = ifelse(direction == "Positive", "red", "blue"), size = 4) +
    geom_text(data = df_plot, aes(x = x, y = y, label = label), vjust = -1.2, size = 4) +
    annotate("text", x = 1.5, y = 1.5, label = "I", fontface = "bold", size = 5, color = "black") +
    annotate("text", x = -1.5, y = 1.5, label = "II", fontface = "bold", size = 5, color = "black") +
    annotate("text", x = -1.5, y = -1.5, label = "III", fontface = "bold", size = 5, color = "black") +
    annotate("text", x = 1.5, y = -1.5, label = "IV", fontface = "bold", size = 5, color = "black") +
    theme_minimal() +
    scale_fill_identity() +
    labs(
      title = paste0("Quadrant Summary: ", dnb_gene, " (", direction, ")"),
      x = "DNB log2FC (M2 / M1)",
      y = paste0("Mean ", direction, " Co-expression log2FC (M2 / M1)")
    ) +
    coord_cartesian(xlim = c(-axis_limit, axis_limit), ylim = c(-axis_limit, axis_limit))
  
  ggsave(
    filename = file.path(outpath, paste0(dnb_gene, "_", direction, "_quadrant_summary.png")),
    plot = p, width = 5.5, height = 5.5
  )
}
plot_quadrant_combined <- function(dnb_gene, expr_mat, signs, outpath = "quadrant_combined_plots", csv_path = "quadrant_combined_data.csv") {
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
  
  M1_clusters <- c("5")
  M2_clusters <- c("12")
  
  # Calculate DNB log2FC
  expr_dnb_m1 <- rowMeans(expr_mat[dnb_gene, M1_clusters, drop = FALSE])
  expr_dnb_m2 <- rowMeans(expr_mat[dnb_gene, M2_clusters, drop = FALSE])
  dnb_log2fc <- log2((expr_dnb_m2 + 1e-6) / (expr_dnb_m1 + 1e-6))
  
  all_genes <- rownames(expr_mat)
  other_genes <- setdiff(all_genes, dnb_gene)
  pos_genes <- other_genes[signs[-1] == "Positive"]
  neg_genes <- other_genes[signs[-1] == "Negative"]
  
  # Quantify aggregate mean log2FC exclusively for positively/negatively correlated sub-modules
  get_avg_fc <- function(genes) {
    if (length(genes) == 0) return(NA)
    vals <- sapply(genes, function(g) {
      m1 <- mean(expr_mat[g, M1_clusters])
      m2 <- mean(expr_mat[g, M2_clusters])
      log2((m2 + 1e-6) / (m1 + 1e-6))
    })
    mean(vals)
  }
  
  pos_fc <- get_avg_fc(pos_genes)
  neg_fc <- get_avg_fc(neg_genes)
  
  # Construct consolidated plot dataframe
  df_plot <- data.frame(
    type = c("Positive", "Negative"),
    x = dnb_log2fc,
    y = c(pos_fc, neg_fc)
  )
  
  # Setup spatial quadrant mapping background
  bg_df <- expand.grid(
    xmin = c(-Inf, 0), xmax = c(0, Inf),
    ymin = c(-Inf, 0), ymax = c(0, Inf),
    fill = c("lightblue", "mistyrose", "lightgray", "lightyellow")
  )
  
  # Project symmetrical axis parameters relative to Cartesian origin
  max_abs_val <- max(abs(c(df_plot$x, df_plot$y)), na.rm = TRUE)
  axis_limit <- ceiling(max_abs_val + 0.5)
  
  # Conditional coordinate classification mapping function
  get_quadrant <- function(x, y) {
    if (is.na(x) || is.na(y)) return("NA")
    if (x > 0 && y > 0) return("I")
    if (x < 0 && y > 0) return("II")
    if (x < 0 && y < 0) return("III")
    if (x > 0 && y < 0) return("IV")
    return("Border")
  }
  
  df_plot$quadrant <- mapply(get_quadrant, df_plot$x, df_plot$y)
  
  # Archive structured results
  df_save <- data.frame(
    gene = dnb_gene,
    direction = df_plot$type,
    dnb_log2fc = round(df_plot$x, 3),
    avg_coexpr_log2fc = round(df_plot$y, 3),
    quadrant = df_plot$quadrant
  )
  
  if (!file.exists(csv_path)) {
    write.csv(df_save, csv_path, row.names = FALSE)
  } else {
    write.table(df_save, csv_path, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  }
  
  # Generate Figure output
  p <- ggplot() +
    geom_rect(data = bg_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(data = df_plot, aes(x = x, y = y, color = type), size = 4) +
    geom_text(data = df_plot, aes(x = x, y = y, label = type), vjust = -1.2, size = 4) +
    scale_color_manual(values = c("Positive" = "red", "Negative" = "blue")) +
    annotate("text", x = axis_limit - 0.5, y = axis_limit - 0.5, label = "I", fontface = "bold", size = 5) +
    annotate("text", x = -axis_limit + 0.5, y = axis_limit - 0.5, label = "II", fontface = "bold", size = 5) +
    annotate("text", x = -axis_limit + 0.5, y = -axis_limit + 0.5, label = "III", fontface = "bold", size = 5) +
    annotate("text", x = axis_limit - 0.5, y = -axis_limit + 0.5, label = "IV", fontface = "bold", size = 5) +
    coord_cartesian(xlim = c(-axis_limit, axis_limit), ylim = c(-axis_limit, axis_limit)) +
    theme_minimal() +
    scale_fill_identity() +
    labs(
      title = paste0("Quadrant Summary (Combined): ", dnb_gene),
      x = "DNB log2FC (Cluster3 / Cluster12)",
      y = "Mean Co-expression log2FC (Cluster3 / Cluster12)"
    )
  
  ggsave(
    filename = file.path(outpath, paste0(dnb_gene, "_combined_quadrant.png")),
    plot = p, width = 5.5, height = 5.5
  )
}

analyze_gene_by_correlation <- function(gene,
                                        seurat_obj,
                                        cluster_id = "6",
                                        trajectory = c("12", "5", "8", "3", "4", "1", "10"),
                                        top_n = 50,
                                        use_deg = TRUE,
                                        deg_genes = NULL,
                                        expr_mat_use = NULL,
                                        log2fc_vector = NULL,
                                        do_permutation_test = TRUE,
                                        n_perm = 1000) {
  if (use_deg && !is.null(deg_genes)) {
    allowed_genes <- deg_genes
  } else {
    allowed_genes <- NULL
  }
  
  result <- get_top_correlated(gene, seurat_obj, cluster_id, top_n, pval_cutoff = 0.05, allowed_genes, method = "spearman",balance_sign = TRUE)
  neighbors <- result$genes
  correlations <- result$correlations
  signs <- result$signs
  
  all_genes <- intersect(c(gene, neighbors), rownames(seurat_obj))
  if (length(all_genes) < 2) return(NULL)
  
  expr_mat <- get_cluster_expr_matrix(all_genes, seurat_obj, trajectory)
  expr_mat <- expr_mat[apply(expr_mat, 1, function(x) sd(x, na.rm = TRUE) > 0), , drop = FALSE]
  all_genes <- rownames(expr_mat)
  neighbors <- setdiff(all_genes, gene)
  if (length(neighbors) < 1) return(NULL)
  
  dir.create("cor_analysis", showWarnings = FALSE)
  plot_expression_heatmap(expr_mat, gene, correlations, signs, paste0("cor_analysis/", gene, "_coexp_heatmap.png"))
  plot_expression_lineplot(expr_mat, gene, paste0("cor_analysis/", gene, "_coexp_lineplot.png"))
  plot_quadrant_combined(gene, expr_mat, signs,
                         outpath = "quadrant_combined_plots",
                         csv_path = "quadrant_combined_data.csv")
  
  if (do_permutation_test && !is.null(expr_mat_use) && !is.null(log2fc_vector)) {
    res1 <- permutation_log2fc_test(
      dnb_gene = gene,
      expr_mat_use = expr_mat_use,
      signs = result$signs,
      direction = "Positive",
      M1_cluster = "5",
      M2_cluster = "12",
      n_perm = n_perm,
      return_full = TRUE,
      log2fc_vector = log2fc_vector
    )
    print(res1$result)
    plot_permutation_density_to_file(res1$permuted_means, res1$observed, res1$result$p_value, gene, "Positive")
    
    res2 <- permutation_log2fc_test(
      dnb_gene = gene,
      expr_mat_use = expr_mat_use,
      signs = result$signs,
      direction = "Negative",
      M1_cluster = "5",
      M2_cluster = "12",
      n_perm = n_perm,
      return_full = TRUE,
      log2fc_vector = log2fc_vector
    )
    print(res2$result)
    plot_permutation_density_to_file(res2$permuted_means, res2$observed, res2$result$p_value, gene, "Negative")
  }
}

dnb_gene <- gene
# Implement statistical validation via permutation framework
permutation_log2fc_test <- function(dnb_gene,
                                    expr_mat_use,
                                    M1_cluster = "5",
                                    M2_cluster = "12",
                                    signs,
                                    direction = c("Positive", "Negative"),
                                    n_perm = 1000,
                                    seed = 42,
                                    min_genes = 1,
                                    enable_test = TRUE,
                                    report_path = "permutation_test_results.csv",
                                    return_full = TRUE,
                                    log2fc_vector = NULL) {
  direction <- match.arg(direction)
  if (!enable_test) return(NULL)
  set.seed(seed)
  #direction =  "Negative"
  coexpr_genes <- setdiff(names(signs)[signs == direction], dnb_gene)
  if (length(coexpr_genes) < min_genes) {
    result <- data.frame(
      gene = dnb_gene,
      direction = direction,
      observed_mean_log2FC = NA,
      p_value = NA,
      n_coexpressed = length(coexpr_genes)
    )
    return(if (return_full) list(result = result, permuted_means = NA, observed = NA) else result)
  }
  
  dnb_log2fc <- log2fc_vector[dnb_gene]
  bg_genes <- setdiff(names(log2fc_vector), dnb_gene)
  rand_dir <- if ((direction == "Positive" && dnb_log2fc >= 0) || (direction == "Negative" && dnb_log2fc < 0)) "up" else "down"
  bg_pool <- if (rand_dir == "up") bg_genes[log2fc_vector[bg_genes] >= 0] else bg_genes[log2fc_vector[bg_genes] < 0]
  bg_pool <- bg_pool[abs(log2fc_vector[bg_pool]) < 7]  # Filter extreme artifact values
  
  observed_mean <- mean(log2fc_vector[coexpr_genes], na.rm = TRUE)
  permuted_means <- replicate(n_perm, {
    mean(log2fc_vector[sample(bg_pool, length(coexpr_genes))], na.rm = TRUE)
  })
  
  p_val <- if ((direction == "Positive" && dnb_log2fc >= 0) || (direction == "Negative" && dnb_log2fc < 0)) {
    mean(permuted_means >= observed_mean)
  } else {
    mean(permuted_means <= observed_mean)
  }
  
  result <- data.frame(
    gene = dnb_gene,
    direction = direction,
    observed_mean_log2FC = round(observed_mean, 3),
    p_value = signif(p_val, 3),
    n_coexpressed = length(coexpr_genes)
  )
  
  if (!file.exists(report_path)) {
    write.csv(result, report_path, row.names = FALSE)
  } else {
    write.table(result, report_path, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  }
  
  if (return_full) {
    return(list(result = result, permuted_means = permuted_means, observed = observed_mean))
  } else {
    return(result)
  }
}


# Global Pre-processing initialization
dev.off()
expr_mat_use <- get_cluster_expr_matrix(rownames(temp_sce), temp_sce, trajectory)
expr_mat_use <- expr_mat_use[rowMeans(expr_mat_use[, c("3", "12")]) != 0, ]
log2fc_vector <- log2((expr_mat_use[, "3"] + 1e-6) / 
                        (expr_mat_use[, "12"] + 1e-6))
names(log2fc_vector) <- rownames(expr_mat_use)
expr_mat_use["Serpinb2",]
log2fc_vector["Serpinb2"]
pheatmap(expr_mat_use[gene_list,],cluster_cols = F,scale = 'row')
# Iterative core execution block (across all defined DNBs)
for (gene in gene_list[44:65]) {
  analyze_gene_by_correlation(gene, temp_sce,
                              cluster_id = "6",
                              trajectory = trajectory,
                              top_n = 50,
                              use_deg = TRUE,
                              deg_genes = deg_genes,
                              expr_mat_use = expr_mat_use,
                              log2fc_vector = log2fc_vector,
                              do_permutation_test = TRUE,
                              n_perm = 1000)
}

gene <- gene_list[40]
gene <- "Abca7"
# Iterate processing across all identified DNB genes
for (gene in gene_list) {
  # analyze_gene_by_ppi(gene, temp_sce, symbol_to_protein, id_to_symbol, ppi_links, trajectory,
  #                     filter_combined_score = 700, use_deg = TRUE, deg_genes = deg_genes)

  analyze_gene_by_correlation(gene, filter_macrophage_cells, cluster_id = "1", trajectory = trajectory, 
                              top_n = 50, use_deg = T, deg_genes = deg_genes, do_permutation_test = TRUE)
}
print(ppi_list)  # Review output mapping of downstream DNB PPI targets
print(coexp_list)  # Review output mapping of downstream DNB co-expressed subsets
plasma_b_cell$seurat_clusters <- plasma_b_cell$anno_clusters
seurat_obj <- temp_sce <- filter_macrophage_cells
seurat_obj <- temp_sce <- plasma_b_cell
seurat_obj <- temp_sce <- filter_sce_mesen
temp_sce$seurat_clusters <- filter_sce_mesen$RNA_snn_res.0.8
seurat_obj$seurat_clusters <- filter_sce_mesen$RNA_snn_res.0.8
table(temp_sce$orig.ident,temp_sce$anno_clusters)
table(seurat_obj$orig.ident,seurat_obj$seurat_clusters)
# Placeholder logic for standard Seurat implementation
degs <- FindMarkers(seurat_obj, ident.1 = "3", ident.2 = "12", group.by = "seurat_clusters",only.pos = FALSE, min.pct =0, logfc.threshold = 0,test.use = "wilcox")
hist(degs[deg_genes,]$avg_log2FC)
degs["Serpinb2",]
# Enforce significance thresholds on log2FC and P-values
deg_genes <- rownames(degs)[abs(degs$avg_log2FC) > 0.5  & degs$p_val_adj < 0.05]
deg_genes <- rownames(degs)[abs(degs$avg_log2FC) > 0.5   & degs$p_val < 0.1]
back_deg_genes <- rownames(degs)[ abs(degs$avg_log2FC) > 1]
grep("^Klf", rownames(degs), value = TRUE)
dev.off()



library(ggplot2)
library(ggrepel)
library(dplyr)

# ==== Structural Parameter Setup ====
quadrant_file <- "quadrant_combined_data.csv"
perm_file <- "permutation_test_results.csv"
direction_to_plot <- "Positive"   # Toggle to "Positive" or "Negative" for respective directions
p_cutoff <- 0.2

# ==== Ingest Source Artifacts ====
df <- read.csv(quadrant_file)
perm_df <- read.csv(perm_file)

# ==== Directional Subset Isolation ====
df_dir <- df[df$direction == direction_to_plot, ]
perm_dir <- perm_df[perm_df$direction == direction_to_plot & perm_df$p_value < p_cutoff, ]

# ==== Map Annotation Labels to Target DNBs ====
sig_genes <- perm_dir$gene
label_df <- df_dir[df_dir$gene %in% sig_genes, ]

# ==== Harmonize Coordinate Bounds ====
max_val <- max(abs(c(df_dir$dnb_log2fc, df_dir$avg_coexpr_log2fc)), na.rm = TRUE)
axis_limit <- ceiling(max_val + 0.5)

# ==== Generate Scatter Visualization ====
ggplot(df_dir, aes(x = dnb_log2fc, y = avg_coexpr_log2fc)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(color = ifelse(direction_to_plot == "Positive", "red", "blue"), size = 3) +
  geom_text_repel(data = label_df, aes(label = gene), size = 3,
                  max.overlaps = 1000, box.padding = 0.01, point.padding = 0.03) +
  theme_minimal() +
  labs(
    title = paste0("DNB vs ", direction_to_plot, " Co-expression Module (p < ", p_cutoff, ")"),
    x = "DNB log2FC (Cluster3 / Cluster12)",
    y = paste0("Mean ", direction_to_plot, " Co-expression log2FC")
  ) 
  #+
  #coord_cartesian(xlim = c(-axis_limit, axis_limit), ylim = c(-axis_limit, axis_limit))


library(ggplot2)
library(ggrepel)
library(dplyr)

# ==== Structural Parameter Setup ====
quadrant_file <- "quadrant_combined_data.csv"
perm_file <- "permutation_test_results.csv"
direction_to_plot <- "Negative"
p_cutoff <- 0.2
top_n <- 10

# ==== Ingest Source Artifacts ====
df <- read.csv(quadrant_file)
df <- df[!duplicated(df), ]
perm_df <- read.csv(perm_file)

df_dir <- df[df$direction == direction_to_plot, ]
perm_dir <- perm_df[perm_df$direction == direction_to_plot & perm_df$p_value < p_cutoff, ]
highlight_color <- ifelse(direction_to_plot == "Positive", "#D7263D", "#377EB8")
# ==== Extract maximum cross-coordinate genes (Top N logic mapping) ====
get_significant_label_df_by_quadrant <- function(
    df, perm_df, direction = "Positive",
    top_n_I = 5, top_n_III = 10,  # Limits for Positive module targets
    top_n_II = 6, top_n_IV = 6    # Limits for Negative module targets
) {
  df_dir <- df[df$direction == direction, ]
  perm_dir <- perm_df[perm_df$direction == direction & perm_df$p_value < p_cutoff, ]
  df_sig <- df_dir[df_dir$gene %in% perm_dir$gene, ]
  
  if (direction == "Positive") {
    # Quadrant I (Activating mechanism)
    df_I <- df_sig %>%
      filter(dnb_log2fc > 0 & avg_coexpr_log2fc > 0) %>%
      arrange(desc(abs(dnb_log2fc) + abs(avg_coexpr_log2fc))) #%>%
      #head(top_n_I)
    
    # Quadrant III (Synergistic downregulation mechanism)
    df_III <- df_sig %>%
      filter(dnb_log2fc < 0 & avg_coexpr_log2fc < 0) %>%
      arrange(desc(abs(dnb_log2fc) + abs(avg_coexpr_log2fc))) #%>%
      #head(top_n_III)
    
    label_df <- bind_rows(df_I, df_III)
    
  } else if (direction == "Negative") {
    # Quadrant II (Negative regulation - Upregulated target)
    df_II <- df_sig %>%
      filter(dnb_log2fc < 0 & avg_coexpr_log2fc > 0) %>%
      arrange(desc(abs(dnb_log2fc) + abs(avg_coexpr_log2fc))) #%>%
      #head(top_n_II)
    
    # Quadrant IV (Negative regulation - Downregulated target)
    df_IV <- df_sig %>%
      filter(dnb_log2fc > 0 & avg_coexpr_log2fc < 0) %>%
      arrange(desc(abs(dnb_log2fc) + abs(avg_coexpr_log2fc))) #%>%
      #head(top_n_IV)
    
    label_df <- bind_rows(df_II, df_IV)
    
  } else {
    stop("direction must be 'Positive' or 'Negative'")
  }
  
  return(label_df)
}

# Positive correlation module: Highlight Quadrant III (downregulated DNBs)
# label_df <- get_significant_label_df_by_quadrant(
#   df, perm_df, direction = "Positive",
#   top_n_I = 2,
#   top_n_III = 10
# )

# Negative correlation module: Highlight repressors/activators in Quadrants II/IV
label_df <- get_significant_label_df_by_quadrant(
  df, perm_df, direction = "Negative",
  top_n_II = 8,
  top_n_IV = 5
)
# ==== Append highlighting logic column for aesthetics ====
df_dir$highlight <- ifelse(df_dir$gene %in% label_df$gene, "Highlighted", "Other")

# ==== Generate Advanced Visualization ====
ggplot(df_dir, aes(x = dnb_log2fc, y = avg_coexpr_log2fc)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  
  # Point size and hue aesthetic mapping
  geom_point(aes(color = highlight, size = highlight), alpha = 0.9) +
  scale_color_manual(values = c("Highlighted" = highlight_color, "Other" = "#B0B0B0")) +
  scale_size_manual(values = c("Highlighted" = 4.5, "Other" = 1.5)) +
  
  # Inject textual annotations
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 3.5,
    max.overlaps = 1000,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey60",
    segment.size = 0.3
  ) +
  
  labs(
    title = paste0("DNB vs ", direction_to_plot, " Co-expression Module (p < ", p_cutoff, ")"),
    x = "DNB log2FC (Cluster12 / Cluster6)",
    y = paste0("Mean ", direction_to_plot, " Co-expression log2FC")
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.line = element_line(size = 0.8),
    axis.ticks = element_line(size = 0.6),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

library(ggplot2)
library(ggrepel)
library(grid)

# Extrapolate absolute coordinate bounds
x_range <- range(df_dir$dnb_log2fc, na.rm = TRUE)
y_range <- range(df_dir$avg_coexpr_log2fc, na.rm = TRUE)

# ==== Refined Visualization Output ====
ggplot(df_dir, aes(x = dnb_log2fc, y = avg_coexpr_log2fc)) +
  # Map structural guides
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  
  # Target nodes
  geom_point(aes(color = highlight, size = highlight), alpha = 0.9) +
  scale_color_manual(values = c("Highlighted" = highlight_color, "Other" = "#B0B0B0")) +
  scale_size_manual(values = c("Highlighted" = 4.5, "Other" = 1.5)) +
  
  # Map functional labels
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 3.5,
    max.overlaps = 1000,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey60",
    segment.size = 0.3
  ) +
  
  # Add custom axis arrows to quadrant plots
  geom_segment(aes(x = x_range[1], xend = x_range[2], y = 0, yend = 0),
               arrow = arrow(length = unit(0.25, "cm")), size = 0.7) +
  geom_segment(aes(x = 0, xend = 0, y = y_range[1], yend = y_range[2]),
               arrow = arrow(length = unit(0.25, "cm")), size = 0.7) +
  
  # Append descriptive context elements
  labs(
    title = paste0("DNB vs ", direction_to_plot, " Co-expression Module (p < ", p_cutoff, ")"),
    x = "DNB log2FC (Cluster3 / Cluster12)",
    y = paste0("Mean ", direction_to_plot, " Co-expression log2FC")
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.line = element_blank(),  # Hide native default axis constraints
    axis.ticks = element_line(size = 0.6),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )










library(ggplot2)
library(ggrepel)
library(dplyr)
library(grid)  # Required for rendering coordinate arrows

# ==== Automatic Range Normalization ====
xlim_max <- max(abs(df_dir$dnb_log2fc), na.rm = TRUE)
ylim_max <- max(abs(df_dir$avg_coexpr_log2fc), na.rm = TRUE)
limit <- ceiling(max(c(xlim_max, ylim_max)) + 0.5)

# ==== DNB vs Co-expression Rendering ====
ggplot(df_dir, aes(x = dnb_log2fc, y = avg_coexpr_log2fc)) +
  
  # Centralized coordinate axes with directional arrows
  geom_segment(aes(x = -limit, xend = limit, y = 0, yend = 0),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               color = "black", linewidth = 0.6) +
  geom_segment(aes(x = 0, xend = 0, y = -limit, yend = limit),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               color = "black", linewidth = 0.6) +
  
  # Scatter Point Mapping
  geom_point(aes(color = highlight, size = highlight), alpha = 0.85) +
  scale_color_manual(values = c("Highlighted" = highlight_color, "Other" = "gray70")) +
  scale_size_manual(values = c("Highlighted" = 4.5, "Other" = 1.5)) +
  
  # Data Feature Labels
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 3.5,
    max.overlaps = 1000,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey60"
  ) +
  
  # Spatial Coordinate Constraints
  coord_cartesian(xlim = c(-limit, limit), ylim = c(-limit, limit)) +
  
  # Annotation & Thematic Styling
  labs(
    title = paste0("DNB vs ", direction_to_plot, " Co-expression Module (p < ", p_cutoff, ")"),
    x = "DNB log2FC (Cluster3 / Cluster12)",
    y = paste0("Positive Co-expression log2FC")
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.line = element_blank(),  # Remove default native axes
    axis.ticks = element_line(size = 0.5),
    axis.text = element_text(size = 13, color = "black", face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )
dev.off()

# Assuming DNB_list and expr_mat are pre-initialized in the global environment
# RPL (Ribosomal Protein Large) and RPS (Ribosomal Protein Small) genes can be extracted via regex, or mapped via predefined lists
rpl_rps_genes <- grep("^il|^Rps", rownames(expression_matrix), value = TRUE)

# Intersect identified DNBs with established ribosomal/inflammatory feature lists
selected_genes <- intersect(DNB_list$X5, rpl_rps_genes)

# Extract relative expression matrices for selected intersecting subsets
DNBs_expression.df <-DNBs_expression(expression_matrix,trajectory,selected_genes)
exps <- t(DNBs_expression.df)
# Generate hierarchically clustered heatmap
library(pheatmap)
pheatmap(exps, 
         scale = "row",        # Standardize rows (Z-score mapping)
         cluster_cols = FALSE, # Preserve topological cluster sequencing
         show_rownames = TRUE, # Display specific gene nomenclatures
         show_colnames = TRUE, # Display originating cluster identities
         main = "RPL and RPS Genes Expression Trend")

target_cluster <- c("1","4","9") # Fibroblast DNB Target: Cluster 4 (Early Critical State)
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster,unique(DNB_list$X4))
exps <- t(DNBs_expression.df)
write.table(exps, "fibroblast_4exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

target_cluster <- c("9","6","7") # Fibroblast DNB Target: Cluster 7 (Late Critical State)
DNBs_expression.df <- DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X6))
exps <- t(DNBs_expression.df)
write.table(exps, "fibroblast_6exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

target_cluster <- c("0","2","1") # Granulocyte DNB Target: Cluster 2 (Early Critical State)
DNBs_expression.df <- DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X2))
exps <- t(DNBs_expression.df)
write.table(exps, "2exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

target_cluster <- c("0","5","7") # Granulocyte DNB Target: Cluster 5 (Early Critical State)
DNBs_expression.df <- DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X5))
exps <- t(DNBs_expression.df)
write.table(exps, "Granulocytes_057exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

target_cluster <- c("3","4","1") # Macrophage DNB Target: Cluster 4 (Late Critical State)
DNBs_expression.df <- DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X4))
exps <- t(DNBs_expression.df)
write.table(exps, "macrophages_341exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

target_cluster <- c("5","8","3") # Macrophage DNB Target: Cluster 8 (Early Critical State)
DNBs_expression.df <- DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X8))
exps <- t(DNBs_expression.df)
write.table(exps, "macrophages_583exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)


target_cluster <- c("12","5","8") # Macrophage DNB Target: Cluster 5 (Early Critical State)
DNBs_expression.df <- DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X5))
exps <- t(DNBs_expression.df)
write.table(exps, "macrophages_1258exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

library(SingleCellExperiment)

# Verify the integration of specific DNB module genes within the global expression matrix
gene_list <- DNB[DNB %in% rownames(subsce_trans)]

# Extract temporal cluster topologies, assuming storage within subsce_trans$cluster
clusters_of_interest <- c("9", "6", "7")

# Initialize result structures
cluster_sd_list <- list()
cluster_cor_list <- list()

for (clust in clusters_of_interest) {
  # Extract cells inherent to the specified cluster parameter
  cells_in_cluster <- colnames(subsce_trans)[subsce_trans$seurat_clusters == clust[1]]

  # Extract expression array (applicable formats: logcounts, counts, normalized expr, etc.)
  expr_mat <- logcounts(subsce_trans)[gene_list, cells_in_cluster, drop = FALSE]

  # Calculate expression standard deviation per feature
  gene_sd <- apply(expr_mat, 1, sd)

  # Calculate Pearson correlation matrices
  gene_cor <- cor(t(as.matrix(expr_mat)), method = "pearson")

  # Store output iteratively
  cluster_sd_list[[as.character(clust)]] <- gene_sd
  cluster_cor_list[[as.character(clust)]] <- gene_cor
}

# Construct structured dataframes from list outputs
df_list <- lapply(names(cluster_sd_list), function(clust) {
  data.frame(
    Gene = names(cluster_sd_list[[clust]]),
    sd = cluster_sd_list[[clust]],
    Cluster = clust,
    stringsAsFactors = FALSE
  )
})

# Aggregate multi-cluster metrics
df_all <- do.call(rbind, df_list)

# Transform matrix topology to wide format, setting headers to "cluster_Corr"
library(tidyr)
df_wide <- df_all %>%
  mutate(Cluster = paste0(Cluster, "_sd")) %>%
  pivot_wider(names_from = Cluster, values_from = sd)

# Display processed header structure
head(df_wide)
write.csv(df_wide, file = "./临界态图/Figure2-2D-cytoscape_node_attributes_wide.csv", row.names = FALSE)
### Isolate specific DNBs for downstream visualization
# Extract Standard Deviation (SD) and Pearson Correlation Coefficients (PCC) for cluster 5
sd_5 <- cluster_sd_list[["6"]]
cor_5 <- cluster_cor_list[["6"]]

# Calculate mean PCC per gene (excluding self-correlations)
mean_pcc_5 <- rowMeans(cor_5[rownames(cor_5) != colnames(cor_5), ], na.rm = TRUE)

# Aggregate SD and PCC metrics
combined_metrics <- data.frame(
  Gene = names(sd_5),
  SD = sd_5,
  PCC = mean_pcc_5[names(sd_5)]
)

# Apply restrictive filtering thresholds (e.g., top 30th percentile)
high_sd_genes <- combined_metrics[order(-combined_metrics$SD), ][1:round(nrow(combined_metrics) * 0.3), ]
high_pcc_genes <- combined_metrics[order(-combined_metrics$PCC), ][1:round(nrow(combined_metrics) * 0.3), ]

# Establish intersecting core network modules
selected_genes <- intersect(high_sd_genes$Gene, high_pcc_genes$Gene)
selected_genes <- c("C1qtnf3","H19","Col11a1")
# Restructure cluster_sd_list into continuous dataframe
sd_df <- lapply(names(cluster_sd_list), function(clust) {
  gene_sd <- cluster_sd_list[[clust]]
  data.frame(
    Gene = names(gene_sd),
    SD = gene_sd,
    Cluster = clust
  )
}) %>% bind_rows()

# Retain exclusively targeted genes of interest
sd_df_selected <- sd_df %>% filter(Gene %in% selected_genes)

# Map Cluster structures as ordered factors
sd_df_selected$Cluster <- factor(sd_df_selected$Cluster, levels = c("9", "6", "7"))

# Generate longitudinal trend line plot
ggplot(sd_df_selected, aes(x = Cluster, y = SD, group = Gene, color = Gene)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  theme_minimal(base_size = 14) +
  labs(title = "Standard Deviation Trend (Selected DNB Genes)", x = "Cluster", y = "SD Value") +
  theme(legend.position = "right") +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),         # Eliminate background grid layers
    plot.background = element_blank(),    # Ensure global plot transparency
    panel.background = element_blank()    # Ensure plot area transparency
  )


library(SingleCellExperiment)
library(dplyr)
library(tidyr)

expr_list <- list()

for (clust in clusters_of_interest) {
  # Extract inherent cluster cells
  cells <- colnames(subsce_trans)[subsce_trans$seurat_clusters == clust]

  # Extract targeted expression matrix
  expr_mat <- logcounts(subsce_trans)[selected_genes, cells, drop = FALSE]

  # Melt to long format and filter non-expressing (0 value) observations
  expr_df <- as.data.frame(as.matrix(expr_mat)) %>%
    tibble::rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "Cell", values_to = "Expression") %>%
    mutate(Cluster = clust) %>%
    filter(Expression != 0)  # <- Critical exclusion boundary

  expr_list[[clust]] <- expr_df
}


# Synthesize mult-cluster global expression matrix
expression_long <- bind_rows(expr_list)
library(ggplot2)
expression_long$Cluster <- factor(expression_long$Cluster, levels = c("9", "6", "7"))
ggplot(expression_long, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.3) +
  facet_wrap(~Gene, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "red", size = 1.5) +
  theme_minimal() +
  labs(title = "Expression of Selected Genes across Clusters",
       y = "Log Expression", x = "Cluster")
ggplot(expression_long, aes(x = Cluster, y = Expression, fill = Cluster)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0, size = 0, alpha = 0) +
  facet_wrap(~Gene, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "red", size = 1.5) +
  theme_minimal()+
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),         # Eliminate background grid layers
    plot.background = element_blank(),    # Ensure global plot transparency
    panel.background = element_blank()    # Ensure plot area transparency
  )
# Execute non-parametric Kruskal-Wallis permutation testing per feature
kruskal_test <- sd_df %>%
  group_by(Gene) %>%
  summarise(p_value = kruskal.test(SD ~ Cluster)$p.value)

# Implement multiple testing significance correction (FDR methodology)
kruskal_test$FDR <- p.adjust(kruskal_test$p_value, method = "fdr")
library(ggplot2)

ggplot(sd_df, aes(x = Cluster, y = SD, fill = Cluster)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3, color = "gray") +
  stat_summary(fun = mean, geom = "point", color = "red", size = 2) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "SD of Gene Expression across Clusters",
       y = "Expression SD", x = "Cluster")
# Extract statistically significant differential modules
sig_genes <- kruskal_test %>% filter(FDR < 0.05)

library(reshape2)
library(ggplot2)

# Synthesize SD lists into unified dataframes
sd_df <- do.call(cbind, cluster_sd_list)
sd_df <- as.data.frame(sd_df)
sd_df$Gene <- rownames(sd_df)

# Matrix pivot to long format
sd_long <- melt(sd_df, id.vars = "Gene", variable.name = "Cluster", value.name = "SD")

# Generate plotting visualization
ggplot(sd_long, aes(x = Cluster, y = SD, group = Gene, color = Gene)) +
  geom_line() +
  geom_point() +
  labs(title = "Expression SD Trend of Genes Across Clusters",
       x = "Cluster", y = "Standard Deviation (SD)") +
  theme_minimal()

# Consolidate standard deviation data arrays
sd_df <- do.call(cbind, cluster_sd_list)
sd_df <- as.data.frame(sd_df)
sd_df$Gene <- rownames(sd_df)

# Map matrix output to long architecture
sd_long <- melt(sd_df, id.vars = "Gene", variable.name = "Cluster", value.name = "SD")

# Assert temporal clustering hierarchies
sd_long$Cluster <- factor(sd_long$Cluster, levels = c("9", "6", "7"))

# Produce visualization renders
ggplot(sd_long, aes(x = Cluster, y = SD, group = Gene)) +
  geom_line(color = "gray", alpha = 0.4) +  # Render underlying trajectories in grey shading
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "red", size = 1.2) +  # Emphasize central tendency with red mappings
  stat_summary(aes(group = 1), fun = mean, geom = "point", color = "red", size = 2) +   # Emphasize mean locus geometry
  labs(title = "SD Expression Trend Across Clusters",
       x = "Cluster", y = "Standard Deviation (SD)") +
  theme_minimal(base_size = 14)

library(ComplexHeatmap)
library(circlize)
# Map custom visual gradients
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
# Generate Complex Heatmap arrays
ht1 <- Heatmap(cluster_cor_list[["12"]],cluster_rows = T, name = "Cor", column_title = "Cluster 12", col = col_fun)
ht2 <- Heatmap(cluster_cor_list[["5"]], cluster_rows = T, name = "Cor", column_title = "Cluster 5", col = col_fun)
ht3 <- Heatmap(cluster_cor_list[["8"]], cluster_rows = T, name = "Cor", column_title = "Cluster 8", col = col_fun)
ht4 <- Heatmap(cluster_cor_list[["3"]], cluster_rows = T, name = "Cor", column_title = "Cluster 3", col = col_fun)
ht5 <- Heatmap(cluster_cor_list[["4"]], cluster_rows = T, name = "Cor", column_title = "Cluster 4", col = col_fun)
ht6 <- Heatmap(cluster_cor_list[["1"]] %>% {.[is.na(.)] <- 0; .},cluster_rows = T,  name = "Cor", column_title = "Cluster 1", col = col_fun)

# Consolidate aligned displays
ht1 + ht2 + ht3 +ht4 + ht5 + ht6

avg_cors <- sapply(cluster_cor_list, function(mat) {
  mean(abs(mat[upper.tri(mat)]), na.rm = TRUE)
})

# Construct dedicated visual mapping frames
cor_df <- data.frame(Cluster = names(avg_cors), AvgCorrelation = avg_cors)
cor_df$Cluster <- factor(cor_df$Cluster, levels = c("0", "2", "1"))

ggplot(cor_df[1:3,], aes(x = Cluster, y = AvgCorrelation, group = 1)) +
  geom_line(color = "red", size = 1.5) +
  geom_point(color = "red", size = 3) +
  theme_minimal(base_size = 14) +
  labs(title = "Variance Dynamic Topology (DNB Candidates)",
       x = "Cluster",
       y = "Average Pearson Correlation") +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),         # Suppress grid formatting constraints
    plot.background = element_blank(),    # Force transparent plot framework mappings
    panel.background = element_blank()    # Force transparent plot panel mappings
  )


# Initialize structured containers to track gene-specific correlation shifts
genewise_cor_list <- list()

for (clust in clusters_of_interest) {
  cells_in_cluster <- colnames(subsce_trans)[subsce_trans$seurat_clusters == clust]
  expr_mat <- logcounts(subsce_trans)[DNB[-1], cells_in_cluster, drop = FALSE]

  # Calculate cross-gene correlative mappings
  gene_cor <- cor(t(as.matrix(expr_mat)), method = "pearson")

  # Define specific gene-to-gene average linkages (extracting magnitude absolutes)
  gene_avg_cor <- apply(gene_cor, 1, function(row) {
    mean(abs(row[-which(is.na(row) | row == 1)]))  # Screen for missing mappings (NA) and self-relationships (value=1)
  })

  # Push logic maps into global environment variable
  genewise_cor_list[[clust]] <- gene_avg_cor
}
# Synthesize bound dataframes based on matrix length
cor_gene_df <- do.call(rbind, lapply(names(genewise_cor_list), function(clust) {
  data.frame(
    Gene = names(genewise_cor_list[[clust]]),
    AvgCorrelation = genewise_cor_list[[clust]],
    Cluster = clust
  )
}))
library(ggplot2)
cor_gene_df$Cluster <- factor(cor_gene_df$Cluster, levels = c("9", "6", "7"))
ggplot(cor_gene_df, aes(x = Cluster, y = AvgCorrelation, group = Gene)) +
  geom_line(color = "gray", alpha = 0.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "red", size = 1.2) +
  theme_minimal(base_size = 14) +
  labs(title = "Variance Dynamic Topology (DNB Candidates)",
       x = "Cluster", y = "Average Correlation Magnitude (Absolute Value)") +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),         # Suppress grid formatting constraints
    plot.background = element_blank(),    # Force transparent plot framework mappings
    panel.background = element_blank()    # Force transparent plot panel mappings
  )

VlnPlot()
## Perform Clustervis hierarchical profiling
#exps <- exps[!duplicated(exps[,1]), ]
#rownames(exps) <- exps[,1]
#exps <- exps[,-1]

exps <- as.matrix(exps)
exps <-na.omit(exps)
head(exps)
# Load ClusterGVis library; Note: Avoid concurrent Monocle 3 implementations to circumvent dependency inheritances & conflicts
library(ClusterGVis)
library(sf)

sce.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

temp_sce2.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
# Supply grouped means as standard reference coordinates, implicitly de-duplicating specified loci
st.data1 <- prepareDataFromscRNA(object = temp_sce,
                                 diffData = top10,
                                 showAverage = TRUE)
ck <- clusterData(exp = exps,
                  #cluster.method = "mfuzz",
                  cluster.method = "kmeans",
                  cluster.num = 4)
# Map supplemental gene identifiers
markGenes = rownames(exps)[sample(1:nrow(exps),32,replace = F)]
pdf('addgene.pdf',height = 10,width = 6,onefile = F)
visCluster(object = ck,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes)
dev.off()
# Configure functional enrichments using standard Mus musculus genomic definitions
pdf('term.pdf',height = 20,width = 12,onefile = F)
library(org.Mm.eg.db)
#library(org.Hs.eg.db)
enrich <- enrichCluster(#object = ck,
                        object = st.data1,
                        OrgDb = org.Mm.eg.db,
                        type = "BP",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 5201314
                        #,
                        #add.gene =TRUE
)

# Erase secondary gene annotations if previously bound for structural mappings
enrich <- subset(enrich, select = -geneID)
head(enrich,3)
pdf('term2.pdf',height = 12,width = 12,onefile = F)
  visCluster(#object = ck,
           object = st.data1,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           #markGenes = markGenes,
           markGenes = top10$gene,
           markGenes.side = "left",
           genes.gp = c('italic',fontsize = 12,col = "black"),
           annoTerm.data = enrich,
           line.side = "left",
           #go.col = rep(ggsci::pal_d3()(15),each = 5),
           go.size = 8,
           cluster.order = c(1:14),
           #mulGroup = c(22,14,37,61,42,4),
           #mline.col = c(ggsci::pal_lancet()(7))
           )
dev.off()
# Render sequential trend representations
visCluster(object = ck,
           plot.type = "line",
           ms.col = c("green","orange","red"), # Customize spectral aesthetics
           add.mline = TRUE  # Overlay centralized medians
)

# Algorithmically classify shifting functional trends per gene expression trace
trends <- apply(exps, 1, function(row_vec) {
  diff_vec <- diff(row_vec)
  if (all(diff_vec > 0)) {
    return("Gradually Increasing")
  } else if (all(diff_vec < 0)) {
    return("Gradually Decreasing")
  } else if (diff_vec[1] > 0 && diff_vec[2] < 0) {
    return("Critical State High (Inverted-U)")
  } else if (diff_vec[1] < 0 && diff_vec[2] > 0) {
    return("Critical State Low (U-shape)")  # Modified classification mapping framework
  } else {
    return("Unclassified")
  }
})

# Isolate corresponding indexes mapping back to classified hierarchical clusters
sorted_idx <- order(trends)

# Synthesize sorted matrices based on ordered trend indexes
exps <- exps[sorted_idx, ]
sorted_trends <- as.data.frame(trends[sorted_idx])
colnames(sorted_trends) <- c("Category")
# Map parameter arguments explicitly tailored to heatmaps
heatmap_params <- list(
  col = colorRampPalette(c("blue", "white", "red"))(100),  # Topographic hue mappings
  scale = "row",  # Apply Z-score constraint standardizations across specific array vectors
  main = "Expression Heatmap"  # Figure title declarations
)

# Produce final Complex Heatmap implementations
exps <- merge(exps, sorted_trends,by = "row.names")
rownames(exps) <- exps[,1]
exps <- exps[,-1]
exps %>%  group_by(Category) %>% pheatmap(exps,
         cluster_rows = F,
         cluster_cols = F,# Disable native row clustering geometries
         clustering_distance_cols = "euclidean",  # Distance algorithmic specification
         clustering_method = "complete",  # Clustering calculation rulesets
         annotation_row = sorted_trends,  # Inject dynamic trend annotations
         col = heatmap_params$col,  # Define custom visualization gradients
         scale = "row",  # Retain standard dimensional normalization boundaries
         main = heatmap_params$main ) # Declare defined title spaces

# Call ComplexHeatmap logic dependencies explicitly
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#4169E1", "white", "#DC143C"))
row_annotation = rowAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#DC143C", "#4169E1","#FF4500","#006400")),
                       labels = c("Gradually Increasing","Gradually Decreasing","Critical State High","Critical State Low"),
                       labels_gp = gpar(col = "white", fontsize = 6)))

Heatmap(t(scale(t(exps))),
            col = col_fun,
            cluster_columns = F,
            cluster_rows = T,
            #left_annotation  = row_annotation,
            row_split = sorted_trends$Category,
            show_heatmap_legend = T,
            border = T,
            show_column_names = T,
            show_row_names = T,
            row_names_gp = gpar(fontsize = 8),
            heatmap_width = unit(1, "npc"),
            heatmap_height = unit(1, "npc")
         )






#################### scRank Methodology: Identification of Optimal Cellular Drug Targets ###################################
##### Package demonstration implementation framework
library(devtools)
devtools::install_github("rikenbit/rTensor")
devtools::install_github("ZJUFanLab/scRank")
help(scRank)
library(rTensor)
library(scRank)
seuratObj <- system.file("extdata", "AML_object.rda", package="scRank")
load(seuratObj)
obj <- CreateScRank(input = seuratObj,
                    species = 'mouse',
                    cell_type = 'labels',
                    target = 'Brd4')
obj <- scRank::Constr_net(obj)







obj <- CreateScRank(input = temp_sce2,  # 'Input' dictates expression signatures utilizing Seurat object metadata to define specified cell phenotypes
                    species = "mouse", # Reference organism parameters
                    #cell_type = 'seurat_clusters',
                    cell_type = 'label',# Specify descriptive column holding cellular identity tags within active metadata frameworks
                    #drug = ,  # Identify defined pharmaceutical agents cataloged explicitly in the default utile_database reference map
                    target = 'F9', # Alternatively pass precise gene identifiers. Bypasses pharmaceutical databases if the intended MoA is directly linked to an explicit target
                    type = "antagonist" , # Designate Mechanisms of Action (MOAs), predominantly parsing definitions as 'antagonist' or 'agonist'. Default config utilizes 'antagonist' parameters
                    if_cluster = FALSE# Boolean logic specifying clustering demands over the targeted scRNA-seq expression dataset. Set natively to FALSE
                    #,var.genes = NULL # var.genes (optional) pass specific vector structures dictating discrete target genes for secondary clustering implementations
              )
utile_database("Aspirin")
# Initialize Gene Regulatory Network (GRN) mappings
obj <- scRank::Constr_net(obj,
                           select_ratio = 0.5, # percentage of cells selected randomly across global totals parameters.
                           n_selection = 10, # Defined thresholds for n_selection generation protocols to construct network arrays tailored across target cell-types.
                           cut_ratio = 0.95, # Identify edge linkages specifying the lowest integrated confidence weight. Sequester weaker values retaining predominant interaction signals
                           keep_ratio = 0.25, # Define ratios to retain within established subset hierarchies
                           min_cells = 25,
                           n.core = NULL,
                           n.core_cp = 4,
                           use_py = F,
                           env = "base")

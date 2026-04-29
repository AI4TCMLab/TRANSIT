###### Validation Set Analysis #####
rm(list= ls())
options(stringsAsFactors = F)
getwd()
setwd("/home/tsh/GSE163129_RAW/")
library(batchelor)
library(cluster)
library(corrplot)
library(dynamicTreeCut)
library(edgeR)
library(gplots)
library(gridExtra)
library(igraph)
library(limma)
library(org.Mm.eg.db)
library(pheatmap)
library(psych)
library(scater)
library(scran)
library(SingleCellExperiment)
library(stringr)
library(BioTIP)
library(progress)
library(tidyverse)
library(SingleCellExperiment)
library(DropletUtils)
library(patchwork)
library(Seurat)
library(slingshot)
library(RColorBrewer)
options(future.globals.maxSize = 1e9)
getwd()
rawdata_path <- "/home/tsh/GSE163129_RAW/"  # Path to single-cell data files for different treatments
filename <- list.files(rawdata_path)
rawdata_path <- paste(rawdata_path,filename,sep = '')
# Use Read10X_h5 to load .H5 files for each timepoint and generate Seurat objects
sceList <- lapply(rawdata_path[1:5], function(x){
  obj <- CreateSeuratObject(counts = Read10X_h5(x),
                            project = str_split(x,'/')[[1]][5] %>% str_split(.,'_')%>%.[[1]]%>%.[2])
})
project <- lapply(rawdata_path[1:5], function(x){
  str_split(x,'/')[[1]][5] %>% str_split(.,'_')%>%.[[1]]%>%.[2]
})
# Merge Seurat objects: use the first element of sceList as the seed and iteratively merge remaining objects.
names(sceList) <- project
sce <- merge(sceList[[1]],sceList[-1],
             add.cell.ids = names(sceList), # add.cell.ids specifies identifiers for merged groupings
             project = 'Bobo'  # project parameter sets the master project name to 'Bobo'
)
head(sce@meta.data, 5)
sce <- PercentageFeatureSet(sce,'^mt',col.name = 'percent_mt')
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3,pt.size = 0,group.by = 'orig.ident')
# Quality Control (QC)
dim(sce)
sce <- subset(sce,subset = nFeature_RNA > 407 & nFeature_RNA < 7825 & percent_mt < 10)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce)
# Cell cycle analysis (Pending/Optional)
# s.genes <- cc.genes.updated.2019$s.genes
# g2m.genes <- cc.genes.updated.2019$g2m.genes
# sce <- CellCycleScoring(sce, 
#                            s.features = s.genes, 
#                            g2m.features = g2m.genes, 
#                            set.ident = TRUE)
# sce$CC.Difference <- sce$S.Score - sce$G2M.Score
# sce <- ScaleData(sce, 
#                      vars.to.regress = "CC.Difference", 
#                      features = rownames(sce))
# sce <- RunPCA(sce, 
#                   features = VariableFeatures(sce), 
#                   nfeatures.print = 10)
# sce <- RunPCA(sce, features = c(s.genes, g2m.genes))
# DimPlot(sce)
#sce <- ScaleData(sce)
#sce <- RunPCA(sce)

# Perform batch correction using "harmony"
sce <- IntegrateLayers(
  object = sce, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:30)
for (i in c(0.01,0.05,0.1,0.2,0.3,0.5,0.8,1,2)){
  sce <- FindClusters(sce,resolution = i)
}
clustree(sce ,prefix = "RNA_snn_res.")
sce <- FindClusters(sce, resolution = 2, cluster.name = "harmony_clusters")
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
sce <- RunTSNE(sce, reduction = "harmony", dims = 1:30, reduction.name = "tsne.harmony")
sce <- JoinLayers(sce)
sce

## Automated annotation via SingleR
mrsd.se <- ImmGenData()
MI.integrated.se <- as.SingleCellExperiment(sce, assay = "RNA")
mrsd.common <- intersect(rownames(MI.integrated.se), rownames(mrsd.se))
mrsd.se <- mrsd.se[mrsd.common,]
MI.integrated.se <- MI.integrated.se[mrsd.common,]
MI.integrated.se <- logNormCounts(MI.integrated.se)
MI.mrsd.pred <- SingleR(test = MI.integrated.se, ref = mrsd.se, method = "single", labels = mrsd.se$label.main)
sce[["single_R"]] <- MI.mrsd.pred$labels
cell.data <- table(sce@meta.data$harmony_clusters,sce@meta.data$single_R)

## Manual annotation based on marker expression
Macrophage_genes <- c("Adgre1", "Cd68", "Csf1r")
Neutrophil_genes <- c("Csf3r", "S100a9", "S100a8")
B_genes <- c("Ms4a1", "Cd79a", "Ly6d")
Monocyte_genes <- c("Ly6c2", "Chil3", "F10")
NK_genes <- c("Nkg7", "Klrb1c", "Gzma")
CD209DC_genes <- c("Cd209a", "Klrd1", "Flt3")
T_genes <- c("Cd3e", "Cd3d", "Lef1")
XcrDC_genes <- c("Xcr1", "Ifi205", "Itgae")
MigratoryDC_genes <- c("Ccr7", "Fscn1", "Cacnb3")
ILC2_genes <- c("Rora", "Cxcr6", "Gata3")
PC_genes <- c("Jchain", "Iglv1", "Mzb1")
Mast_genes <- c("Cma1", "Kit", "Rab27b")
Idents(sce) <- "orig.ident"
custom_order <- c("Steady-state","Day1","Day3","Day5","Day7")
sce@meta.data$orig.ident <- factor(sce@meta.data$orig.ident, levels = custom_order)
DimPlot(sce, reduction = "umap.harmony", group.by = c("orig.ident","harmony_clusters"),label = TRUE)
DimPlot(sce, reduction = "tsne.harmony", group.by = c("orig.ident","harmony_clusters"),label = TRUE)

knowmarkers <- factor(c(
  "Adgre1", "Cd68", "Csf1r", # Macrophage markers
  "Csf3r", "S100a9", "S100a8", # Neutrophil markers
  "Ms4a1", "Cd79a", "Ly6d", # B cell markers
  "Ly6c2", "Chil3", "F10", # Monocyte markers
  "Nkg7", "Klrb1c", "Gzma", # NK cell markers
  "Cd209a", "Klrd1", "Flt3", # CD209DC markers
  "Cd3e", "Cd3d", "Lef1",  # T cell markers
  "Xcr1", "Ifi205", "Itgae",  # XcrDC markers
  "Ccr7", "Fscn1", "Cacnb3", # MigratoryDC markers
  "Rora", "Cxcr6", "Gata3", # ILC2 markers
  "Jchain", "Iglv1", "Mzb1", # PC markers
  "Cma1", "Kit", "Rab27b" # Mast cell markers
))
Idents(sce) <- "harmony_clusters"
print(knowmarkers)
DotPlot(object = sce, features =knowmarkers)+RotatedAxis()+scale_x_discrete("")+scale_y_discrete("")
FeaturePlot(sce, features = knowmarkers,reduction = "umap.harmony")
VlnPlot(sce,features = Macrophage_genes,group.by = c("seurat_clusters"),pt.size = 0)
table(sce$orig.ident,Idents(sce))
anno_human <- c("Macrophage","Macrophage", "Macrophage", "Macrophage", "Neutrophil", 
                "Macrophage", "Macrophage", "B", "Macrophage", "B", 
                "Monocyte", "Neutrophil", "Monocyte", "Macrophage", 
                "Neutrophil", "Macrophage", "Macrophage", "Macrophage", 
                "Macrophage", "Neutrophil", "Monocyte", "CD209DC", "T", 
                "CD209DC", "NK", "Macrophage", "T", "Monocyte", "Macrophage", 
                "Macrophage", "Macrophage", "Macrophage", "NK", "XcrDC", 
                "MigratoryDC", "MigratoryDC", "MigratoryDC", "ILC2", "B", 
                "Macrophage", "Macrophage", "B", "PC", "Neutrophil")
names(anno_human) <- levels(sce)
# Update cluster identities (Idents)
sce <- RenameIdents(sce, anno_human)
# Append CellType information to metadata
sce$anno_human <- Idents(sce)
head(sce@meta.data)
DimPlot(sce, reduction = "tsne.harmony", group.by = c("seurat_clusters","orig.ident","anno_human"),label = T)
DimPlot(sce, reduction = "umap.harmony", group.by = c("seurat_clusters","orig.ident","anno_human"),label = T)


# Rename specific sub-clusters
sce$anno_human[match(colnames(sce[,Idents(sce)== "12"]),colnames(sce))] =  "Macrophage"

# Filter out mixed/contaminant clusters (12, 27, 34, 36, 43)
Idents(sce) <- "harmony_clusters"
sce  <- subset(sce,subset= harmony_clusters != "12")
sce  <- subset(sce,subset= harmony_clusters != "27")
sce  <- subset(sce,subset= harmony_clusters != "34")
sce  <- subset(sce,subset= harmony_clusters != "36")
sce  <- subset(sce,subset= harmony_clusters != "43")
Idents(sce) <- "anno_human"
DimPlot(sce, reduction = "tsne.harmony", group.by = c("seurat_clusters","orig.ident","anno_human"),label = T)
DimPlot(sce, reduction = "umap.harmony", group.by = c("seurat_clusters","orig.ident","anno_human"),label = T)
DimPlot(sce, reduction = "umap.harmony", group.by = c("anno_human"),label = T)
DimPlot(sce, reduction = "umap.harmony", split.by =  c("orig.ident"),label = T)
#saveRDS(sce, file = "seurat_obj.rds")
sce <- readRDS("../GSE163129_RAW/seurat_obj.rds")


 # Visualize shifts in cell proportions across timepoints
Idents(sce) <-"anno_human"
table(sce$orig.ident)
table(Idents(sce),sce$orig.ident)
cell_data <- as.data.frame(prop.table(table(Idents(sce),sce$orig.ident),margin = 2))
colnames(cell_data) <- c("细胞类型", "时间", "细胞丰度")
# Sort by chronological order
cell_data <- arrange(cell_data, rev(时间))
# Define custom cell type order
cell_data <- mutate(cell_data, 细胞类型 = factor(细胞类型, levels = c("Macrophage",
                                                              "Neutrophil","B",
                                                              "Monocyte","NK",
                                                              "CD209DC","T",
                                                              "XcrDC",
                                                              "MigratoryDC",
                                                              "ILC2",
                                                              "PC")))

## Optional: sort by cell abundance
## cell_data <- mutate(cell_data, 细胞类型 = reorder(细胞类型, 细胞丰度))
# Create stacked bar charts
plot <- ggplot(cell_data, aes(x = 细胞丰度, y = fct_rev(时间), #fill = 细胞类型, 
                              fill = fct_rev(细胞类型),
)) +
  geom_col() +
  labs(x = "细胞丰度", y = "时间", fill = "细胞类型") +
  theme_minimal()
# Define custom colors for cell types
color_palette <-  c("Neutrophil" = "#FF9900","Macrophage" = "#FF6666",
                    "B" = "#CC9933","Monocyte" = "#99CC33",
                    "NK"="#99CCFF","CD209DC"="#33CC33",
                    "T"="#99CCCC","XcrDC"="#4169E1",
                    "MigratoryDC"="#9966CC","ILC2"="#FF99CC",
                    "PC"="#CC3399")
plot <- plot + scale_fill_manual(values = color_palette)
plot <- plot + guides(colour = guide_legend(reverse = TRUE))
# Adjust axis labels and legend formatting for better readability of temporal data
plot <- plot +theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    axis.text = element_text(size = 12, color = "black"),
                    axis.title = element_text(size = 14, color = "black"),
                    legend.title= element_text(size = 16, color = "black"),
                    legend.text= element_text(size = 14, color = "black"))
# Render the plot with border styling
plot <- plot + geom_col(color = "black", size = 0.5)
print(plot)

### Marker Gene DotPlot
Idents(sce) <- "anno_human"
cell_order<-c("Macrophage",
              "Neutrophil","B",
              "Monocyte","NK",
              "CD209DC","T",
              "XcrDC",
              "MigratoryDC",
              "ILC2",
              "PC")
sce@meta.data$anno_human <- factor(sce@meta.data$anno_human, levels = cell_order)
DotPlot(object = sce, features =as_factor(knowmarkers[1:33]),dot.scale = 6)+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+ theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_color_continuous(low="#498EA4",high =  "#E54924")+ 
  theme(legend.position="right")  + 
  labs(title = "Cluster markers", y = "", x="")

FeaturePlot(sce, features = c("Adgre1","S100a9","Cd79a","Ly6c2","Nkg7","Cd209a","Cd3e","Xcr1","Fscn1"),cols = c("gray", "red"), pt.size=0.1)

FeaturePlot(sce, features = c("S100a9"),label = TRUE,cols = c("gray", "red"), pt.size=0.1)


markers_sce <- FindAllMarkers(sce, only.pos = TRUE, min.pct =0.25, logfc.threshold = 0.25)







#### Subset Extraction ####
Idents(sce) <-"anno_human"
levels(Idents(sce))
wanted_clusters <- c("Macrophage","Monocyte")
wanted_clusters <- c("Neutrophil")
subsce <- subset(sce,subset = anno_human %in% wanted_clusters)
## Filter ribosomal genes if necessary
rb.genes <- rownames(subsce)[grep("^Rp[sl]",rownames(subsce))]
Ct <- GetAssayData(object = subsce, layer = "counts")
percent.ribo <- Matrix::colSums(Ct[rb.genes,])/Matrix::colSums(Ct)*100
subsce <- AddMetaData(subsce, percent.ribo, col.name = "percent.ribo")
VlnPlot(subsce, features = "percent.ribo", pt.size = 0.1 ) + NoLegend()
subsce <- subset(subsce,subset = percent.ribo<10)
Idents(subsce) <-"orig.ident"
Idents(subsce)
# Re-process sub-lineages using "harmony"
library(harmony)
subsce <- RunHarmony(subsce,c("orig.ident"))
subsce <- FindVariableFeatures(subsce)
subsce <- ScaleData(subsce,features = rownames(subsce))
subsce <- RunPCA(subsce)
subsce <- FindNeighbors(subsce, dims = 1:30, features = VariableFeatures(object = subsce))
subsce <- FindClusters(subsce, resolution = 1.2) # Macrophages
subsce <- FindClusters(subsce, resolution = 0.25) # Granulocytes
subsce <- RunUMAP(subsce, dims = 1:30, reduction = "pca")
subsce <- RunUMAP(subsce, dims = 1:15, reduction = "pca")
DimPlot(subsce, reduction = "umap", group.by = c("anno_human", "seurat_clusters"),label = T)
DimPlot(subsce, reduction = "umap", group.by = c("seurat_clusters"),label = T)
DimPlot(subsce, reduction = "umap", split.by = c("orig.ident"),label = T)
clustree(subsce ,prefix = "RNA_snn_res.")
markers_subsce <- FindAllMarkers(subsce, only.pos = TRUE, min.pct =0.25, logfc.threshold = 0.25)
write.csv(markers_subsce,file = "test_Neutrophil_markers.csv")
table(subsce$orig.ident,Idents(subsce))

DotPlot(object = subsce, features =knowmarkers,scale = FALSE)+RotatedAxis()+scale_x_discrete("")+scale_y_discrete("")
### Sub-lineage Annotation ####
# Monocyte/Macrophage populations
knownmarkers <-c( "Ly6c2","Plac8",  # Monocytes
                  "Arg1","Il1b","Ccl6","F10", # Pro_inflammation
                  "Ifit1","Ifit3","Isg15", "Rsad2",  # IFN response
                  "Stmn1","Cenpa","Pclaf", # Proliferation
                  "Pf4", "Fn1","Slc7a2","Sdc3", # Angiogenesis
                  "Trem2","Spp1","Gpnmb","Gdf15","Fabp5","Apoe","Ms4a7", # Reparative
                  "Ccr2","H2-Aa","H2-Eb1","H2-DMa",  # Antigen-presenting
                  "Col5a2","Col1a1","Col3a1","Dcn", # Collagens
                  "Mgl2","Folr2","Adgre1" # RCM
                  
)
## Distinguish Monocyte, M1, and M2 subtypes
## Monocytes (Clusters 5, 7, 12, 13, 18)
VlnPlot(subsce,features =c("Ly6c2","Plac8","Ccr2"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Ly6c2","Plac8","Ccr2"),reduction = "umap" )
## M1 Subtype (iNOS, Il1b, Il6, Tnfa)
VlnPlot(subsce,features =c("Nos2","Il1b","Il6","Il23","Mmp8","Mmp9","Tnf"),group.by = c("seurat_clusters"),pt.size = 0)

## M2 Subtype (Arg1, Tgfb, Il10, Vegfa/b/c/d)
VlnPlot(subsce,features =c("Arg1","Tgfb1","Tgfb2","Tgfb3","Il10","Il12","Vegfa"),group.by = c("seurat_clusters"),pt.size = 0)



## Leukocyte chemotaxis & Interferon-stimulated populations (Clusters 7, 8)
VlnPlot(subsce,features =c("S100a9","S100a8","Ccl7","Rsad2","Isg15","Il18","Irf7","Cxcl10"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Ifit1","Ifit3","Isg15", "Rsad2"),reduction = "umap" )
## Anti-inflammatory potential (Clusters 1, 2, 3, 4)
VlnPlot(subsce,features =c("Arg1","Trem2","Spp1","Gpnmb","Ctsd","Cd63","Tgfb1","Tgfb2","Il10","Il12","Hmox1","Fabp4","Fabp5"),group.by = c("seurat_clusters"),pt.size = 0)
## Anti-inflammatory potential (Cluster 20)
VlnPlot(subsce,features =c("Arg1","Trem2","Prdx1","Tgfb1","Tgfb2","Tgfb3","Il10","Il12","Gpnmb","Hmox1","Fabp4","Fabp5"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Trem2","Spp1","Gpnmb","Gdf15","Fabp5","Apoe","Ms4a7"),reduction = "umap" )
## Monocyte-derived antigen presentation (Clusters 0, 7, 10; cluster 7 optional)
VlnPlot(subsce,features =c("Ccr2","H2-Eb1","H2-DMa","Il1b","Tnfsf9","Tnip3"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Ccr2","H2-Aa","H2-Eb1","H2-DMa"),reduction = "umap" )

## Pro-angiogenic (Clusters 6, 22)
VlnPlot(subsce,features =c("Pf4", "Fn1","Slc7a2","Sdc3","Pdgfra","Hif1a","Vegfd","Vegfa","Vegfc","Kdr"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c( "Pf4", "Fn1","Saa3","Fcna"),reduction = "umap" )
## Proliferative phenotypes (Clusters 14, 11, 21)
VlnPlot(subsce,features =c( "Stmn1","Birc5","Top2a","Ube2c" ),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Stmn1","Cenpa","Pclaf"),reduction = "umap" )


## Tissue-resident macrophages (Cluster 8)
VlnPlot(subsce,features =c("Cd14","Mmp12","Cxc3r1","Il1b","Ccr2","H2-Eb1"),group.by = c("seurat_clusters"),pt.size = 0)
## Tissue-resident macrophages (Cluster 8)
VlnPlot(subsce,features =c( "Timd4","Lyve1","Folr2","Mrc1","Cd163","Ccr2","Igf1"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Mgl2","Folr2","Adgre1"),reduction = "umap" )

## Clusters 15/9 may represent polarization transition states
FeaturePlot(subsce, features = c("Arg1","Il1b","Ccl6","F10"),reduction = "umap" )
## Cluster 17 represents the collagen-secreting sub-population
FeaturePlot(subsce, features = c("Col5a2","Col1a1","Col3a1","Dcn"),reduction = "umap" )

VlnPlot(subsce, features = c("Col5a2","Col1a1","Col3a1","Dcn"),group.by = c("seurat_clusters"),pt.size = 0)

# Supplementary Figure 2 visual validation
VlnPlot(subsce, features = c("Lyve1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Cd163"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Timd4"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
FeaturePlot(subsce, features = c("Lyve1","Cd163","Timd4"),reduction = "umap" )



VlnPlot(subsce, features = c("H2-Eb1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("H2-Ab1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Cd74"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
FeaturePlot(subsce, features = c("H2-Eb1","H2-Ab1","Cd74"),reduction = "umap" )



VlnPlot(subsce, features = c("Ccr2"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Ly6c2"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Plac8"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
FeaturePlot(subsce, features = c("Ccr2","Ly6c2","Plac8"),reduction = "umap" )

VlnPlot(subsce, features = c("Cd68"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Fcgr1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Itgam"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
FeaturePlot(subsce, features = c("Cd68","Fcgr1","Itgam"),reduction = "umap" )

VlnPlot(subsce, features = c("Ace"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Ear2"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Itgal"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
FeaturePlot(subsce, features = c("Ace","Ear2","Itgal"),reduction = "umap" )

FeaturePlot(subsce, features = c("S100a8","S100a9","Plac8"),reduction = "umap" )

VlnPlot(subsce, features = c("Fcrls"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Rgs10"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Adgre1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
FeaturePlot(subsce, features = c("Fcrls","Rgs10","Adgre1"),reduction = "umap" )




VlnPlot(subsce, features = c("Trem2"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Gpnmb"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
VlnPlot(subsce, features = c("Spp1"), pt.size = 0, assay = "RNA", cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7", "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1")) + NoLegend() + geom_boxplot(width=0.1, fill="white")
FeaturePlot(subsce, features = c("Trem2","Gpnmb","Spp1"),reduction = "umap" )




inflammation_human  <- c("Antigen_presenting", "Reparative", "RCM", "Reparative", "RCM", "Pro_inflammation", 
  "IFN", "Neut_like", "Pro_inflammation", "Proliferation", "ECM_like", 
  "Pro_inflammation", "Proliferation", "Reparative", "ECM_like", "Pro_inflammation", 
  "RCM", "Reparative", "Proliferation", "Proliferation", "B_like", "Reparative", 
  "Proliferation")
names(inflammation_human) <- levels(subsce)
# Update cluster identities (RenameIdents)
subsce <- RenameIdents(subsce, inflammation_human)
# Append CellType information to metadata
subsce$inflammation_human <- Idents(subsce)

head(subsce@meta.data)
DimPlot(subsce, reduction = "umap", group.by = c("seurat_clusters","orig.ident","inflammation_human"),label = T)
DimPlot(subsce, reduction = "umap", group.by = c("seurat_clusters","inflammation_human"),label = T)
DimPlot(subsce, reduction = "umap", split.by = c("orig.ident"),label = T)

#saveRDS(subsce, file = "seurat_subsce.rds")
subsce <- readRDS("./seurat_subsce.rds")





######### Immune Compartment: Target Lineage Extraction #################################
# Based on Monocle 3 trajectory results and prior dataset insights, exclude proliferative/antigen-presenting clusters and their immediate progenitors
Idents(temp_sce) <- "seurat_clusters"
Idents(temp_sce) <- "inflammation_human"
wanted_clusters <- c("8","11","15","5","3","21","17","3","13","1")
temp_sce <- subset(subsce,subset = seurat_clusters %in% wanted_clusters)
#temp_sce <- subset(temp_sce,subset = seurat_clusters != c("Proliferation"))
temp_sce$"label" <- Idents(temp_sce)
library(harmony)
## Filter ribosomal genes if necessary
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
temp_sce <- FindClusters(temp_sce, resolution = 1)
temp_sce <- RunUMAP(temp_sce, dims = 1:30, reduction = "pca")

DimPlot(temp_sce, reduction = "umap", group.by = c("anno_human", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("inflammation_human", "seurat_clusters","label"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("inflammation_human", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", split.by = c("orig.ident"),label = T)
FeaturePlot(temp_sce, features = c("Ly6c2","Ccr2","Arg1","Trem2"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Angptl4","Ccnb2","Acta2","Sfrp2"),reduction = "umap" )


table(Idents(temp_sce),temp_sce$orig.ident)
cell_data <-as.data.frame(prop.table(table(Idents(temp_sce),temp_sce$orig.ident),margin = 2))
cell_data <-as.data.frame(table(Idents(temp_sce),temp_sce$orig.ident),margin = 2)

table(Idents(subsce),subsce$orig.ident)
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

markers_temp_sce <- FindAllMarkers(temp_sce, only.pos = TRUE, min.pct =0.25, logfc.threshold = 0.25)
write.csv(markers_temp_sce,file = "test_macrophages_markers.csv")



### Targeted Sub-lineage Extraction (Optional) ####
Idents(temp_sce) <- "seurat_clusters"
wanted_clusters <- c("1","4","9") # DNB candidates for FIBROBLAST (Cluster 4)
wanted_clusters <- c("9","6","7") # DNB candidates for FIBROBLAST (Cluster 6)

wanted_clusters <- c("0","1","3") # DNB candidates for Granulocyte (Cluster 1)
wanted_clusters <- c("2","0","1") # DNB candidates for Granulocyte (Cluster 0)


wanted_clusters <- c("5","4","7") # DNB candidates for Monocyte/Macrophage (Cluster 4)
wanted_clusters <- c("1","3","0") # DNB candidates for Monocyte/Macrophage (Cluster 3)
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
temp_sce2 <- RunTSNE(temp_sce2, dims = 1:30, reduction = "harmony")
temp_sce2 <- RunUMAP(temp_sce2, dims = 1:30, reduction = "harmony")
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
write.csv(temp_sce2.markers,file = "fibroblast_967_DEGs.csv")
write.csv(temp_sce2.markers,file = "temp_sce2markers.csv")

temp_sce2.markers <- FindAllMarkers(temp_sce, only.pos = FALSE, min.pct =0, logfc.threshold = 0)
write.csv(temp_sce2.markers,file = "all_unfiltered_markers.csv")





##################### dyno Multi-method Pseudotime Analysis #############################
devtools::install_github("dynverse/dyno")
devtools::install_github("dynverse/dynmethods")
library(dyno)
library(tidyverse)
library(Matrix)
library(Seurat)
# Add raw counts and normalised expression
# Seurat matrices require transposition so that rows are cells and columns are genes
dataset <- wrap_expression(
  counts = t(temp_sce[["RNA"]]$counts),
  expression = t(temp_sce[["RNA"]]$data)
)

dataset <- wrap_expression(
  counts = t(subsce[["RNA"]]$counts),
  expression = t(subsce[["RNA"]]$data)
)


# Add prior information; the starting "cell id" is added here. 
# Visualizations can be adjusted later based on specific trajectory inference results.
dataset <- add_prior_information(
  dataset,
  start_id = "0day_AAACCTGAGTTATCGC-1"
)  # Fibroblast sub-population
dataset <- add_prior_information(
  dataset,
  start_id = "Steady-state_AAAGGTATCTAGCATG-1"
) # Monocyte/Macrophage sub-population
dataset <- add_prior_information(
  dataset,
  start_id = "0day_ACGATACAGGGTCGAT-1"
) # Neutrophil
# Add cluster information; "seurat_clusters" is used directly here
dataset <- add_grouping(
  dataset,
  temp_sce$seurat_clusters
)
dataset <- add_grouping(
  dataset,
  subsce$seurat_clusters
)
## dynwrap::test_docker_installation(detailed = TRUE) Check if Docker can be called
## Select the best parameters and methods for the dataset; it is recommended to try several methods
guidelines <- guidelines_shiny(dataset)
methods_selected <- guidelines$methods_selected
## paga_tree is selected here
model_paga <- infer_trajectory(dataset, methods_selected[1])
model_paga <- model_paga %>% add_dimred(dyndimred::dimred_mds, expression_source = dataset$expression)
## Preliminary plotting with many adjustable parameters
plot_dimred(
  model_paga, 
  expression_source = dataset$expression, 
  grouping = dataset$grouping
)
system()
## Add rooting gene list to define the biological differentiation direction
model <- model_paga %>% 
  add_root_using_expression(c("Angptl4","Ccnb2","Acta2","Sfrp2"), dataset$expression)

model <- model_paga %>% 
  add_root_using_expression(c("Ly6c2","Il1b","Stmn1","Trem2","Ccr2","Arg1"), dataset$expression)
### Mark milestones with marker genes
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
    Pro_inflammation = c("Ccr2"),
    #Proliferation = c("Stmn1"),
    Reparative = c("Trem2"),
    Collagens =c("Arg1")
    
  ),
  dataset$expression
)
## Visualization
model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = dataset$expression)
plot_dimred(
  model, 
  color_cells = "pseudotime",
  #color_cells = "grouping",
  expression_source = dataset$expression, 
  grouping = dataset$grouping,
  label_milestones = TRUE
)


### Display gene expression

plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "Col3a1",
  #color_density = "grouping",
  grouping = dataset$grouping,
  label_milestones = TRUE
)

plot_dimred(
  model, 
  expression_source = dataset$expression, 
  #color_cells = "grouping",
  #feature_oi = "Ly6c2",
  color_density = "grouping",
  grouping = dataset$grouping,
  label_milestones = TRUE
)


# Global overview of the most predictive genes
plot_heatmap(
  model,
  expression_source = dataset$expression,
  grouping = dataset$grouping,
  features_oi = 50
)


## Important genes at branching points

branching_milestone <- model$milestone_network %>% group_by(from) %>% filter(n() > 1) %>% pull(from) %>% first()

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







################### BioTIP Analysis ################################################
# BioTIP Main Functions Walk Through
# samplesL
subsce_trans <- as.SingleCellExperiment(temp_sce)
subsce_trans <- as.SingleCellExperiment(subsce)
samplesL <- split(rownames(colData(subsce_trans)),f = colData(subsce_trans)$seurat_clusters)
# Transcript Pre-selection
dec.pois <- modelGeneVarByPoisson(subsce_trans)
hvg <- getTopHVGs(dec.pois, n=4000)
hvg <- intersect(hvg, rownames(subsce_trans))
dat <- subsce_trans[hvg,]
logmat <- as.matrix(logcounts(dat))
global <- as.matrix(logcounts(subsce_trans)) 
## Draw scatter plot
gene_mean <- apply(global, 1, mean)
gene_sd <- apply(global, 1, sd)
blue_gene_indices <- which(names(gene_mean) %in% hvg)
red_gene_indices <- which(names(gene_mean) %in% cluster_specificgenes)
red_gene_indices <- which(gene_sd >1)
points(gene_mean[red_gene_indices], gene_sd[red_gene_indices], col = "red", pch = 10)
points(gene_mean[blue_gene_indices], gene_sd[blue_gene_indices], col = "lightblue", pch = 10)
# Add legend
legend("topright", legend = c("Red Genes"), col = "pink", pch = 16)
plot(gene_mean, gene_sd, xlab = "Mean Expression", ylab = "Standard Deviation")

ggplot(data, aes(x = gene_mean, y = gene_sd)) +
  geom_point(size = data$size)

# Draw smooth scatter plot
smoothScatter(x = gene_mean, y = gene_sd, nbin = 128,
              colramp = colorRampPalette(c("white", blues9)),
              nrpoints = 100, ret.selection = FALSE,
              transformation = function(x) x^0.25,
              postPlotHook = box,
              xlab = "Mean Expression", ylab = "Standard Deviation",
              xaxs = par("xaxs"), yaxs = par("yaxs"))




# Select highly variable genes
cut.preselect = 0.10
cut.preselect = 0.15
cut.preselect = 0.20
# cut.preselect = 0.01 threshold increased as it runs too slowly
set.seed(2020)
pb <- txtProgressBar(min = 0, max = 100, style = 3)
testres <- optimize.sd_selection(logmat[,unlist(samplesL)], samplesL, B=100, cutoff=cut.preselect, times=.75, percent=0.8)
#save(testres, file=paste0("Fibroblasts_biotip.RData"), compress=TRUE)
#save(testres, file=paste0("macrophages_biotip.RData"), compress=TRUE)
#save(testres, file=paste0("test_Granulocytes_biotip.RData"), compress=TRUE)
#save(testres, file=paste0("test_macrophages_biotip.RData"), compress=TRUE)
testres <- load("test_Granulocytes_biotip.RData")
testres <- load
testres_list<- list()
#### Test
nested_list<- list()
#### Test
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
cut.fdr = 0.10
cut.fdr = 0.12
cut.fdr = 0.14
cut.fdr = 0.20
igraphL <- getNetwork(testres, fdr = cut.fdr)
cluster <- getCluster_methods(igraphL)

## Network partition using random walk
dev.off()
par(mfrow=c(3,3))
cluster = getCluster_methods(igraphL)
#i = 1
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

plot(tmp,vertex.label = V(tmp)$name,vertex.color=colrs[V(tmp)$community],vertex.size =5,mark.groups=cluster[[i]])
nodes_in_community <- V(tmp)$name[V(tmp)$community == 10]
# Create a subgraph containing only target community nodes
subgraph <- induced_subgraph(tmp, nodes_in_community)
plot(subgraph,vertex.label = V(subgraph)$name,vertex.color=colrs[V(subgraph)$community],vertex.size = 6)
# Calculate node degree
node_degrees <- degree(subgraph)
# Determine minimum and maximum node degrees
min_degree <- min(node_degrees)
max_degree <- max(node_degrees)
# Map node degree to cold-warm color spectrum
color_palette <- colorRampPalette(c("#8FD2E6", "white","#ED884C" ))
color_palette <- colorRampPalette(c( "white","#ED884C" ))
node_colors <- color_palette(100)[as.integer(100 * (node_degrees - min_degree) / (max_degree - min_degree)) + 1]

# Map expression values to cold-warm color spectrum
local_sce <- subset(temp_sce ,subset = seurat_clusters == "0")
#local_sce <- subset(subsce ,subset = seurat_clusters == "0")
local_sce <- as.SingleCellExperiment(local_sce)
local_sce <- as.matrix(logcounts(local_sce)) 
gene_sd <- apply(local_sce, 1, sd)
gene_sd <- gene_sd[V(subgraph)$name] 
gene_exp <- apply(local_sce, 1, mean)
gene_exp <- gene_exp[V(subgraph)$name] 
# Determine minimum and maximum SD values
min_exp <- min(gene_exp)
max_exp <- max(gene_exp)
node_size <- gene_sd*10
node_colors <- color_palette(100)[as.integer(100 * (gene_exp - min_exp) / (max_exp - min_exp)) ]
# Use label propagation algorithm for layout
layout <-  layout_with_kk(subgraph)
# Draw the plot, adjusting node colors according to node degree
plot(subgraph,
     vertex.label = V(subgraph)$name,
     vertex.color = node_colors,
     vertex.size = node_size*1.5,
     #vertex.label.dist = 1.5,
     vertex.label.cex = 0.66,
     layout = layout,
     vertex.label.color = "black",
     edge.alpha = 0.25)

plot(subgraph,
     vertex.label = V(subgraph)$name,
     vertex.color = node_colors,
     vertex.size = node_size*1.1,
     vertex.label.cex = 0.5,
     layout = layout,
     vertex.label.color = "black",
     edge.color = "gray50",
     edge.alpha = 0.5)

# Save edge data
edges <- as_edgelist(subgraph)
# Fibroblast cluster
write.table(edges, "7edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "7gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
write.table(edges, "4edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "4gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
# Fibroblast cluster - new 4 and 6 
write.table(edges, "fibroblast_4edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "fibroblast_4gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
write.table(edges, "fibroblast_6edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "fibroblast_6gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)

# Neutrophil cluster
write.table(edges, "Granulocytes_1edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "Granulocytes_1gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
write.table(edges, "Granulocytes_0edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "Granulocytes_0gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)

# Monocyte/Macrophage cluster
write.table(edges, "macrophages_4edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "macrophages_4gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
# Monocyte/Macrophage cluster
write.table(edges, "macrophages_3edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "macrophages_3gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)




# Map PCC correlation to edge thickness

# Identify Putative Critical Transition Signals (CTSs) by MCI score
membersL <- getMCI(cluster, testres, adjust.size = FALSE, fun='BioTIP')
# Plot DNB scores for modules in each cell cluster
dev.off()
cut.minsize = 60
pdf('./验证图/Figure2-3.pdf',height = 12,width = 20,onefile = F)
par(oma = c(0, 0, 0, 0))  # Adjust outer margins to 0
par(mar = c(1, 1, 1, 1))  # Adjust inner margins to 1
plotBar_MCI(membersL, ylim=c(0,3), minsize = cut.minsize)
dev.off()
# Identify the highest-scoring modules
topMCI = getTopMCI(membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]], min=cut.minsize, n=3)
maxMCIms <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], min=cut.minsize, n=3)
maxMCI = getMaxStats(membersL[['MCI']], maxMCIms[['idx']])
CTS.Lib = getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])
CTS.Lib.Symbol <- CTS.Lib
maxMCI <- sort(maxMCI, decreasing = TRUE)
maxMCI <- head(maxMCI, 5)
# Output DNB
# Determine the maximum vector length
max_length <- max(lengths(CTS.Lib))
df <- data.frame(matrix(NA, nrow = max_length, ncol = length(CTS.Lib)))
for (i in 1:length(CTS.Lib)) {
  col_name <- names(CTS.Lib)[i]
  col_data <- CTS.Lib[[i]]
  df[, i] <- c(col_data, rep(NA, max_length - length(col_data)))
}
colnames(df) <- names(CTS.Lib)
write.csv(df, file = "fibrobalst_DNB.csv", row.names = TRUE)
write.csv(df, file = "./验证图/Granulocytes_DNB2.0.csv", row.names = TRUE)
write.csv(df, file = "test_macropahges_DNB.csv", row.names = TRUE)
write.csv(df, file = "test_Granulocytes_DNB.csv.csv", row.names = TRUE)
names(CTS.Lib.Symbol)
# Shrink the matrix and calculate IC.shrink.score
M <- cor.shrink(logmat[,unlist(samplesL)], Y = NULL, MARGIN = 1, shrink = TRUE)
# Analyze if the DNB score is significant
C = 200
simuMCI = list()
set.seed(2020)
for (i in 1:length(CTS.Lib)){
  n <- length(CTS.Lib[[i]])
  simuMCI[[i]] <- simulationMCI(n, samplesL, logmat,  B=C, fun="BioTIP", M=M)
}
dev.off()
pdf('./验证图/Figure2-4.pdf',height = 16,width = 12,onefile = F)
par(mfrow=c(5,1))
for (i in 1:length(CTS.Lib)){
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las=2,
                      main=paste0("Cluster ", names(maxMCI)[i], "; ",
                                  length(CTS.Lib[[i]]), " genes", "\n","vs. ",
                                  "100 times of gene-permutation"),
                      which2point=names(maxMCI)[i])
}

dev.off()
par(mfrow=c(3,1))
for (i in 1:3){
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las=2,
                      main=paste0("Cluster ", names(maxMCI)[i], "; ",
                                  length(CTS.Lib[[i]]), " genes", "\n","vs. ",
                                  "100 times of gene-permutation"),
                      which2point=names(maxMCI)[i])
}

# Confirm candidate tipping points using Ic index and Delta-Ic
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
par(mfrow=c(3,2))
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
par(mfrow=c(3,2))
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

pdf('./验证图/Figure2-5.pdf',height = 16,width = 12,onefile = F)

ylim = 0.3
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
ylim = 0.3
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





dev.off()












############# monocle3 Pseudotime Analysis #############################
library(devtools)
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(Seurat)
library(leidenbase)
DefaultAssay(temp_sce) <- "RNA"
DefaultAssay(subsce) <- "RNA"
# Prepare data for Monocle3
# expression_matrix: Gene-by-cell expression matrix
expression_matrix = temp_sce[["RNA"]]$counts
expression_matrix = subsce[["RNA"]]$counts
# cell_metadata: Metadata corresponding to the Seurat object
cell_metadata = data.frame(temp_sce@meta.data)
cell_metadata = data.frame(subsce@meta.data)
# gene_metadata: Gene annotation info, gene symbols can be used
gene_annotation = data.frame(expression_matrix[,1])
gene_annotation[,1] = row.names(gene_annotation)
colnames(gene_annotation)=c("gene_short_name")
## Construct Monocle3 CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
# Data preprocessing
cds <- preprocess_cds(cds, num_dim = 50,norm_method = c("log"))

## Dimensionality reduction, default is "Umap"
cds <- reduce_dimension(cds,reduction_method="UMAP",cores=5)
plot_cells(cds,color_cells_by = "seurat_clusters")
## Cluster cells; low resolution ensures cells are grouped as one for better visualization
cds <- cluster_cells(cds,resolution = 0.0000001)
cds <- cluster_cells(cds,cluster_method = "louvain", reduction_method ="UMAP")
## Pseudotime
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster=FALSE, label_leaves=FALSE,label_branch_points=FALSE)

## Select specific cells as the starting point
## Here we assume HEpiD cells are the starting point
myselect <- function(cds,select.classify,my_select){
  cell_ids <- which(colData(cds)[,select.classify] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                             (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes}


cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "8"))
cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "10"))
cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "2"))
## Synchronize Seurat UMAP information to maintain cell distribution consistency
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(temp_sce , reduction = "umap")
int.embed <- Embeddings(subsce , reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

#root_group = colnames(cds)[clusters(cds) == 1]
#cds = order_cells(cds, root_cells = root_group)

## Pseudotime values for different cell types
## Higher values represent higher differentiation; this is for demonstration, not actual differentiation
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
                                                 graph_label_size=1)+ scale_color_manual(values = pal_jco("default", alpha = 0.6)(23))
# Display differential genes
Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=6)
# Select TOP genes by Moran's I
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by = "seurat_clusters",min_expr=0.5, ncol= 2,cell_size=1.5) + scale_color_manual(values = pal_jco("default", alpha = 0.6)(11)) 
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by = "pseudotime",min_expr=0.5, ncol= 2,cell_size=1.5) + scale_color_manual(values = pal_jco("default", alpha = 0.6)(11)) 
# Select specific genes for visualization
plot_cells(cds, genes=c("Mt2","Mt1","Cxcl3"),show_trajectory_graph=FALSE)#IR
pData(cds) 
# Pseudotime scatter plot
plot_genes_in_pseudotime(cds[c("Mt2","Mt1","Cxcl3"),],
                         color_cells_by="seurat_clusters",
                         min_expr=0.5, ncol= 2,cell_size=1.5)+ scale_color_manual(values = pal_jco("default", alpha = 0.6)(7)) 

############################ monocle2 Pseudotime #####################################
# Convert Seurat object counts to a sparse matrix
library(monocle)
Idents(temp_sce2) <- "label"
sparse_data <-  as(as.matrix(temp_sce2[["RNA"]]$counts),"sparseMatrix") 
# Create the original matrix containing cell data
mdata <- new("AnnotatedDataFrame",data = temp_sce2@meta.data)
# Create a matrix containing gene_short_name
fData <- data.frame(gene_short_name=row.names(sparse_data),row.names =row.names(sparse_data))
fd <- new("AnnotatedDataFrame",data = fData)
# Construct Monocle CDS object
monocle_cds <- newCellDataSet(cellData = sparse_data,
                              phenoData = mdata,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
# Calculate size factors and dispersions
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
# Filter lowly expressed cells
monocle_cds <- detectGenes(monocle_cds,min_expr = 3)

# Identify highly variable genes
# Use variable features selected by Seurat
expressed_genes <- VariableFeatures(temp_sce2)
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],
                                      fullModelFormulaStr = "~label")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) ## qval < 0.01 threshold.
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)

# Dimensionality reduction
monocle_cds <- reduceDimension(monocle_cds,max_components = 2,reduction_method = "DDRTree")
# Order cells in pseudotime
monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds,root_state = "3")
monocle_cds <- orderCells(monocle_cds,root_state = "5")
# Plotting
# Coloring can also be set to 'Pseudotime'
p <- plot_cell_trajectory(monocle_cds,color_by = "Pseudotime",size = 1,show_backbone = TRUE) + facet_wrap("~orig.ident",nrow=1)
p0 <- plot_cell_trajectory(monocle_cds,color_by = "label",size = 1,show_backbone = TRUE) + facet_wrap("~orig.ident",nrow=1)
p1 <- plot_cell_trajectory(monocle_cds,color_by = "label",size = 1,show_backbone = TRUE)
theme(legend.position = 'none',panel.border = element_blank())
#+scale_color_manual(values= color) # custom colors
# Add tree structure plot
p2 <- plot_complex_cell_trajectory(monocle_cds,x=1,y=2,color_by = "label")+
  theme(legend.title = element_blank())
#+scale_color_manual(values= color) # custom colors
plot_cell_trajectory(monocle_cds,color_by = "State",size = 1,show_backbone = TRUE) 
p0|p|p1|p2


## Identify pseudotime-associated genes
Time_diff <- differentialGeneTest(monocle_cds[ordering_genes,],cores =1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
# Filter genes expressed in at least 10 cells
Time_diff <- Time_diff %>%
  filter(num_cells_expressed > 10)
# Filter top 100 significantly varying genes
Time_genes <- top_n(Time_diff,n=100,desc(qval)) %>%  pull(gene_short_name) %>% as.character()
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p4 <- plot_pseudotime_heatmap(monocle_cds[Time_genes,],num_clusters = 4,show_rownames = T,return_heatmap = T)
## Extract pseudotime-associated genes by cluster
clusters <- cutree(p4$tree_row,k=4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_clusters"
table(clustering)




##### Jaccard Similarity Calculation ####
ceshi_markers <- read.csv(file = "../E-MTAB-7895/macrophagesmarkers.csv",row.names = 1)
ceshi_markers <- temp_sce.markers 
ceshi_markers <- read.csv(file = "../E-MTAB-7895/Granulocytesmarkers.csv",row.names = 1)
ceshi_markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> ceshi_top100


valid_markers <- markers_temp_sce
valid_markers <- markers_subsce
valid_markers <- read.csv(file = "../GSE163129_RAW/test_Neutrophil_markers.csv",row.names = 1)
valid_markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> valid_top100

# Extract genes from ceshi_markers with absolute log2FC > 0.5
ceshi_top100 <- ceshi_markers %>%
  filter(abs(avg_log2FC) > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)))

# Extract genes from valid_markers with absolute log2FC > 0.5  
valid_top100 <- valid_markers %>%
  filter(abs(avg_log2FC) > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)))

# Extract unique clusters from each dataframe
ceshi_clusters <- unique(ceshi_top100$cluster)
valid_clusters <- unique(valid_top100$cluster)

# Create Jaccard similarity matrix
jaccard_matrix <- matrix(nrow = length(ceshi_clusters), ncol = length(valid_clusters))
colnames(jaccard_matrix) <- valid_clusters
rownames(jaccard_matrix) <- ceshi_clusters
# Calculate Jaccard similarity
for (i in seq_along(ceshi_clusters)) {
  ceshi_cluster <- ceshi_top100$gene[ceshi_top100$cluster == ceshi_clusters[i]]
  for (j in seq_along(valid_clusters)) {
    valid_cluster <- valid_top100$gene[valid_top100$cluster == valid_clusters[j]]
    intersection <- length(intersect(ceshi_cluster, valid_cluster))
    union <- length(union(ceshi_cluster, valid_cluster))
    jaccard_matrix[i, j] <- intersection / union
    #jaccard_matrix[i, j] <- intersection
  }
}


# Create heatmap
dev.off()
pheatmap(jaccard_matrix,
         main = "Jaccard Similarity Heatmap",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("lightblue", "white", "red"))(100),
         fontsize_row = 20,
         fontsize_col = 20,
         display_numbers = TRUE,
         number_color = "black")












fisher_test <- function(x, y) {
  a <- length(intersect(x, y))
  b <- length(x) - a
  c <- length(y) - a
  d <- max(0, 100 - a - b - c)
  p_value <- tryCatch({
    fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
  }, error = function(e) {
    NA
  })
  p_value
}

# Calculate Fisher's Exact Test p-value matrix
fisher_matrix <- matrix(nrow = length(ceshi_clusters), ncol = length(valid_clusters))
colnames(jaccard_matrix) <- valid_clusters
rownames(jaccard_matrix) <- ceshi_clusters


# Calculate Jaccard similarity
for (i in seq_along(ceshi_clusters)) {
  ceshi_cluster <- ceshi_top100$gene[ceshi_top100$cluster == ceshi_clusters[i]]
  for (j in seq_along(valid_clusters)) {
    valid_cluster <- valid_top100$gene[valid_top100$cluster == valid_clusters[j]]
    fisher_matrix[i, j] <- fisher_test(ceshi_cluster, valid_cluster)
  }
}




















set.seed(123)

# Generate the first set of data
cluster1 <- lapply(1:6, function(i) sample(letters, 10, replace = TRUE))

# Generate the second set of data  
cluster2 <- lapply(1:5, function(i) sample(letters, 10, replace = TRUE))



jaccard_sim <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}


# Calculate Jaccard similarity matrix
jaccard_matrix <- matrix(NA, nrow = 6, ncol = 5)
rownames(jaccard_matrix) <- paste0("Cluster", 1:6)
colnames(jaccard_matrix) <- paste0("Cluster", 1:5)

for (i in 1:6) {
  for (j in 1:5) {
    jaccard_matrix[i, j] <- jaccard_sim(cluster1[[i]], cluster2[[j]])
  }
}

fisher_test <- function(x, y) {
  a <- length(intersect(x, y))
  b <- length(x) - a
  c <- length(y) - a
  d <- max(0, 100 - a - b - c)
  p_value <- tryCatch({
    fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
  }, error = function(e) {
    NA
  })
  p_value
}

# Calculate Fisher's Exact Test p-value matrix
fisher_matrix <- matrix(NA, nrow = 6, ncol = 5)
rownames(fisher_matrix) <- paste0("Cluster", 1:6)
colnames(fisher_matrix) <- paste0("Cluster", 1:5)

for (i in 1:6) {
  for (j in 1:5) {
    fisher_matrix[i, j] <- fisher_test(cluster1[[i]], cluster2[[j]])
  }
}


data <- expand.grid(row = rownames(jaccard_matrix), 
                    col = colnames(jaccard_matrix))
data$jaccard <- as.vector(jaccard_matrix)
data$p_value <- as.vector(fisher_matrix)

library(ggplot2)

ggplot(data, aes(x = col, y = row, fill = jaccard)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p_value)), color = "white", size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

############################### DNB Expression Trend and Clustering Analysis ################################
# Custom expression matrix function
### DNB average expression value function
i=2
h=1
Idents(temp_sce2) <- "label"
DNBs_expression <- function(expression_matrix,target_cluster,DNB_names){
  failed_genes <- c()
  # Remove duplicate transcription factor names
  DNB_names <- unique(na.omit(DNB_names))
  #tryCatch({
  average_expression <- c()
  df <- data.frame(matrix(nrow = length(target_cluster), ncol = length(DNB_names)))
  rownames(df) <- target_cluster
  colnames(df) <- DNB_names
  # Extract cell indices for the target cluster
  for (i in seq_along(target_cluster)){
    cell_indices <- which(temp_sce2$label == target_cluster[i])
    for (h in 1:length(DNB_names)){
      # First, check if the gene exists in the single-cell matrix
      if (DNB_names[h] %in% rownames(expression_matrix)) {
        average_expression <- mean(expression_matrix[DNB_names[h], cell_indices])
        
      }else {
        failed_genes <<- c(failed_genes, DNB_names[h])
        average_expression <- NA # Calculation failed, return NA
      }
      df[i, DNB_names[h]] <- average_expression
    }
    if (length(failed_genes) > 0) {
      
      print(paste("Failed DNBs list:", paste(unique(failed_genes), collapse = ", ")))
      cat("\n")}
  }
  return(df)
  #})
}
DNB_list <- as.list(read.csv("./fibrobalst_DNB.csv"))
DNB_list <- as.list(read.csv("./Granulocytes_DNB.csv"))
DNB_list <- as.list(read.csv("./test_macropahges_DNB.csv"))
DNB_list <- as.list(read.csv("./test_Granulocytes_DNB.csv.csv"))
Idents(temp_sce2) # check the single-cell object, confirm label as current identity
Idents(subsce) 

expression_matrix = temp_sce2[["RNA"]]$scale.data #normalized_scale

target_cluster <- c("1","4","9") # Fibroblast dnb 4 early
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster,unique(DNB_list$X4))
exps <- t(DNBs_expression.df)
write.table(exps, "fibroblast_4exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)





target_cluster <- c("9","6","7") # Fibroblast dnb 7 late
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X6))
exps <- t(DNBs_expression.df)
write.table(exps, "fibroblast_6exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)


target_cluster <- c("0","1","3") # Granulocytes dnb 1 early
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X1))
exps <- t(DNBs_expression.df)
write.table(exps, "Granulocytes_013exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

target_cluster <- c("2","0","1") # Granulocytes dnb 2 early
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X0))
exps <- t(DNBs_expression.df)
write.table(exps, "Granulocytes_201exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)


target_cluster <- c("5","4","7") # macrophages dnb 4 early
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X4))
exps <- t(DNBs_expression.df)
write.table(exps, "macrophages_547exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

target_cluster <- c("1","3","0") # macrophages dnb 8 late
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X3))
exps <- t(DNBs_expression.df)
write.table(exps, "macrophages_130exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)





?VlnPlot()
## Clustervis clustering analysis
#exps <- exps[!duplicated(exps[,1]), ]
#rownames(exps) <- exps[,1]
#exps <- exps[,-1]

exps <- as.matrix(exps)
exps <-na.omit(exps)
head(exps)
# ClusterGVis: Note conflict with Monocle3; if inherited function errors occur, unload Monocle3
library(ClusterGVis)
library(sf)

markers_subsce %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

markers_temp_sce %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
# Use cluster means as input; genes are deduplicated by default
st.data1 <- prepareDataFromscRNA(object = temp_sce,
                                  #object = subsce,
                                  diffData = top10,
                                  showAverage = TRUE)
ck <- clusterData(exp = exps,
                  #cluster.method = "mfuzz",
                  cluster.method = "kmeans",
                  cluster.num = 4)
# add gene name
markGenes = rownames(exps)[sample(1:nrow(exps),32,replace = F)]
pdf('addgene.pdf',height = 10,width = 6,onefile = F)
visCluster(object = ck,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes)
dev.off()
# Enrichment analysis, mouse genome
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

# If genes were added in the previous step, they need to be removed here for plotting
enrich <- subset(enrich, select = -geneID)
head(enrich,3)
pdf('term4.pdf',height = 12,width = 12,onefile = F)
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
  go.size = 10,
  cluster.order = c(1:10),
  #mulGroup = c(22,14,37,61,42,4),
  #mline.col = c(ggsci::pal_lancet()(7))
)
dev.off()
# Trend line display
visCluster(object = ck,
           plot.type = "line",
           ms.col = c("green","orange","red"), # change color
           add.mline = TRUE  # remove median line
)

# Calculate expression trends for each row of genes
trends <- apply(exps, 1, function(row_vec) {
  diff_vec <- diff(row_vec)
  if (all(diff_vec > 0)) {
    return("Gradually increasing")
  } else if (all(diff_vec < 0)) {
    return("Gradually decreasing")
  } else if (diff_vec[1] > 0 && diff_vec[2] < 0) {
    return("Critically high")
  } else if (diff_vec[1] < 0 && diff_vec[2] > 0) {
    return("Critically low")  # Modified to "central high"
  } else {
    return("Unclassified")
  }
})

# Sort clusters based on trends
sorted_idx <- order(trends)

# Reorder expression matrix and trend labels based on sorted index
exps <- exps[sorted_idx, ]
sorted_trends <- as.data.frame(trends[sorted_idx])
colnames(sorted_trends) <- c("Category")

# Load ComplexHeatmap package
library(ComplexHeatmap)
library(circlize)
pdf('heatmap4.pdf',height = 16,width = 20,onefile = F)
col_fun = colorRamp2(c(-2, 0, 2), c("#4169E1", "white", "#DC143C"))
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
dev.off()


############# Monocle3 Pseudotime Analysis #############################
library(devtools)
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(Seurat)
library(leidenbase)
DefaultAssay(temp_sce) <- "RNA"
DefaultAssay(subsce) <- "RNA"

# Construct data required for Monocle3
# expression_matrix: Row represents gene, Column represents cell expression matrix
expression_matrix = temp_sce[["RNA"]]$counts
expression_matrix = subsce[["RNA"]]$counts

# cell_metadata: Corresponds to the Seurat object meta.data
cell_metadata = data.frame(temp_sce@meta.data)
cell_metadata = data.frame(subsce@meta.data)

# gene_metadata: Gene annotation information, gene symbols are applicable
gene_annotation = data.frame(expression_matrix[,1])
gene_annotation[,1] = row.names(gene_annotation)
colnames(gene_annotation)=c("gene_short_name")

## Construct Monocle3 CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# Data preprocessing
cds <- preprocess_cds(cds, num_dim = 50,norm_method = c("log"))

## Dimensionality reduction, defaults to "UMAP"
cds <- reduce_dimension(cds,reduction_method="UMAP",cores=5)
plot_cells(cds,color_cells_by = "seurat_clusters")

## Clustering: Low resolution is used to group cells into a single entity for optimal visualization
cds <- cluster_cells(cds,resolution = 0.0000001)
cds <- cluster_cells(cds,cluster_method = "louvain", reduction_method ="UMAP")

## Pseudotime inference
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster=FALSE, label_leaves=FALSE,label_branch_points=FALSE)

## Select specific cells as the starting point
## Here we assume HEpiD cells represent the developmental origin
myselect <- function(cds,select.classify,my_select){
  cell_ids <- which(colData(cds)[,select.classify] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                             (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes}


cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "8"))
cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "10"))
cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "2"))

## Integrate Seurat UMAP embeddings to ensure spatial consistency across objects
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(temp_sce , reduction = "umap")
int.embed <- Embeddings(subsce , reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

# root_group = colnames(cds)[clusters(cds) == 1]
# cds = order_cells(cds, root_cells = root_group)

## Pseudotime values across cell types
## Higher values indicate higher differentiation states; provided for demonstration purposes
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
                                                 graph_label_size=1)+ scale_color_manual(values = pal_jco("default", alpha = 0.6)(23))

# Visualization of differential genes along the trajectory
Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=6)

# Identify TOP genes based on Moran's I
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by = "seurat_clusters",min_expr=0.5, ncol= 2,cell_size=1.5) + scale_color_manual(values = pal_jco("default", alpha = 0.6)(11)) 
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by = "pseudotime",min_expr=0.5, ncol= 2,cell_size=1.5) + scale_color_manual(values = pal_jco("default", alpha = 0.6)(11)) 

# Selective visualization of genes of interest
plot_cells(cds, genes=c("Mt2","Mt1","Cxcl3"),show_trajectory_graph=FALSE)#IR
pData(cds) 

# Pseudotime scatter plots
plot_genes_in_pseudotime(cds[c("Mt2","Mt1","Cxcl3"),],
                         color_cells_by="seurat_clusters",
                         min_expr=0.5, ncol= 2,cell_size=1.5)+ scale_color_manual(values = pal_jco("default", alpha = 0.6)(7)) 

############################ Monocle2 Pseudotime Analysis #####################################
# Convert Seurat counts to a sparse matrix for Monocle2 ingestion
library(monocle)
Idents(temp_sce2) <- "label"
sparse_data <-  as(as.matrix(temp_sce2[["RNA"]]$counts),"sparseMatrix") 

# Construct the metadata object for cell data
mdata <- new("AnnotatedDataFrame",data = temp_sce2@meta.data)

# Construct the gene annotation object (gene_short_name)
fData <- data.frame(gene_short_name=row.names(sparse_data),row.names =row.names(sparse_data))
fd <- new("AnnotatedDataFrame",data = fData)

# Instantiate the Monocle CDS object
monocle_cds <- newCellDataSet(cellData = sparse_data,
                              phenoData = mdata,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

# Estimate size factors and dispersions
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

# Filter lowly expressed genes
monocle_cds <- detectGenes(monocle_cds,min_expr = 3)

# Feature selection: Identify highly variable genes
# Leverage Seurat-calculated variable features
expressed_genes <- VariableFeatures(temp_sce2)
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],
                                      fullModelFormulaStr = "~label")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) ## Enforce strict qval < 0.01 threshold.
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)

# Manifold learning and dimensionality reduction
monocle_cds <- reduceDimension(monocle_cds,max_components = 2,reduction_method = "DDRTree")

# Pseudotime ordering
monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds,root_state = "3")
monocle_cds <- orderCells(monocle_cds,root_state = "5")

# Trajectory Visualization
# Colors can also represent 'Pseudotime'
p <- plot_cell_trajectory(monocle_cds,color_by = "Pseudotime",size = 1,show_backbone = TRUE) + facet_wrap("~orig.ident",nrow=1)
p0 <- plot_cell_trajectory(monocle_cds,color_by = "label",size = 1,show_backbone = TRUE) + facet_wrap("~orig.ident",nrow=1)
p1 <- plot_cell_trajectory(monocle_cds,color_by = "label",size = 1,show_backbone = TRUE)
theme(legend.position = 'none',panel.border = element_blank())
#+scale_color_manual(values= color) # Custom palette

# Plot complex trajectory trees
p2 <- plot_complex_cell_trajectory(monocle_cds,x=1,y=2,color_by = "label")+
  theme(legend.title = element_blank())
#+scale_color_manual(values= color) # Custom palette
plot_cell_trajectory(monocle_cds,color_by = "State",size = 1,show_backbone = TRUE) 
p0|p|p1|p2


## Identification of pseudotime-dependent genes
Time_diff <- differentialGeneTest(monocle_cds[ordering_genes,],cores =1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

# Filter for genes expressed in at least 10 cells
Time_diff <- Time_diff %>%
  filter(num_cells_expressed > 10)

# Select top 100 significantly varying genes along the trajectory
Time_genes <- top_n(Time_diff,n=100,desc(qval)) %>%  pull(gene_short_name) %>% as.character()
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p4 <- plot_pseudotime_heatmap(monocle_cds[Time_genes,],num_clusters = 4,show_rownames = T,return_heatmap = T)

## Extract pseudotime-associated genes categorized by cluster
clusters <- cutree(p4$tree_row,k=4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_clusters"
table(clustering)




##### Jaccard Similarity Calculation ####
ceshi_markers <- read.csv(file = "../E-MTAB-7895/macrophagesmarkers.csv",row.names = 1)
ceshi_markers <- temp_sce.markers 
ceshi_markers <- read.csv(file = "../E-MTAB-7895/Granulocytesmarkers.csv",row.names = 1)
ceshi_markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> ceshi_top100


valid_markers <- markers_temp_sce
valid_markers <- markers_subsce
valid_markers <- read.csv(file = "../GSE163129_RAW/test_Neutrophil_markers.csv",row.names = 1)
valid_markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> valid_top100

# Extract genes from discovery set (ceshi) with absolute log2FC > 0.5
ceshi_top100 <- ceshi_markers %>%
  filter(abs(avg_log2FC) > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)))

# Extract genes from validation set (valid) with absolute log2FC > 0.5  
valid_top100 <- valid_markers %>%
  filter(abs(avg_log2FC) > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)))

# Extract unique cluster identities from both dataframes
ceshi_clusters <- unique(ceshi_top100$cluster)
valid_clusters <- unique(valid_top100$cluster)

# Initialize Jaccard similarity matrix
jaccard_matrix <- matrix(nrow = length(ceshi_clusters), ncol = length(valid_clusters))
colnames(jaccard_matrix) <- valid_clusters
rownames(jaccard_matrix) <- ceshi_clusters

# Compute pairwise Jaccard similarity indices
for (i in seq_along(ceshi_clusters)) {
  ceshi_cluster <- ceshi_top100$gene[ceshi_top100$cluster == ceshi_clusters[i]]
  for (j in seq_along(valid_clusters)) {
    valid_cluster <- valid_top100$gene[valid_top100$cluster == valid_clusters[j]]
    intersection <- length(intersect(ceshi_cluster, valid_cluster))
    union <- length(union(ceshi_cluster, valid_cluster))
    jaccard_matrix[i, j] <- intersection / union
    #jaccard_matrix[i, j] <- intersection
  }
}


# Inter-study Similarity Heatmap
dev.off()
pheatmap(jaccard_matrix,
         main = "Jaccard Similarity Heatmap",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("lightblue", "white", "red"))(100),
         fontsize_row = 20,
         fontsize_col = 20,
         display_numbers = TRUE,
         number_color = "black")


fisher_test <- function(x, y) {
  a <- length(intersect(x, y))
  b <- length(x) - a
  c <- length(y) - a
  d <- max(0, 100 - a - b - c)
  p_value <- tryCatch({
    fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
  }, error = function(e) {
    NA
  })
  p_value
}

# Calculate Fisher's Exact Test p-value matrix for overlap significance
fisher_matrix <- matrix(nrow = length(ceshi_clusters), ncol = length(valid_clusters))
colnames(jaccard_matrix) <- valid_clusters
rownames(jaccard_matrix) <- ceshi_clusters


# Execute similarity indexing
for (i in seq_along(ceshi_clusters)) {
  ceshi_cluster <- ceshi_top100$gene[ceshi_top100$cluster == ceshi_clusters[i]]
  for (j in seq_along(valid_clusters)) {
    valid_cluster <- valid_top100$gene[valid_top100$cluster == valid_clusters[j]]
    fisher_matrix[i, j] <- fisher_test(ceshi_cluster, valid_cluster)
  }
}



set.seed(123)

# Generate first dataset (simulated)
cluster1 <- lapply(1:6, function(i) sample(letters, 10, replace = TRUE))

# Generate second dataset (simulated)  
cluster2 <- lapply(1:5, function(i) sample(letters, 10, replace = TRUE))



jaccard_sim <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}


# Compute Jaccard similarity matrix
jaccard_matrix <- matrix(NA, nrow = 6, ncol = 5)
rownames(jaccard_matrix) <- paste0("Cluster", 1:6)
colnames(jaccard_matrix) <- paste0("Cluster", 1:5)

for (i in 1:6) {
  for (j in 1:5) {
    jaccard_matrix[i, j] <- jaccard_sim(cluster1[[i]], cluster2[[j]])
  }
}

fisher_test <- function(x, y) {
  a <- length(intersect(x, y))
  b <- length(x) - a
  c <- length(y) - a
  d <- max(0, 100 - a - b - c)
  p_value <- tryCatch({
    fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
  }, error = function(e) {
    NA
  })
  p_value
}

# Compute Fisher's Exact Test significance matrix
fisher_matrix <- matrix(NA, nrow = 6, ncol = 5)
rownames(fisher_matrix) <- paste0("Cluster", 1:6)
colnames(fisher_matrix) <- paste0("Cluster", 1:5)

for (i in 1:6) {
  for (j in 1:5) {
    fisher_matrix[i, j] <- fisher_test(cluster1[[i]], cluster2[[j]])
  }
}


data <- expand.grid(row = rownames(jaccard_matrix), 
                    col = colnames(jaccard_matrix))
data$jaccard <- as.vector(jaccard_matrix)
data$p_value <- as.vector(fisher_matrix)

library(ggplot2)

# Quantitative Overlap Visualization
ggplot(data, aes(x = col, y = row, fill = jaccard)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p_value)), color = "white", size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

############################### DNB Expression Dynamics & Clustering ################################
# Custom expression matrix utilities
### Function to calculate DNB mean expression values
i=2
h=1
Idents(temp_sce2) <- "label"
DNBs_expression <- function(expression_matrix,target_cluster,DNB_names){
  failed_genes <- c()
  # Remove redundant transcription factor identifiers
  DNB_names <- unique(na.omit(DNB_names))
  #tryCatch({
  average_expression <- c()
  df <- data.frame(matrix(nrow = length(target_cluster), ncol = length(DNB_names)))
  rownames(df) <- target_cluster
  colnames(df) <- DNB_names
  # Retrieve cell indices for the target cluster
  for (i in seq_along(target_cluster)){
    cell_indices <- which(temp_sce2$label == target_cluster[i])
    for (h in 1:length(DNB_names)){
      # Validate gene existence in the single-cell matrix
      if (DNB_names[h] %in% rownames(expression_matrix)) {
        average_expression <- mean(expression_matrix[DNB_names[h], cell_indices])
        
      }else {
        failed_genes <<- c(failed_genes, DNB_names[h])
        average_expression <- NA # Calculation failed, return null value
      }
      df[i, DNB_names[h]] <- average_expression
    }
    if (length(failed_genes) > 0) {
      
      print(paste("List of failed DNBs:", paste(unique(failed_genes), collapse = ", ")))
      cat("\n")}
  }
  return(df)
  #})
}
DNB_list <- as.list(read.csv("./fibrobalst_DNB.csv"))
DNB_list <- as.list(read.csv("./Granulocytes_DNB.csv"))
DNB_list <- as.list(read.csv("./test_macropahges_DNB.csv"))
DNB_list <- as.list(read.csv("./test_Granulocytes_DNB.csv.csv"))
Idents(temp_sce2) # Verification: Ensure 'label' is set as the current identity
Idents(subsce) 

expression_matrix = temp_sce2[["RNA"]]$scale.data # Normalized and scaled data

target_cluster <- c("1","4","9") # Fibroblast DNB (Cluster 4, early state)
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster,unique(DNB_list$X4))
exps <- t(DNBs_expression.df)
write.table(exps, "fibroblast_4exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)





target_cluster <- c("9","6","7") # Fibroblast DNB (Cluster 6, late state)
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X6))
exps <- t(DNBs_expression.df)
write.table(exps, "fibroblast_6exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)


target_cluster <- c("0","1","3") # Granulocytes DNB (Cluster 1, early state)
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X1))
exps <- t(DNBs_expression.df)
write.table(exps, "Granulocytes_013exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

target_cluster <- c("2","0","1") # Granulocytes DNB (Cluster 0, early state)
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X0))
exps <- t(DNBs_expression.df)
write.table(exps, "Granulocytes_201exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)


target_cluster <- c("5","4","7") # Macrophages DNB (Cluster 4, early state)
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X4))
exps <- t(DNBs_expression.df)
write.table(exps, "macrophages_547exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

target_cluster <- c("1","3","0") # Macrophages DNB (Cluster 3, late state)
DNBs_expression.df <-DNBs_expression(expression_matrix,target_cluster, unique(DNB_list$X3))
exps <- t(DNBs_expression.df)
write.table(exps, "macrophages_130exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)





?VlnPlot()
## ClusterGVis clustering analysis
#exps <- exps[!duplicated(exps[,1]), ]
#rownames(exps) <- exps[,1]
#exps <- exps[,-1]

exps <- as.matrix(exps)
exps <-na.omit(exps)
head(exps)
# ClusterGVis: Note potential conflict with Monocle3. Unload Monocle3 if inherited function errors occur.
library(ClusterGVis)
library(sf)

markers_subsce %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

markers_temp_sce %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

# Utilize cluster means as input; gene deduplication is active by default
st.data1 <- prepareDataFromscRNA(object = temp_sce,
                                  #object = subsce,
                                  diffData = top10,
                                  showAverage = TRUE)
ck <- clusterData(exp = exps,
                  #cluster.method = "mfuzz",
                  cluster.method = "kmeans",
                  cluster.num = 4)

# Gene label annotation
markGenes = rownames(exps)[sample(1:nrow(exps),32,replace = F)]
pdf('addgene.pdf',height = 10,width = 6,onefile = F)
visCluster(object = ck,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes)
dev.off()

# GO Enrichment Analysis (Mus musculus genome)
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

# Remove geneID column if added in the previous step to facilitate plotting
enrich <- subset(enrich, select = -geneID)
head(enrich,3)

pdf('term4.pdf',height = 12,width = 12,onefile = F)
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
  go.size = 10,
  cluster.order = c(1:10),
  #mulGroup = c(22,14,37,61,42,4),
  #mline.col = c(ggsci::pal_lancet()(7))
)
dev.off()

# Trend line visualization
visCluster(object = ck,
           plot.type = "line",
           ms.col = c("green","orange","red"), # Color customization
           add.mline = TRUE  # Remove median line
)

# Quantification of temporal expression trends for each gene row
trends <- apply(exps, 1, function(row_vec) {
  diff_vec <- diff(row_vec)
  if (all(diff_vec > 0)) {
    return("Monotonic Increase")
  } else if (all(diff_vec < 0)) {
    return("Monotonic Decrease")
  } else if (diff_vec[1] > 0 && diff_vec[2] < 0) {
    return("Critical State High")
  } else if (diff_vec[1] < 0 && diff_vec[2] > 0) {
    return("Critical State Low")  # Representing "Intermediate High" logic
  } else {
    return("Unclassified")
  }
})

# Execute trend-based cluster sorting
sorted_idx <- order(trends)

# Reorder expression matrix and trend metadata based on sorted indices
exps <- exps[sorted_idx, ]
sorted_trends <- as.data.frame(trends[sorted_idx])
colnames(sorted_trends) <- c("Category")

# Sophisticated Heatmap Rendering
library(ComplexHeatmap)
library(circlize)
pdf('heatmap4.pdf',height = 16,width = 20,onefile = F)
col_fun = colorRamp2(c(-2, 0, 2), c("#4169E1", "white", "#DC143C"))
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
dev.off()

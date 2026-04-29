###### Validation Set #####
rm(list= ls())
options(stringsAsFactors = F)
getwd()
setwd("/home/tsh/valid_E-MTAB-9816/")
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
rawdata_path <- "/home/tsh/valid_E-MTAB-9816/"  # Obtain paths for single-cell data files across different treatments
filename <- list.files(rawdata_path)
rawdata_path <- paste(rawdata_path,filename,sep = '')
# Read10X loads barcodes and related files for each timepoint to generate Seurat objects
sceList <- lapply(rawdata_path, function(x){
  obj <- CreateSeuratObject(counts = Read10X(x),
                            project = str_split(x,'/')[[1]][5])
})
# The merge function consolidates Seurat objects, using the first element of sceList as the seed and iteratively merging others.
names(sceList) <- filename
sce <- merge(sceList[[1]],sceList[-1],
             add.cell.ids = names(sceList), # add.cell.ids specifies identifiers for merged groupings
             project = 'Bobo'  # project parameter sets the master project name to 'Bobo'
)
head(sce@meta.data, 5)

# Calculate the percentage of mitochondrial genes
sce <- PercentageFeatureSet(sce,'^mt',col.name = 'percent_mt')
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3,pt.size = 0,group.by = 'orig.ident')
# Quality Control (QC)
dim(sce)
# Filter cells based on criteria provided in literature
sce <- subset(sce,
              subset = nCount_RNA <20000 & nFeature_RNA>200 & nFeature_RNA<4000 & percent_mt<10)
dim(sce)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce)
sce <- FindNeighbors(sce, dims = 1:30, reduction = "pca")
for (i in c(0.01,0.05,0.1,0.2,0.3,0.5,0.8,1,2)){
  sce <- FindClusters(sce,resolution = i)
}
clustree(sce ,prefix = "RNA_snn_res.")
# Batch correction using "Harmony"
sce <- IntegrateLayers(
  object = sce, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:30)
sce <- FindClusters(sce, resolution = 2, cluster.name = "harmony_clusters")
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
sce <- RunTSNE(sce, reduction = "harmony", dims = 1:30, reduction.name = "tsne.harmony")



# Batch correction using "CCA"
sce <- IntegrateLayers(
  object = sce, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

sce <- FindNeighbors(sce, reduction = "integrated.cca", dims = 1:30)
sce <- FindClusters(sce, resolution = 0.6, cluster.name = "cca_clusters")
sce <- RunUMAP(sce, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
sce <- RunTSNE(sce, reduction = "integrated.cca", dims = 1:30, reduction.name = "tsne.cca")
DimPlot(sce, reduction = "umap.cca", group.by = c("orig.ident","seurat_clusters"),label = TRUE)
DimPlot(sce, reduction = "tsne.cca", group.by = c("orig.ident","seurat_clusters"),label = TRUE)
sce <- JoinLayers(sce)
sce
Idents(sce) <- "orig.ident"
custom_order <- c("0day","1day","3day","7day","14day","28day")
sce@meta.data$orig.ident <- factor(sce@meta.data$orig.ident, levels = custom_order)
DimPlot(sce, reduction = "umap.cca", group.by = c("orig.ident","seurat_clusters"),label = TRUE)
DimPlot(sce, reduction = "tsne.cca", group.by = c("orig.ident","seurat_clusters"),label = TRUE)
saveRDS(sce, file = "sce.rds")
sce <- readRDS("./sce.rds")
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
.rs.restartR()
BiocManager::install("gypsum")
remotes::install_github("LTLA/celldex")
devtools::install_github('thomasp85/farver')
BiocManager::install("celldex")
library(celldex)
## Automated annotation via SingleR
mrsd.se <- ImmGenData()
mrsd.se <- MouseRNAseqData()
MI.integrated.se <- as.SingleCellExperiment(sce, assay = "RNA")
mrsd.common <- intersect(rownames(MI.integrated.se), rownames(mrsd.se))
mrsd.se <- mrsd.se[mrsd.common,]
MI.integrated.se <- MI.integrated.se[mrsd.common,]
MI.integrated.se <- logNormCounts(MI.integrated.se)
MI.mrsd.pred <- SingleR(test = MI.integrated.se, ref = mrsd.se, method = "single", labels = mrsd.se$label.main)
sce[["single_R"]] <- MI.mrsd.pred$labels
cell.data <- table(sce@meta.data$single_R,sce@meta.data$seurat_clusters)

## Manual annotation based on canonical markers
knowmarkers <- factor(c(
  "Col1a1","Gsn","Dcn", # Fibroblast I
  "Wif1","Dkk3", # Fibroblast II
  "Mt1","Mt2", "Timp1", # Fibroblast III
  "Cthrc1","Postn","Acta2", # Myofibroblasts
  "Pecam1","Cdh5","Fabp4", #Endothelium
  "Chil3","Plac8","Ly6c2","Cd14", # Monocytes I
  "C1qa","Cd68","Ms4a7", #Macrophages
  "H2-Ab1","Cd74",  # Dendritic cells
  "Cd79a","H2-DMb2", #B cells
  "S100a9","S100a8", # Granulocytes
  "Nkg7","Cd3d","Ms4a4b","Gzma", # NK cell
  "Rgs5","Kcnj8","P2ry14", #SMCs_Pericytes
  "Plp1","Kcna1" #Schwann cells
))
print(knowmarkers)
Idents(sce) <- "seurat_clusters"
DotPlot(object = sce, features =knowmarkers)+RotatedAxis()+scale_x_discrete("")+scale_y_discrete("")
FeaturePlot(sce, features = knowmarkers,reduction = "umap.cca")
VlnPlot(sce,features = Macrophage_genes,group.by = c("seurat_clusters"),pt.size = 0)
table(sce$orig.ident,Idents(sce))
anno_human <- c("Fibroblasts", "Fibroblasts", "Fibroblasts", "Macrophages", "Endothelial", "B", "Fibroblasts", "Fibroblasts", "NK", "Fibroblasts", "Endothelial", "Granulocytes", "SMCs_Pericytes", "Fibroblasts", "NK", "Monocytes", "Fibroblasts", "Monocytes", "Myofibroblasts", "Endothelial", "Myofibroblasts", "T", "Schwann")
names(anno_human) <- levels(sce)
# Update Identities (Idents)
sce <- RenameIdents(sce, anno_human)
# Append CellType information to metadata
sce$anno_human <- Idents(sce)
head(sce@meta.data)
DimPlot(sce, reduction = "tsne.cca", group.by = c("seurat_clusters","orig.ident","anno_human"),label = T)
DimPlot(sce, reduction = "umap.cca", group.by = c("seurat_clusters","orig.ident","anno_human"),label = T)
DimPlot(sce, reduction = "umap.cca", group.by = c("anno_human"),label = T)
DimPlot(sce, reduction = "umap.cca", split.by =  c("orig.ident"),label = T)
FeaturePlot(sce, features = c("Col1a1","Ly6c2","Cd68","S100a9"),label = TRUE,cols = c("gray", "red"), pt.size=0.1)
#saveRDS(sce, file = "seurat_obj.rds")
sce <- readRDS("./seurat_obj.rds")


# Visualize changes in cell proportions across timepoints
Idents(sce) <-"anno_human"
table(sce$orig.ident)
table(Idents(sce),sce$orig.ident)
cell_data <- as.data.frame(prop.table(table(Idents(sce),sce$orig.ident),margin = 2))
colnames(cell_data) <- c("细胞类型", "时间", "细胞丰度")
# Sort by chronological order
cell_data <- arrange(cell_data, rev(时间))
# Define custom cell type order
cell_data <- mutate(cell_data, 细胞类型 = factor(细胞类型, levels = c("Fibroblasts",  "Myofibroblasts",  "SMCs_Pericytes", "Endothelial","Monocytes","Macrophages","Granulocytes", "B", "NK",  "T", "Schwann")))

## Optional: sort by abundance
## cell_data <- mutate(cell_data, 细胞类型 = reorder(细胞类型, 细胞丰度))
# Create stacked bar charts
plot <- ggplot(cell_data, aes(x = 细胞丰度, y = fct_rev(时间), #fill = 细胞类型, 
                              fill = fct_rev(细胞类型),
)) +
  geom_col() +
  labs(x = "细胞丰度", y = "时间", fill = "细胞类型") +
  theme_minimal()
# Define custom color palette for cell types
color_palette <-  c("Macrophages" = "#FF9900","Fibroblasts" = "#FF6666",
                    "Endothelial" = "#CC9933","Monocyte" = "#99CC33","NK" = "#69CC33",
                    "SMCs_Pericytes"="#99CCFF","B"="#33CC33","Granulocytes" = "#09CCFF",
                    "T"="#99CCCC","Monocytes"="#4169E1",
                    "Myofibroblasts"="#9966CC","Schwann"="#FF99CC",
                    "T"="#CC3399")
plot <- plot + scale_fill_manual(values = color_palette)
plot <- plot + guides(colour = guide_legend(reverse = TRUE))
# Adjust axis labels and legend for better readability of temporal data
plot <- plot +theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    axis.text = element_text(size = 12, color = "black"),
                    axis.title = element_text(size = 14, color = "black"),
                    legend.title= element_text(size = 16, color = "black"),
                    legend.text= element_text(size = 14, color = "black"))
# Display the plot with border styling
plot <- plot + geom_col(color = "black", size = 0.5)
print(plot)

### Marker DotPlot
Idents(sce) <- "anno_human"
cell_order<-c("Fibroblasts",  "Myofibroblasts",  "SMCs_Pericytes", "Endothelial","Monocytes","Macrophages","Granulocytes", "B", "NK",  "T", "Schwann")
sce@meta.data$anno_human <- factor(sce@meta.data$anno_human, levels = cell_order)
knowmarkers <- factor(c(
  "Col1a1","Gsn","Dcn",  "Mt1","Mt2",  # Fibroblast I
  "C1qa","Cd68","Ms4a7", #Macrophages
  "Pecam1","Cdh5","Fabp4", #Endothelium
  "H2-Ab1","Cd74",  # Dendritic cells
  "Cd79a","H2-DMb2", #B cells
  "Nkg7","Ms4a4b","Gzma", # NK cell

  "S100a9","S100a8", # Granulocytes
  "Rgs5","Kcnj8","P2ry14", #SMCs_Pericytes
  "Chil3","Plac8","Ly6c2","Cd14", # Monocytes I
  "Cthrc1","Postn","Acta2", # Myofibroblasts
  "Cd3d",
  "Plp1","Kcna1" #Schwann cells
))
DotPlot(object = sce, features =as_factor(knowmarkers),dot.scale = 6)+
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
wanted_clusters <- c("Fibroblasts","Myofibroblasts")
wanted_clusters <- c("Macrophage","Monocyte")
wanted_clusters <- c("Neutrophil")
subsce <- subset(sce,subset = anno_human %in% wanted_clusters)
## Filter ribosomal genes if necessary
rb.genes <- rownames(subsce)[grep("^Rp[sl]",rownames(subsce))]
Ct <- GetAssayData(object = subsce, layer = "counts")
percent.ribo <- Matrix::colSums(Ct[rb.genes,])/Matrix::colSums(Ct)*100
subsce <- AddMetaData(subsce, percent.ribo, col.name = "percent.ribo")
VlnPlot(subsce, features = "percent.ribo", pt.size = 0.1 ) + NoLegend()
subsce <- subset(subsce,subset = percent.ribo<30)
Idents(subsce) <-"orig.ident"
Idents(subsce)
# Re-process subsets using "harmony"
library(harmony)
subsce <- RunHarmony(subsce,c("orig.ident"))
subsce <- FindVariableFeatures(subsce)
subsce <- ScaleData(subsce,features = rownames(subsce))
subsce <- RunPCA(subsce)
subsce <- FindNeighbors(subsce, dims = 1:30, features = VariableFeatures(object = subsce))
subsce <- FindClusters(subsce, resolution = 0.6) # Fibroblast Macrophages
subsce <- FindClusters(subsce, resolution = 0.25) # Granulocytes
subsce <- RunUMAP(subsce, dims = 1:30, reduction = "pca")
subsce <- RunUMAP(subsce, dims = 1:15, reduction = "pca")
DimPlot(subsce, reduction = "umap", group.by = c("anno_human"),label = T)
DimPlot(subsce, reduction = "umap", group.by = c("anno_human", "seurat_clusters"),label = T)
DimPlot(subsce, reduction = "umap", group.by = c("seurat_clusters"),label = T)
DimPlot(subsce, reduction = "umap", split.by = c("orig.ident"),label = T)
clustree(subsce ,prefix = "RNA_snn_res.")
markers_subsce <- FindAllMarkers(subsce, only.pos = TRUE, min.pct =0.25, logfc.threshold = 0.25)
write.csv(markers_subsce,file = "test_Neutrophil_markers.csv")
table(subsce$orig.ident,Idents(subsce))

DotPlot(object = subsce, features =anti_inflammatory,scale = FALSE)+RotatedAxis()+scale_x_discrete("")+scale_y_discrete("")
## Establish trajectory logic; markers may not fully align with all requirements
Pro_inflammatory <- c("Il1b", "Il6","Cxcl1", "Cxcl2", "Cxcl5", "Cxcl8",'Mt2', 'Mt1','Angptl4',"Gsn","Pdpn")
Angiogenesis <- c("Vegfb","Vegfa","Lrg1","Pdgfra","Il10","Kdr","Hif1a")
Proliferation <- c('Stmn1', 'H2afz', 'Cks2', 'Cenpa', 'Hmgb2')
anti_inflammatory <- c('Col1a3','Col1a1', 'Postn',"Acta2","Cthrc1","Comp","Sfrp2")
# Fibroblast population markers
knownmarkers <- c( "Il6","Cxcl1",'Mt2', 'Mt1',"Gsn","Pdpn",
                   'Angptl4',"Vegfb","Vegfa","Pdgfra","Hif1a","Fgfr1","Igfr","Egfr",
                   'Stmn1', 'H2afz', 'Cks2', 'Cenpa', 'Hmgb2',
                   'Col1a1', 'Postn',"Acta2","Cthrc1","Comp","Sfrp2")
# Monocyte/Macrophage population markers
knownmarkers <-c( "Ly6c2","Plac8",  # Monocytes
                  "Arg1","Il1b","Ccl6","F10", # Transition state
                  "Ifit1","Ifit3","Isg15", "Rsad2",  # Leukocyte chemotaxis & IFN-stimulated
                  "Stmn1","Cenpa","Pclaf", # Proliferative sub-population
                  "Pf4", "Fn1","Slc7a2","Sdc3", # Pro-angiogenic
                  "Trem2","Spp1","Gpnmb","Gdf15","Fabp5","Apoe","Ms4a7", # Reparative phenotype
                  "Ccr2","H2-Aa","H2-Eb1","H2-DMa",  # Monocyte-derived antigen presentation 
                  "Col5a2","Col1a1","Col3a1","Dcn", # Collagen-like
                  "Mgl2","Folr2","Adgre1" # Resident macrophages
                  
)

# Neutrophil population markers
knownmarkers <-c( "Cxcl3", "Ccl4", "Ccl3",  "Cxcl2", "Basp1", "Cish", "Il1rn","Xbp1","Hif1", # Neutro1-leukocyte and pro-inflammatory
                  "Arg1","Il1b","Ccl6","F10", # Transition state
                  "Ifit1","Ifit3","Isg15", "Rsad2",  # Leukocyte chemotaxis & IFN-stimulated
                  "Stmn1","Cenpa","Pclaf", # Proliferative sub-population
                  "Pf4", "Fn1","Slc7a2","Sdc3", # Pro-angiogenic
                  "Trem2","Spp1","Gpnmb","Gdf15","Fabp5","Apoe","Ms4a7", # Reparative phenotype
                  "Ccr2","H2-Aa","H2-Eb1","H2-DMa",  # Monocyte-derived antigen presentation 
                  "Col5a2","Col1a1","Col3a1","Dcn", # Collagen-like
                  "Mgl2","Folr2","Adgre1" # Resident macrophages
                  
)




Idents(subsce) <-"seurat_clusters"
subsce <- JoinLayers(subsce)
sce.markers <- FindAllMarkers(subsce, only.pos = TRUE, min.pct =0.25, logfc.threshold = 0.25)
write.csv(sce.markers,file = "allmarkers.csv")
write.csv(sce.markers,file = "Granulocytesmarkers.csv")
DotPlot(object = subsce, features =Pro_inflammatory)+RotatedAxis()+scale_x_discrete("")+scale_y_discrete("")
DimPlot(subsce, reduction = "umap", group.by = c("anno_human", "seurat_clusters"),label = T)
FeaturePlot(subsce, features = knownmarkers,reduction = "umap" )
table(Idents(subsce),subsce$orig.ident)
cell_data <- as.data.frame(prop.table(table(Idents(subsce),subsce$orig.ident),margin = 2))




FeaturePlot(subsce, features = c("Ly6a","Mt1","Mt2","Angptl4","Ccl2","Prg4","Csf1","Tgfb1","Cxcl1","Cxcl5"),reduction = "umap" )
FeaturePlot(subsce, features = c("Stmn1","Tubb3","Acta2","Tpm2","Cthrc1","H2afz","Ckap4"),reduction = "umap" )
FeaturePlot(subsce, features = c("Pdgfra","Hif1a","Vegfd","Mmp14"),reduction = "umap" )
FeaturePlot(subsce, features = c("Comp","Sfrp2","Col1a1","Dkk3","Postn","Ly6a","Acta2","Ccnb2"),reduction = "umap" )
FeaturePlot(subsce, features = c("Pdgfra","Ly6a","Postn","Ccl2","Angptl4","Cilp","Acta2","Ccnb2","Cdk1","Col1a1","Comp","Sfrp2","Sfrp2","Tgfb1","Scx","Thbs4","Ccn4"),reduction = "umap" )
VlnPlot(subsce,features =c("Pdgfra","Ly6a","Postn","Ccl2","Angptl4","Cilp","Acta2","Ccnb2","Cdk1","Col1a1","Comp","Sfrp2","Tgfb1","Scx","Thbs4","Ccn4"),group.by = c("seurat_clusters"),pt.size = 0)


### Cell Annotation ####
## Clusters 7, 24
VlnPlot(subsce,features =c("Ly6a","Mt1","Mt2","Angptl4","Ccl2","Prg4","Csf1","Tgfb1","Cxcl1","Cxcl5"),group.by = c("seurat_clusters"),pt.size = 0)
## Clusters 18, 20
VlnPlot(subsce,features =c("Postn","Cilp","Acta2","Ccnb2","Cdk1","Mki67"),group.by = c("seurat_clusters"),pt.size = 0)
## Clusters 6, 9, 10, 11, 22, 28  
VlnPlot(subsce,features =c("Postn","Cilp","Acta2","Cthrc1","Col1a1","Lyz2","Notch2"),group.by = c("seurat_clusters"),pt.size = 0)
## Cluster 11
VlnPlot(subsce,features =c("Postn","Comp","Sfrp2","Tgfb1","Scx","Thbs4","Ccn4"),group.by = c("seurat_clusters"),pt.size = 0)
## Cluster 24
VlnPlot(subsce,features =c("Postn","Pdgfra","Hif1a","Vegfd","Mmp14","Angptl4","Angpt2"),group.by = c("seurat_clusters"),pt.size = 0)

#### Save Current Clustering Results ####
## saveRDS(subsce, file = "subsce222.rds")
subsce <- readRDS("./subsce.rds")


### Cell Annotation ####

## Primary identification of Monocyte, M1, and M2 subtypes
## Monocytes (Clusters 5, 7, 12, 13, 18)
VlnPlot(subsce,features =c("Ly6c2","Plac8","Ccr2"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Ly6c2","Plac8","Ccr2"),reduction = "umap" )
## M1 subtype (iNOS, Il1b, Il6, Tnfa)
VlnPlot(subsce,features =c("Nos2","Il1b","Il6","Il23","Mmp8","Mmp9","Tnf"),group.by = c("seurat_clusters"),pt.size = 0)

## M2 subtype (Arg1, Tgfb, Il10, Vegfa/b/c/d)
VlnPlot(subsce,features =c("Arg1","Tgfb1","Tgfb2","Tgfb3","Il10","Il12","Vegfa"),group.by = c("seurat_clusters"),pt.size = 0)



## Leukocyte chemotaxis & Interferon-stimulated population (Clusters 7, 22; recommend merging 22 and 6)
VlnPlot(subsce,features =c("S100a9","S100a8","Ccl7","Rsad2","Isg15","Il18","Irf7","Cxcl10"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Ifit1","Ifit3","Isg15", "Rsad2"),reduction = "umap" )
## Anti-inflammatory potential (Clusters 1, 2, 3, 4)
VlnPlot(subsce,features =c("Arg1","Trem2","Spp1","Gpnmb","Ctsd","Cd63","Tgfb1","Tgfb2","Il10","Il12","Hmox1","Fabp4","Fabp5"),group.by = c("seurat_clusters"),pt.size = 0)
## Anti-inflammatory potential (Clusters 1, 2, 3, 4)
VlnPlot(subsce,features =c("Arg1","Trem2","Gdf15","Psap","Pld3","Nceh1","Tgfb1","Tgfb2","Tgfb3","Il10","Il12","Gpnmb","Hmox1","Fabp4","Fabp5"),group.by = c("seurat_clusters"),pt.size = 0)
## Anti-inflammatory potential (Cluster 20)
VlnPlot(subsce,features =c("Arg1","Trem2","Prdx1","Tgfb1","Tgfb2","Tgfb3","Il10","Il12","Gpnmb","Hmox1","Fabp4","Fabp5"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Trem2","Spp1","Gpnmb","Gdf15","Fabp5","Apoe","Ms4a7"),reduction = "umap" )
## Monocyte-derived antigen presentation (Clusters 0, 7, 10; Cluster 7 optional)
VlnPlot(subsce,features =c("Ccr2","H2-Eb1","H2-DMa","Il1b","Tnfsf9","Tnip3"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Ccr2","H2-Aa","H2-Eb1","H2-DMa"),reduction = "umap" )

## Pro-angiogenic (Clusters 6, 22)
VlnPlot(subsce,features =c("Pf4", "Fn1","Slc7a2","Sdc3","Pdgfra","Hif1a","Vegfd","Vegfa","Vegfc","Kdr"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c( "Pf4", "Fn1","Saa3","Fcna"),reduction = "umap" )
## Proliferative phenotype (Clusters 14, 11, 21)
VlnPlot(subsce,features =c( "Stmn1","Birc5","Top2a","Ube2c" ),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Stmn1","Cenpa","Pclaf"),reduction = "umap" )


## Tissue-resident macrophages (Cluster 8)
VlnPlot(subsce,features =c("Cd14","Mmp12","Cxc3r1","Il1b","Ccr2","H2-Eb1"),group.by = c("seurat_clusters"),pt.size = 0)
## Tissue-resident macrophages (Cluster 8)
VlnPlot(subsce,features =c( "Timd4","Lyve1","Folr2","Mrc1","Cd163","Ccr2","Igf1"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(subsce, features = c("Mgl2","Folr2","Adgre1"),reduction = "umap" )

## Clusters 15/9 may represent polarization transition states
FeaturePlot(subsce, features = c("Arg1","Il1b","Ccl6","F10"),reduction = "umap" )
## Cluster 17 is the collagen-expressing sub-population
FeaturePlot(subsce, features = c("Col5a2","Col1a1","Col3a1","Dcn"),reduction = "umap" )

Idents(subsce) <- "seurat_clusters"
inflammation_human <- c("Antigen-presenting", "Reparative", "Reparative", "Reparative", "Reparative", "Monocytes", "Neut", "IFN", "RCM", "Pro_inflammation", "Antigen-presenting", "Proliferation", "Monocytes", "Monocytes", "Proliferation", "Pro_inflammation", "Lymphatic", "Collagens", "Monocytes", "Monocytes", "Reparative", "Proliferation", "Neut", "Bcell")
names(inflammation_human) <- levels(subsce)
# Update Identifiers
subsce <- RenameIdents(subsce, inflammation_human)
# Append refined CellType information to metadata
subsce$inflammation_human <- Idents(subsce)
head(subsce@meta.data)
Idents(subsce) <- "inflammation_human"
DimPlot(subsce, reduction = "umap", split.by = c("orig.ident"),label = T)
DimPlot(subsce, reduction = "umap", group.by = c("inflammation_human", "seurat_clusters"),label = T)







### Targeted Subset Extraction ####
Idents(subsce) <- "anno_human"
wanted_clusters <- c("Fibroblasts") # Fibroblast lineage
wanted_clusters <- c("0","10","3","5","1","4","2","7","6","8") # Fibroblast lineage indices
wanted_clusters <- c("5","12","13","18","15","9","1","2","3","4","20","0","10","14","11","21","17") # Monocyte/Macrophage indices
temp_sce <- subset(subsce,subset = seurat_clusters %in% wanted_clusters)
temp_sce <- RunHarmony(temp_sce,c( "orig.ident" ))
temp_sce <- NormalizeData(temp_sce)
temp_sce <- FindVariableFeatures(temp_sce)
temp_sce <- ScaleData(temp_sce,features = rownames(temp_sce))
temp_sce <- RunPCA(temp_sce)
temp_sce <- FindNeighbors(temp_sce, dims = 1:30, features = VariableFeatures(object = temp_sce))
temp_sce <- FindClusters(temp_sce, resolution = 0.8)
temp_sce <- RunUMAP(temp_sce, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
temp_sce <- RunUMAP(temp_sce, dims = 1:30, reduction = "pca")
DimPlot(temp_sce, reduction = "umap", group.by = c("anno_human", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap.harmony", group.by = c("anno_human", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", split.by = c("orig.ident"),label = T)
#saveRDS(temp_sce, file = "valid_temp_sce.rds")
temp_sce <- readRDS("./valid_temp_sce.rds")
temp_sce <- subset(temp_sce,subset = seurat_clusters != "14")
table(Idents(temp_sce),temp_sce$orig.ident)
cell_data <-as.data.frame(prop.table(table(Idents(temp_sce),temp_sce$orig.ident),margin = 2))
cell_data <-as.data.frame(table(Idents(temp_sce),temp_sce$orig.ident),margin = 2)
colnames(cell_data) <- c("细胞类型", "时间", "细胞丰度")
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



FeaturePlot(temp_sce, features = c("Ly6a","Mt1","Mt2","Angptl4","Ccl2","Prg4","Csf1","Tgfb1","Cxcl1","Cxcl5"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Postn","Cilp","Acta2","Ccnb2","Cdk1","Mki67"),reduction = "umap" )
VlnPlot(temp_sce,features =c("Postn","Cilp","Acta2","Ccnb2","Cdk1","Mki67"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(temp_sce, features = c("Postn","Pdgfra","Hif1a","Vegfd","Vegfc","Mmp14","Angptl4","Dpep1","Dcn","Mmp2"),reduction = "umap" )
VlnPlot(temp_sce,features =c("Postn","Pdgfra","Hif1a","Vegfd","Vegfc","Mmp14","Angptl4","Dpep1","Dcn","Mmp2"),group.by = c("seurat_clusters"),pt.size = 0)
FeaturePlot(temp_sce, features = c("Postn","Cilp","Tnc","Fn1","Acta2","Cthrc1","Col1a1"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Postn","Comp","Sfrp2","Tgfb1","Scx","Thbs4","Ccn4","Ccn5","Col8a2"),reduction = "umap" )
DotPlot(object = temp_sce, features =knownmarkers)+RotatedAxis()+scale_x_discrete("")+scale_y_discrete("")

temp_sce2.markers <- FindAllMarkers(temp_sce, only.pos = TRUE, min.pct =0.25, logfc.threshold = 0.25)
write.csv(temp_sce2.markers,file = "valid_fibroblastmarkers2.0.csv")
write.csv(temp_sce2.markers,file = "macrophagesmarkers.csv")
FeaturePlot(temp_sce, features = c("Angptl4","Cks2","Sfrp2","Comp"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Ly6c2","Plac8","Ccr2"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Arg1","Il1b","Ccl6","F10"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Stmn1","Cenpa","Pclaf"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Trem2","Spp1","Gpnmb","Fabp5"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Ccr2","H2-Aa","H2-Eb1","H2-DMa"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Col5a2","Col1a1","Col3a1","Dcn"),reduction = "umap" )

## Save the final result
## saveRDS(temp_sce, file = "temp_sce.rds")


temp_sce <- readRDS("./temp_sce.rds") # Fibroblasts




##saveRDS(temp_sce, file = "inflammation_temp_sce.rds") # Monocytes/Macrophages
##saveRDS(temp_sce, file = "filtered_inflammation_temp_sce.rds")

temp_sce <- readRDS("../E-MTAB-7895/inflammation_temp_sce.rds")
temp_sce <- readRDS("../E-MTAB-7895/filtered_inflammation_temp_sce.rds")

##saveRDS(subsce, file = "Granulocytes_temp_sce.rds") # Neutrophils

temp_sce <- readRDS("../E-MTAB-7895/Granulocytes_temp_sce.rds")


######### Immune Compartment: Target Lineage Extraction #################################
# Based on Monocle 3 trajectory results, exclude proliferative and antigen-presenting sub-populations along with their respective progenitor sub-clusters.
Idents(temp_sce) <- "seurat_clusters"
Idents(temp_sce) <- "inflammation_human"
wanted_clusters <- c("6","12","13","3","5","2","0","14","10","4","9")
temp_sce <- subset(temp_sce,subset = seurat_clusters %in% wanted_clusters)
#temp_sce <- subset(temp_sce,subset = seurat_clusters != c("Proliferation"))
temp_sce$"label" <- Idents(temp_sce)

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
temp_sce <- FindClusters(temp_sce, resolution = 1.2)
temp_sce <- FindClusters(temp_sce, resolution = 1)
temp_sce <- FindClusters(temp_sce, resolution = 0.36)
temp_sce <- RunTSNE(temp_sce, dims = 1:30, reduction = "harmony")
temp_sce <- RunUMAP(temp_sce, dims = 1:30, reduction = "harmony")
temp_sce <- RunTSNE(temp_sce, dims = 1:30, reduction = "pca")
temp_sce <- RunUMAP(temp_sce, dims = 1:30, reduction = "pca")
DimPlot(temp_sce, reduction = "tsne", group.by = c("anno_human", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("anno_human", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("inflammation_human", "seurat_clusters","label"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("inflammation_human", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"),label = T)
DimPlot(temp_sce, reduction = "umap", split.by = c("orig.ident"),label = T)
FeaturePlot(temp_sce, features = c("Ly6c2","Arg1","Trem2","Col3a1"),reduction = "umap" )
FeaturePlot(temp_sce, features = c("Angptl4","Ccnb2","Acta2","Sfrp2"),reduction = "umap" )

markers_tempsce <- FindAllMarkers(temp_sce, only.pos = TRUE, min.pct =0.25, logfc.threshold = 0.25)

table(Idents(temp_sce),temp_sce$orig.ident)
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




# Export count matrix for all cellular expressed genes
write.table(temp_sce[["RNA"]]$counts, "counts.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Export data matrix for all cellular expressed genes
write.table(temp_sce[["RNA"]]$data, "data.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Export scale.data matrix for all cellular expressed genes
write.table(temp_sce[["RNA"]]$scale.data, "scale_data.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Export meta.data attribute matrix for all cells
write.table(temp_sce@meta.data, "meta.data.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)



### Secondary Targeted Sub-cluster Extraction (Optional) ####
Idents(temp_sce) <- "seurat_clusters"
wanted_clusters <- c("1","4","9") # DNB candidates for FIBROBLAST (Cluster 4)
wanted_clusters <- c("9","6","7") # DNB candidates for FIBROBLAST (Cluster 6)
wanted_clusters <- c("3","10","4") # DNB candidates for FIBROBLAST (Cluster 10)
wanted_clusters <- c("0","2","1") # DNB candidates for Granulocyte (Cluster 2)
wanted_clusters <- c("3","4","1") # DNB candidates for Monocyte/Macrophage (Cluster 4)
wanted_clusters <- c("5","8","3") # DNB candidates for Monocyte/Macrophage (Cluster 8)
wanted_clusters <- c("12","5","8") # DNB candidates for Monocyte/Macrophage (Cluster 8)
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
write.csv(temp_sce2.markers,file = "成纤维967差异基因.csv")
write.csv(temp_sce2.markers,file = "temp_sce2markers.csv")

temp_sce2.markers <- ?FindAllMarkers(temp_sce, only.pos = FALSE, min.pct =0, logfc.threshold = 0)
write.csv(temp_sce2.markers,file = "all_unfiltered_markers.csv")
############# Slingshot Pseudotime Inference #############################
# Convert to SingleCellExperiment format
filter_sce_trans <- as.SingleCellExperiment(temp_sce)
Slingshot <- slingshot(filter_sce_trans, 
                       clusterLabels = 'seurat_clusters',  # Select column name for cell annotation in colData
                       reducedDim = 'UMAP', 
                       reweight = FALSE, # Non-overlapping lineages
                       start.clus= "10",  # Specify starting point
                       end.clus = NULL     # Specify ending point
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
# Gradient colors
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) # Re-encode colors into 100 gradients to simulate a continuous scale
plotcol <- colors[cut(Slingshot$slingPseudotime_1, breaks=100)] # Segment lineage 1 into 100 intervals, assigning colors based on interval
plotcol <- colors[cut(Slingshot$slingPseudotime_2, breaks=100)] 


plotcol[is.na(plotcol)] <- "lightgrey" # Assign grey to points not belonging to lineage 1
Slingshot$plotcol <- plotcol
plot.new()
plot(reducedDims(Slingshot)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(Slingshot), lwd=1, col=brewer.pal(9,"Set1"))
legend("right",
       legend = paste0("lineage",1:1),
       col = unique(brewer.pal(3,"Set1")),
       inset=0.8,
       pch = 16)

#################### Slingshot Sub-trajectory Extraction Analysis ############################################
# Reconstruct colData as a data.frame (as standard SCE colData is not a native data.frame)
coldata <- colData(Slingshot)
coldata <- data.frame(celltype = coldata@listData$RNA_snn_res.0.36, 
                      sampleId = coldata@listData$orig.ident,
                      plotcol = coldata@listData$plotcol)
rownames(coldata) = Slingshot@colData@rownames

# Filter for lineage 1 information
# Cell barcodes
filter_cell <- dplyr::filter(coldata, plotcol != "lightgrey")
filter_cell <- rownames(filter_cell)
head(filter_cell)
# Count matrix
counts <- Slingshot@assays@data@listData$counts
filter_counts <- counts[,filter_cell]
filter_counts[1:5,1:5]
# Optional: Random sampling of 2000 cells to optimize computational time
# set.seed(111)
# scell <- sample(colnames(filter_counts), size = 2000)

# Apply sampled matrix
filter_counts = filter_counts[, scell]
dim(filter_counts)

# Instantiate a new SingleCellExperiment object with filtered counts
filter_sim <- SingleCellExperiment(assays = List(counts = filter_counts))

# Align colData
filter_coldata = colData(Slingshot)[colnames(filter_counts), 1:3]
filter_sim@colData = filter_coldata

# Align dimensionality reduction coordinates
rd = reducedDim(Slingshot)
filter_rd <- rd[colnames(filter_counts),]
reducedDims(filter_sim) <- SimpleList(UMAP = filter_rd)

# Execute K-Means clustering on UMAP coordinates
set.seed(111)
cl <- kmeans(filter_rd, centers = 6)$cluster # Specify number of clusters
head(cl)
colData(filter_sim)$kmeans <- cl

## Visualize K-Means clustering results
library(RColorBrewer)
mycolors = brewer.pal(6,"Set1") # Color palette

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
  guides(color = guide_legend(override.aes = list(size = 5)) # Adjust point size in legend
  )  


# Trajectory inference using K-Means defined sub-populations
filter_sim <- slingshot(filter_sim,
                        clusterLabels = 'kmeans',
                        reducedDim = 'UMAP',
                        start.clus= "2", # Specify starting cluster
                        end.clus = '5' # Specify terminal cluster
)
head(colnames(colData(filter_sim)))

plot(reducedDims(filter_sim)$UMAP, pch=16, asp = 1)
lines(SlingshotDataSet(filter_sim), lwd=2, col=mycolors)


#################### Pseudotime-associated Gene Analysis ##########################
BiocManager::install("tradeSeq")
library(tradeSeq)
# Fit negative binomial generalized additive model (NB-GAM)
counts <- filter_sim@assays@data$counts
crv <- SlingshotDataSet(filter_sim)
# Determine the number of knots (basis functions) required for NB-GAM modeling.
# This function evaluates optimal knot numbers and may take significant computational time.
# ~16 min for 2k cells depending on trajectory complexity.
set.seed(111)
icMat <- evaluateK(counts = counts, 
                   sds = crv, 
                   k = 3:10,    # Maximum of 12
                   nGenes = 200, # Number of genes included; defaults to 500
                   verbose = T)



# Selected nknots = 6 for final modeling.
set.seed(111)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
# Fit negative binomial GAM
# ~13 min for 2k cells
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




############# Monocle3 Pseudotime Analysis #############################
library(devtools)
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(Seurat)
library(leidenbase)
DefaultAssay(temp_sce) <- "RNA"
DefaultAssay(subsce) <- "RNA"

# Prepare data for Monocle3
# expression_matrix: Row=gene, Column=cell expression matrix
expression_matrix = temp_sce[["RNA"]]$counts
expression_matrix = subsce[["RNA"]]$counts

# cell_metadata: Metadata from Seurat object
cell_metadata = data.frame(temp_sce@meta.data)
cell_metadata = data.frame(subsce@meta.data)

# gene_metadata: Gene annotation (gene symbols)
gene_annotation = data.frame(expression_matrix[,1])
gene_annotation[,1] = row.names(gene_annotation)
colnames(gene_annotation)=c("gene_short_name")

## Construct Monocle3 CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# Preprocessing
cds <- preprocess_cds(cds, num_dim = 50,norm_method = c("log"))

## Dimensionality reduction (UMAP)
cds <- reduce_dimension(cds,reduction_method="UMAP",cores=5)
plot_cells(cds,color_cells_by = "seurat_clusters")

## Clustering: Minimal resolution used to ensure grouping as a single entity for better visualization
cds <- cluster_cells(cds,resolution = 0.0000001)
cds <- cluster_cells(cds,cluster_method = "louvain", reduction_method ="UMAP")

## Pseudotime learning
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster=FALSE, label_leaves=FALSE,label_branch_points=FALSE)

## Select developmental origin
## Assuming HEpiD cells represent the starting point
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
cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "0"))
cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "4"))

## Synchronize UMAP coordinates with Seurat objects to maintain spatial consistency
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(temp_sce , reduction = "umap.harmony")
int.embed <- Embeddings(temp_sce , reduction = "umap")
int.embed <- Embeddings(subsce , reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

# root_group = colnames(cds)[clusters(cds) == 1]
# cds = order_cells(cds, root_cells = root_group)

## Plot pseudotime values across lineages
## Higher values represent higher differentiation states; provided for demonstration
dev.off()
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

# Feature display for trajectory-dependent genes
Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=6)

# Select top genes based on Moran's I coefficient
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by="subanno_human",min_expr=0.5, ncol= 2,cell_size=1.5) + scale_color_manual(values = pal_jco("default", alpha = 0.6)(11)) 

# Visualization of specific genes
plot_cells(cds, genes=c("Mt2","Mt1","Cxcl3"),show_trajectory_graph=FALSE)#IR

# Pseudotime scatter visualization
plot_genes_in_pseudotime(cds[c("Mt2","Mt1","Cxcl3"),],
                         color_cells_by="seurat_clusters",
                         min_expr=0.5, ncol= 2,cell_size=1.5)+ scale_color_manual(values = pal_jco("default", alpha = 0.6)(7)) 


############################ Monocle 2 Pseudotime Analysis #####################################
# Convert Seurat count data to a sparse matrix for Monocle 2 ingestion
library(monocle)
Idents(temp_sce2) <- "label"
sparse_data <-  as(as.matrix(temp_sce2[["RNA"]]$counts),"sparseMatrix") 

# Construct original matrix containing cell metadata
mdata <- new("AnnotatedDataFrame",data = temp_sce2@meta.data)

# Construct matrix for gene_short_name
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

# Feature selection: identify highly variable genes
# Leverage Seurat-selected variable features
expressed_genes <- VariableFeatures(temp_sce2)
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],
                                      fullModelFormulaStr = "~label")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) ## Threshold set to q < 0.01.
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)

# Dimensionality reduction (DDRTree)
monocle_cds <- reduceDimension(monocle_cds,max_components = 2,reduction_method = "DDRTree")

# Pseudotime ordering
monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds,root_state = "3")
monocle_cds <- orderCells(monocle_cds,root_state = "1")

# Trajectory plotting
# Alternative coloring: 'Pseudotime'
p <- plot_cell_trajectory(monocle_cds,color_by = "Pseudotime",size = 1,show_backbone = TRUE) + facet_wrap("~orig.ident",nrow=1)
p0 <- plot_cell_trajectory(monocle_cds,color_by = "label",size = 1,show_backbone = TRUE) + facet_wrap("~orig.ident",nrow=1)
p1 <- plot_cell_trajectory(monocle_cds,color_by = "label",size = 1,show_backbone = TRUE)
theme(legend.position = 'none',panel.border = element_blank())

# Complex trajectory tree visualization
p2 <- plot_complex_cell_trajectory(monocle_cds,x=1,y=2,color_by = "label")+
  theme(legend.title = element_blank())
plot_cell_trajectory(monocle_cds,color_by = "State",size = 1,show_backbone = TRUE) 
p0|p|p1|p2


## Identification of pseudotime-associated genes
Time_diff <- differentialGeneTest(monocle_cds[ordering_genes,],cores =1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

# Filter genes expressed in at least 10 cells
Time_diff <- Time_diff %>%
  filter(num_cells_expressed > 10)

# Select top 100 significantly varying genes
Time_genes <- top_n(Time_diff,n=100,desc(qval)) %>%  pull(gene_short_name) %>% as.character()
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p4 <- plot_pseudotime_heatmap(monocle_cds[Time_genes,],num_clusters = 4,show_rownames = T,return_heatmap = T)

## Categorize pseudotime-associated genes by cluster
clusters <- cutree(p4$tree_row,k=4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_clusters"
table(clustering)

##################### DESTINY Pseudotime Analysis #############################
BiocManager::install("destiny")
library(destiny) 
library(Biobase) 
data(guo_norm)

### Automated Biobase data structure construction from Seurat object
library(Biobase)
ct <-GetAssayData(object = temp_sce)
ct<-ct[VariableFeatures(temp_sce),]
ct <- as.ExpressionSet(as.data.frame(t(ct)))

# Add metadata annotations
# Annotations accessible via ct$column or ct[['column']]
ct$celltype <- DPT@meta.data[,c("integrated_merge_cluster")]
dm <- DiffusionMap(ct,k = 10)
palette(cube_helix(4)) # Configure color palette
plot(dm, pch = 20, # Stylized points
     col_by = "celltype")

## Manual construction of ExpressionSet object
ct <- ExpressionSet(assayData=as.matrix(temp_sce[["RNA"]]$counts),
                    phenoData= mdata,featureData = fd)
dm <- DiffusionMap(ct,k = 3,n_pcs = 50)
plot(dm)
palette(cube_helix(6)) # Continuous scale using cube_helix
# palette(hue_pal()(6)) # Default ggplot2 colors alternative
plot(dm, pch = 20, 
     col_by = 'num_cells', # Color by vector or discrete value
     legend_main = 'Cell stage')

# 2D visualization
plot(dm, 1:2, pch = 20, col_by = 'num_cells',
     legend_main = 'Cell stage')

# 3D visualization using rgl
library(rgl)
plot3d(eigenvectors(dm)[, 1:3],
       col = log2(guo_norm$num_cells),
       type = 's', radius = .01)
view3d(theta = 10, phi = 30, zoom = .8)

# Interaction: Use mouse to rotate the 3D plot
rgl.close()

##################### dyno Multi-method Pseudotime Analysis #############################
devtools::install_github("dynverse/dyno")
devtools::install_github("dynverse/dynmethods")
library(dyno)
library(tidyverse)
library(Matrix)
library(Seurat)

# Add raw counts and normalized expression
# Seurat matrices require transposition so that rows are cells and columns are genes
dataset <- wrap_expression(
  counts = t(temp_sce[["RNA"]]$counts),
  expression = t(temp_sce[["RNA"]]$data)
)

dataset <- wrap_expression(
  counts = t(subsce[["RNA"]]$counts),
  expression = t(subsce[["RNA"]]$data)
)


# Add prior information; the starting "cell ids" are specified here. 
# Visualizations can be adjusted later based on specific trajectory inference results.
dataset <- add_prior_information(
  dataset,
  start_id = "0day_AAACCTGAGTTATCGC-1"
)  # Fibroblast sub-population
dataset <- add_prior_information(
  dataset,
  start_id = "0day_AAAGATGCATTAGCCA-1"
) # Monocyte/Macrophage sub-population
dataset <- add_prior_information(
  dataset,
  start_id = "0day_ACGATACAGGGTCGAT-1"
) # Neutrophil

# Add cluster information to the dataset, utilizing "seurat_clusters" directly
dataset <- add_grouping(
  dataset,
  temp_sce$seurat_clusters
)
dataset <- add_grouping(
  dataset,
  subsce$seurat_clusters
)

## dynwrap::test_docker_installation(detailed = TRUE) Verify if Docker installation is functional
## Determine optimal parameters and methods for the dataset; exploring multiple methods is recommended
guidelines <- guidelines_shiny(dataset)
methods_selected <- guidelines$methods_selected

## Selected "paga_tree" for trajectory inference
model_paga <- infer_trajectory(dataset, methods_selected[1])
model_paga <- model_paga %>% add_dimred(dyndimred::dimred_mds, expression_source = dataset$expression)

## Preliminary visualization with various tunable parameters
plot_dimred(
  model_paga, 
  expression_source = dataset$expression, 
  grouping = dataset$grouping
)
system()

## Add rooting gene lists to define the biological differentiation directionality
model <- model_paga %>% 
  add_root_using_expression(c("Angptl4","Ccnb2","Acta2","Sfrp2"), dataset$expression)

model <- model_paga %>% 
  add_root_using_expression(c("Ly6c2","Il1b","Stmn1","Trem2","Col5a2"), dataset$expression)

### Annotate trajectory milestones using canonical marker genes
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

## Trajectory Visualization
model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = dataset$expression)
plot_dimred(
  model, 
  color_cells = "pseudotime",
  #color_cells = "grouping",
  expression_source = dataset$expression, 
  grouping = dataset$grouping,
  label_milestones = TRUE
)


### Display gene expression dynamics along the trajectory
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
  color_cells = "grouping",
  #feature_oi = "Ly6c2",
  color_density = "grouping",
  grouping = dataset$grouping,
  label_milestones = TRUE
)


# Global heatmap overview of the top predictive genes
plot_heatmap(
  model,
  expression_source = dataset$expression,
  grouping = dataset$grouping,
  features_oi = 50
)


## Identification of critical genes driving trajectory bifurcation points
branching_milestone <- model$milestone_network %>% group_by(from) %>% filter(n() > 1) %>% pull(from) %>% first()

branch_feature_importance <- calculate_branching_point_feature_importance(model, expression_source=dataset$expression, milestones_oi = branching_milestone)

branching_point_features <- branch_feature_importance %>% top_n(20, importance) %>% pull(feature_id)

plot_heatmap(
  model,
  expression_source = dataset$expression,
  features_oi = branching_point_features
)

# Spatial mapping of the top 20 genes significant to branching points
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
## Plot scatter plot
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

# Plot smooth scatter plot
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
# cut.preselect = 0.01: Runs too slowly, increase threshold
set.seed(2020)
pb <- txtProgressBar(min = 0, max = 100, style = 3)
testres <- optimize.sd_selection(logmat[,unlist(samplesL)], samplesL, B=100, cutoff=cut.preselect, times=.75, percent=0.8)
#save(testres, file=paste0("Fibroblasts_biotip.RData"), compress=TRUE)
#save(testres, file=paste0("macrophages_biotip.RData"), compress=TRUE)
#save(testres, file=paste0("Granulocytes_biotip.RData"), compress=TRUE)
#save(testres, file=paste0("macrophages_biotip2.0.RData"), compress=TRUE)
#save(testres, file=paste0("valid_Fibroblasts_biotip.RData"), compress=TRUE)
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
cut.fdr = 0.05
cut.fdr = 0.10
cut.fdr = 0.12
cut.fdr = 0.20
igraphL <- getNetwork(testres, fdr = cut.fdr)
cluster <- getCluster_methods(igraphL)

## Network partition using random walk
dev.off()
par(mfrow=c(4,4))
cluster = getCluster_methods(igraphL)
#i = 11
#i = 12
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
nodes_in_community <- V(tmp)$name[V(tmp)$community ==2]
# Create a subgraph containing only target community nodes
subgraph <- induced_subgraph(tmp, nodes_in_community)
plot(subgraph,vertex.label = V(subgraph)$name,vertex.color=colrs[V(subgraph)$community],vertex.size = 6)
# Calculate node degrees
node_degrees <- degree(subgraph)
# Determine minimum and maximum node degrees
min_degree <- min(node_degrees)
max_degree <- max(node_degrees)
# Map node degrees to cold-warm color spectrum
color_palette <- colorRampPalette(c("#8FD2E6", "white","#ED884C" ))
color_palette <- colorRampPalette(c( "white","#ED884C" ))
node_colors <- color_palette(100)[as.integer(100 * (node_degrees - min_degree) / (max_degree - min_degree)) + 1]

# Map expression values to cold-warm color spectrum
local_sce <- subset(temp_sce ,subset = seurat_clusters == "10")
#local_sce <- subset(subsce ,subset = seurat_clusters == "5")
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
# Use Fruchterman-Reingold algorithm for layout
layout <- layout_with_fr(subgraph)
# Plot graph; adjust node colors based on degree
plot(subgraph,
     vertex.label = V(subgraph)$name,
     vertex.color = node_colors,
     vertex.size = node_size*1.2,
     #vertex.label.dist = 1.5,
     vertex.label.cex = 0.66,
     layout = layout,
     vertex.label.color = "black")
# Save edge data
edges <- as_edgelist(subgraph)
# Fibroblast population
write.table(edges, "7edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "7gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
write.table(edges, "4edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "4gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
# Fibroblast population - new 4 and 6 
write.table(edges, "fibroblast_4edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "fibroblast_4gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
write.table(edges, "fibroblast_6edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "fibroblast_6gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)

# Validate Fibroblast population - 10 
write.table(edges, "fibroblast_10edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "fibroblast_10gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)



# Neutrophil population
write.table(edges, "2edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "2gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)

# Monocyte/Macrophage population
write.table(edges, "macrophages_4edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "macrophages_4gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)
# Monocyte/Macrophage population
write.table(edges, "macrophages_8edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "macrophages_8gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)

write.table(edges, "macrophages_5edges.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)
write.table(as.data.frame(gene_sd), "macrophages_5gene_sd.tsv", sep = "\t", quote = T, row.names = T, col.names = T)






# Map PCC correlation to edge thickness

# Putative Critical Transition Signals (CTSs) Identification by MCI score
membersL <- getMCI(cluster, testres, adjust.size = FALSE, fun='BioTIP')
# Plot DNB scores for modules in each cell cluster
dev.off()
cut.minsize = 60
pdf('../GSE163129_RAW/验证图/Figure 3-5.pdf',height = 12,width = 20,onefile = F)
par(oma = c(0, 0, 0, 0))  # Adjust outer margins to 0
par(mar = c(1, 1, 1, 1))  # Adjust inner margins to 1
plotBar_MCI(membersL, ylim=c(0,6), minsize = cut.minsize)
dev.off()
# Identify the highest-scoring modules
topMCI = getTopMCI(membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]], min=cut.minsize, n=5)
maxMCIms <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], min=cut.minsize, n=5)
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
write.csv(df, file = "Granulocytes_DNB.csv", row.names = TRUE)
write.csv(df, file = "macropahges_DNB.csv", row.names = TRUE)
write.csv(df, file = "macropahges_DNB2.0.csv", row.names = TRUE)
write.csv(df, file = "fibroblast_DNB3.0.csv", row.names = TRUE)
names(CTS.Lib.Symbol)
# Shrink matrix and calculate IC.shrink.score
M <- cor.shrink(logmat[,unlist(samplesL)], Y = NULL, MARGIN = 1, shrink = TRUE)
# Analyze if DNB score is significant
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
pdf('../GSE163129_RAW/验证图/Figure3-6.pdf',height = 12,width = 10,onefile = F)
dev.off()
par(mfrow=c(5,1))
for (i in 1:5){
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las=2,
                      main=paste0("Cluster ", names(maxMCI)[i], "; ",
                                  length(CTS.Lib[[i]]), " genes", "\n","vs. ",
                                  "100 times of gene-permutation"),
                      which2point=names(maxMCI)[i])
}

# Identify candidate tipping points via Ic index and Delta-Ic
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
ylim = 1.5
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
ylim = 1.5
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

pdf('../GSE163129_RAW/验证图/Figure3-4.pdf',height = 12,width = 10,onefile = F)
dev.off()
ylim = 0.5
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
ylim = 0.5
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
############################### DNB Expression Dynamics and Clustering Analysis ################################

### Function: Calculate Mean Expression for DNB Modules
i=2
h=1
Idents(temp_sce2) <- "label"

DNBs_expression <- function(expression_matrix, target_cluster, DNB_names) {
  failed_genes <- c()
  # Remove redundant or missing identifiers
  DNB_names <- unique(na.omit(DNB_names))
  average_expression <- c()
  df <- data.frame(matrix(nrow = length(target_cluster), ncol = length(DNB_names)))
  rownames(df) <- target_cluster
  colnames(df) <- DNB_names

  # Iterate through targeted clusters and gene sets
  for (i in seq_along(target_cluster)) {
    cell_indices <- which(temp_sce2$label == target_cluster[i])
    for (h in 1:length(DNB_names)) {
      # Verify gene presence in the current expression matrix
      if (DNB_names[h] %in% rownames(expression_matrix)) {
        average_expression <- mean(expression_matrix[DNB_names[h], cell_indices])
      } else {
        failed_genes <<- c(failed_genes, DNB_names[h])
        average_expression <- NA # Return NA if calculation fails
      }
      df[i, DNB_names[h]] <- average_expression
    }
    if (length(failed_genes) > 0) {
      print(paste("Failed to map genes:", paste(unique(failed_genes), collapse = ", ")))
      cat("\n")
    }
  }
  return(df)
}

# Ingest DNB lists from previous BioTIP results
DNB_list <- as.list(read.csv("./fibrobalst_DNB.csv"))
DNB_list <- as.list(read.csv("./Granulocytes_DNB.csv"))
DNB_list <- as.list(read.csv("./macropahges_DNB.csv"))
DNB_list <- as.list(read.csv("./macropahges_DNB2.0.csv"))
DNB_list <- as.list(read.csv("./fibroblast_DNB2.0.csv"))

# Verify object identity and metadata slots
Idents(temp_sce2) 
Idents(subsce) 

# Coordinate matrix extraction: scaled and normalized slots
expression_matrix = temp_sce2[["RNA"]]$scale.data
expression_matrix = t(scale(t(temp_sce2[["RNA"]]$counts)))
expression_matrix = t(scale(t(temp_sce2[["RNA"]]$data)))

# Targeted Extraction: Fibroblast DNB (Early State)
target_cluster <- c("1", "4", "9") 
DNBs_expression.df <- DNBs_expression(expression_matrix, target_cluster, unique(DNB_list$X4))
exps <- t(DNBs_expression.df)
write.table(exps, "fibroblast_4exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Targeted Extraction: Fibroblast DNB (Late State)
target_cluster <- c("9", "6", "7") 
DNBs_expression.df <- DNBs_expression(expression_matrix, target_cluster, unique(DNB_list$X6))
exps <- t(DNBs_expression.df)
write.table(exps, "fibroblast_6exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Targeted Extraction: Fibroblast DNB (Complex Lineage)
target_cluster <- c("3", "10", "4") 
DNBs_expression.df <- DNBs_expression(expression_matrix, target_cluster, unique(DNB_list$X10))
exps <- t(DNBs_expression.df)
write.table(exps, "fibroblast_3104exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Targeted Extraction: Granulocyte DNB (Early Stage)
target_cluster <- c("0", "2", "1") 
DNBs_expression.df <- DNBs_expression(expression_matrix, target_cluster, unique(DNB_list$X2))
exps <- t(DNBs_expression.df)
write.table(exps, "2exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Targeted Extraction: Macrophage DNB (Late Stage)
target_cluster <- c("3", "4", "1") 
DNBs_expression.df <- DNBs_expression(expression_matrix, target_cluster, unique(DNB_list$X4))
exps <- t(DNBs_expression.df)
write.table(exps, "macrophages_341exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Targeted Extraction: Macrophage DNB (Early Stage Axis 5-8-3)
target_cluster <- c("5", "8", "3") 
DNBs_expression.df <- DNBs_expression(expression_matrix, target_cluster, unique(DNB_list$X8))
exps <- t(DNBs_expression.df)
write.table(exps, "macrophages_583exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

# Targeted Extraction: Macrophage DNB (Early Stage Axis 12-5-8)
target_cluster <- c("12", "5", "8") 
DNBs_expression.df <- DNBs_expression(expression_matrix, target_cluster, unique(DNB_list$X5))
exps <- t(DNBs_expression.df)
write.table(exps, "macrophages_1258exps.tsv", sep = "\t", quote = FALSE, row.names = T, col.names = T)

## Functional Clustering using ClusterGVis
library(ClusterGVis)
library(sf)

exps <- na.omit(as.matrix(exps))

# Filter top biomarkers across temporal clusters
sce.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top10
markers_tempsce %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top10
temp_sce2.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top10

# Aggregate cluster-level means for high-resolution plotting
st.data1 <- prepareDataFromscRNA(object = temp_sce, diffData = top10, showAverage = TRUE)
ck <- clusterData(exp = exps, cluster.method = "kmeans", cluster.num = 4)

# Plot heatmap with annotated representative genes
markGenes = rownames(exps)[sample(1:nrow(exps), 32, replace = FALSE)]
pdf('addgene.pdf', height = 10, width = 6, onefile = FALSE)
visCluster(object = ck, plot.type = "heatmap", column_names_rot = 45, markGenes = markGenes)
dev.off()

# Biological Process (BP) Gene Ontology Enrichment
library(org.Mm.eg.db)
enrich <- enrichCluster(object = st.data1, OrgDb = org.Mm.eg.db, type = "BP", 
                        pvalueCutoff = 0.05, topn = 5, seed = 5201314)

# Clean result object for dual-visualization (Heatmap + Enrichment)
enrich <- subset(enrich, select = -geneID)
pdf('term.pdf', height = 14, width = 12, onefile = FALSE)
visCluster(object = st.data1, plot.type = "both", column_names_rot = 45, show_row_dend = FALSE, 
           markGenes = top10$gene, markGenes.side = "left", genes.gp = c('italic', fontsize = 12, col = "black"), 
           annoTerm.data = enrich, line.side = "left", go.size = 8, cluster.order = c(1:15))
dev.off()

# Trajectory Trend Analysis
visCluster(object = ck, plot.type = "line", ms.col = c("green", "orange", "red"), add.mline = TRUE)

# Categorize temporal expression patterns per gene
trends <- apply(exps, 1, function(row_vec) {
  diff_vec <- diff(row_vec)
  if (all(diff_vec > 0)) {
    return("Gradually Increasing")
  } else if (all(diff_vec < 0)) {
    return("Gradually Decreasing")
  } else if (diff_vec[1] > 0 && diff_vec[2] < 0) {
    return("Critically High")
  } else if (diff_vec[1] < 0 && diff_vec[2] > 0) {
    return("Critically Low")
  } else {
    return("Unclassified")
  }
})

# Reorder matrix and labels based on trajectory classification
sorted_idx <- order(trends)
exps <- exps[sorted_idx, ]
sorted_trends <- as.data.frame(trends[sorted_idx])
colnames(sorted_trends) <- c("Category")

# Pheatmap: Integrated Trend Visualization
heatmap_params <- list(
  col = colorRampPalette(c("blue", "white", "red"))(100),
  scale = "row",
  main = "DNB Expression Heatmap"
)

exps_merged <- merge(exps, sorted_trends, by = "row.names")
rownames(exps_merged) <- exps_merged[,1]
exps_merged <- exps_merged[,-1]
exps_merged %>% group_by(Category) %>% 
  pheatmap(exps, cluster_rows = FALSE, cluster_cols = FALSE, clustering_distance_cols = "euclidean", 
           clustering_method = "complete", annotation_row = sorted_trends, 
           col = heatmap_params$col, scale = "row", main = heatmap_params$main)

# ComplexHeatmap: Advanced Visualization with Splitting
library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(-2, 0, 2), c("#4169E1", "white", "#DC143C"))
row_annotation = rowAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#DC143C", "#4169E1", "#FF4500", "#006400")),
  labels = c("Increasing", "Decreasing", "Critical-High", "Critical-Low"),
  labels_gp = gpar(col = "white", fontsize = 6)))

pdf('heatmap3.pdf', height = 16, width = 20, onefile = FALSE)
Heatmap(t(scale(t(exps))), col = col_fun, cluster_columns = FALSE, cluster_rows = TRUE, 
        row_split = sorted_trends$Category, show_heatmap_legend = TRUE, border = TRUE, 
        show_column_names = TRUE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), 
        heatmap_width = unit(1, "npc"), heatmap_height = unit(1, "npc"))
dev.off()

############################### Cross-Dataset Conservation (Jaccard Index) ################################

# Ingest marker lists from discovery (ceshi) and validation sets
ceshi_markers <- read.csv(file = "../E-MTAB-7895/macrophagesmarkers.csv", row.names = 1)
ceshi_markers <- read.csv(file = "../E-MTAB-7895/Granulocytesmarkers.csv", row.names = 1)
ceshi_markers <- read.csv(file = "../E-MTAB-7895/newfibroblastmarkers.csv", row.names = 1)

ceshi_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> ceshi_top100

valid_markers <- read.csv(file = "../GSE163129_RAW/test_Neutrophil_markers.csv", row.names = 1)
valid_markers <- read.csv(file = "./valid_fibroblastmarkers.csv", row.names = 1)

valid_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> valid_top100

# Stringent filtering: Extract active features (|log2FC| > 0.5)
ceshi_top100 <- ceshi_markers %>% group_by(cluster) %>% filter(abs(avg_log2FC) > 0.5) %>% arrange(desc(abs(avg_log2FC)))
valid_top100 <- valid_markers %>% group_by(cluster) %>% filter(abs(avg_log2FC) > 0.5) %>% arrange(desc(abs(avg_log2FC)))

ceshi_clusters <- unique(ceshi_top100$cluster)
valid_clusters <- unique(valid_top100$cluster)

# Initialize Jaccard Similarity Matrix
jaccard_matrix <- matrix(nrow = length(ceshi_clusters), ncol = length(valid_clusters))
colnames(jaccard_matrix) <- valid_clusters
rownames(jaccard_matrix) <- ceshi_clusters

# Compute intersection over union for cluster marker conservation
for (i in seq_along(ceshi_clusters)) {
  ceshi_cluster_genes <- ceshi_top100$gene[ceshi_top100$cluster == ceshi_clusters[i]]
  for (j in seq_along(valid_clusters)) {
    valid_cluster_genes <- valid_top100$gene[valid_top100$cluster == valid_clusters[j]]
    intersection <- length(intersect(ceshi_cluster_genes, valid_cluster_genes))
    union_set <- length(union(ceshi_cluster_genes, valid_cluster_genes))
    jaccard_matrix[i, j] <- intersection / union_set
  }
}

# Visualize cross-dataset conservation metrics
dev.off()
pheatmap(jaccard_matrix, main = "Jaccard Similarity Heatmap: Marker Conservation", 
         cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("lightblue", "white", "red"))(100), 
         fontsize_row = 20, fontsize_col = 20, display_numbers = TRUE, number_color = "black")

# Significance Testing: Fisher's Exact Test for Overlap
fisher_test <- function(x, y) {
  a <- length(intersect(x, y))
  b <- length(x) - a
  c <- length(y) - a
  d <- max(0, 100 - a - b - c) # Assume a local universe of top 100 markers
  p_value <- tryCatch({ fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value }, error = function(e) { NA })
  return(p_value)
}

fisher_matrix <- matrix(nrow = length(ceshi_clusters), ncol = length(valid_clusters), dimnames = list(ceshi_clusters, valid_clusters))

for (i in seq_along(ceshi_clusters)) {
  ceshi_cluster_genes <- ceshi_top100$gene[ceshi_top100$cluster == ceshi_clusters[i]]
  for (j in seq_along(valid_clusters)) {
    valid_cluster_genes <- valid_top100$gene[valid_top100$cluster == valid_clusters[j]]
    fisher_matrix[i, j] <- fisher_test(ceshi_cluster_genes, valid_cluster_genes)
  }
}

# Meta-analysis: Merging Jaccard indices and P-values
data_grid <- expand.grid(row = rownames(jaccard_matrix), col = colnames(jaccard_matrix))
data_grid$jaccard <- as.vector(jaccard_matrix)
data_grid$p_value <- as.vector(fisher_matrix)

# Render Quantitative Marker Conservation Grid
library(ggplot2)
ggplot(data_grid, aes(x = col, y = row, fill = jaccard)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p_value)), color = "white", size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
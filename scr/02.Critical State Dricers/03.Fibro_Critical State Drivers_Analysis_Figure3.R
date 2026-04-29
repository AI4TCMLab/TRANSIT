# ==================================================================================================
# MODULE: Figure 3 Global Trajectory, Tipping Point Networks & Doublet Filtering
# 
# DESCRIPTION: 
#   Integrates DoubletFinder screening, global cell type annotations, Monocle 2 
#   pseudotime trajectories, and complex dynamic network biomarker (DNB) topologies.
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# #### Figure 3.A: Environment Initialization & DoubletFinder Workflow ####
# --------------------------------------------------------------------------------------------------
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
options("repos"     = c(CRAN = "https://mirrors.westlake.edu.cn/CRAN/"))

if (!require("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("mengxu98/scop")
library(scop)

# Initialize and subset fibroblast specific populations
temp_sce <- readRDS("/home/tsh/E-MTAB-7895/temp_sce.rds") 
temp_sce <- subset(temp_sce, subset = orig.ident %in% c("0day", "1day", "3day", "5day", "7day"))
table(temp_sce$orig.ident)

# [Historical exploratory commands retained]
object$seurat_clusters
object$

# --------------------------------------------------------------------------------------------------
# #### DoubletFinder Implementation ####
# --------------------------------------------------------------------------------------------------
library(DoubletFinder)
library(Seurat)
library(tidyverse)
library(dplyr)

# Dynamic object targeting
temp_sce <- readRDS("/home/tsh/E-MTAB-7895/temp_sce.rds")                            # Fibroblasts
temp_sce <- readRDS("/home/tsh/E-MTAB-7895/filtered_inflammation_temp_sce.rds")      # Macrophages
object   <- subset(temp_sce, subset = orig.ident %in% c("0day", "1day", "3day", "5day", "7day"))

# 1. Partition object by sample source
list_sce <- SplitObject(object, split.by = "orig.ident")

# 2. Iterative DoubletFinder processing across temporal samples
processed_list <- lapply(list_sce, function(seu_temp) {
  
  s_name  <- unique(seu_temp$orig.ident)
  n_cells <- ncol(seu_temp)
  
  # --- Strategy: Bypass samples with insufficient cellular coverage ---
  if (n_cells < 50) {
    message("Skipping sample ", s_name, ": Cell count is only ", n_cells)
    return(NULL)
  }
  
  # --- Foundational Preprocessing ---
  seu_temp <- NormalizeData(seu_temp, verbose = FALSE) %>% 
              FindVariableFeatures(verbose = FALSE) %>% 
              ScaleData(verbose = FALSE)
  
  # Dynamically allocate Principal Components to circumvent boundary errors
  current_pcs <- min(15, n_cells - 1)
  seu_temp    <- RunPCA(seu_temp, npcs = current_pcs, verbose = FALSE)
  
  # --- Core DoubletFinder Hyperparameter Optimization ---
  sweep.res   <- paramSweep(seu_temp, PCs = 1:current_pcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn       <- find.pK(sweep.stats)
  opt_pK      <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # Establish homotypic proportion corrections if major_type constraints exist
  if ("seurat_clusters" %in% colnames(seu_temp@meta.data)) {
    homo_prop <- modelHomotypic(seu_temp$seurat_clusters)
  } else {
    homo_prop <- 0 # Default to uncorrected if cluster identities are missing
  }
  
  # Estimate expected doublet rates (baseline set to 5.0%, scalable based on loading density)
  nExp     <- round(0.050 * n_cells) 
  nExp_adj <- round(nExp * (1 - homo_prop))
  
  # Execute doublet classification manifold
  seu_temp <- doubletFinder(seu_temp, 
                            PCs  = 1:current_pcs, 
                            pN   = 0.25, 
                            pK   = opt_pK, 
                            nExp = nExp_adj, 
                            sct  = FALSE)
  
  # --- Standardize Output Nomenclatures ---
  # Harmonize dynamic column suffixes to enable downstream object merging
  meta <- seu_temp@meta.data
  meta$pANN_final <- meta[[grep("^pANN_", colnames(meta), value = TRUE)]]
  meta$DF_final   <- meta[[grep("^DF.classifications_", colnames(meta), value = TRUE)]]
  
  seu_temp@meta.data <- meta
  return(seu_temp)
})

# 3. Reintegrate processed metrics into the global object framework
processed_list   <- Filter(Negate(is.null), processed_list)
all_doublet_meta <- bind_rows(lapply(processed_list, function(x) x@meta.data))

# Append harmonized metadata
object <- AddMetaData(object, all_doublet_meta[, c("pANN_final", "DF_final")])

# Extract plotting coordinate frames
plot_df <- object@meta.data %>%
           select(seurat_clusters, pANN = pANN_final, DF = DF_final)

# Derive summary statistics (Mean and Standard Error of the Mean [SEM]) for Panel D
summary_d <- plot_df %>%
             group_by(seurat_clusters) %>%
             summarise(
               mean_pANN = mean(as.numeric(pANN), na.rm = TRUE),
               sd_pANN   = sd(as.numeric(pANN), na.rm = TRUE),
               n         = n()
             ) %>%
             mutate(se_pANN = sd_pANN / sqrt(n))

# Derive compositional frequencies for Panel E
summary_e <- plot_df %>%
             filter(!is.na(DF)) %>%  # Core exclusion: Remove unclassified cells from bypassed sparse samples
             group_by(seurat_clusters, DF) %>%
             tally() %>%
             mutate(freq = n / sum(n))

p1 <- ggplot(summary_d, aes(x = seurat_clusters, y = mean_pANN)) +
      geom_bar(stat = "identity", fill = "grey40", color = "black", width = 0.7) +
      geom_errorbar(aes(ymin = mean_pANN - se_pANN, ymax = mean_pANN + se_pANN), width = 0.2) +
      theme_classic() +
      labs(y = "Average pANN", x = "", title = "Average pANN by Cell Type") +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title   = element_text(hjust = 0.5)
      )

# Enforce factor hierarchy for consistent stacked presentation
summary_e$DF <- factor(summary_e$DF, levels = c("Doublet", "Singlet"))

p2 <- ggplot(summary_e, aes(x = seurat_clusters, y = freq, fill = DF)) +
      geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7) +
      scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
      scale_fill_manual(values = c("Doublet" = "#E41A1C", "Singlet" = "#D9D9D9")) +
      theme_classic() +
      labs(y = "Frequency (%)", x = "", fill = "Classification") +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")
      )

print(p1 + p2)

# --------------------------------------------------------------------------------------------------
# #### Cluster Identification & Marker Projection ####
# --------------------------------------------------------------------------------------------------
DimPlot(temp_sce, reduction = "umap", split.by = c("orig.ident"), label = TRUE)
DimPlot(temp_sce, reduction = "umap", group.by = c("seurat_clusters"), label = TRUE)
factor(temp_sce@meta.data$seurat_clusters)

# Assign functional overarching identities (major_type) based on transcriptomic profiles
temp_sce$major_type <- dplyr::case_when(
  temp_sce$seurat_clusters %in% c("1", "4", "14", "12", "11") ~ "IR",
  temp_sce$seurat_clusters %in% c("5", "8", "10", "13")       ~ "F-cyc",
  temp_sce$seurat_clusters %in% c("9")                        ~ "F-CI",
  temp_sce$seurat_clusters %in% c("3", "6", "7", "0", "2")    ~ "MYO",
  TRUE                                                        ~ "Other"
) |> factor(levels = c("IR", "F-CI", "F-cyc", "MYO"))

# Establish standardized cell lineage color arrays
cell_colors <- c(
  "IR"    = "#CC79A7",   # Sky blue (Note: Hex represents a pink/purple tone historically)
  "F-cyc" = "#E69F00",   # Orange
  "MYO"   = "#009E73",   # Bluish green
  "F-CI"  = "lightblue"  # Blue
)

CellDimPlot(
  srt       = temp_sce,
  group.by  = c("seurat_clusters"),
  reduction = "umap",
  theme_use = "theme_blank"
)

FeatureDimPlot(
  srt              = temp_sce,
  features         = c("Atf4", "Ccne2", "Mki67", "Cdh2"),
  compare_features = TRUE,
  label            = TRUE,
  label_insitu     = TRUE,
  reduction        = "umap",
  theme_use        = "theme_blank",
  pt.size          = 2
)

# Extract differential signatures
macro_markers <- FindAllMarkers(temp_sce, group.by = "seurat_clusters")
macro_markers <- FindAllMarkers(temp_sce, group.by = "major_type")
write.csv(macro_markers, "./大论文图和数据/Figure 2/Fibro_major_type.csv")

# Isolate top driving markers
markers <- macro_markers %>%
           dplyr::group_by(cluster) %>%
           dplyr::top_n(n = 3, wt = avg_log2FC)

marers_genes <- unique(markers$gene)
marers_genes <- c(
  "Dbndd1", "Alx3", "Tnnc2",    # Cluster 0
  "Pi16", "Atf4", "Angptl4",    # Cluster 1
  "Cilp", "C1qtnf3", "Wnt2",    # Cluster 2
  "Shisa6", "Wif1", "Cxcl12",   # Cluster 3
  "Cxcl3", "Mt2", "Nfkb1",      # Cluster 4
  "E2f2", "Stmn1", "Fbxo5",     # Cluster 5
  "Runx1", "Acta2", "Cdh2",     # Cluster 6
  "Eln", "Sfrp2", "Vegfd",      # Cluster 7
  "Serpinb2", "Mki67", "Ccnb1", # Cluster 8
  "Ccne2", "Serpinb2", "Dkk2",  # Cluster 9
  "Ccdc39", "Cenpa", "Hmmr",    # Cluster 10
  "Muc13", "Egr4", "Lsmem1",    # Cluster 11
  "Egln3", "Hmox1", "Tafa2",    # Cluster 12
  "Gipc3", "Vwf", "Depp1",      # Cluster 13
  "Arl5c", "Rpgrip1", "Gsdme"   # Cluster 14
)

# Render customized GroupHeatmap via scop
ht <- GroupHeatmap(
  srt                     = temp_sce,
  features                = marers_genes,
  group.by                = c("seurat_clusters"),
  heatmap_palette         = "RdBu",
  # cell_annotation       = c("Phase", "G2M_score", "Cdh2"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names          = TRUE, 
  row_names_side          = "left",
  add_dot                 = FALSE, 
  add_reticle             = FALSE
)
print(ht$plot)

# --------------------------------------------------------------------------------------------------
# #### ClusterGVis Enrichment Heatmaps ####
# --------------------------------------------------------------------------------------------------
macro_markers <- FindAllMarkers(temp_sce, group.by = "seurat_clusters")
markers       <- macro_markers %>%
                 dplyr::group_by(cluster) %>%
                 dplyr::top_n(n = 20, wt = avg_log2FC)

head(markers)
markGenes <- unique(markers$gene)[sample(1:length(unique(markers$gene)), 45, replace = FALSE)]

markGenes <- c(
  "Prr18", "4921513I03Rik", "Gm13067", "Qprt", "Cxcl5", "Hcn1", "Lhx6",
  "9330175M20Rik", "Gm12059", "Lvrn", "Cdca2", "Gm11638", "Gm38973",
  "Il11", "Thbs4", "E2f2", "Gm49492", "Scube2", "Wfdc21", "Hist1h2ab",
  "Adh1", "Gm16153", "Ahsg", "Hmx3", "2810429I04Rik", "Wnt10b", "Frmd3",
  "Gm42664", "Moxd1", "Zfp185", "Zbtb16", "Zfp536", "F2rl3", "Emx2",
  "Esco2", "Gm21123", "Gipc3", "Nusap1", "Alx3", "Nr1i3", "Hist1h4j",
  "Kif14", "Igfbp3", "Clec1a", "Galnt15"
)

library(ClusterGVis)
Idents(temp_sce) <- "major_type"

data <- prepareDataFromscRNA(object        = temp_sce,
                             diffData      = markers,
                             keep.uniqGene = FALSE,  
                             showAverage   = TRUE)

enrich <- enrichCluster(object       = data,
                        OrgDb        = org.Mm.eg.db,
                        type         = "BP",
                        organism     = "mus",
                        pvalueCutoff = 0.05,
                        topn         = 5,
                        seed         = 123)

head(enrich)

pdf("FigureS1.Fibro_markers_heatmap.pdf", width = 16, height = 10)
visCluster(
  object           = data, 
  plot.type        = "both",
  column_names_rot = 45,
  show_row_dend    = FALSE,
  markGenes        = markGenes,
  markGenes.side   = "left",
  annoTerm.data    = enrich,
  line.side        = "left",
  # cluster.order  = c("IR", "Pro_inflammation", "Proliferation", "Reparative", "Antigen-presenting", "Collagens"),
  # Assuming standard categorical lengths
  go.col           = rep(jjAnno::useMyCol("stallion", n = 4), each = 5)[1:20],
  add.bar          = TRUE
)
dev.off()

# Secondary exploratory mappings
Idents(temp_sce) <- "seurat_clusters"
FeaturePlot(temp_sce, features = c("Ccn5", "Sfrp2", "Gucy1a1", "Fbln1", "Serpine2", "Pi16"), reduction = "umap")
FeaturePlot(temp_sce, features = c("Tgfb1", "Thbs4", "Crlf1", "Col15a1", "Cthrc1"), reduction = "umap")
VlnPlot(temp_sce, features = c("Ccn5", "Sfrp2", "Gucy1a1", "Fbln1", "Postn"))
DotPlot(temp_sce, features = c("Ccn5", "Sfrp2", "Gucy1a1", "Fbln1"))
markers <- FindAllMarkers(temp_sce)

# --------------------------------------------------------------------------------------------------
# #### Monocle 2 Trajectory Inference Initiation ####
# --------------------------------------------------------------------------------------------------
temp_sce <- subset(temp_sce, subset = seurat_clusters %in% c("1", "4", "9", "6", "7", "0"))
factor(temp_sce@meta.data$seurat_clusters)
stopifnot("seurat_clusters" %in% colnames(temp_sce@meta.data))

library(monocle)
Idents(temp_sce) <- "seurat_clusters"

# Convert standard sparse matrix formats
sparse_data <- as(as.matrix(temp_sce[["RNA"]]$counts), "sparseMatrix")
mdata       <- new("AnnotatedDataFrame", data = temp_sce@meta.data)
fData       <- data.frame(gene_short_name = row.names(sparse_data), row.names = row.names(sparse_data))
fd          <- new("AnnotatedDataFrame", data = fData)

monocle_cds <- newCellDataSet(cellData            = sparse_data,
                              phenoData           = mdata,
                              featureData         = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily    = negbinomial.size())

# Parameter estimation and thresholding
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 3)

# Integrate native Seurat variable features
expressed_genes <- VariableFeatures(temp_sce)
diff_test_res   <- differentialGeneTest(monocle_cds[expressed_genes, ],
                                        fullModelFormulaStr = "~seurat_clusters")

# Enforce strict significance thresholding
ordering_genes  <- row.names(subset(diff_test_res, qval < 0.01)) 
monocle_cds     <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)

# save(monocle_cds, file = "./大论文图和数据/Figure 2/Fibro_monocle2.Rdata", compress = TRUE)

library(igraph)
# Dimensionality reduction via discriminative dimensionality reduction with trees (DDRTree)
monocle_cds <- monocle::reduceDimension(monocle_cds, max_components = 2, norm_method = "log", reduction_method = "DDRTree")
monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds, root_state = "2")
monocle_cds <- orderCells(monocle_cds, root_state = "1")

# [Historical syntax artifact retained]
??

# Render trajectory structures
p <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 0.6, show_backbone = TRUE, show_branch_points = FALSE) + 
     facet_wrap("~orig.ident", nrow = 1)

p0 <- plot_cell_trajectory(monocle_cds, color_by = "major_type", cell_size = 0.5, show_backbone = TRUE, show_branch_points = TRUE) + 
      facet_wrap("~orig.ident", nrow = 1)

p1 <- plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters", cell_size = 0.6, show_backbone = FALSE, show_branch_points = FALSE) +
      theme(legend.position = 'none', panel.border = element_blank()) +
      scale_color_manual(values = color)

plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 0.6, show_backbone = TRUE) +
  theme(legend.position = 'none', panel.border = element_blank())

p2 <- plot_complex_cell_trajectory(monocle_cds, x = 1, y = 2, color_by = "label") +
      theme(legend.title = element_blank())

plot_cell_trajectory(monocle_cds, color_by = "State", size = 1, show_backbone = TRUE)
p0 | p | p1 | p2

plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters", size = 1, show_backbone = TRUE) + facet_wrap("~orig.ident", nrow = 1)

plot_complex_cell_trajectory(monocle_cds, x = 1, y = 2, color_by = "seurat_clusters") +
  theme(legend.title = element_blank())


# ==================================================================================================
# #### Figure 3.B: Temporal Area Proportions ####
# ==================================================================================================
Idents(temp_sce) <- "seurat_clusters"

# Compute proportionate and absolute abundances
cell_data <- as.data.frame(prop.table(table(Idents(temp_sce), temp_sce$orig.ident), margin = 2))
cell_data <- as.data.frame(table(Idents(temp_sce), temp_sce$orig.ident), margin = 2)
colnames(cell_data) <- c("Cell_Type", "Timepoint", "Abundance")

plot_abundance <- ggplot(cell_data, aes(x = Timepoint, y = Abundance, fill = Cell_Type)) +
                  geom_col() +
                  labs(x = "Timepoint", y = "Cell Abundance", fill = "Cell Type") +
                  theme_minimal() +
                  theme(
                    panel.background = element_rect(fill = "transparent", color = NA),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line        = element_line(colour = "grey")
                  ) +
                  facet_wrap(~ Cell_Type, scales = "free_x", ncol = 6) +
                  geom_area()

print(plot_abundance)


# ==================================================================================================
# #### Figure 3.C & 3.D: BioTIP Network Partitioning and Criticality ####
# ==================================================================================================
# save(testres, file=paste0("Fibroblasts_biotip.RData"), compress=TRUE)

# 1. Network Partition Initialization
subsce_trans <- as.SingleCellExperiment(temp_sce)
samplesL     <- split(rownames(colData(subsce_trans)), f = colData(subsce_trans)$seurat_clusters)

# Model High Variance Genes (Poisson assumption)
dec.pois <- modelGeneVarByPoisson(subsce_trans)
hvg      <- getTopHVGs(dec.pois, n = 4000)
hvg      <- intersect(hvg, rownames(subsce_trans))
dat      <- subsce_trans[hvg, ]
logmat   <- as.matrix(logcounts(dat))
global   <- as.matrix(logcounts(subsce_trans))

cut.fdr  <- 0.12
igraphL  <- getNetwork(testres, fdr = cut.fdr)
cluster  <- getCluster_methods(igraphL)

# 2. Subgraph Extraction and Rendering (Focus on Sub-network 2)
i   = 2
tmp = igraphL[[i]]
E(tmp)$width     <- E(tmp)$weight * 3
V(tmp)$community <- cluster[[i]]$membership
mark.groups      <- table(cluster[[i]]$membership)
colrs            <- rainbow(length(mark.groups), alpha = 0.3)
V(tmp)$label     <- NA

plot(tmp, vertex.color = colrs[V(tmp)$community], vertex.size = 5, mark.groups = cluster[[i]])
table(V(tmp)$community)
which(V(tmp)$community == 3)

plot(tmp, vertex.label = V(tmp)$name, vertex.color = colrs[V(tmp)$community], vertex.size = 5, mark.groups = cluster[[i]])

# Isolate nodes native to community 3
nodes_in_community <- V(tmp)$name[V(tmp)$community == 3]
subgraph           <- induced_subgraph(tmp, nodes_in_community)

plot(subgraph, vertex.label = V(subgraph)$name, vertex.color = colrs[V(subgraph)$community], vertex.size = 6)

# Map dynamic metrics to graph geometries (Degrees and Expressions)
node_degrees  <- degree(subgraph)
min_degree    <- min(node_degrees)
max_degree    <- max(node_degrees)

color_palette <- colorRampPalette(c("#8FD2E6", "white", "#ED884C"))
color_palette <- colorRampPalette(c("white", "#ED884C"))
node_colors   <- color_palette(100)[as.integer(100 * (node_degrees - min_degree) / (max_degree - min_degree)) + 1]

# Contextualize expression arrays against the subgraph
local_sce <- subset(temp_sce, subset = seurat_clusters == "4")
local_sce <- as.SingleCellExperiment(local_sce)
local_sce <- as.matrix(logcounts(local_sce))

gene_sd  <- apply(local_sce, 1, sd)[V(subgraph)$name]
gene_exp <- apply(local_sce, 1, mean)[V(subgraph)$name]

min_exp <- min(gene_exp)
max_exp <- max(gene_exp)

node_size   <- gene_sd * 10
node_colors <- color_palette(100)[as.integer(100 * (gene_exp - min_exp) / (max_exp - min_exp))]

# Define Force-Directed Layout
layout <- layout_with_fr(subgraph)

plot(subgraph,
     vertex.label       = V(subgraph)$name,
     vertex.color       = node_colors,
     vertex.size        = node_size * 3,
     # vertex.label.dist = 1.5,
     vertex.label.cex   = 0.8,
     layout             = layout,
     vertex.label.color = "black")

# 3. Graph Topology Exportation
edges <- as_edgelist(subgraph)

# Fibroblast Sub-lineage Architectures (Clusters 4 & 6)
write.table(edges, "fibroblast_4edges.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(as.data.frame(gene_sd), "fibroblast_4gene_sd.tsv", sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
write.table(edges, "fibroblast_6edges.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(as.data.frame(gene_sd), "fibroblast_6gene_sd.tsv", sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)


# --------------------------------------------------------------------------------------------------
# #### Critical Transition Signals (CTSs) and DNB Module Inference ####
# --------------------------------------------------------------------------------------------------
membersL <- getMCI(cluster, testres, adjust.size = FALSE, fun = 'BioTIP')

cut.minsize = 60
par(oma = c(0, 0, 0, 0)) 
par(mar = c(1, 1, 1, 1)) 
plotBar_MCI(membersL, ylim = c(0, 20), minsize = cut.minsize)
dev.off()

# Extrapolate highest scoring DNB topological modules
target_nodes <- CTS.Lib.Symbol[["4"]]
g            <- igraphL[["4"]]
subg         <- induced_subgraph(g, vids = V(g)[name %in% target_nodes])

edgelist_sub   <- as_edgelist(subg, names = TRUE)
edge_attrs_sub <- edge_attr(subg)
result_internal <- data.frame(edgelist_sub, edge_attrs_sub, check.names = FALSE)

head(result_internal)
# write.csv(result_internal,"./大论文图和数据/C1_Fibro_DNB_early_network.csv")
membersL$sd$`4`[["3"]]

# Execute robust extraction of DNB nodes
topMCI         <- getTopMCI(membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]], min = cut.minsize, n = 13)
maxMCIms       <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], min = cut.minsize, n = 13)
maxMCI         <- getMaxStats(membersL[['MCI']], maxMCIms[['idx']])
CTS.Lib        <- getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])
CTS.Lib.Symbol <- CTS.Lib

maxMCI <- sort(maxMCI, decreasing = TRUE)
maxMCI <- head(maxMCI, 13)

# Structure and export DNB targets
max_length <- max(lengths(CTS.Lib))
df <- data.frame(matrix(NA, nrow = max_length, ncol = length(CTS.Lib)))

for (i in 1:length(CTS.Lib)) {
  col_name <- names(CTS.Lib)[i]
  col_data <- CTS.Lib[[i]]
  df[, i]  <- c(col_data, rep(NA, max_length - length(col_data)))
}
colnames(df) <- names(CTS.Lib)
write.csv(df, file = "fibrobalst_DNB.csv", row.names = TRUE)

names(CTS.Lib.Symbol)
table(CTS.Lib.Symbol)

# Shrinkage estimation of covariance network
M <- cor.shrink(logmat[, unlist(samplesL)], Y = NULL, MARGIN = 1, shrink = TRUE)

# Statistical validation of MCI magnitudes via bootstrapping
C       = 1000
simuMCI = list()
set.seed(2020)

for (i in 1:length(CTS.Lib)){
  n          <- length(CTS.Lib[[i]])
  simuMCI[[i]] <- simulationMCI(n, samplesL, logmat, B = C, fun = "BioTIP", M = M)
}

dev.off()
par(mfrow = c(5, 1))
for (i in 1:length(CTS.Lib)){
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las = 2,
                      main = paste0("Cluster ", names(maxMCI)[i], "; ",
                                    length(CTS.Lib[[i]]), " genes", "\n", "vs. ",
                                    "100 times of gene-permutation"),
                      which2point = names(maxMCI)[i])
}

dev.off()
par(mfrow = c(1, 1))
for (i in 1:10){
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las = 2,
                      main = paste0("Cluster ", names(maxMCI)[i], "; ",
                                    length(CTS.Lib[[i]]), " genes", "\n", "vs. ",
                                    "1000 times of gene-permutation"),
                      which2point = names(maxMCI)[i])
}

# Sub-visualization of targeted clusters
i <- 5
sel_clusters <- c("1", "4", "9", "6", "7")
maxMCI_sub   <- maxMCI[i][sel_clusters]
simuMCI_sub  <- simuMCI[[i]][sel_clusters, ]

plot_MCI_Simulation(maxMCI[i], simuMCI_sub, las = 2,
                    main = paste0("Cluster ", names(maxMCI)[i], "; ",
                                  length(CTS.Lib[[i]]), " genes", "\n", "vs. ",
                                  "1000 times of gene-permutation"),
                    which2point = names(maxMCI)[i])
dev.off()

# Confirmatory extraction using Delta-IC dynamics
C = 1000
BioTIP_scores <- SimResults_g <- list()
set.seed(101010)

for(i in 1:length(CTS.Lib)){
  CTS                 <- CTS.Lib.Symbol[[i]]
  n                   <- length(CTS)
  BioTIP_scores[[i]]  <- getIc(logmat[, unlist(samplesL)], samplesL, CTS, fun = "BioTIP", shrink = TRUE, PCC_sample.target = 'none')
  SimResults_g[[i]]   <- simulation_Ic(n, samplesL, logmat, B = C, fun = "BioTIP", shrink = TRUE, PCC_sample.target = 'none')
}
names(BioTIP_scores) <- names(SimResults_g) <- names(CTS.Lib)

ylim = 0.4
i    = 5
par(mfrow = c(1, 2))

for(i in 1:length(BioTIP_scores)){
  n           <- length(CTS.Lib[[i]])
  interesting <- which(names(samplesL) == names(BioTIP_scores[i]))
  
  plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], las = 2, ylab = "Ic.shrink", ylim = c(0, ylim),
                     main        = paste("Cluster ", names(CTS.Lib)[i], "_", n, "genes", "\n", "vs. ", "1000 gene-permutations"),
                     fun         = "matplot", 
                     which2point = interesting)
                     
  plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]],
                     main = paste("Delta Ic*", n, "genes"), ylab = NULL,
                     xlim = range(c(BioTIP_scores[[i]][names(BioTIP_scores)[i]], SimResults_g[[i]])))
}

n <- length(CTS.Lib[[i]])

# Explicit filter for functional subpopulations of interest
sel_clusters     <- c("1", "4", "9", "6", "7")
score_sub        <- BioTIP_scores[[i]][sel_clusters]
SimResults_g_sub <- SimResults_g[[i]][sel_clusters, ]

plot_Ic_Simulation(score_sub, SimResults_g_sub,
                   las         = 2,
                   ylab        = "Ic.shrink",
                   ylim        = c(0, ylim),
                   main        = paste("Cluster ", names(CTS.Lib)[i], "_", n, "genes", "\n", "vs. ", "1000 gene-permutations"),
                   fun         = "matplot",
                   which2point = interesting)

plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]],
                   main = paste("Delta Ic*", n, "genes"), ylab = NULL,
                   xlim = range(c(BioTIP_scores[[i]][names(BioTIP_scores)[i]], SimResults_g[[i]])))
dev.off()

# Trajectory sequence plotting
x         <- c(1, 4, 9, 6, 7)
DNB_count <- c(68, 76, 121, 81, 92)
y         <- c(1.2935080, 1.5503377, 1.3433640, 1.6148139, 1.5758125)

x <- factor(x, levels = c(1, 4, 9, 6, 7))

plot(x, y, type = "line", col = "black", xlab = "Cluster", ylab = "DNB_Score", main = "DNB_Score")
text(x, y, labels = DNB_count, pos = 2, cex = 1, col = "red")

save(testres, igraphL, file = paste0("Fibroblasts_biotip.RData"), compress = TRUE)
save(list = c("testres", "igraphL", "cluster", "cut.fdr", "membersL", "topMCI", 
              "maxMCIms", "maxMCI", "CTS.Lib", "CTS.Lib.Symbol", "df", "simuMCI", 
              "BioTIP_scores", "SimResults_g"), 
     file = "ALL_Fibroblasts_biotip20260303.RData", compress = TRUE)


# --------------------------------------------------------------------------------------------------
# #### Figure 3.D: Quantitative Export of PCC and SD Profiles ####
# --------------------------------------------------------------------------------------------------
DNB_vec       <- CTS.Lib.Symbol[["4"]]
exprs_matrix  <- global
cluster_info  <- temp_sce$seurat_clusters
selected_clusters <- c(1, 4, 9, 6, 7, 0)

# Secure vector arrays against matrix dimensionality
exprs_matrix_subset <- exprs_matrix[DNB_vec, , drop = FALSE]

# Initialize storage objects
pcc_results <- list()
sd_results  <- list()

# Calculate correlative vectors (PCC) and deviations (SD)
for (cluster in selected_clusters) {
  
  cluster_cells <- which(cluster_info == cluster)
  cluster_exprs <- exprs_matrix_subset[, cluster_cells, drop = FALSE]
  
  # Pearson correlation map construction
  pcc_matrix <- cor(t(cluster_exprs), method = "pearson")
  
  # Translate adjacency matrix to EdgeList dataframes
  pcc_df <- as.data.frame(as.table(pcc_matrix))
  colnames(pcc_df) <- c("Gene1", "Gene2", "PCC")
  
  # Omit recursive identity lines (Matrix Diagonal)
  pcc_df <- pcc_df[pcc_df$Gene1 != pcc_df$Gene2, ]
  
  pcc_results[[as.character(cluster)]] <- pcc_df
  
  # Absolute variation tracking
  gene_sd <- apply(cluster_exprs, 1, sd)
  sd_results[[as.character(cluster)]] <- data.frame(gene = names(gene_sd), sd = gene_sd)
}

# Automated directory generation and structural batch exportation
output_dir <- "pcc_sd_results" 
dir.create(output_dir, showWarnings = FALSE)

for (cluster in selected_clusters) {
  # Write Correlative Networks
  pcc_file <- file.path(output_dir, paste0("cluster_", cluster, "_PCC.csv"))
  write.csv(pcc_results[[as.character(cluster)]], pcc_file, row.names = FALSE)
  
  # Write Deviations
  sd_file <- file.path(output_dir, paste0("cluster_", cluster, "_SD.csv"))
  write.csv(sd_results[[as.character(cluster)]], sd_file, row.names = FALSE)
}

cat("PCC and SD data matrices successfully exported!\n")
# ==================================================================================================
# MODULE: Figure 3.E-G Functional Co-expression Topology & Network Visualizations
# 
# DESCRIPTION: 
#   Generates 2D quadrant scatter plots for DNB co-expression targets, mfuzz 
#   temporal clustering (ClusterGVis), and tri-layered Sankey networks mapping 
#   DNBs to functional subsets and co-expression directionality.
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# #### Figure 3.E: 2D Quadrant Scatter Plot (DNB Target Gene Regulation) ####
# --------------------------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(readr)

# Load aggregated DNB co-expression data
# Historical files referenced: TableS13.Macro, TableS14.Fibro, TableS15.Fibro
res_df <- read_csv("./大论文图和数据/TableS12.Fibro_dnb_coexp_plot.csv")

# Transform p-values and assign aesthetic mapping parameters
res_df <- res_df %>%
  mutate(
    log_positive_p_value = -log10(positive_p_value + 1e-6),
    log_negative_p_value = -log10(negative_p_value + 1e-6),
    # Scale node size based on the most significant (minimum) -log10(p-value)
    size_factor = pmin(log_positive_p_value, log_negative_p_value),
    # Map color attributes based on significance thresholds (p < 0.2)
    color_factor = case_when(
      positive_p_value < 0.2 & negative_p_value < 0.2 ~ "red",   # Dual significance
      negative_p_value < 0.2                          ~ "green", # Negative significance only
      positive_p_value < 0.2                          ~ "blue",  # Positive significance only
      TRUE                                            ~ "grey"   # Non-significant
    )
  )

# Introduce plot coordinate scalars to prevent dense cluster overlapping
res_df2 <- res_df %>%
  mutate(
    x_plot = ifelse(negative_mean_log2FC > 0, negative_mean_log2FC * 1.4, negative_mean_log2FC),
    y_plot = ifelse(positive_mean_log2FC < 0, positive_mean_log2FC * 1.6, positive_mean_log2FC)
  )

# Isolate significant genes for textual annotation
label_df <- res_df2 %>% filter(color_factor != "grey")

# 1. Main Quadrant Plot Configuration
main_plot <- ggplot(res_df2, aes(x = x_plot, y = y_plot)) +
  geom_point(aes(color = color_factor, size = size_factor), alpha = 0.6) +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 3.6, box.padding = 0.4, point.padding = 0.6,
    max.overlaps = Inf, color = "black", segment.color = "black",
    segment.size = 0.4, min.segment.length = 0, force = 2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  scale_size_continuous(range = c(3, 10), name = "-log10(P-value)") +
  scale_color_manual(
    values = c("red" = "red", "green" = "green", "blue" = "blue", "grey" = "grey"),
    name = "Significance",
    labels = c("Positive Sig.", "Negative Sig.", "Not Sig.", "Dual Sig.")
  ) +
  # Map visual coordinates back to biological log2FC representations
  scale_x_continuous(
    name = "Inflammation Mean Log2FC",
    breaks = c(-1, 0, 1, 2), labels = c(-1, 0, 1, 2), sec.axis = dup_axis()
  ) +
  scale_y_continuous(
    name = "Repair Mean Log2FC",
    breaks = c(-1.6, 0, 1, 2, 3), labels = c(-1, 0, 1, 2, 3)
  ) +
  labs(title = "2D Quadrant Evaluation of Target Gene Regulation") +
  theme_minimal() +
  theme(
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 14),
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 10),
    plot.title       = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(20, 20, 20, 20)
  ) +
  coord_cartesian(xlim = c(-1.5, 3.5), ylim = c(-2.4, 4), expand = FALSE)

# 2. Sub-plot: Isolate and emphasize genes strictly within the Fourth Quadrant (x > 0, y < 0)
subset_df <- res_df %>%
  filter(negative_mean_log2FC > 0 & positive_mean_log2FC < 0) %>%
  filter(color_factor != "grey") %>%
  na.omit()

sub_plot <- ggplot(subset_df, aes(x = negative_mean_log2FC, y = positive_mean_log2FC)) +
  geom_point(aes(color = color_factor, size = size_factor), alpha = 0.6) +
  geom_text_repel(aes(label = gene), size = 3.6, box.padding = 0.35, point.padding = 0.5, max.overlaps = 10, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  scale_size_continuous(range = c(3, 10), name = "-log10(P-value)") +
  scale_color_manual(values = c("red" = "red", "green" = "green", "blue" = "blue", "grey" = "grey")) +
  labs(
    title = "Fourth Quadrant Genes",
    x = "Inflammation Mean Log2FC",
    y = "Repair Mean Log2FC"
  ) +
  theme_minimal() +
  theme(
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 14),
    plot.title       = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(20, 20, 20, 20),
    legend.position  = "none"
  ) +
  scale_x_continuous(breaks = seq(0, 3, by = 1)) +
  scale_y_continuous(breaks = seq(-3, 0, by = 1)) +
  coord_cartesian(expand = TRUE)

# 3. Render Combined Panel
grid.arrange(main_plot, sub_plot, ncol = 2)


# --------------------------------------------------------------------------------------------------
# #### Figure 3.F: ClusterGVis Temporal Clustering of DNB Co-expression Networks ####
# --------------------------------------------------------------------------------------------------
load(file = ".Fibro_analysis_DNB_data.RData")
library(ClusterGVis)

# Synthesize unique co-expression gene universe
all_coexpr_genes <- unique(unlist(dnb_coexpr_genes_list))
all_coexpr_genes <- unique(result_table$TARGET)
cat("Total number of unique coexpressed genes across all DNBs:", length(all_coexpr_genes), "\n")

# Filter expression matrix for active clustering components
expr_mat <- expr_mat_use[c(all_coexpr_genes, "Atf4", "Mki67"), ]
expr_mat <- expr_mat[apply(expr_mat, 1, function(x) sd(x, na.rm = TRUE) > 0), , drop = FALSE]

# Execute mfuzz soft clustering framework
ck <- clusterData(exp = expr_mat, cluster.method = "mfuzz", cluster.num = 5)

# Perform Biological Process (BP) Enrichment on defined clusters
enrich <- enrichCluster(
  object       = ck,
  OrgDb        = org.Mm.eg.db,
  type         = "BP",
  pvalueCutoff = 0.05,
  topn         = 20,
  seed         = 5201314
)

# Core markers highlighted alongside topological clusters
markGenes <- c("Stmn1", "Cdk2", "Mki67", "Mcm6", "Ccne1", "Col4a4", "Eln", "Adamtsl2", "Pi16", "Dkk3",
               "Cxcl5", "Pcna", "Hif1a", "Cks1b", "Dnmt1", "Col3a1", "Cdh2", "Runx1", "Pdgfa", "Bhlhe41",
               "Angptl4", "Atf4", "Jun", "Pdgfra", "Mt1", "Nfkb1", "Acta2")

# Generate structured multi-panel visualization
pdf('Figure3.E-1.pdf', height = 10, width = 14, onefile = FALSE)
visCluster(
  object           = ck,
  plot.type        = "both",
  column_names_rot = 45,
  show_row_dend    = FALSE,
  markGenes        = markGenes,
  markGenes.side   = "left",
  genes.gp         = c('italic', fontsize = 12, col = "black"),
  annoTerm.data    = enrich,
  line.side        = "left",
  go.col           = rep(ggsci::pal_d3()(5), each = 5),
  go.size          = "pval",
  add.bar          = TRUE, 
  add.line         = FALSE,
  add.box          = TRUE,
  textbar.pos      = c(0.8, 0.2)
)
dev.off()

# Export annotation dictionaries
write.csv(enrich, "./TableS_Fibroblasts_Coexp_enrich.csv")
write.csv(ck$long.res, "./TableS_Fibroblasts_Coexp_Clusters2.csv")

# Classify Co-expression Targets by Functional Modules
inflammation_genes <- ck$wide.res %>% filter(cluster == 5) %>% .[["gene"]]
ECM_genes          <- ck$wide.res %>% filter(cluster == 2) %>% .[["gene"]]
Prolif_genes       <- ck$wide.res %>% filter(cluster %in% c(1, 3)) %>% .[["gene"]]
Vascular_genes     <- ck$wide.res %>% filter(cluster == 4) %>% .[["gene"]]

result_table <- result_table %>%
  mutate(TYPE = case_when(
    TARGET %in% inflammation_genes ~ "Inflam",  
    TARGET %in% ECM_genes          ~ "Fibrosis",     
    TARGET %in% Prolif_genes       ~ "Prolif",      
    TARGET %in% Vascular_genes     ~ "Vascular",     
    TRUE                           ~ "Other"                           
  ))


# --------------------------------------------------------------------------------------------------
# #### Figure 3.G: Tri-Layered Sankey Network (DNB -> Cluster -> Directionality) ####
# --------------------------------------------------------------------------------------------------
library(networkD3)
library(dplyr)
library(readxl)
library(tibble)
library(pheatmap)

# Ingest DNB definitions
result_table <- read.csv("./Fibro_analysis_DNB_data.csv")
unique_dnbs  <- DNB_vec

# Active DNB targets for focused Sankey topology
active_dnbs <- c("Cxcl1", "Nfkb1", "Nr4a1", "Tnfaip6", "Adam12", "Snai1")

for (dnb in active_dnbs) {
  
  # Acquire interaction topologies
  dnb_data <- read_xlsx("./大论文图和数据/ASII_Rd_Sankeyplot.xlsx")
  colnames(dnb_data)[2] <- "TARGET" 
  colnames(dnb_data)[5] <- "COR" 
  
  # Filter active network interactions
  dnb_data <- result_table[result_table$DNB %in% dnb, ] %>% filter(abs(Cor_values) > 0.1)
  
  # Establish hierarchical node orders
  dnb_gene_order <- dnb_data %>%
    count(TYPE, TARGET) %>%
    arrange(TYPE, desc(n)) %>%
    pull(TARGET)
  
  cat(paste0("\n=== Structure & Connectivity Map for DNB: ", dnb, " ===\n"))
  
  # Construct hierarchical nodes lists
  dnb_nodes      <- dnb
  cluster_nodes  <- paste0("Cluster_", sort(unique(dnb_data$TYPE)))
  coexgene_nodes <- dnb_gene_order
  
  nodes <- data.frame(name = c(dnb_nodes, cluster_nodes, coexgene_nodes))
  
  # Build Layer 1: Target Gene to Functional Cluster linkages
  layer1_links_dnb <- dnb_data %>%
    group_by(TARGET, TYPE) %>%
    summarise(Value = n(), .groups = 'drop') %>%
    mutate(
      source = match(DNB, nodes$name) - 1,
      target = match(paste0("Cluster_", TYPE), nodes$name) - 1,
      group  = as.factor(TYPE)  
    ) %>%
    dplyr::select(source, target, value = Value, group)
  
  # Build Layer 2: Functional Cluster to Co-expression directionality mappings
  layer2_links_dnb <- dnb_data %>%
    mutate(
      source = match(paste0("Cluster_", TYPE), nodes$name) - 1,
      target = match(TARGET, nodes$name) - 1,
      value  = 1,
      group  = as.factor(COR)  # Correlation vector alignment (Positive/Negative)
    ) %>%
    dplyr::select(source, target, value, group)
  
  # Consolidate Multi-Layered Sankey Topology
  links_dnb <- bind_rows(layer1_links_dnb, layer2_links_dnb)
  
  # Define continuous color mappings via JavaScript injection for D3.js engine
  dnb_all_groups     <- unique(c(as.character(layer1_links_dnb$group), as.character(layer2_links_dnb$group)))
  dnb_cluster_groups <- grep("^Cluster_", dnb_all_groups, value = TRUE)
  dnb_cluster_colors <- rainbow(length(dnb_cluster_groups))
  
  dnb_direction_colors <- c()
  if (any(c("Positive", "positive") %in% dnb_all_groups)) dnb_direction_colors[c("Positive", "positive")] <- "red"
  if (any(c("Negative", "negative") %in% dnb_all_groups)) dnb_direction_colors[c("Negative", "negative")] <- "green"
  if ("unknown" %in% dnb_all_groups)                      dnb_direction_colors["unknown"] <- "gray"
  
  dnb_all_colors  <- c(dnb_cluster_colors, dnb_direction_colors)
  dnb_all_domains <- c(dnb_cluster_groups, names(dnb_direction_colors))
  
  dnb_color_scale_js <- paste0('d3.scaleOrdinal() .domain(["', 
                               paste(dnb_all_domains, collapse = '", "'), 
                               '"]) .range(["', 
                               paste(dnb_all_colors, collapse = '", "'), '"])')
  
  # Render networkD3 Sankey Object
  sankey_plot <- sankeyNetwork(
    Links       = links_dnb,
    Nodes       = nodes,
    Source      = "source",
    Target      = "target",
    Value       = "value",
    NodeID      = "name",
    LinkGroup   = "group",
    units       = "",
    fontSize    = 10,
    nodeWidth   = 20,
    height      = 900,
    width       = 900,
    nodePadding = 8,
    sinksRight  = FALSE,
    iterations  = 0,                      # Restrict physical node permutations
    colourScale = JS(dnb_color_scale_js)  # Inject mapped aesthetic JS constraints
  )
  
  # HTML Object Export
  output_dir <- "/home/tsh/E-MTAB-7895/大论文图和数据/"
  html_file  <- file.path(output_dir, paste0("Figure3.G_Sankey_", dnb, ".html"))
  pdf_file   <- file.path(output_dir, paste0("Figure3.G_Sankey_", dnb, ".pdf"))
  
  library(htmlwidgets)
  saveWidget(sankey_plot, html_file)
  
  # Heatmap Construction representing Sankey Target Expressions
  existing_genes <- dnb_data$TARGET[dnb_data$TARGET %in% rownames(expr_mat_use)]
  
  if (length(existing_genes) > 0) {
    heatmap_data <- dnb_data[, c("TARGET", "logFC", "Cor_values")] %>% 
                    column_to_rownames("TARGET") %>% 
                    .[dnb_gene_order, , drop = FALSE]
    
    annotation_row <- data.frame(
      Correlation = dnb_data$COR[match(rownames(heatmap_data), dnb_data$TARGET)],
      Cluster     = as.factor(dnb_data$TYPE[match(rownames(heatmap_data), dnb_gene_order)])
    )
    
    # Secure row names via unique suffixes if necessary
    rownames(annotation_row) <- rownames(heatmap_data) <- make.unique(rownames(heatmap_data), sep = "-")
    
    # Synchronize hierarchical clustering arrays
    ord            <- order(annotation_row$Cluster)
    annotation_row <- annotation_row[ord, , drop = FALSE]
    heatmap_data   <- heatmap_data[ord, , drop = FALSE]
    
    annotation_colors <- list(
      Correlation = c(Positive = "red", Negative = "green", unknown = "gray"),
      Cluster     = setNames(rainbow(length(unique(annotation_row$Cluster))), unique(annotation_row$Cluster))
    )
    
    pheatmap(heatmap_data,
             cluster_rows      = FALSE,
             cluster_cols      = FALSE,
             scale             = "row",
             color             = colorRampPalette(c("#2E75B6", "white", "#C00000"))(100),
             show_rownames     = FALSE,
             show_colnames     = TRUE,
             annotation_row    = annotation_row,
             annotation_colors = annotation_colors)
  }
  
  # PDF Rendering utilizing webshot virtual engine
  if (requireNamespace("webshot", quietly = TRUE)) {
    tryCatch({
      webshot::webshot(html_file, pdf_file)
    }, error = function(e) {
      print(paste("Failed to generate PDF for Sankey:", e$message))
    })
  }
}
# ==================================================================================================
# MODULE: Figure 3.D Cross-Cluster Correlation Topology & Heatmap Projection
# ==================================================================================================
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ========================
# 1. Deplete overlapping loci to circumvent auto-correlation artifacts
# ========================
DNB_vec <- CTS.Lib.Symbol[["4"]]
DNB_vec <- unique(DNB_vec)
all_coexpr_genes <- unique(all_coexpr_genes)
CTS.Lib.Symbol[["5"]]
# Optional toggle: Bypass exclusion of DNB-overlapping genes within the co-expression module
all_coexpr_genes_use <- setdiff(all_coexpr_genes, DNB_vec)
all_coexpr_genes_use <- all_coexpr_genes
# ========================
# 2. Assign clustering identities native to the Seurat object
# ========================
Idents(temp_sce) <- temp_sce$seurat_clusters

clusters_to_plot <- c("4", "1", "9")

# ========================
# 3. Define local functional module: Extract cluster-specific correlation matrices
# ========================
get_cluster_cor_mat <- function(
    seu,
    cluster_id,
    DNB_genes,
    coexpr_genes,
    assay = "RNA",
    slot = "data",
    method = "pearson",
    min_cells = 10
) {
  cells_use <- WhichCells(seu, idents = cluster_id)
  
  if (length(cells_use) < min_cells) {
    stop(paste0("Cluster ", cluster_id, " cell count falls below minimum threshold: ", min_cells))
  }
  
  # Retrieve expression matrix while strictly preserving row and column naming attributes
  expr_mat <- GetAssayData(seu, assay = assay, slot = slot)[, cells_use, drop = FALSE]
  rownames(expr_mat) <- rownames(GetAssayData(seu, assay = assay, slot = slot))  # Enforce row names
  colnames(expr_mat) <- cells_use  # Enforce column names via cell identifiers
  
  expr_mat <- as.matrix(expr_mat)
  
  # Secondary stringency filter against the derived expression matrix row names
  DNB_use <- intersect(DNB_genes, rownames(expr_mat))
  coexpr_use <- intersect(coexpr_genes, rownames(expr_mat))
  
  if (length(DNB_use) == 0) {
    stop(paste0("Cluster ", cluster_id, ": No valid DNB genes mapped to the expression matrix"))
  }
  if (length(coexpr_use) == 0) {
    stop(paste0("Cluster ", cluster_id, ": No valid co-expression targets mapped to the expression matrix"))
  }
  
  mat1 <- expr_mat[DNB_use, , drop = FALSE]
  mat2 <- expr_mat[coexpr_use, , drop = FALSE]
  
  # Initialize core matrices for correlation coefficients and corresponding p-values
  cor_mat <- matrix(NA, nrow = length(DNB_use), ncol = length(coexpr_use))
  pval_mat <- matrix(NA, nrow = length(DNB_use), ncol = length(coexpr_use))
  
  # Assign definitive dimensional descriptors to the initialized matrices
  rownames(cor_mat) <- DNB_use
  colnames(cor_mat) <- coexpr_use
  rownames(pval_mat) <- DNB_use
  colnames(pval_mat) <- coexpr_use
  
  # Execute base cor.test to iteratively compute coefficients and significance metrics
  for (i in 1:length(DNB_use)) {
    for (j in 1:length(coexpr_use)) {
      test <- cor.test(mat1[i, , drop = FALSE], mat2[j, , drop = FALSE], method = method)
      cor_mat[i, j] <- test$estimate  # Correlation coefficient
      pval_mat[i, j] <- test$p.value  # Statistical significance
    }
  }
  
  # Historical toggle: Constrain output to strictly significant interactions (p < 0.05)
  #cor_mat[pval_mat >= 0.05] <- 0
  
  return(cor_mat)
}

clean_cor_mat <- function(mat) {
  # Retain original coordinate nomenclatures
  rownames_mat <- rownames(mat)
  colnames_mat <- colnames(mat)
  
  mat <- as.matrix(mat)
  
  # Purge artifacts: Convert NA, NaN, and Infinite values to baseline (0)
  mat[is.na(mat)] <- 0
  mat[is.nan(mat)] <- 0
  mat[is.infinite(mat)] <- 0
  
  # Re-apply structural coordinate names
  rownames(mat) <- rownames_mat
  colnames(mat) <- colnames_mat
  
  return(mat)
}

# ========================
# 4. Generate discrete sub-cluster arrays and apply significance gating
# ========================
cor_list <- lapply(clusters_to_plot, function(clu) {
  mat <- get_cluster_cor_mat(
    seu = temp_sce,
    cluster_id = clu,
    DNB_genes = DNB_vec,
    coexpr_genes = all_coexpr_genes_use,
    assay = "RNA",
    slot = "data",
    method = "pearson"
  )
  print(rownames(mat))  # Output structural row names for verification
  print(colnames(mat))  # Output structural column names for verification
  mat <- clean_cor_mat(mat)
  mat
})
names(cor_list) <- paste0("Cluster ", clusters_to_plot)

# Configure custom continuous color mapping for Pearson correlations
col_fun <- colorRamp2(
  c(-0.4, -0.2, 0, 0.2, 0.4),  # Injected high-resolution discrete breakpoints to amplify gradient contrasts
  c("#2166AC", "#1D68B5", "white", "#E31A1C", "#B2182B")  # High-intensity terminal hues
)

ht_list <- NULL

head(mat[1:10,1:10])



# Iterate ComplexHeatmap generation per sub-cluster
nm <- "Cluster 4"
for (nm in names(cor_list)) {
  mat <- cor_list[[nm]]
  
  ht <- Heatmap(
    mat,
    name = "Pearson",                                                              
    col = col_fun,
    cluster_rows = T,
    cluster_columns = T,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    show_row_names = TRUE,  # Render row annotations
    show_column_names = T,  # Render column annotations
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 1),  # Optional font size modulation for high-density columns
    column_title = nm,
    column_title_gp = gpar(fontsize = 1, fontface = "bold"),
    border = TRUE,
    heatmap_legend_param = list(
      title = "Pearson\ncorrelation",
      at = c(-0.4, 0, 0.4)
    )
  )
  
  if (is.null(ht_list)) {
    ht_list <- ht
  } else {
    ht_list <- ht_list + ht
  }
}

# Render composite concatenated heatmap array
ht_list_drawn <- draw(
  ht_list,
  merge_legends = TRUE,
  heatmap_legend_side = "right"
)





library(ComplexHeatmap)

# 1. Isolate the baseline matrix from the target anchor plot
first_nm <- names(cor_list)[1]
first_mat <- cor_list[[first_nm]]
first_nm <- "Cluster 4"
identical(first_mat, mat)
# 2. Manually compute the hierarchical dendrogram topology for the anchor plot
# This enforces topological congruence across all subsequent visual projections
row_dend <- as.dendrogram(hclust(dist(first_mat, method = "euclidean"), method = "ward.D2"))
col_dend <- as.dendrogram(hclust(dist(t(first_mat), method = "euclidean"), method = "ward.D2"))

ht_list <- NULL
nm <- "Cluster 4"
for (nm in names(cor_list)) {
  mat <- cor_list[[nm]]
  
  # Determine logical hierarchy (Render dendrograms exclusively on the anchor plot to eliminate redundancy)
  is_first <- (nm == first_nm)
  
  ht <- Heatmap(
    mat,
    name = "Pearson",
    col = col_fun,
    
    # --- Critical: Enforce row topological consistency ---
    cluster_rows = row_dend,             # Apply the anchor plot's hierarchical row tree
    row_dend_side = "left",
    show_row_dend = is_first,            # Restrict visual rendering to the anchor plot
    
    # --- Critical: Enforce column topological consistency ---
    cluster_columns = col_dend,          # Apply the anchor plot's hierarchical column tree
    show_column_dend = is_first,         # Restrict visual rendering to the anchor plot
    
    # Secondary aesthetic configurations
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 2),
    column_title = nm,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    border = TRUE,
    heatmap_legend_param = list(
      title = "Pearson\ncorrelation",
      at = c(-0.4, 0, 0.4)
    )
  )
  
  if (is.null(ht_list)) {
    ht_list <- ht
  } else {
    ht_list <- ht_list + ht
  }
}

# Final execution of topologically synchronized arrays
draw(ht_list, merge_legends = TRUE, heatmap_legend_side = "right")


library(ComplexHeatmap)
ht = draw(ht)

ord <-  ComplexHeatmap::column_order(ht)
colnames(mat)
Heatmap(
  mat[,ord],
  name = "Pearson",                                                              
  col = col_fun,
  cluster_rows = row_dend,  
  cluster_columns = F,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = TRUE,  # Render row annotations
  show_column_names = T,  # Render column annotations
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 12),  # Modulate font aesthetics for variable column densities
  column_title = nm,
  column_title_gp = gpar(fontsize = 6, fontface = "bold"),
  border = TRUE,
  heatmap_legend_param = list(
    title = "Pearson\ncorrelation",
    at = c(-0.4, 0, 0.4)
  )
)
heatmap_c4_order <- expr_mat_use[colnames(mat),]
heatmap_c4_order <- heatmap_c4_order[ord,]

pheatmap(heatmap_c4_order,
         #main = paste("Expression Heatmap for DNB:", dnb, "\n(Red: Positive correlation, Green: Negative correlation)"),
         cluster_rows = F,
         cluster_cols = FALSE,
         scale = "row",
         #border_color = "white",
         color = colorRampPalette(c("#2E75B6", "white", "#C00000"))(100),
         show_rownames = T,
         show_colnames = TRUE,
         fontsize_row = 1,
         # fontsize_col = 6,
         # angle_col = "45",
         # cellwidth = 10,
         # cellheight = max(8, 200 / length(rownames(heatmap_data))),
         #annotation_row = annotation_row,
         #annotation_colors = annotation_colors
)


head(cor_list[["Cluster 1"]][1:5,1:5])
head(cor_list[["Cluster 4"]][1:5,1:5])
head(cor_list[["Cluster 9"]][1:5,1:5])
# ==================================================================================================
# MODULE: DNB vs Stage 1 / Stage 2 Marker Correlation Quadrant Evaluation
# 
# DESCRIPTION: 
#   Computes Spearman correlations between Dynamic Network Biomarker (DNB) genes 
#   and stage-specific markers. Aggregates effect sizes via Fisher z-transformation 
#   and evaluates significance utilizing Wilcoxon signed-rank tests.
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# 1. Environment Setup & User Inputs
# --------------------------------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(ggrepel)

# Input Parameters
seurat_obj   <- readRDS("/home/tsh/E-MTAB-7895/temp_sce.rds")
DNB_vec      <- CTS.Lib.Symbol[["4"]]
clusters_use <- c("1", "4", "9", "6", "7", "0")
assay_use    <- "RNA"
slot_use     <- "data"
cor_method   <- "pearson"
min_cells    <- 10

# --------------------------------------------------------------------------------------------------
# 2. Optimized Helper Functions
# --------------------------------------------------------------------------------------------------

# Safely compute correlation, bypassing invariant vectors or sparse observations
safe_cor <- function(x, y, method = "pearson", min_n = 10) {
  valid_idx <- which(!is.na(x) & !is.na(y))
  if (length(valid_idx) < min_n || sd(x[valid_idx]) == 0 || sd(y[valid_idx]) == 0) return(NA_real_)
  suppressWarnings(cor(x[valid_idx], y[valid_idx], method = method))
}

# Aggregate correlation coefficients via Fisher z-transformation for robust central tendency
fisher_mean_cor <- function(r_vec) {
  r_vec <- na.omit(r_vec)
  if (length(r_vec) == 0) return(NA_real_)
  r_vec <- pmin(pmax(r_vec, -0.999999), 0.999999) # Clamp to prevent Inf at limits
  tanh(mean(atanh(r_vec)))
}

# Vectorized correlation computation between a single DNB and a marker set
compute_gs_cor <- function(dnb_expr, marker_genes, expr_mat, method = "pearson", min_n = 10) {
  marker_genes <- intersect(marker_genes, rownames(expr_mat))
  if (length(marker_genes) == 0) return(tibble(marker = character(0), cor = numeric(0)))
  
  cor_vec <- vapply(marker_genes, function(g) {
    safe_cor(dnb_expr, expr_mat[g, ], method = method, min_n = min_n)
  }, numeric(1))
  
  tibble(marker = marker_genes, cor = cor_vec)
}

# Summarize marker-wise correlations and perform Wilcoxon signed-rank test
summarize_cor <- function(cor_tbl) {
  cor_vec <- na.omit(cor_tbl$cor)
  if (length(cor_vec) == 0) return(tibble(mean_cor = NA, median_cor = NA, n_markers = 0L, p_value = NA))
  
  p_val <- tryCatch(wilcox.test(cor_vec, mu = 0, alternative = "two.sided")$p.value, error = function(e) NA_real_)
  
  tibble(
    mean_cor   = fisher_mean_cor(cor_vec),
    median_cor = median(cor_vec),
    n_markers  = length(cor_vec),
    p_value    = p_val
  )
}

# --------------------------------------------------------------------------------------------------
# 3. Stage-Specific Marker Extraction & Data Preparation
# --------------------------------------------------------------------------------------------------
Idents(seurat_obj) <- "seurat_clusters"

# Retrieve predefined markers (Assuming 'result_table' is loaded in environment)
S1_markers <- result_table %>% filter(TYPE %in% c("Fibrosis", "Inflam")) %>% pull(TARGET) %>% unique()
S2_markers <- result_table %>% filter(TYPE %in% c("Prolif", "Vascular")) %>% pull(TARGET) %>% unique()

# Deplete overlapping loci to ensure orthogonal axis representation
marker_overlap   <- intersect(S1_markers, S2_markers)
S1_markers_clean <- setdiff(S1_markers, marker_overlap)
S2_markers_clean <- setdiff(S2_markers, marker_overlap)

cat("Cleaned S1 markers:", length(S1_markers_clean), "| Cleaned S2 markers:", length(S2_markers_clean), "\n")

# Extract bounded expression matrix
cells_use   <- colnames(seurat_obj)[seurat_obj$seurat_clusters %in% clusters_use]
expr_mat    <- as.matrix(GetAssayData(seurat_obj, assay = assay_use, slot = slot_use)[, cells_use, drop = FALSE])
DNB_vec_use <- intersect(DNB_vec, rownames(expr_mat))

cat("DNB genes validated in matrix:", length(DNB_vec_use), "\n")

# --------------------------------------------------------------------------------------------------
# 4. Core Analytical Execution (DNB evaluation)
# --------------------------------------------------------------------------------------------------
result_df <- map_dfr(DNB_vec_use, function(dnb) {
  dnb_expr   <- expr_mat[dnb, ]
  cor_s1_tbl <- compute_gs_cor(dnb_expr, S1_markers_clean, expr_mat, cor_method, min_cells)
  cor_s2_tbl <- compute_gs_cor(dnb_expr, S2_markers_clean, expr_mat, cor_method, min_cells)
  
  bind_cols(tibble(DNB = dnb), 
            summarize_cor(cor_s1_tbl) %>% rename_with(~paste0(., "_s1")),
            summarize_cor(cor_s2_tbl) %>% rename_with(~paste0(., "_s2")))
})

# Vectorized multiple testing correction, quadrant assignment, and visual scaling
result_df <- result_df %>%
  mutate(
    padj_s1 = p.adjust(p_value_s1, method = "BH"),
    padj_s2 = p.adjust(p_value_s2, method = "BH"),
    
    # Vectorized Quadrant Assignment (Replacing slow rowwise() loop)
    quadrant = case_when(
      is.na(mean_cor_s1) | is.na(mean_cor_s2) ~ NA_character_,
      mean_cor_s1 >= 0 & mean_cor_s2 >= 0     ~ "Q1: S1+ / S2+",
      mean_cor_s1 <  0 & mean_cor_s2 >= 0     ~ "Q2: S1- / S2+",
      mean_cor_s1 <  0 & mean_cor_s2 <  0     ~ "Q3: S1- / S2-",
      mean_cor_s1 >= 0 & mean_cor_s2 <  0     ~ "Q4: S1+ / S2-"
    ),
    
    sig_class = case_when(
      padj_s1 < 0.05 & padj_s2 < 0.05 ~ "Both significant",
      padj_s1 < 0.05                  ~ "S1 significant",
      padj_s2 < 0.05                  ~ "S2 significant",
      TRUE                            ~ "Not significant"
    ),
    
    size_value = pmax(-log10(padj_s1 + 1e-300), -log10(padj_s2 + 1e-300), na.rm = TRUE)
  )

write.csv(result_df, "DNB_S1_S2_quadrant_summary.csv", row.names = FALSE)

# --------------------------------------------------------------------------------------------------
# 5. Export Detailed Pairwise Interaction Matrices
# --------------------------------------------------------------------------------------------------
# Streamlined binding of iteration loops
detail_df <- bind_rows(
  map_dfr(DNB_vec_use, ~ compute_gs_cor(expr_mat[.x, ], S1_markers_clean, expr_mat, cor_method, min_cells) %>% mutate(DNB = .x, marker_set = "S1")),
  map_dfr(DNB_vec_use, ~ compute_gs_cor(expr_mat[.x, ], S2_markers_clean, expr_mat, cor_method, min_cells) %>% mutate(DNB = .x, marker_set = "S2"))
)
write.csv(detail_df, "DNB_marker_pairwise_correlations.csv", row.names = FALSE)

# --------------------------------------------------------------------------------------------------
# 6. Quadrant Scatter Visualization
# --------------------------------------------------------------------------------------------------
label_df <- result_df %>% filter(sig_class != "Not significant")

p <- ggplot(result_df, aes(x = mean_cor_s1, y = mean_cor_s2)) +
  geom_point(aes(color = sig_class, size = size_value), alpha = 0.75) +
  geom_text_repel(
    data = label_df, aes(label = DNB),
    size = 3.8, box.padding = 0.4, point.padding = 0.3, max.overlaps = 100
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c(
      "Both significant" = "#D73027",
      "S1 significant"   = "#4575B4",
      "S2 significant"   = "#1A9850",
      "Not significant"  = "grey70"
    )
  ) +
  scale_size_continuous(range = c(2.5, 8)) +
  labs(
    title = "Quadrant Evaluation: DNB Correlations with S1/S2 Markers",
    x     = "Mean correlation with S1 markers",
    y     = "Mean correlation with S2 markers",
    color = "Significance",
    size  = "-log10(p-adj)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)

ggsave("DNB_S1_S2_quadrant_plot.pdf", p, width = 8, height = 6)
ggsave("DNB_S1_S2_quadrant_plot.png", p, width = 8, height = 6, dpi = 300)
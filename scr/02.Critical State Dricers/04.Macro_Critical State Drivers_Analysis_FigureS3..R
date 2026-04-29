# ==================================================================================================
# MODULE: Supplementary Figure 3 (Macrophages) - Trajectory, Enrichment & DNB Analysis
# 
# DESCRIPTION: 
#   Comprehensive analysis of macrophage subpopulations including UMAP projections, 
#   marker identification, Monocle 2 pseudotime inference, functional enrichment 
#   (ClusterGVis), and Dynamic Network Biomarker (BioTIP) critical state validation.
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# #### Figure 3S.A: Macrophage Subpopulation Initialization & Visualizations ####
# --------------------------------------------------------------------------------------------------
temp_sce <- readRDS("/home/tsh/E-MTAB-7895/filtered_inflammation_temp_sce.rds") # Macrophages
temp_sce <- subset(temp_sce, subset = orig.ident %in% c("0day", "1day", "3day", "5day", "7day"))
table(temp_sce$orig.ident)

DimPlot(temp_sce, reduction = "umap", group.by = c("orig.ident", "seurat_clusters", "inflammation_human"), label = TRUE)
DimPlot(temp_sce, reduction = "umap", group.by = c("seurat_clusters"), label = TRUE)

FeatureDimPlot(
  srt              = temp_sce,
  features         = c("Ly6c2", "Cxcl3", "Trem2"),
  compare_features = TRUE,
  label            = TRUE,
  label_insitu     = TRUE,
  reduction        = "umap",
  theme_use        = "theme_blank",
  pt.size          = 1
)

# Extract cluster-specific differential markers
macro_markers <- FindAllMarkers(temp_sce, group.by = "seurat_clusters")
markers       <- macro_markers %>%
                 dplyr::group_by(cluster) %>%
                 dplyr::top_n(n = 3, wt = avg_log2FC)

marers_genes <- unique(markers$gene)
marers_genes <- c(
  "Apoc4", "Cd83", "Adra1a",    # Cluster 0
  "Fat3", "Spata7", "Cd4",      # Cluster 1
  "Dmkn", "Lrg1", "Ly6c2",      # Cluster 2
  "Pde4d", "F7", "Tubb2a",      # Cluster 3
  "Fabp5", "Lpl", "Trem2",      # Cluster 4
  "Nos2", "Egln3", "Cxcl3",     # Cluster 5
  "Adgre4", "Vnn3", "Plac8",    # Cluster 6
  "H2-Ab1", "H2-Aa", "H2-Eb1",  # Cluster 7
  "Ppbp", "Il6", "Arg1",        # Cluster 8
  "Card11", "Dusp9", "Mmp12",   # Cluster 9
  "Col11a1", "Col4a5", "Bmp5",  # Cluster 10
  "Nlrp12", "Lcn2", "Sell",     # Cluster 11
  "Lrg1", "Snai1", "Ccr2",      # Cluster 12
  "Wnt5a", "C1qtnf1", "Nrap"    # Cluster 13
  # "Sox9", "Anxa2", "Anxa2"    # Cluster 14 (Commented out historically)
)
head(markers)

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

factor(temp_sce@meta.data$seurat_clusters)

# --------------------------------------------------------------------------------------------------
# #### Monocle 2 Trajectory Inference (Macrophage Lineage) ####
# --------------------------------------------------------------------------------------------------
# Target specific developmental trajectory clusters
temp_sce <- subset(temp_sce, subset = seurat_clusters %in% c("12", "5", "8", "3", "4", "1"))
# Active override based on subsequent exploratory logic
temp_sce <- subset(temp_sce, subset = seurat_clusters %in% c("10", "13", "4", "9", "7", "1", "0", "2", "5", "14", "3", "6", "12", "11", "8"))

factor(temp_sce@meta.data$seurat_clusters)
stopifnot("seurat_clusters" %in% colnames(temp_sce@meta.data))

library(monocle)
Idents(temp_sce) <- "seurat_clusters"

# Construct Monocle 2 compatible CellDataSet (CDS)
sparse_data <- as(as.matrix(temp_sce[["RNA"]]$counts), "sparseMatrix")
mdata       <- new("AnnotatedDataFrame", data = temp_sce@meta.data)
fData       <- data.frame(gene_short_name = row.names(sparse_data), row.names = row.names(sparse_data))
fd          <- new("AnnotatedDataFrame", data = fData)

monocle_cds <- newCellDataSet(cellData            = sparse_data,
                              phenoData           = mdata,
                              featureData         = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily    = negbinomial.size())

# Parameter estimation
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 10)

# Identify highly variable genes governing the trajectory
expressed_genes <- VariableFeatures(temp_sce)
diff_test_res   <- differentialGeneTest(monocle_cds[expressed_genes, ],
                                        fullModelFormulaStr = "~inflammation_human")
# Enforce strict q-value thresholding
ordering_genes  <- row.names(subset(diff_test_res, qval < 0.01)) 
monocle_cds     <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)

# Dimensionality reduction via DDRTree
library(igraph)
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, norm_method = "log", reduction_method = "DDRTree")

# Order cells along the pseudotime manifold
monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds, root_state = "2")
monocle_cds <- orderCells(monocle_cds, root_state = "1")



# Trajectory Visualizations
p  <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 0.6, show_backbone = TRUE, show_branch_points = FALSE) + facet_wrap("~orig.ident", nrow = 1)
p0 <- plot_cell_trajectory(monocle_cds, color_by = "inflammation_human", cell_size = 0.5, show_backbone = TRUE, show_branch_points = FALSE) + facet_wrap("~orig.ident", nrow = 1)
p1 <- plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters", cell_size = 0.6, show_backbone = FALSE, show_branch_points = FALSE) +
      theme(legend.position = 'none', panel.border = element_blank()) +
      scale_color_manual(values = color) # Assuming 'color' array is defined in global environment

# --------------------------------------------------------------------------------------------------
# #### ClusterGVis: Trajectory-Dependent Gene Ontology (GO) Enrichment ####
# --------------------------------------------------------------------------------------------------
# Note: Re-inject pseudotime coordinates back to Seurat metadata (Requires Monocle 3 `cds` object from prior analysis)
# temp_sce$pseudotime <- pseudotime(cds, reduction_method = "UMAP")

macro_markers <- FindAllMarkers(temp_sce, group.by = "inflammation_human")
markers       <- macro_markers %>%
                 dplyr::group_by(cluster) %>%
                 dplyr::top_n(n = 50, wt = avg_log2FC)

head(markers)
markGenes <- c("Ly6c2", "Cd177", "Dmkn", "Cd209a", "Chil3",
               "Nos2", "Cxcl3", "Il6", "Ccr2", "F10", "Cd14",
               "C1qc", "Gpnmb", "Fabp5", "Bhlhe41", "Trem2",
               "H2-DMa", "H2-Aa", "H2-Ab1", "H2-DMb1", "Ighm",
               "Ccna2", "Hmmr", "Birc5", "Depdc1b", "Ccne2",
               "Col11a1", "Col4a6", "Myoc")

library(ClusterGVis)
Idents(temp_sce) <- "inflammation_human"

# Aggregate multi-cluster expression matrices
data <- prepareDataFromscRNA(object        = temp_sce,
                             diffData      = markers,
                             keep.uniqGene = FALSE,  
                             showAverage   = TRUE)

# Execute Biological Process (BP) enrichment mapping
enrich <- enrichCluster(object       = data,
                        OrgDb        = org.Mm.eg.db,
                        type         = "BP",
                        organism     = "mus",
                        pvalueCutoff = 0.05,
                        topn         = 5,
                        seed         = 123)

head(enrich)



# Render comprehensive GO mapping heatmap
pdf("FigureS1.Macro_markers_heatmap.pdf", width = 16, height = 10)
visCluster(
  object           = data, 
  plot.type        = "both",
  column_names_rot = 45,
  show_row_dend    = FALSE,
  markGenes        = markGenes,
  markGenes.side   = "left",
  annoTerm.data    = enrich,
  line.side        = "left",
  cluster.order    = c("Monocytes", "Pro_inflammation", "Proliferation", "Reparative", "Antigen-presenting", "Collagens"),
  go.col           = rep(jjAnno::useMyCol("stallion", n = 6), each = 5)[1:30],
  add.bar          = TRUE
)
dev.off()

write.csv(macro_markers, "./macro_markers.csv")

# --------------------------------------------------------------------------------------------------
# #### Pseudotime-Dependent Gene Expression Smoothing ####
# --------------------------------------------------------------------------------------------------
selected_cells <- temp_sce[, temp_sce$seurat_clusters %in% c(13, 6, 12, 3, 5, 2, 0, 1)]
# Active override to refine selected sub-populations
selected_cells <- temp_sce[, temp_sce$seurat_clusters %in% c(13, 6, 12, 3, 11, 8)]

# Predefine functional gene signatures
inflammation_genes <- c("Nos2", "Cxcl3", "Il6", "Ccr2", "F10")
repair_genes       <- c("C1qc", "Gpnmb", "Fabp5", "Bhlhe41", "Trem2") 
MHC_genes          <- c("H2-DMa", "H2-Aa", "H2-Ab1", "H2-DMb1", "Ighm") 
Prolif_genes       <- c("Ccna2", "Hmmr", "Birc5", "Depdc1b", "Ccne2")

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

genes_use <- c("Ccna2", "Hmmr", "Birc5", "Depdc1b", "Ccne2")

# Extract expression matrix strictly for features of interest
expr_mat <- GetAssayData(selected_cells, slot = "data")[genes_use, , drop = FALSE]

# Align assigned pseudotime coordinates
pt <- selected_cells$pseudotime
pt <- pt[colnames(expr_mat)]   

# Pivot matrix to long format for ggplot logic compatibility
plot_df <- as.data.frame(t(as.matrix(expr_mat)))
plot_df$cell       <- rownames(plot_df)
plot_df$pseudotime <- pt[plot_df$cell]

plot_long <- plot_df %>%
  pivot_longer(cols = all_of(genes_use), names_to = "gene", values_to = "expression") %>%
  filter(!is.na(pseudotime)) %>%
  mutate(log_expr = log1p(expression))

# Generate LOESS smoothed trajectory profiles
ggplot(plot_long, aes(x = pseudotime, y = log_expr, color = gene)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.75, linewidth = 1.4) +
  theme_classic(base_size = 14) +
  labs(x = "Pseudotime", y = "Log(expression + 1)", color = NULL) +
  theme(legend.position = "right", axis.line = element_line(color = "black"))


# ==================================================================================================
# #### Figure 3S.B: Area Proportions (Excluding Late Timepoints) ####
# ==================================================================================================
`%!in%` <- function(x, table) { !(x %in% table) }

Idents(temp_sce) <- "seurat_clusters"
cell_data <- as.data.frame(prop.table(table(Idents(temp_sce), temp_sce$orig.ident), margin = 2))
cell_data <- as.data.frame(table(Idents(temp_sce), temp_sce$orig.ident), margin = 2)
colnames(cell_data) <- c("Cell_Type", "Timepoint", "Abundance")

# Restrict compositional analysis to early-stage phenotypes (0-7 days)
cell_data           <- subset(cell_data, Timepoint %!in% c("14day", "28day"))
cell_data$Cell_Type <- factor(cell_data$Cell_Type, levels = c(12, 5, 8, 3, 4, 1))

plot_comp <- ggplot(cell_data, aes(x = Timepoint, y = Abundance, fill = Cell_Type)) +
  geom_col() +
  labs(x = "Day", y = "Count", fill = "Clusters") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "grey")
  ) +
  scale_fill_manual(values = color) +  
  facet_wrap(~ Cell_Type, scales = "free_x", ncol = 6) +
  geom_area()

print(plot_comp)


# ==================================================================================================
# #### Figure 3S.C: BioTIP Network Partitioning & Macrophage Criticality ####
# ==================================================================================================
subsce_trans <- as.SingleCellExperiment(temp_sce)
samplesL     <- split(rownames(colData(subsce_trans)), f = colData(subsce_trans)$seurat_clusters)

# Model transcripts utilizing Poisson distribution for HVG extraction
dec.pois <- modelGeneVarByPoisson(subsce_trans)
hvg      <- getTopHVGs(dec.pois, n = 4000)
hvg      <- intersect(hvg, rownames(subsce_trans))
dat      <- subsce_trans[hvg, ]
logmat   <- as.matrix(logcounts(dat))
global   <- as.matrix(logcounts(subsce_trans))

cut.fdr <- 0.10
# Execute BioTIP network partition logic
igraphL <- getNetwork(testres, fdr = cut.fdr)
cluster <- getCluster_methods(igraphL)

# Subgraph Extraction (Community 3 specific logic)
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

nodes_in_community <- V(tmp)$name[V(tmp)$community == 3]
subgraph           <- induced_subgraph(tmp, nodes_in_community)

plot(subgraph, vertex.label = V(subgraph)$name, vertex.color = colrs[V(subgraph)$community], vertex.size = 6)

# Graph Topology Mapping (Nodes & Degrees)
node_degrees  <- degree(subgraph)
min_degree    <- min(node_degrees)
max_degree    <- max(node_degrees)

color_palette <- colorRampPalette(c("white", "#ED884C"))
node_colors   <- color_palette(100)[as.integer(100 * (node_degrees - min_degree) / (max_degree - min_degree)) + 1]

# Re-map absolute expression values to geometric gradients
local_sce <- subset(temp_sce, subset = seurat_clusters == "4")
local_sce <- as.SingleCellExperiment(local_sce)
local_sce <- as.matrix(logcounts(local_sce))

gene_sd  <- apply(local_sce, 1, sd)[V(subgraph)$name]
gene_exp <- apply(local_sce, 1, mean)[V(subgraph)$name]

min_exp <- min(gene_exp)
max_exp <- max(gene_exp)

node_size   <- gene_sd * 10
node_colors <- color_palette(100)[as.integer(100 * (gene_exp - min_exp) / (max_exp - min_exp))]

# Implement Fruchterman-Reingold Force-Directed algorithm
layout <- layout_with_fr(subgraph)

plot(subgraph,
     vertex.label       = V(subgraph)$name,
     vertex.color       = node_colors,
     vertex.size        = node_size * 3,
     vertex.label.cex   = 0.8,
     layout             = layout,
     vertex.label.color = "black")

# Export validated edge topological features
edges <- as_edgelist(subgraph)

write.table(edges, "fibroblast_4edges.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(as.data.frame(gene_sd), "fibroblast_4gene_sd.tsv", sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)
write.table(edges, "fibroblast_6edges.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(as.data.frame(gene_sd), "fibroblast_6gene_sd.tsv", sep = "\t", quote = TRUE, row.names = TRUE, col.names = TRUE)

# --------------------------------------------------------------------------------------------------
# #### Critical Transition Signals (CTSs) and Permutation Constraints ####
# --------------------------------------------------------------------------------------------------
# Putative Critical Transition Signals (CTSs) Identification by MCI score
membersL <- getMCI(cluster, testres, adjust.size = FALSE, fun = 'BioTIP')



cut.minsize = 60
par(oma = c(0, 0, 0, 0)) 
par(mar = c(1, 1, 1, 1)) 
plotBar_MCI(membersL, ylim = c(0, 20), minsize = cut.minsize)
dev.off()

# Extract top scoring DNB modular architectures
subg            <- induced_subgraph(subgraph, vids = V(subgraph)[name %in% V(subgraph)])
edgelist_sub    <- as_edgelist(subgraph, names = TRUE)
edge_attrs_sub  <- edge_attr(subgraph)
result_internal <- data.frame(edgelist_sub, edge_attrs_sub, check.names = FALSE)
head(result_internal) 

membersL$sd$`5`[["3"]]

topMCI         <- getTopMCI(membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]], min = cut.minsize, n = 12)
maxMCIms       <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], min = cut.minsize, n = 12)
maxMCI         <- getMaxStats(membersL[['MCI']], maxMCIms[['idx']])
CTS.Lib        <- getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])
CTS.Lib.Symbol <- CTS.Lib

maxMCI <- sort(maxMCI, decreasing = TRUE)
maxMCI <- head(maxMCI, 12)

# Vector normalization and array exportation
max_length <- max(lengths(CTS.Lib))
df <- data.frame(matrix(NA, nrow = max_length, ncol = length(CTS.Lib)))
for (i in 1:length(CTS.Lib)) {
  col_name <- names(CTS.Lib)[i]
  col_data <- CTS.Lib[[i]]
  df[, i]  <- c(col_data, rep(NA, max_length - length(col_data)))
}
colnames(df) <- names(CTS.Lib)
write.csv(df, file = "Macro_DNB.csv", row.names = TRUE)

names(CTS.Lib.Symbol)

# Empirical Shrinkage estimation
M <- cor.shrink(logmat[, unlist(samplesL)], Y = NULL, MARGIN = 1, shrink = TRUE)

# Execute robust background permutation testing (Bootstrapping C=500 passes)
C       = 500
simuMCI = list()
set.seed(2020)

for (i in 1:length(CTS.Lib)){
  n            <- length(CTS.Lib[[i]])
  simuMCI[[i]] <- simulationMCI(n, samplesL, logmat, B = C, fun = "BioTIP", M = M)
}

dev.off()
par(mfrow = c(5, 1))
for (i in 1:length(CTS.Lib)){
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las = 2,
                      main = paste0("Cluster ", names(maxMCI)[i], "; ",
                                    length(CTS.Lib[[i]]), " genes", "\n", "vs. ",
                                    "1000 times of gene-permutation"),
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

i = 1
sel_clusters <- c("12", "5", "8", "3", "4", "1")
maxMCI_sub   <- maxMCI[i][sel_clusters]
simuMCI_sub  <- simuMCI[[i]][sel_clusters, ]

plot_MCI_Simulation(maxMCI[i], simuMCI_sub, las = 2,
                    main = paste0("Cluster ", names(maxMCI)[i], "; ",
                                  length(CTS.Lib[[i]]), " genes", "\n", "vs. ",
                                  "1000 times of gene-permutation"),
                    which2point = names(maxMCI)[i])

dev.off()



# Confirm analytical tipping points employing Information Criterion (IC) and Delta-IC matrices
C             = 500
BioTIP_scores <- SimResults_g <- list()
set.seed(101010)

for(i in 1:length(CTS.Lib)){
  CTS                 <- CTS.Lib.Symbol[[i]]
  n                   <- length(CTS)
  BioTIP_scores[[i]]  <- getIc(logmat[, unlist(samplesL)], samplesL, CTS, fun = "BioTIP", shrink = TRUE, PCC_sample.target = 'none')
  SimResults_g[[i]]   <- simulation_Ic(n, samplesL, logmat, B = C, fun = "BioTIP", shrink = TRUE, PCC_sample.target = 'none')
}
names(BioTIP_scores) <- names(SimResults_g) <- names(CTS.Lib)

ylim = 1
i    = 1
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

# Mask visualizations strictly to target sub-clusters
sel_clusters     <- c("12", "5", "8", "3", "4", "1")
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

# Final Trajectory Assessment Render
x <- c(12, 5, 8, 3, 4, 1)
y <- c(1.108058, 2.497589, 1.725561, 1.959227, 2.312956, 2.388201)

# Enforce strict discrete numerical ordering to prevent plotting artifacts
x <- factor(x, levels = c(12, 5, 8, 3, 4, 1))

plot(x, y, type = "line", col = "black", xlab = "Cluster", ylab = "DNB_Score", main = "DNB_Score")
plot(x, y, type = "l", col = "black", xlab = "Cluster", ylab = "DNB_Score", main = "Selected Cell Scores")

save(list = c("testres", "igraphL", "cluster", "cut.fdr", "membersL", "topMCI", 
              "maxMCIms", "maxMCI", "CTS.Lib", "CTS.Lib.Symbol", "df", "simuMCI", 
              "BioTIP_scores", "SimResults_g"), 
     file = "./大论文图和数据/FigureS3.巨噬附图/ALL_Macrophages_biotip20260316.RData", compress = TRUE)

# ==================================================================================================
# MODULE: Figure 3.D-G Macrophage Lineage - Correlation Topologies & Network Dynamics
# 
# DESCRIPTION: 
#   Quantifies Pearson Correlation Coefficients (PCC) and Standard Deviations (SD), 
#   generates 2D quadrant assessments, and constructs topologically synchronized 
#   ComplexHeatmap arrays to visualize DNB-to-Coexpression shifts across 
#   macrophage critical transitions (Clusters 5, 12, 8).
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# #### Figure 3.D: Quantitative Export (PCC & SD Matrices) ####
# --------------------------------------------------------------------------------------------------
# Extract expression matrix and cluster metadata
DNB_vec           <- CTS.Lib.Symbol[["5"]]
exprs_matrix      <- global
cluster_info      <- temp_sce$seurat_clusters
selected_clusters <- c(12, 5, 8, 3, 4, 1)

# Secure vector arrays against matrix dimensionality
exprs_matrix_subset <- exprs_matrix[DNB_vec, , drop = FALSE]

# Initialize result containers
pcc_results <- list()
sd_results  <- list()

# Calculate correlative vectors (PCC) and deviations (SD) per cluster
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
  pcc_file <- file.path(output_dir, paste0("macro_cluster_", cluster, "_PCC.csv"))
  write.csv(pcc_results[[as.character(cluster)]], pcc_file, row.names = FALSE)
  
  sd_file  <- file.path(output_dir, paste0("macro_cluster_", cluster, "_SD.csv"))
  write.csv(sd_results[[as.character(cluster)]], sd_file, row.names = FALSE)
}

cat("PCC and SD data matrices for Macrophage clusters successfully exported!\n")

# --------------------------------------------------------------------------------------------------
# #### Figure 3.E: 2D Quadrant Evaluation (Macro-specific) ####
# --------------------------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)

# Load aggregated DNB co-expression data for Macrophages
res_df <- read_csv("./大论文图和数据/TableS10.Macro_dnb_coexp_plot.csv")

# Transform significance metrics and assign aesthetic factors
res_df <- res_df %>%
  mutate(
    log_positive_p_value = -log10(positive_p_value + 1e-6),
    log_negative_p_value = -log10(negative_p_value + 1e-6),
    size_factor          = pmin(log_positive_p_value, log_negative_p_value),
    color_factor         = case_when(
      positive_p_value < 0.2 & negative_p_value < 0.2 ~ "red",  
      negative_p_value < 0.2                          ~ "green",
      positive_p_value < 0.2                          ~ "blue", 
      TRUE                                            ~ "grey"  
    )
  )

# Introduce plot coordinate scalars for visualization clarity
res_df2 <- res_df %>%
  mutate(
    x_plot = ifelse(negative_mean_log2FC > 0, negative_mean_log2FC * 1.4, negative_mean_log2FC),
    y_plot = ifelse(positive_mean_log2FC < 0, positive_mean_log2FC * 1.6, positive_mean_log2FC)
  )

label_df <- res_df2 %>% filter(color_factor != "grey")

# Render Quadrant evaluation plot
main_plot <- ggplot(res_df2, aes(x = x_plot, y = y_plot)) +
  geom_point(aes(color = color_factor, size = size_factor), alpha = 0.6) +
  geom_text_repel(
    data = label_df, aes(label = gene),
    size = 3.6, box.padding = 0.4, point.padding = 0.6, max.overlaps = Inf,
    color = "black", segment.color = "black", segment.size = 0.4, min.segment.length = 0, force = 2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  scale_size_continuous(range = c(3, 10), name = "-log10(P-value)") +
  scale_color_manual(
    values = c("red" = "red", "green" = "green", "blue" = "blue", "grey" = "grey"),
    name   = "Significance",
    labels = c("Positive Sig.", "Negative Sig.", "Not Sig.", "Dual Sig.")
  ) +
  scale_x_continuous(name = "Inflammation Mean Log2FC", breaks = c(-1, 0, 1, 2), labels = c(-1, 0, 1, 2), sec.axis = dup_axis()) +
  scale_y_continuous(name = "Repair Mean Log2FC", breaks = c(-5, -2.5, 0, 2), labels = c(-5, -2.5, 0, 2)) +
  labs(title = "2D Quadrant Evaluation of Target Gene Regulation") +
  theme_minimal() +
  theme(
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 14),
    plot.title       = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(20, 20, 20, 20)
  ) +
  coord_cartesian(xlim = c(-2, 2.5), ylim = c(-5, 2), expand = FALSE)

print(main_plot)

# --------------------------------------------------------------------------------------------------
# #### ClusterGVis & Mfuzz Classification ####
# --------------------------------------------------------------------------------------------------
library(ClusterGVis)

# Synthesize unique co-expression gene universe
all_coexpr_genes <- unique(result_table$TARGET)
expr_mat         <- expr_mat_use[all_coexpr_genes, ]
expr_mat         <- expr_mat[apply(expr_mat, 1, function(x) sd(x, na.rm = TRUE) > 0), , drop = FALSE]

# Execute mfuzz soft clustering
ck <- clusterData(exp = expr_mat, cluster.method = "mfuzz", cluster.num = 5)

# Extract functional module labels
inflammation_genes <- ck$wide.res %>% filter(cluster %in% c(1, 3, 4)) %>% .[["gene"]]
Repair_genes       <- ck$wide.res %>% filter(cluster %in% c(2, 5))    %>% .[["gene"]]

result_table <- result_table %>%
  mutate(TYPE = case_when(
    TARGET %in% inflammation_genes ~ "Inflam",  
    TARGET %in% Repair_genes       ~ "Repair",      
    TRUE                           ~ "Other"                           
  ))

save(list = c("expr_mat", "ck", "result_table", "res_df", "expr_mat_use"), 
     file = "./大论文图和数据/FigureS3.巨噬附图/ALL_Macrophages_DNB_exp_20260317.RData", compress = TRUE)

# --------------------------------------------------------------------------------------------------
# #### Figure 3.G: Multilayer Sankey Network ####
# --------------------------------------------------------------------------------------------------
library(networkD3)
library(readxl)

# Target DNB genes for visual profiling
active_dnbs <- c("F10", "Tes", "Higd1a", "Prdx6", "Lyz2", "Rpl3", "Lyar")

for (dnb in unique_dnbs) {
  # Load Sankey topology metadata
  dnb_data <- read_xlsx("./大论文图和数据/ASII_Rd_Sankeyplot.xlsx")
  colnames(dnb_data)[2] <- "TARGET"; colnames(dnb_data)[5] <- "COR"
  
  # Filter interactions and establish hierarchical node sequence
  dnb_data       <- result_table[result_table$DNB %in% dnb, ] %>% filter(abs(Cor_values) > 0.1)
  dnb_gene_order <- dnb_data %>% count(TYPE, TARGET) %>% arrange(TYPE, desc(n)) %>% pull(TARGET)
  
  nodes <- data.frame(name = c(dnb, paste0("Cluster_", sort(unique(dnb_data$TYPE))), dnb_gene_order))
  
  # Layer 1 & 2 construction (Source -> Target -> Correlation Group)
  layer1 <- dnb_data %>% group_by(DNB, TYPE) %>% summarise(value = n(), .groups = 'drop') %>%
            mutate(source = match(DNB, nodes$name)-1, target = match(paste0("Cluster_", TYPE), nodes$name)-1, group = as.factor(paste0("Cluster_", TYPE)))
            
  layer2 <- dnb_data %>% mutate(source = match(paste0("Cluster_", TYPE), nodes$name)-1, target = match(TARGET, nodes$name)-1, value = 1, group = as.factor(COR))
  
  links_dnb <- bind_rows(layer1, layer2)
  
  # JavaScript color mapping for D3.js
  dnb_color_scale_js <- 'd3.scaleOrdinal().domain(["Positive", "positive", "Negative", "negative", "unknown"]).range(["red", "red", "green", "green", "gray"])'
  
  sankey_plot <- sankeyNetwork(Links = links_dnb, Nodes = nodes, Source = "source", Target = "target", 
                               Value = "value", NodeID = "name", LinkGroup = "group", colourScale = JS(dnb_color_scale_js))
  
  # Export HTML widget
  saveWidget(sankey_plot, file.path("/home/tsh/E-MTAB-7895/大论文图和数据/FigureS3.巨噬附图/", paste0("Figure3.F-2_Sankey_", dnb, ".html")))
}

# --------------------------------------------------------------------------------------------------
# #### Figure 3.D: Topologically Synchronized Correlation Heatmaps ####
# --------------------------------------------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)

DNB_vec          <- CTS.Lib.Symbol[["5"]]
clusters_to_plot <- c("5", "12", "8")
Idents(temp_sce) <- temp_sce$seurat_clusters

# Generate sub-cluster correlation matrices
cor_list <- lapply(clusters_to_plot, function(clu) {
  mat <- get_cluster_cor_mat(temp_sce, clu, DNB_vec, all_coexpr_genes_use, assay = "RNA", slot = "data")
  clean_cor_mat(mat)
})
names(cor_list) <- paste0("Cluster ", clusters_to_plot)

# Fix topology from Cluster 5 (Anchor Cluster)
first_nm  <- names(cor_list)[1]
first_mat <- cor_list[[first_nm]]
row_dend  <- as.dendrogram(hclust(dist(first_mat, method = "euclidean"), method = "ward.D2"))
col_dend  <- as.dendrogram(hclust(dist(t(first_mat), method = "euclidean"), method = "ward.D2"))

col_fun   <- colorRamp2(c(-0.4, -0.2, 0, 0.2, 0.4), c("#2166AC", "#1D68B5", "white", "#E31A1C", "#B2182B"))

ht_list <- NULL
for (nm in names(cor_list)) {
  mat      <- cor_list[[nm]]
  is_first <- (nm == first_nm)
  
  ht <- Heatmap(mat, name = "Pearson", col = col_fun,
                cluster_rows    = row_dend, 
                show_row_dend   = is_first,
                cluster_columns = col_dend, 
                show_column_dend = is_first,
                column_title    = nm, column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                show_row_names  = TRUE, row_names_gp = gpar(fontsize = 4),
                show_column_names = TRUE, column_names_gp = gpar(fontsize = 2),
                border          = TRUE)
  
  ht_list <- if (is.null(ht_list)) ht else ht_list + ht
}

# Final synchronized render
draw(ht_list, merge_legends = TRUE, heatmap_legend_side = "right")

# Output historical comparison check
head(cor_list[["Cluster 5"]][1:5, 1:5])
head(cor_list[["Cluster 12"]][1:5, 1:5])
head(cor_list[["Cluster 8"]][1:5, 1:5])


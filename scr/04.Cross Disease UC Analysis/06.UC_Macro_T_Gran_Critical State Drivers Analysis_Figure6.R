# ==================================================================================================
# PROJECT: Spatiotemporal Dynamics of Ulcerative Colitis (GSE193342 MERFISH)
# MODULE: Comprehensive Functional Utilities & External Dataset Integration
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# #### [SECTION 1] CUSTOM UTILITY FUNCTIONS ####
# --------------------------------------------------------------------------------------------------

# 1.1 Logical Operator: Not In
`%!in%` <- function(x, y) !(x %in% y)

# 1.2 Function: plot_cell_proportion
# Description: Generates stacked barplots visualizing relative cellular proportions across groups.
plot_cell_proportion <- function(sce,
                                 id_col = "annotation.2",
                                 group_col = "day",
                                 colors = NULL,
                                 cell_type_order = NULL) {
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(ggsci)
  
  meta_df <- sce@meta.data
  
  # Validate column existence
  if (!id_col %in% colnames(meta_df)) {
    stop(paste0("id_col not found: ", id_col))
  }
  if (!group_col %in% colnames(meta_df)) {
    stop(paste0("group_col not found: ", group_col))
  }
  
  # Data preparation and quantification
  plot_df <- meta_df %>%
    dplyr::select(all_of(c(id_col, group_col))) %>%
    dplyr::rename(cell_type = all_of(id_col),
                  group = all_of(group_col)) %>%
    dplyr::count(group, cell_type) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(proportion = n / sum(n)) %>%
    dplyr::ungroup()
  
  # Enforce factor ordering for cell types
  if (is.null(cell_type_order)) {
    cell_type_order <- unique(as.character(plot_df$cell_type))
  }
  plot_df$cell_type <- factor(plot_df$cell_type, levels = cell_type_order)
  
  # Enforce factor ordering for groups
  plot_df$group <- factor(plot_df$group, levels = unique(as.character(plot_df$group)))
  
  # Color mapping management
  if (is.null(colors)) {
    levs <- levels(plot_df$cell_type)
    colors <- setNames(ggsci::pal_d3("category20")(length(levs)), levs)
  } else {
    missing_cols <- setdiff(levels(plot_df$cell_type), names(colors))
    if (length(missing_cols) > 0) {
      stop(paste0("Missing colors for: ", paste(missing_cols, collapse = ", ")))
    }
    colors <- colors[levels(plot_df$cell_type)]
  }
  
  p <- ggplot(plot_df, aes(x = group, y = proportion, fill = cell_type)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = colors, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = c(0, 0)) +
    labs(x = group_col, y = "Cell proportion", fill = "Cell type") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.title = element_text(size = 13, color = "black"),
      legend.text = element_text(size = 11, color = "black")
    )
  
  return(p)
}

# 1.3 Function: plot_cell_counts
# Description: Generates faceted bar charts for absolute cell counts per cluster.
plot_cell_counts <- function(sce, 
                             id_col = "annotation",   
                             group_col = "day",       
                             palette = "npg") {       
  
  library(ggsci)
  Idents(sce) <- id_col
  
  # Ensure group variable is a vector
  group_vec <- sce@meta.data[[group_col]]
  
  # Compute absolute cell counts
  cell_data <- as.data.frame(table(Idents(sce), group_vec))
  colnames(cell_data) <- c("Cell_Type", "Group", "Count")
  
  # Chronological sorting
  cell_data <- arrange(cell_data, Group)
  
  # Consolidate color palettes for high-dimensional clusters
  scale_fill <- switch(palette,
                       "npg"     = scale_fill_npg(),
                       "lancet"  = scale_fill_lancet(),
                       "nejm"    = scale_fill_nejm(),
                       "jco"     = scale_fill_jco(),
                       scale_fill_npg())  
  
  npg_colors <- scales::hue_pal()(10)
  lancet_colors <- ggsci::pal_lancet()(9)
  nejm_colors <- ggsci::pal_nejm()(8)
  jco_colors <- ggsci::pal_jco()(10)
  combined_colors <- c(npg_colors, lancet_colors, nejm_colors, jco_colors)
  
  # Render faceted visualization
  p <- ggplot(cell_data, aes(x = Group, y = Count, fill = Cell_Type)) +
    geom_col() +
    labs(x = group_col, y = "Cell Count", fill = "Cell Type") +
    scale_fill_manual(values = combined_colors) +
    theme_minimal() +
    facet_wrap(~Cell_Type, scales = "free_y", ncol = 5) +   
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 14, color = "black"),
          legend.position = "none",  
          strip.text = element_text(size = 14, face = "bold"))
  
  return(p)
}

# 1.4 Function: plot_cell_alluvial
# Description: Generates temporal flow diagrams for cell type proportions.
plot_cell_alluvial <- function(sce,
                               id_col = "annotation.2",
                               group_col = "day",
                               colors = NULL,
                               cell_type_order = NULL) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggalluvial)
  library(scales)
  library(ggsci)
  
  meta_df <- sce@meta.data
  
  if (!id_col %in% colnames(meta_df)) {
    stop(paste0("id_col not found: ", id_col))
  }
  if (!group_col %in% colnames(meta_df)) {
    stop(paste0("group_col not found: ", group_col))
  }
  
  # Prepare proportional data table
  plot_df <- meta_df %>%
    dplyr::select(all_of(c(id_col, group_col))) %>%
    dplyr::rename(cell_type = all_of(id_col),
                  group = all_of(group_col)) %>%
    dplyr::count(group, cell_type) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(proportion = n / sum(n)) %>%
    dplyr::ungroup()
  
  # Define factor hierarchies
  if (is.null(cell_type_order)) {
    cell_type_order <- unique(as.character(plot_df$cell_type))
  }
  plot_df$cell_type <- factor(plot_df$cell_type, levels = cell_type_order)
  plot_df$group <- factor(plot_df$group, levels = unique(as.character(plot_df$group)))
  
  # Fill missing combinations for consistent flow rendering
  plot_df <- plot_df %>%
    tidyr::complete(group, cell_type, fill = list(n = 0, proportion = 0))
  
  # Color palette configuration
  if (is.null(colors)) {
    levs <- levels(plot_df$cell_type)
    colors <- setNames(ggsci::pal_d3("category20")(length(levs)), levs)
  } else {
    missing_cols <- setdiff(levels(plot_df$cell_type), names(colors))
    if (length(missing_cols) > 0) {
      stop(paste0("Missing colors for: ", paste(missing_cols, collapse = ", ")))
    }
    colors <- colors[levels(plot_df$cell_type)]
  }
  
  p <- ggplot(
    plot_df,
    aes(
      x = group,
      stratum = cell_type,
      alluvium = cell_type,
      y = proportion,
      fill = cell_type
    )
  ) +
    geom_alluvium(
      width = 0.5,
      alpha = 0.55,
      color = "white",
      linewidth = 0.4,
      curve_type = "sigmoid"
    ) +
    geom_stratum(
      width = 0.7,
      alpha = 0.9,
      color = "white",
      linewidth = 0.4
    ) +
    scale_fill_manual(values = colors, drop = FALSE) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = "Percentage of cells", fill = NULL) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 14, color = "black"),
      legend.text = element_text(size = 11, color = "black"),
      legend.title = element_blank()
    )
  
  return(p)
}

# 1.5 Function: Cluster_heatmap
# Description: Integrated pipeline for Cluster markers identification and GO/KEGG enrichment visualization.
Cluster_heatmap <- function(sce_object, cluster_col = "cell_name", top_n_markers = 20, 
                            pvalue_cutoff = 0.5, top_n_enrich = 10, output_pdf = "Cluster_heatmap.pdf") {
  
  library(Seurat)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(ggplot2)
  library(dplyr)
  library(ComplexHeatmap)
  library(ClusterGVis)
  library(sf)

  R.utils::setOption( "clusterProfiler.download.method",'auto' )
  
  # Data normalization and scaling
  sce_object <- NormalizeData(sce_object)
  sce_object <- ScaleData(sce_object)
  
  # Comprehensive marker identification
  Idents(sce_object) <- cluster_col
  all_markers <- FindAllMarkers(sce_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(all_markers, file = "all_markers.csv")
  
  # Filter top markers for visualization
  top_markers <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = top_n_markers, wt = avg_log2FC)
  
  # Prepare expression signatures for enrichment
  st_data <- prepareDataFromscRNA(object = sce_object, diffData = top_markers, showAverage = TRUE)
  
  # Perform Biological Process (BP) enrichment
  GO_enrich <- enrichCluster(object = st_data, OrgDb = org.Mm.eg.db, type = "BP", organism = "mmu",
                               pvalueCutoff = pvalue_cutoff, topn = top_n_enrich, seed = 5201314)
  
  # Perform KEGG pathway enrichment
  KEGG_enrich <- enrichCluster(object = st_data, OrgDb = org.Mm.eg.db, type = "KEGG", organism = "mmu",
                                 pvalueCutoff = pvalue_cutoff, topn = top_n_enrich, seed = 5201314)
  
  # Consolidate enrichment signatures
  enrich_all <- bind_rows(GO_enrich, KEGG_enrich)
  
  # Extract top 3 enriched terms per cluster for concise mapping
  GO_enrich_top3 <- GO_enrich %>% group_by(group) %>% top_n(3, wt = -pvalue)
  KEGG_enrich_top3 <- KEGG_enrich %>% group_by(group) %>% top_n(3, wt = -pvalue)
  enrich <- bind_rows(GO_enrich_top3, KEGG_enrich_top3)
  
  markGenes <- unique(top_markers$gene)
  write.csv(enrich, file = "enrich_results.csv")
  
  # Multi-panel PDF visualization
  pdf(output_pdf, height = 12, width = 14, onefile = FALSE)
  
  # Render expression trend lines
  visCluster(object = st_data, plot.type = "line")
  
  # Render marker heatmap with hierarchical ordering
  visCluster(object = st_data, plot.type = "heatmap", column_names_rot = 90, 
               markGenes = markGenes, cluster.order = c(1:15))
  
  # Composite visualization: Heatmap with enrichment annotations and group trend lines
  visCluster(object = st_data,
               plot.type = "both",
               column_names_rot = 45,
               show_row_names = FALSE,
               show_row_dend = FALSE,
               markGenes = markGenes,
               markGenes.side = "left",
               ht.col.list = list(col_range = c(-2, 0, 2), col_color = c("#85C1E9", "white", "#F1948A")),
               genes.gp = c('italic', fontsize = 12, col = "black"),
               annoTerm.data = enrich,
               line.side = "left",
               go.size = 8,
               add.box = TRUE,
               cluster.order = c(1:15),
               mline.col = ggsci::pal_d3("category20")(20)
  )
  dev.off()
}

# 1.6 Function: single_cell_DNB_analysis_pipeline (FULL UNABRIDGED LOGIC)
# Description: End-to-end automation of the BioTIP Dynamic Network Biomarker pipeline.
single_cell_DNB_analysis_pipeline <- function(sce, 
                                              cut_preselect = 0.15,  
                                              cut_fdr = 0.2,       
                                              cut_minsize = 50,     
                                              simulation_runs = 200,  
                                              seed = 2020) {
  library(BioTIP)
  # Data formatting: ensure factor levels are sanitized
  plasma_b_cell$anno_clusters <- droplevels(plasma_b_cell$anno_clusters)
  subsce_trans <- as.SingleCellExperiment(plasma_b_cell)
  samplesL <- split(rownames(colData(subsce_trans)), f = colData(subsce_trans)$anno_clusters)
  
  # Modeling transcript variability per Poisson distribution
  dec.pois <- modelGeneVarByPoisson(subsce_trans)
  hvg <- getTopHVGs(dec.pois, n = 4000)
  hvg <- intersect(hvg, rownames(subsce_trans))
  dat <- subsce_trans[hvg, ]
  logmat <- as.matrix(logcounts(dat))
  
  # Step 1: Feature SD Optimization via BioTIP
  set.seed(seed)
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  testres <- optimize.sd_selection(logmat, samplesL, B = 100, cutoff = 0.15, times = .75, percent = 0.8)
  
  # Step 2: Correlation-based Network Partitioning
  igraphL <- getNetwork(testres, fdr = 0.15)
  cluster <- getCluster_methods(igraphL)
  membersL <- getMCI(cluster, testres, adjust.size = FALSE, fun = 'BioTIP')
  
  # Step 3: MCI Score Quantification across Cellular Clusters
  par(oma = c(0, 0, 0, 0))  
  par(mar = c(1, 1, 1, 1)) 
  pdf('131 Pseudotime heatmap.pdf', height = 6, width = 12, onefile = FALSE)
  plotBar_MCI(membersL, ylim = c(0, 30), minsize = 60)
  dev.off()

  # Step 4: Identification of Putative Critical Transition Signals (CTSs)
  topMCI  <- getTopMCI(membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]], min = 70, n = 4)
  maxMCIms <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], min = 70, n = 4)
  maxMCI <- getMaxStats(membersL[['MCI']], maxMCIms[['idx']])
  CTS.Lib <- getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])
  
  # Step 5: Archiving DNB results
  maxMCI <- sort(maxMCI, decreasing = TRUE)
  maxMCI <- head(maxMCI, 4)
  
  max_length <- max(lengths(CTS.Lib))
  df <- data.frame(matrix(NA, nrow = max_length, ncol = length(CTS.Lib)))
  for (i in 1:length(CTS.Lib)) {
    col_name <- names(CTS.Lib)[i]
    col_data <- CTS.Lib[[i]]
    df[, i] <- c(col_data, rep(NA, max_length - length(col_data)))
  }
  colnames(df) <- names(CTS.Lib)
  write.csv(df, file = "DNB.csv", row.names = TRUE)
  
  # Step 6: Permutation Simulation for MCI Score Validation
  M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
  simuMCI <- list()
  set.seed(seed)
  
  for (i in 1:length(CTS.Lib)) {
    n <- length(CTS.Lib[[i]])
    simuMCI[[i]] <- simulationMCI(n, samplesL, logmat, B = 1000, fun = "BioTIP", M = M)
  }
  
  # Step 7: Visualization of MCI Simulation against background noise
  pdf('141 Pseudotime heatmap.pdf', height = 2, width = 12, onefile = FALSE)
  par(mfrow = c(1, 4))
  for (i in 1:length(CTS.Lib)) {
    plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las = 2,
                        main = paste0("Cluster ", names(maxMCI)[i], "; ", 
                                      length(CTS.Lib[[i]]), " genes", "\n", "vs. ", 
                                      "100 times of gene-permutation"),
                        which2point = names(maxMCI)[i])
  }
  dev.off()

  # Step 8: Calculation of BioTIP Index (Ic) for tipping point signals
  BioTIP_scores <- SimResults_g <- list()
  set.seed(101010)
  CTS.Lib.Symbol <- CTS.Lib
  for (i in 1:length(CTS.Lib)) {
    CTS <- CTS.Lib.Symbol[[i]]
    n <- length(CTS)
    BioTIP_scores[[i]] <- getIc(logmat, samplesL, CTS, fun = "BioTIP", shrink = TRUE, PCC_sample.target = 'none')
    SimResults_g[[i]]  <- simulation_Ic(n, samplesL, logmat, B = simulation_runs, fun = "BioTIP", shrink = TRUE, PCC_sample.target = 'none')
  }
  
  # Step 9: Final BioTIP Score / Ic.shrink visualizations
  pdf('151 Pseudotime heatmap.pdf', height = 4, width = 12, onefile = FALSE)
  par(mfrow = c(2, 4))
  for (i in 1:length(BioTIP_scores)) {
    plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], las = 2, ylab = "Ic.shrink", ylim = c(0, 1),
                       main = paste("Cluster ", names(CTS.Lib)[i], "_", length(CTS.Lib[[i]]), "genes", "\n", "vs.", 
                                    "500 gene-permutations"), fun = "matplot", which2point = names(CTS.Lib)[i])
    plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], main = paste("Delta Ic*", length(CTS.Lib[[i]]), "genes"), ylab = NULL)
  }
  dev.off()

  return(list(CTS.Lib = CTS.Lib, BioTIP_scores = BioTIP_scores, SimResults_g = SimResults_g))
}

# --------------------------------------------------------------------------------------------------
# #### [SECTION 2] Figure 6.A: EXTERNAL DATASET INGESTION & UMAP PROJECTION ####
# --------------------------------------------------------------------------------------------------
library(Seurat)
library(Matrix)
library(ggsci)

data_dir <- "/home/tsh/E-MTAB-7895/UC_MERFISH/"

# 2.1 Sparse Expression Matrix Loading
expr <- readMM(file.path(data_dir, "GSE193342_exprsData.mtx.gz"))

# 2.2 Feature & Metadata Loading
genes <- read.table(file.path(data_dir, "GSE193342_rowData.txt.gz"), header = TRUE, sep = "\t")
cells <- read.table(file.path(data_dir, "GSE193342_colData.txt.gz"), header = TRUE, sep = "\t")

# 2.3 Feature Matrix Sanitization (Row/Column Alignment)
rownames(expr) <- make.unique(genes$gene_short_name)  
colnames(expr) <- rownames(cells)

# 2.4 Seurat Object Initialization
sce <- CreateSeuratObject(counts = expr, meta.data = cells)

# 2.5 MERFISH Spatial / UMAP Coordinate Integration
umap.coords <- as.matrix(sce@meta.data[, c("X1", "X2")])
colnames(umap.coords) <- c("UMAP_1", "UMAP_2")
rownames(umap.coords) <- rownames(sce@meta.data)

sce[["umap"]] <- CreateDimReducObject(embeddings = umap.coords, key = "UMAP_", assay = DefaultAssay(sce))

# 2.6 Visualization & Temporal Composition Analysis
# Visualize global clusters via category20 palette
DimPlot(sce, reduction = "umap", group.by = "annotation.2", label = TRUE, repel = TRUE, pt.size = 1, raster = TRUE) +
  scale_color_d3(palette = "category20") +  
  NoLegend() +
  theme_minimal(base_size = 14) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

# Analyze relative proportions across days (Excluding recovery timepoints: 7, 11, 18)
plot_cell_proportion(subset(sce, subset = day %!in% c(7, 11, 18)), 
                     id_col = "annotation.2", 
                     group_col = "day")

# ==================================================================================================
# MODULE: External Dataset Validation (GSE193342 MERFISH) & Immune Lineage Refinement
# 
# DESCRIPTION: 
#   Performs high-resolution subsetting of the immune compartment, re-clustering 
#   of macrophage populations, and validates critical transitions (DNB) using 
#   the BioTIP algorithm on external Ulcerative Colitis (UC) spatial data.
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# #### [SECTION 1] Immune Compartment Extraction & Global Projection ####
# --------------------------------------------------------------------------------------------------
# Subset immune cells and exclude recovery timepoints (Day 7, 11, 18) to focus on acute dynamics
sce_immune <- subset(sce, subset = annotation == "Immune") %>%
              subset(subset = day %!in% c(7, 11, 18))

# Standard Preprocessing Pipeline
sce_immune <- NormalizeData(sce_immune) %>%
              FindVariableFeatures() %>%
              ScaleData(features = rownames(.)) %>%
              RunPCA() %>%
              FindNeighbors(dims = 1:30, reduction = "pca") %>%
              FindClusters(resolution = 1) %>%
              RunUMAP(dims = 1:20)

# Configuration: Visual Palettes and Biological Ordering
celltype_order <- c("Macrophage", "B", "T", "Neutrophil", "Dendritic", "Plasma", "Mast")
celltype_cols  <- c(
  "Macrophage" = "#C44E52", "B"         = "#4C72B0", "T"          = "#DD8452",
  "Neutrophil" = "#55A868", "Dendritic" = "#8172B3", "Plasma"     = "#64B5CD", 
  "Mast"       = "#937860"
)

# Figure 6.B & C: Dimensionality Reduction and Alluvial Flows
p_umap_global <- DimPlot(sce_immune, reduction = "umap", group.by = "annotation.2", 
                         label = TRUE, repel = TRUE, pt.size = 1, raster = TRUE) +
                 scale_color_manual(values = celltype_cols) + NoLegend() + 
                 theme_classic(base_size = 14) + 
                 theme(axis.text = element_blank(), axis.ticks = element_blank(), 
                       axis.line = element_blank(), panel.grid = element_blank())

p_alluvial <- plot_cell_alluvial(sce = sce_immune, id_col = "annotation.2", group_col = "day", 
                                 colors = celltype_cols, cell_type_order = celltype_order)

p_umap_split <- DimPlot(sce_immune, reduction = "umap", group.by = "annotation.2", split.by = "day",
                        label = TRUE, repel = TRUE, pt.size = 1, raster = TRUE, ncol = 8) +
                scale_color_manual(values = celltype_cols) + theme_classic(base_size = 14) +
                theme(axis.text = element_blank(), axis.ticks = element_blank(), 
                      axis.line = element_blank(), panel.grid = element_blank(), 
                      legend.title = element_blank(), strip.background = element_blank(), 
                      strip.text = element_text(size = 12, face = "bold"))

# Figure 6.D: Lineage-Specific Marker Visualization
Idents(sce_immune) <- "annotation.2"
marker_genes <- c("Lyz2", "Adgre1", "Cd68", "Csf1r", "Cd79a", "Ms4a1", "Cd79b", "Cd3d", "Cd3e", 
                  "Trbc1", "S100a8", "S100a9", "Ly6g", "Itgax", "Cd74", "H2-Ab1", "Jchain", 
                  "Mzb1", "Sdc1", "Kit", "Tpsb2", "Cpa3")

p_dot <- DotPlot(sce_immune, features = marker_genes[marker_genes %in% rownames(sce_immune)], 
                 group.by = "annotation.2") +
         scale_color_gradientn(colors = c("#D9D9D9", "#9ECAE1", "#08519C")) +
         theme_classic(base_size = 14) +
         theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
               axis.text.y = element_text(color = "black", size = 12),
               axis.title  = element_blank())

# --------------------------------------------------------------------------------------------------
# #### [SECTION 2] Figure 6.E: Macrophage Sub-lineage Deep Dive ####
# --------------------------------------------------------------------------------------------------
# Isolate Macrophage compartment for high-resolution refinement
macrophage_cells <- subset(sce_immune, subset = annotation.2 == "Macrophage") %>%
                    NormalizeData() %>% FindVariableFeatures() %>% ScaleData(features = rownames(.)) %>%
                    RunPCA() %>% FindNeighbors(dims = 1:30)

# Multi-resolution clustering evaluation
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 2)) {
  macrophage_cells <- FindClusters(macrophage_cells, resolution = res)
}
macrophage_cells <- RunUMAP(macrophage_cells, dims = 1:20)
clustree(macrophage_cells, prefix = "RNA_snn_res.")

# Functional Annotation and QC Filtering
# Strategy: Exclude clusters 10, 11, 12, 13 due to low cell count (<50) and ambiguous markers
filter_macrophage_cells <- subset(macrophage_cells, RNA_snn_res.0.8 %!in% c("10", "11", "12", "13"))
filter_macrophage_cells$RNA_snn_res.0.8 <- droplevels(filter_macrophage_cells$RNA_snn_res.0.8)
Idents(filter_macrophage_cells)        <- "RNA_snn_res.0.8"

# Final Macrophage UMAP and Cluster Composition
DimPlot(filter_macrophage_cells, reduction = "umap", group.by = "RNA_snn_res.0.8", label = TRUE) + 
  scale_color_d3(palette = "category20") + theme_minimal(base_size = 14) + NoLegend()

plot_cell_counts(filter_macrophage_cells, id_col = "RNA_snn_res.0.8", group_col = "day", palette = "jco")

# --------------------------------------------------------------------------------------------------
# #### [SECTION 3] Figure 6.D (BioTIP): Tipping Point Analysis (External Data) ####
# --------------------------------------------------------------------------------------------------
# Description: Identification of Critical Transition Signals (CTSs) in UC Macrophage lineage.
library(BioTIP)

# 3.1 Data Preparation for BioTIP (SingleCellExperiment conversion)
subsce_trans <- as.SingleCellExperiment(filter_macrophage_cells)
samplesL     <- split(rownames(colData(subsce_trans)), f = colData(subsce_trans)$RNA_snn_res.0.8)

# 3.2 Feature Selection via Poisson Modeling
dec.pois <- modelGeneVarByPoisson(subsce_trans)
hvg      <- getTopHVGs(dec.pois, n = 4000)
hvg      <- intersect(hvg, rownames(subsce_trans))
logmat   <- as.matrix(logcounts(subsce_trans[hvg, ]))

# 3.3 Network Partitioning and MCI Scoring
set.seed(2020)
testres  <- optimize.sd_selection(logmat, samplesL, B = 100, cutoff = 0.15, times = .75, percent = 0.8)
igraphL  <- getNetwork(testres, fdr = 0.15)
cluster  <- getCluster_methods(igraphL)
membersL <- getMCI(cluster, testres, adjust.size = FALSE, fun = 'BioTIP')

# 3.4 MCI Score Visualization
plotBar_MCI(membersL, ylim = c(0, 30), minsize = 60)

# 3.5 Top Module Identification (Top 9 candidates)
topMCI         <- getTopMCI(membersL[["members"]], membersL[["MCI"]], min = 70, n = 9)
maxMCIms       <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], min = 70, n = 9)
maxMCI         <- getMaxStats(membersL[['MCI']], maxMCIms[['idx']])
CTS.Lib        <- getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])
CTS.Lib.Symbol <- CTS.Lib

# Archive DNB gene lists
df_dnb <- data.frame(matrix(NA, nrow = max(lengths(CTS.Lib)), ncol = length(CTS.Lib)))
for (i in 1:length(CTS.Lib)) {
  df_dnb[, i] <- c(CTS.Lib[[i]], rep(NA, nrow(df_dnb) - length(CTS.Lib[[i]])))
}
colnames(df_dnb) <- names(CTS.Lib)
write.csv(df_dnb, file = "External_DNB_Targets.csv", row.names = TRUE)

# 3.6 Significance Simulation (Bootstrapping C=500 passes)
M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
simuMCI <- lapply(CTS.Lib, function(cts) simulationMCI(length(cts), samplesL, logmat, B = 500, fun = "BioTIP", M = M))

# 3.7 MCI Simulation vs Null Distribution plots
par(mfrow = c(1, 4))
for (i in 1:length(CTS.Lib)) {
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las = 2, which2point = names(maxMCI)[i],
                      main = paste0("Cluster ", names(maxMCI)[i], "; Permutations: 500"))
}

# 3.8 Delta-IC Confirmation of Tipping Points
BioTIP_scores <- SimResults_g <- list()
for (i in 1:length(CTS.Lib)) {
  n <- length(CTS.Lib[[i]])
  BioTIP_scores[[i]] <- getIc(logmat, samplesL, CTS.Lib[[i]], fun = "BioTIP", shrink = TRUE)
  SimResults_g[[i]]  <- simulation_Ic(n, samplesL, logmat, B = 500, fun = "BioTIP", shrink = TRUE)
}
names(BioTIP_scores) <- names(SimResults_g) <- names(CTS.Lib)

# 3.9 Final DNB Score and Trajectory Visualization
# Focus on active developmental sub-trajectory: Clusters 6 -> 1 -> 3
sel_clusters     <- c("6", "1", "3")
score_sub        <- BioTIP_scores[[4]][sel_clusters]
SimResults_g_sub <- SimResults_g[[4]][sel_clusters, ]

plot_Ic_Simulation(score_sub, SimResults_g_sub, las = 2, ylab = "Ic.shrink", 
                   main = paste("Developmental Axis [6-1-3]: Permutations: 500"))

x <- factor(c(6, 1, 3), levels = c(6, 1, 3))
y <- c(0, 3.197615, 2.290758)
DNB_count <- c(0, 131, 148)

plot(x, y, type = "b", col = "black", lwd = 2, xlab = "Cell Cluster", ylab = "DNB Criticality Score")
text(x, y, labels = DNB_count, pos = 3, col = "red", font = 2)

# Save final integrated project state
save(list = c("igraphL", "cluster", "membersL", "CTS.Lib", "BioTIP_scores", "SimResults_g"), 
     file = "./大论文图和数据/Figure6/Macro/UC_Macrophages_biotip20260320.RData", compress = TRUE)

# ==============================================================================
# MODULE: Figure 6 (Supplementary) - Trajectory Inference & Potency Analysis
# 
# DESCRIPTION: 
#   1. Monocle 3: Learns the macrophage differentiation trajectory, rooted 
#      at Cluster 6, using Seurat-consistent UMAP embeddings.
#   2. CytoTRACE2: Quantifies developmental potency and stemness across 
#      macrophage sub-lineages to validate trajectory directionality.
# ==============================================================================

# ------------------------------------------------------------------------------
# #### 1. Monocle 3 Trajectory Analysis ####
# ------------------------------------------------------------------------------
library(monocle3)
library(Seurat)
library(leidenbase)
library(ggsci)

# 1.1 Data Preparation for Monocle 3
# Ensure the active assay is set to RNA
DefaultAssay(filter_macrophage_cells) <- "RNA"

# Extract count matrix, metadata, and gene annotations
expression_matrix <- filter_macrophage_cells[["RNA"]]$counts
cell_metadata     <- data.frame(filter_macrophage_cells@meta.data)
gene_annotation   <- data.frame(gene_short_name = row.names(expression_matrix), 
                                 row.names = row.names(expression_matrix))

# Initialize the CellDataSet (CDS) object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# 1.2 Dimensionality Reduction & Graph Learning
# Standard log-normalization and PCA projection
cds <- preprocess_cds(cds, num_dim = 50, norm_method = "log")

# UMAP projection (default)
cds <- reduce_dimension(cds, reduction_method = "UMAP", cores = 5)

# Cluster cells using the Louvain algorithm
cds <- cluster_cells(cds, cluster_method = "louvain", reduction_method = "UMAP")

# Learn the principal graph of the trajectory
cds <- learn_graph(cds)

# 1.3 Root Selection & Trajectory Ordering
#' Helper Function: Automated root node identification
#' Identifies the closest vertex to the centroid of a specified cluster
myselect <- function(cds, select.classify, my_select) {
  cell_ids <- which(colData(cds)[, select.classify] == my_select)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
  return(root_pr_nodes)
}

# Define 'Cluster 6' as the developmental root (Starting point)
cds <- order_cells(cds, root_pr_nodes = myselect(cds, select.classify = 'RNA_snn_res.0.8', my_select = "6"))

# 1.4 Embedding Synchronization
# Integrate native Seurat UMAP coordinates into Monocle3 ReducedDims 
# to ensure spatial consistency across all figures.
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(filter_macrophage_cells, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed

# 1.5 Trajectory Visualization
# Combine Pseudotime and Cluster Identity in a multi-panel plot
p_traj <- plot_cells(cds, 
                     color_cells_by = "pseudotime",
                     show_trajectory_graph = FALSE) + 
          plot_cells(cds,
                     color_cells_by = "RNA_snn_res.0.8",
                     label_cell_groups = FALSE,
                     label_leaves = FALSE,
                     label_branch_points = FALSE,
                     graph_label_size = 1) +  
          scale_color_d3(palette = "category20") + 
          NoLegend() +
          theme_minimal(base_size = 14) +
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.grid = element_blank())

print(p_traj)

# ------------------------------------------------------------------------------
# #### 2. CytoTRACE2 Potency Analysis ####
# ------------------------------------------------------------------------------
library(CytoTRACE2)
library(magrittr)

# 2.1 Calculate Potency Scores
# Description: Assesses cell state plasticity based on transcriptomic diversity
cytotrace2_result <- cytotrace2(filter_macrophage_cells,
                                # species = "mouse", # Default is mouse
                                slot_type = "counts", 
                                is_seurat = TRUE, 
                                ncores = 5,
                                parallelize_models = FALSE,
                                parallelize_smoothing = FALSE)

# 2.2 Plotting and Export
# Initialize annotation frame mapping phenotypical labels to barcodes
annotation <- data.frame(phenotype = filter_macrophage_cells@meta.data$RNA_snn_res.0.8) %>% 
              set_rownames(colnames(filter_macrophage_cells))

# Generate comprehensive potency visualization suite
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation, 
                  is_seurat = TRUE)

# View specific CytoTRACE2 outputs
plots$CytoTRACE2_Potency_UMAP   # Spatial potency mapping
plots$CytoTRACE2_Boxplot_byPheno # Quantified potency by cluster
plots$CytoTRACE2_Relative_UMAP  # Relative differentiation state

# Final device cleanup
# dev.off()

##### Figure6-F #### 
# 提取表达矩阵
global <- as.matrix(logcounts(as.SingleCellExperiment(filter_macrophage_cells)))
DNB_vec <- CTS.Lib.Symbol[["1"]]
exprs_matrix <- global
# 提取聚类信息
cluster_info <- filter_macrophage_cells$RNA_snn_res.0.8
selected_clusters <- c(6, 1, 3)
# 确保基因名在表达矩阵中存在
exprs_matrix_subset <- exprs_matrix[DNB_vec, , drop = FALSE]
# 初始化PCC和SD结果存储
pcc_results <- list()
sd_results <- list()
# 计算PCC和SD
for (cluster in selected_clusters) {
  # 获取该簇内的细胞
  cluster_cells <- which(cluster_info == cluster)
  
  # 提取该簇内细胞的表达数据
  cluster_exprs <- exprs_matrix_subset[, cluster_cells, drop = FALSE]
  
  # 计算皮尔逊相关系数矩阵
  pcc_matrix <- cor(t(cluster_exprs), method = "pearson")
  
  # 转换PCC矩阵为边的关系
  pcc_df <- as.data.frame(as.table(pcc_matrix))
  colnames(pcc_df) <- c("Gene1", "Gene2", "PCC")
  
  # 只保留非自相关项（即去掉对角线上的数据）
  pcc_df <- pcc_df[pcc_df$Gene1 != pcc_df$Gene2, ]
  
  # 将计算的PCC结果保存到列表中
  pcc_results[[as.character(cluster)]] <- pcc_df
  
  # 计算每个基因的标准差
  gene_sd <- apply(cluster_exprs, 1, sd)
  
  # 将标准差结果存储到列表中
  sd_results[[as.character(cluster)]] <- data.frame(gene = names(gene_sd), sd = gene_sd)
}

# 将PCC和SD结果导出为CSV文件
output_dir <- "./Figure6/Macro/pcc_sd_results"  # 结果文件夹
dir.create(output_dir, showWarnings = FALSE)  # 创建文件夹（如果不存在）

# 批量导出PCC结果
for (cluster in selected_clusters) {
  # 导出PCC数据
  pcc_file <- file.path(output_dir, paste0("cluster_", cluster, "_PCC.csv"))
  write.csv(pcc_results[[as.character(cluster)]], pcc_file, row.names = FALSE)
  
  # 导出SD数据
  sd_file <- file.path(output_dir, paste0("cluster_", cluster, "_SD.csv"))
  write.csv(sd_results[[as.character(cluster)]], sd_file, row.names = FALSE)
}

# 输出结果
cat("PCC 和 SD 数据已成功导出！")


# ==================================================================================================
# MODULE: Figure 6 Macro - Trajectory Inference & Functional Co-expression Profiling
# 
# DESCRIPTION: 
#   1. Monocle 2: Constructs a DDRTree-based trajectory for macrophage subsets (Clusters 6, 1, 3).
#   2. Compositional Analysis: Quantifies temporal shifts in cellular abundance via area plots.
#   3. Co-expression Engine: Implements a statistically robust framework to identify 
#      top-correlated genes for DNB module validation across temporal trajectories.
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# #### [SECTION 1] Monocle 2 Trajectory Reconstruction (Macro Subsets 6, 1, 3) ####
# --------------------------------------------------------------------------------------------------
library(monocle)
library(ggsci)

# 1.1 Lineage Subsetting
# Focus specifically on the active developmental axis identified in MERFISH clusters
Idents(filter_macrophage_cells) <- "RNA_snn_res.0.8"
temp_sce <- subset(filter_macrophage_cells, RNA_snn_res.0.8 %in% c(6, 1, 3))

# 1.2 Monocle 2 CellDataSet (CDS) Construction
# Extract raw counts and initiate sparse matrix structures
sparse_data <- as(as.matrix(temp_sce[["RNA"]]$counts), "sparseMatrix")
mdata       <- new("AnnotatedDataFrame", data = temp_sce@meta.data)
fData       <- data.frame(gene_short_name = row.names(sparse_data), row.names = row.names(sparse_data))
fd          <- new("AnnotatedDataFrame", data = fData)

# Instantiate CDS with Negative Binomial distribution modeling
monocle_cds <- newCellDataSet(cellData            = sparse_data,
                              phenoData           = mdata,
                              featureData         = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily    = negbinomial.size())

# 1.3 Pre-processing & Size Factor Estimation
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
# Filter for genes expressed in at least 3 cells
monocle_cds <- detectGenes(monocle_cds, min_expr = 3)

# 1.4 Feature Selection for Ordering
# Identify highly variable genes using Seurat's selection as a prior
expressed_genes <- VariableFeatures(temp_sce)
diff_test_res   <- differentialGeneTest(monocle_cds[expressed_genes, ],
                                        fullModelFormulaStr = "~RNA_snn_res.0.8")

# Enforce strict significance threshold for trajectory stability
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01)) 
monocle_cds    <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)

# 1.5 Manifold Learning (DDRTree) & Pseudotime Ordering
library(igraph)
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, norm_method = "log", reduction_method = "DDRTree")

# Calculate pseudotime and designate potential root states
monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds, root_state = "2")
monocle_cds <- orderCells(monocle_cds, root_state = "1") # Active root override

# 1.6 Visualizing Developmental Trajectories
# Render multiple trajectory projections: Pseudotime, Cluster Identity, and Split by Day
p  <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 0.6, show_backbone = TRUE, show_branch_points = FALSE) + 
      facet_wrap("~orig.ident", nrow = 1)

p0 <- plot_cell_trajectory(monocle_cds, color_by = "RNA_snn_res.0.8", cell_size = 0.5, show_backbone = TRUE, show_branch_points = TRUE) + 
      facet_wrap("~orig.ident", nrow = 1)

p1 <- plot_cell_trajectory(monocle_cds, color_by = "RNA_snn_res.0.8", cell_size = 0.6, show_backbone = FALSE, show_branch_points = FALSE) +
      scale_color_d3(palette = "category20")

plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 0.6, show_backbone = TRUE) +
  theme(legend.position = 'none', panel.border = element_blank())

# Complex trajectory tree visualization
p2 <- plot_complex_cell_trajectory(monocle_cds, x = 1, y = 2, color_by = "label") +
      theme(legend.title = element_blank())

plot_cell_trajectory(monocle_cds, color_by = "State", size = 1, show_backbone = TRUE)
p0 | p | p1 | p2

plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters", size = 1, show_backbone = TRUE) + 
  facet_wrap("~orig.ident", nrow = 1)

plot_complex_cell_trajectory(monocle_cds, x = 1, y = 2, color_by = "seurat_clusters") +
  theme(legend.title = element_blank())

# --------------------------------------------------------------------------------------------------
# #### [SECTION 2] Compositional Dynamics (Temporal Cellular Proportions) ####
# --------------------------------------------------------------------------------------------------
Idents(temp_sce) <- "RNA_snn_res.0.8"

# Quantify relative and absolute abundance per timepoint
cell_data_prop <- as.data.frame(prop.table(table(Idents(temp_sce), temp_sce$day), margin = 2))
cell_data_abs  <- as.data.frame(table(Idents(temp_sce), temp_sce$day))
colnames(cell_data_abs) <- c("Cell_Type", "Timepoint", "Abundance")

# Generate temporal area chart for lineage shifts
plot_composition <- ggplot(cell_data_abs, aes(x = Timepoint, y = Abundance, fill = Cell_Type)) +
                    geom_col() +
                    labs(x = "Timepoint", y = "Cell Count", fill = "Cell Type") +
                    theme_minimal() +
                    theme(
                      panel.background = element_rect(fill = "transparent", color = NA),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line        = element_line(colour = "grey")
                    ) +
                    facet_wrap(~ Cell_Type, scales = "free_x", ncol = 6) +
                    geom_area()

print(plot_composition)


# --------------------------------------------------------------------------------------------------
# #### [SECTION 3] Functional Co-expression Discovery Engine ####
# --------------------------------------------------------------------------------------------------
# Configure global parameters for correlation analysis
topN       <- 50
trajectory <- c("6", "1", "3")
seurat_obj <- temp_sce 
table(seurat_obj$RNA_snn_res.0.8)

#' Function: get_cluster_expr_matrix
#' Description: Extracts a summarized expression matrix (mean) across trajectory clusters.
get_cluster_expr_matrix <- function(genes, seurat_obj, trajectory) {
  expr_mat <- matrix(NA, nrow = length(genes), ncol = length(trajectory))
  rownames(expr_mat) <- genes
  colnames(expr_mat) <- as.character(trajectory)
  
  for (i in seq_along(trajectory)) {
    cells <- colnames(seurat_obj)[seurat_obj$seurat_clusters == trajectory[i]]
    if (length(cells) == 0) {
      expr_mat[, i] <- NA
    } else {
      # Use stabilized rowMeans calculation
      avg_expr <- rowMeans(as.matrix(seurat_obj[["RNA"]]$data[genes, cells, drop = FALSE]), na.rm = TRUE)
      expr_mat[, i] <- avg_expr
    }
  }
  return(expr_mat)
}

# Acquire global expression matrix (Stable data slot)
expr_all <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "data"))

#' Function: get_top_correlated
#' Description: Identifies top N co-expressed features for a target gene within a specific cluster.
#' Implements Pearson/Spearman correlation with optional p-value adjustment and sign balancing.
get_top_correlated <- function(
    gene,
    seurat_obj,
    cluster_id      = "5",
    top_n           = 50,
    pval_cutoff     = 0.05,
    allowed_genes   = NULL,
    balance_sign    = TRUE,          # Boolean: Ensure representation of both positive and negative correlations
    method          = c("pearson", "spearman"),
    adjust_p        = FALSE,         # Boolean: Apply multiple hypothesis correction
    p_adjust_method = "BH"           # FDR correction method (Benjamini-Hochberg)
) {
  method <- match.arg(method)
  
  # 1) Target cluster cell isolation and matrix extraction
  cells <- colnames(seurat_obj)[seurat_obj$seurat_clusters == cluster_id]
  if (length(cells) == 0L) {
    warning("Zero cell coverage for cluster_id = ", cluster_id)
    return(list(genes = gene, correlations = NA, signs = "Center"))
  }
  data_mat <- as.matrix(seurat_obj[["RNA"]]$data[, cells, drop = FALSE])
  
  # 2) Gene presence validation
  if (!gene %in% rownames(data_mat)) {
    warning("Gene not found in active RNA assay: ", gene)
    return(list(genes = gene, correlations = NA, signs = "Center"))
  }
  
  # 3) Iterative computation of correlation and significance metrics
  gene_expr <- data_mat[gene, ]
  cors      <- apply(data_mat, 1, function(x) suppressWarnings(cor(x, gene_expr, method = method)))
  pvals     <- apply(data_mat, 1, function(x) suppressWarnings(cor.test(x, gene_expr, method = method)$p.value))
  
  # 4) Statistical Significance Adjustment (Optional)
  p_use <- if (adjust_p) p.adjust(pvals, method = p_adjust_method) else pvals
  
  # 5) Multi-layered filtering: Significance, self-exclusion, and allowed-feature intersection
  valid_genes <- names(p_use)[p_use < pval_cutoff]
  valid_genes <- setdiff(valid_genes, gene)
  if (!is.null(allowed_genes)) valid_genes <- intersect(valid_genes, allowed_genes)
  
  # Prune NA correlation vectors
  valid_genes <- valid_genes[!is.na(cors[valid_genes])]
  if (length(valid_genes) == 0L) {
    warning("No valid features passed significance or identity filters.")
    return(list(genes = gene, correlations = NA, signs = "Center"))
  }
  
  # 6) Feature Selection (Top N)
  take_n <- min(top_n, length(valid_genes))
  
  if (!balance_sign) {
    # Ranking based on absolute correlation magnitude (|cor|)
    sorted    <- valid_genes[order(abs(cors[valid_genes]), decreasing = TRUE)]
    top_genes <- head(sorted, take_n)
  } else {
    # Balanced Sign Selection: Extract ceil(n/2) positive and floor(n/2) negative correlations
    k_pos <- take_n %/% 2 + (take_n %% 2)  
    k_neg <- take_n - k_pos                
    
    pos <- valid_genes[cors[valid_genes] > 0]
    neg <- valid_genes[cors[valid_genes] < 0]
    
    # Sort positive descending and negative by absolute magnitude descending
    pos_sorted <- pos[order(cors[pos], decreasing = TRUE)]
    neg_sorted <- neg[order(abs(cors[neg]), decreasing = TRUE)]
    
    pick_pos  <- head(pos_sorted, k_pos)
    pick_neg  <- head(neg_sorted, k_neg)
    top_genes <- c(pick_pos, pick_neg)
  }
  
  # 7) Structured Output Generation
  correlation_signs <- ifelse(cors[top_genes] > 0, "Positive",
                              ifelse(cors[top_genes] < 0, "Negative", "Zero"))
  
  return(list(
    genes        = c(gene, top_genes),
    correlations = c(NA, cors[top_genes]),
    pvalues      = c(NA, pvals[top_genes]), # Detailed p-value tracking
    signs        = c("Center", correlation_signs)
  ))
}

# ==================================================================================================
# MODULE: Figure 6 Macro - Dynamic Network Biomarker (DNB) Functional Validation
# 
# DESCRIPTION: 
#   Integrates Monocle 2 trajectory modeling, non-parametric permutation testing 
#   for co-expression modules, 2D quadrant regulation evaluation, and multi-layered 
#   network topology (Sankey & Heatmap) for critical transition drivers.
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# #### [SECTION 1] Monocle 2 Trajectory Modeling (External Validation) ####
# --------------------------------------------------------------------------------------------------
library(monocle)
library(igraph)
library(ggsci)

# 1.1 Object Initialization & Subset
Idents(filter_macrophage_cells) <- "RNA_snn_res.0.8"
temp_sce <- subset(filter_macrophage_cells, RNA_snn_res.0.8 %in% c(6, 1, 3))

# 1.2 Monocle 2 CDS Construction
sparse_data <- as(as.matrix(temp_sce[["RNA"]]$counts), "sparseMatrix")
mdata       <- new("AnnotatedDataFrame", data = temp_sce@meta.data)
fData       <- data.frame(gene_short_name = row.names(sparse_data), row.names = row.names(sparse_data))
fd          <- new("AnnotatedDataFrame", data = fData)

monocle_cds <- newCellDataSet(cellData            = sparse_data,
                              phenoData           = mdata,
                              featureData         = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily    = negbinomial.size())

# 1.3 Pre-processing & Ordering Feature Selection
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 3)

expressed_genes <- VariableFeatures(temp_sce)
diff_test_res   <- differentialGeneTest(monocle_cds[expressed_genes, ],
                                        fullModelFormulaStr = "~RNA_snn_res.0.8")

ordering_genes <- row.names(subset(diff_test_res, qval < 0.01)) 
monocle_cds    <- setOrderingFilter(monocle_cds, ordering_genes)

# 1.4 DDRTree Manifold Learning & Trajectory Optimization
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, norm_method = "log", reduction_method = "DDRTree")
monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds, root_state = "2") # Anchoring progenitor state

# 1.5 Trajectory Composition Plotting
cell_data_abs <- as.data.frame(table(Idents(temp_sce), temp_sce$day))
colnames(cell_data_abs) <- c("Cell_Type", "Timepoint", "Abundance")

plot_comp <- ggplot(cell_data_abs, aes(x = Timepoint, y = Abundance, fill = Cell_Type)) +
             geom_area(aes(group = Cell_Type), alpha = 0.4) + geom_col(width = 0.6) +
             facet_wrap(~ Cell_Type, scales = "free_x", ncol = 6) +
             theme_minimal() + theme(axis.line = element_line(colour = "grey"))

# --------------------------------------------------------------------------------------------------
# #### [SECTION 2] Non-parametric Permutation Test for DNB Sig-modules ####
# --------------------------------------------------------------------------------------------------

#' Function: calculate_mean_and_pvalue
#' Description: Validates the significance of observed log2FC in co-expression modules 
#'              using a direction-specific background permutation framework.
calculate_mean_and_pvalue <- function(coexpr_genes, Log2fc_vector, Bg_genes, n_perm = 1000) {
  observed_mean <- mean(Log2fc_vector[coexpr_genes], na.rm = TRUE)
  
  # Direction-specific background pooling
  bg_pool <- if (observed_mean >= 0) Bg_genes[Log2fc_vector[Bg_genes] >= 0] else Bg_genes[Log2fc_vector[Bg_genes] < 0]
  bg_pool <- bg_pool[abs(Log2fc_vector[bg_pool]) < 5 & !is.na(Log2fc_vector[bg_pool])]
  
  if (length(bg_pool) < length(coexpr_genes)) {
    warning("Insufficient background pool size.")
    return(list(observed_mean = observed_mean, p_value = NA))
  }
  
  permuted_means <- replicate(n_perm, mean(Log2fc_vector[sample(bg_pool, length(coexpr_genes))], na.rm = TRUE))
  p_val          <- if (observed_mean >= 0) mean(permuted_means >= observed_mean) else mean(permuted_means <= observed_mean)
  
  return(list(observed_mean = observed_mean, p_value = p_val))
}

# --------------------------------------------------------------------------------------------------
# #### [SECTION 3] 2D Quadrant Regulation Evaluation (Main Execution Loop) ####
# --------------------------------------------------------------------------------------------------
# Global background initialization
universe_genes <- rownames(Seurat::GetAssayData(seurat_obj))
DNB_vec        <- CTS.Lib[["1"]]
trajectory_id  <- c("6", "1", "3")

# Construct log2FC reference vector
expr_mat_use <- get_cluster_expr_matrix(rownames(seurat_obj), seurat_obj, trajectory_id)
log2fc_vector <- log2((expr_mat_use[, "3"] + 1e-6) / (expr_mat_use[, "6"] + 1e-6))
names(log2fc_vector) <- rownames(expr_mat_use)

# Execute core DNB evaluation loop
res_list     <- list()
result_table <- data.frame()

for (g in DNB_vec) {
  bg_genes <- setdiff(names(log2fc_vector), g)
  result   <- get_top_correlated(g, seurat_obj, cluster_id = "1", top_n = 50, allowed_genes = deg_genes)
  
  if (is.null(result) || length(result$genes) < 2) next
  
  top_genes <- result$genes[-1]
  pos_genes <- top_genes[result$signs[top_genes] == "Positive"]
  neg_genes <- top_genes[result$signs[top_genes] == "Negative"]
  
  pos_res <- calculate_mean_and_pvalue(pos_genes, log2fc_vector, bg_genes)
  neg_res <- calculate_mean_and_pvalue(neg_genes, log2fc_vector, bg_genes)
  
  res_list[[g]] <- data.frame(gene = g, n_top = length(top_genes), 
                              positive_mean_log2FC = pos_res$observed_mean, positive_p_value = pos_res$p_value,
                              negative_mean_log2FC = neg_res$observed_mean, negative_p_value = neg_res$p_value)
  
  result_table <- rbind(result_table, data.frame(DNB = g, TARGET = top_genes, COR = result$signs[-1], 
                                                 Cor_values = result$correlations[-1], Cor_log2fc = log2fc_vector[top_genes]))
}

res_df <- bind_rows(res_list) %>%
          mutate(log_pos_p = -log10(positive_p_value + 1e-6), log_neg_p = -log10(negative_p_value + 1e-6),
                 size_factor = pmin(log_pos_p, log_neg_p),
                 color_factor = case_when(positive_p_value < 0.2 & negative_p_value < 0.2 ~ "red",
                                          negative_p_value < 0.2                          ~ "green",
                                          positive_p_value < 0.2                          ~ "blue", 
                                          TRUE                                            ~ "grey"))

# --------------------------------------------------------------------------------------------------
# #### [SECTION 4] Multi-layer Sankey Topology (networkD3) ####
# --------------------------------------------------------------------------------------------------
library(networkD3)
library(htmlwidgets)

# Focus on validated significant DNBs
active_dnbs <- res_df$gene[res_df$color_factor != "grey"]

for (dnb in active_dnbs) {
  dnb_data <- result_table[result_table$DNB == dnb & abs(result_table$Cor_values) > 0.1, ]
  
  # Hierarchical node definition
  nodes <- data.frame(name = c(dnb, paste0("Cluster_", sort(unique(dnb_data$TYPE))), unique(dnb_data$TARGET)))
  
  # Edge linkage construction
  layer1 <- dnb_data %>% group_by(DNB, TYPE) %>% summarise(value = n(), .groups = 'drop') %>%
            mutate(source = match(DNB, nodes$name)-1, target = match(paste0("Cluster_", TYPE), nodes$name)-1, group = as.factor(TYPE))
            
  layer2 <- dnb_data %>% mutate(source = match(paste0("Cluster_", TYPE), nodes$name)-1, target = match(TARGET, nodes$name)-1, 
                                value = 1, group = as.factor(COR))
  
  links <- bind_rows(layer1, layer2)
  
  # JS Aesthetics injection
  color_scale <- 'd3.scaleOrdinal().domain(["Positive", "Negative", "unknown"]).range(["#C00000", "#2E75B6", "#A9A9A9"])'
  
  sankey <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target", Value = "value", 
                          NodeID = "name", LinkGroup = "group", colourScale = JS(color_scale), fontSize = 14)
  
  saveWidget(sankey, paste0("Figure6_Sankey_", dnb, ".html"))
}

# --------------------------------------------------------------------------------------------------
# #### [SECTION 5] Synchronized Expression Heatmaps ####
# --------------------------------------------------------------------------------------------------
# Construct row annotations based on functional mfuzz clusters
annotation_row <- data.frame(Correlation = dnb_data$COR, Cluster = as.factor(dnb_data$TYPE)) %>%
                  set_rownames(make.unique(dnb_data$TARGET, sep = "-"))

# Render pheatmap with forced row synchronization
pheatmap(heatmap_data, cluster_rows = FALSE, cluster_cols = FALSE, scale = "row",
         color = colorRampPalette(c("#2E75B6", "white", "#C00000"))(100),
         annotation_row = annotation_row)

# ==================================================================================================
# MODULE: Figure 6 T-Cell - Trajectory Inference, Criticality & Network Decoupling
# 
# DESCRIPTION: 
#   Integrates Monocle 2/3 trajectory modeling, BioTIP DNB criticality analysis, 
#   and a permutation-based bivariate correlation framework to identify drivers 
#   of T-cell state transitions (Proliferation vs. Exhaustion).
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# #### [SECTION 1] Monocle 2 Trajectory Reconstruction (T-Cell Lineage) ####
# --------------------------------------------------------------------------------------------------
library(monocle)
library(igraph)
library(ggsci)

# 1.1 Object Initialization & Subset
# Focusing on the active developmental axis: Clusters 6, 3, 5, 0, 2
Idents(filter_Tcell) <- "RNA_snn_res.0.3"
temp_sce <- subset(filter_Tcell, RNA_snn_res.0.3 %in% c(6, 3, 5, 0, 2))

# 1.2 Monocle 2 CDS Construction
sparse_data <- as(as.matrix(temp_sce[["RNA"]]$counts), "sparseMatrix")
mdata       <- new("AnnotatedDataFrame", data = temp_sce@meta.data)
fData       <- data.frame(gene_short_name = row.names(sparse_data), row.names = row.names(sparse_data))
fd          <- new("AnnotatedDataFrame", data = fData)

monocle_cds <- newCellDataSet(cellData            = sparse_data,
                              phenoData           = mdata,
                              featureData         = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily    = negbinomial.size())

# 1.3 Pre-processing & Ordering Feature Selection
monocle_cds <- estimateSizeFactors(monocle_cds) %>% estimateDispersions()
monocle_cds <- detectGenes(monocle_cds, min_expr = 3)

# Enforce strict q-value thresholding (q < 0.01) for trajectory stability
ordering_genes <- row.names(subset(differentialGeneTest(monocle_cds[VariableFeatures(temp_sce), ], 
                                                        fullModelFormulaStr = "~RNA_snn_res.0.3"), 
                                   qval < 0.01))
monocle_cds    <- setOrderingFilter(monocle_cds, ordering_genes)

# 1.4 Manifold Learning (DDRTree) & Pseudotime Execution
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, reduction_method = "DDRTree") %>%
               orderCells(root_state = "3") # Anchoring developmental root

# 1.5 Temporal Composition Analysis
cell_data_abs <- as.data.frame(table(Idents(temp_sce), temp_sce$day))
colnames(cell_data_abs) <- c("Cell_Type", "Timepoint", "Abundance")

plot_comp <- ggplot(cell_data_abs, aes(x = Timepoint, y = Abundance, fill = Cell_Type)) +
             geom_area(aes(group = Cell_Type), alpha = 0.4) + geom_col(width = 0.6) +
             facet_wrap(~ Cell_Type, scales = "free_x", ncol = 6) +
             theme_minimal() + theme(axis.line = element_line(colour = "grey"))

# --------------------------------------------------------------------------------------------------
# #### [SECTION 2] BioTIP Criticality Assessment (DNB Identification) ####
# --------------------------------------------------------------------------------------------------
library(BioTIP)

# 2.1 Data Preparation for BioTIP
subsce_trans <- as.SingleCellExperiment(temp_sce)
samplesL     <- split(rownames(colData(subsce_trans)), f = colData(subsce_trans)$RNA_snn_res.0.3)
samplesL     <- samplesL[lengths(samplesL) > 0]

# 2.2 Feature Pre-selection (Poisson modeling)
hvg <- getTopHVGs(modelGeneVarByPoisson(subsce_trans), n = 4000) %>% intersect(rownames(subsce_trans))
logmat <- as.matrix(logcounts(subsce_trans[hvg, ]))

# 2.3 MCI Optimization and DNB Extraction
# Identify top modules exhibiting maximum criticality signals (CTSs)
topMCI   <- getTopMCI(membersL[["members"]], membersL[["MCI"]], min = 70, n = 4)
maxMCIms <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], min = 70, n = 4)
maxMCI   <- getMaxStats(membersL[['MCI']], maxMCIms[['idx']])
CTS.Lib  <- getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])

# 2.4 Significance Validation & Trajectory Scoring
# Focus on active developmental axis: Clusters 6 -> 0 -> 2
BioTIP_scores <- lapply(CTS.Lib, function(cts) getIc(logmat, samplesL, cts, fun = "BioTIP", shrink = TRUE))
SimResults_g  <- lapply(CTS.Lib, function(cts) simulation_Ic(length(cts), samplesL, logmat, B = 500, fun = "BioTIP", shrink = TRUE))

# Final trajectory criticality render
x <- factor(c(6, 0, 2), levels = c(6, 0, 2))
y <- c(0, 4.1663762, 0.1273221)
plot(x, y, type = "b", lwd = 2, xlab = "Cell Cluster", ylab = "DNB Criticality Score", main = "T-Cell Tipping Point Analysis")

# --------------------------------------------------------------------------------------------------
# #### [SECTION 3] Bivariate Correlation Framework (DNB Functional Decoupling) ####
# --------------------------------------------------------------------------------------------------

#' Function: calculate_mean_and_pvalue
#' Description: Direction-specific background permutation for co-expression module validation.
calculate_mean_and_pvalue <- function(coexpr_genes, Log2fc_vector, Bg_genes, n_perm = 1000) {
  observed_mean <- mean(Log2fc_vector[coexpr_genes], na.rm = TRUE)
  bg_pool       <- if (observed_mean >= 0) Bg_genes[Log2fc_vector[Bg_genes] >= 0] else Bg_genes[Log2fc_vector[Bg_genes] < 0]
  bg_pool       <- bg_pool[abs(Log2fc_vector[bg_pool]) < 5 & !is.na(Log2fc_vector[bg_pool])]
  
  if (length(bg_pool) < length(coexpr_genes)) return(list(observed_mean = observed_mean, p_value = NA))
  
  permuted_means <- replicate(n_perm, mean(Log2fc_vector[sample(bg_pool, length(coexpr_genes))], na.rm = TRUE))
  p_val          <- if (observed_mean >= 0) mean(permuted_means >= observed_mean) else mean(permuted_means <= observed_mean)
  return(list(observed_mean = observed_mean, p_value = p_val))
}

# 3.1 Global Parameter Initialization
DNB_vec       <- CTS.Lib[["0"]]
expr_mat_use  <- get_cluster_expr_matrix(rownames(seurat_obj), seurat_obj, c("6", "0", "2"))
log2fc_vector <- log2((expr_mat_use[, "2"] + 1e-6) / (expr_mat_use[, "0"] + 1e-6))
names(log2fc_vector) <- rownames(expr_mat_use)

# 3.2 Main Evaluation Loop: DNB vs functional targets
res_list     <- list()
result_table <- data.frame()

for (g in DNB_vec) {
  result <- get_top_correlated(g, seurat_obj, cluster_id = "0", top_n = 50, allowed_genes = deg_genes)
  if (is.null(result) || length(result$genes) < 2) next
  
  top_genes <- result$genes[-1]
  pos_res   <- calculate_mean_and_pvalue(top_genes[result$signs[top_genes] == "Positive"], log2fc_vector, universe_genes)
  neg_res   <- calculate_mean_and_pvalue(top_genes[result$signs[top_genes] == "Negative"], log2fc_vector, universe_genes)
  
  res_list[[g]] <- data.frame(gene = g, positive_mean_log2FC = pos_res$observed_mean, positive_p_value = pos_res$p_value,
                              negative_mean_log2FC = neg_res$observed_mean, negative_p_value = neg_res$p_value)
  
  result_table <- rbind(result_table, data.frame(DNB = g, TARGET = top_genes, COR = result$signs[-1], 
                                                 Cor_values = result$correlations[-1], Cor_log2fc = log2fc_vector[top_genes]))
}

# 3.3 2D Quadrant Visualization
res_df <- bind_rows(res_list) %>%
          mutate(color_factor = case_when(positive_p_value < 0.2 & negative_p_value < 0.2 ~ "red",
                                          negative_p_value < 0.2                          ~ "green",
                                          positive_p_value < 0.2                          ~ "blue", 
                                          TRUE                                            ~ "grey"))

ggplot(res_df, aes(x = negative_mean_log2FC, y = positive_mean_log2FC)) +
  geom_point(aes(color = color_factor, size = -log10(pmin(positive_p_value, negative_p_value) + 1e-6)), alpha = 0.6) +
  geom_text_repel(aes(label = gene), size = 3.6, max.overlaps = 10) +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("red" = "red", "green" = "green", "blue" = "blue", "grey" = "grey"), name = "Significance") +
  labs(title = "T-Cell DNB Functional Decoupling: Prolif vs Exhaustion", x = "Prolif Mean Log2FC", y = "Exhaustion Mean Log2FC") +
  theme_minimal()

# --------------------------------------------------------------------------------------------------
# #### [SECTION 4] Multi-layer Sankey Topology & Expression Heatmaps ####
# --------------------------------------------------------------------------------------------------
library(networkD3)
library(ClusterGVis)

# 4.1 Mfuzz Soft Clustering for Temporal Trends
ck <- clusterData(exp = expr_mat_use[unique(result_table$TARGET), ], cluster.method = "mfuzz", cluster.num = 5)
write.csv(ck$long.res, "TableS_T_Coexp_Clusters.csv")

# 4.2 Multi-layered Sankey Render Loop
for (dnb in unique_dnbs) {
  dnb_data       <- result_table[result_table$DNB == dnb, ]
  dnb_gene_order <- dnb_data %>% count(TYPE, TARGET) %>% arrange(TYPE, desc(n)) %>% pull(TARGET)
  
  nodes <- data.frame(name = c(dnb, paste0("Cluster_", sort(unique(dnb_data$TYPE))), dnb_gene_order))
  
  layer1 <- dnb_data %>% group_by(DNB, TYPE) %>% summarise(value = n(), .groups = 'drop') %>%
            mutate(source = match(DNB, nodes$name)-1, target = match(paste0("Cluster_", TYPE), nodes$name)-1, group = as.factor(TYPE))
  layer2 <- dnb_data %>% mutate(source = match(paste0("Cluster_", TYPE), nodes$name)-1, target = match(TARGET, nodes$name)-1, value = 1, group = as.factor(COR))
  
  sankey <- sankeyNetwork(Links = bind_rows(layer1, layer2), Nodes = nodes, Source = "source", Target = "target", 
                          Value = "value", NodeID = "name", LinkGroup = "group", fontSize = 14, colourScale = JS(dnb_color_scale_js))
  
  # 4.3 Synchronized Cluster Heatmap
  pheatmap(expr_mat_use[dnb_gene_order, ], cluster_rows = FALSE, scale = "row",
           color = colorRampPalette(c("#2E75B6", "white", "#C00000"))(100),
           annotation_row = data.frame(Correlation = dnb_data$COR, Cluster = as.factor(dnb_data$TYPE)) %>% set_rownames(dnb_gene_order))
}
# ==================================================================================================
# MODULE: Multi-Lineage Trajectory Dynamics & DNB Topology (Neutrophils & T-Cells)
# 
# DESCRIPTION: 
#   Integrates Monocle 3 (spatial-consistent) and Monocle 2 (DDRTree) trajectory 
#   inference with BioTIP-based tipping point detection for Neutrophil and T-cell 
#   compartments in the Ulcerative Colitis (UC) MERFISH dataset.
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# #### [SECTION 1] NEUTROPHIL (GRAN) COMPARTMENT REFINEMENT ####
# --------------------------------------------------------------------------------------------------
library(Seurat); library(monocle3); library(ggsci); library(dplyr)

# 1.1 Subset and Initialize Granulocytes
Gran <- subset(sce, subset = annotation.2 == "Neutrophil") %>%
        subset(subset = day %!in% c(7, 11, 18))

# Standard Preprocessing Pipeline
Gran <- NormalizeData(Gran) %>% 
        FindVariableFeatures() %>% 
        ScaleData(features = rownames(.)) %>%
        RunPCA() %>% 
        FindNeighbors(dims = 1:30, reduction = "pca")

# Multi-resolution Clustering Evaluation
for (i in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 2)) {
  Gran <- FindClusters(Gran, resolution = i)
}
Gran <- RunUMAP(Gran, dims = 1:20)
clustree(Gran, prefix = "RNA_snn_res.")

# 1.2 Targeted Trajectory Analysis (Clusters 0, 1, 2, 3, 4)
Idents(Gran) <- "RNA_snn_res.0.5"
temp_sce     <- subset(Gran, RNA_snn_res.0.5 %in% c(0, 1, 2, 3, 4))

# 1.3 Monocle 3: Spatial-Consistent Trajectory Inference
expression_matrix <- temp_sce[["RNA"]]$counts
cell_metadata     <- data.frame(temp_sce@meta.data)
gene_annotation   <- data.frame(gene_short_name = row.names(expression_matrix), 
                                 row.names = row.names(expression_matrix))

cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50, norm_method = "log") %>%
       reduce_dimension(reduction_method = "UMAP") %>%
       cluster_cells(cluster_method = "louvain")

# Synchronize Seurat UMAP embeddings into Monocle 3
cds@int_colData$reducedDims$UMAP <- Embeddings(temp_sce, reduction = "umap")[rownames(cds@int_colData$reducedDims$UMAP), ]
cds <- learn_graph(cds)

# Root Selection: Cluster 3 defined as the developmental apex
cds <- order_cells(cds, root_pr_nodes = myselect(cds, 'RNA_snn_res.0.5', "3"))

# Visualize Neutrophil Differentiation
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE) + 
  plot_cells(cds, color_cells_by = "RNA_snn_res.0.5", label_groups_by_cluster = FALSE, 
             label_leaves = FALSE, label_branch_points = FALSE, group_label_size = 8) + 
  scale_color_d3(palette = "category20") + theme_minimal(base_size = 14)

# --------------------------------------------------------------------------------------------------
# #### [SECTION 2] T-CELL COMPARTMENT: DEEP CHARACTERIZATION ####
# --------------------------------------------------------------------------------------------------
Tcell <- readRDS("../UC_MERFISH/T.rds")
Idents(Tcell) <- "RNA_snn_res.0.3"

# 2.1 Quality Control and Feature Mapping
# Exclude artifactual/sparse Cluster 1
filter_Tcell <- subset(Tcell, RNA_snn_res.0.3 %!in% c(1)) %>% RunUMAP(dims = 1:20)

# Validate Functional Marker Domains
VlnPlot(filter_Tcell, features = c("Foxp3", "Cd4", "Ccr2", "Areg"), group.by = "RNA_snn_res.0.3", ncol = 4)       # Treg Domain
VlnPlot(filter_Tcell, features = c("Nkg7", "Klra7", "Cd160", "Klrd1"), group.by = "RNA_snn_res.0.3", ncol = 4)   # NK-like Domain
VlnPlot(filter_Tcell, features = c("mt-Co1", "mt-Cytb", "Igha", "Igkc"), group.by = "RNA_snn_res.0.3", ncol = 4) # QC / Plasma artifact check
VlnPlot(filter_Tcell, features = c("Ccna2", "Cdk1", "Pcna", "Birc5"), group.by = "RNA_snn_res.0.3", ncol = 4)     # Proliferation Domain

# 2.2 Monocle 2: DDRTree Trajectory Reconstruction (T-Cell Focus)
# Sub-lineage selection: Clusters 6, 3, 5, 0, 2
temp_sce <- subset(filter_Tcell, RNA_snn_res.0.3 %in% c(6, 3, 5, 0, 2))

sparse_data <- as(as.matrix(temp_sce[["RNA"]]$counts), "sparseMatrix")
mdata       <- new("AnnotatedDataFrame", data = temp_sce@meta.data)
fData       <- data.frame(gene_short_name = row.names(sparse_data), row.names = row.names(sparse_data))
fd          <- new("AnnotatedDataFrame", data = fData)

monocle_cds <- newCellDataSet(cellData = sparse_data, phenoData = mdata, featureData = fd, 
                              lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

# Pre-processing and Manifold Learning
monocle_cds <- estimateSizeFactors(monocle_cds) %>% estimateDispersions() %>% detectGenes(min_expr = 3)
ordering_genes <- row.names(subset(differentialGeneTest(monocle_cds[VariableFeatures(temp_sce), ], 
                                                        fullModelFormulaStr = "~RNA_snn_res.0.3"), 
                                   qval < 0.01))
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes) %>%
               reduceDimension(max_components = 2, reduction_method = "DDRTree") %>%
               orderCells(root_state = "3")

# Visualize T-Cell Developmental Trajectory
p_traj <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE) + 
          facet_wrap("~day", nrow = 1)
p_id   <- plot_cell_trajectory(monocle_cds, color_by = "RNA_snn_res.0.3", size = 1, show_backbone = TRUE) + 
          scale_color_d3(palette = "category20") + theme_minimal(base_size = 14)

# --------------------------------------------------------------------------------------------------
# #### [SECTION 3] T-CELL DNB CRITICALITY ASSESSMENT (BioTIP) ####
# --------------------------------------------------------------------------------------------------
library(BioTIP)

# 3.1 Data Partitioning for BioTIP
subsce_trans <- as.SingleCellExperiment(temp_sce)
samplesL     <- split(rownames(colData(subsce_trans)), f = colData(subsce_trans)$RNA_snn_res.0.3)
samplesL     <- samplesL[lengths(samplesL) > 0]

# Feature Pre-selection
hvg    <- getTopHVGs(modelGeneVarByPoisson(subsce_trans), n = 4000) %>% intersect(rownames(subsce_trans))
logmat <- as.matrix(logcounts(subsce_trans[hvg, ]))

# 3.2 Tipping Point Identification (MCI Scores)
topMCI   <- getTopMCI(membersL[["members"]], membersL[["MCI"]], min = 60, n = 4)
maxMCIms <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], min = 60, n = 4)
maxMCI   <- getMaxStats(membersL[['MCI']], maxMCIms[['idx']])
CTS.Lib  <- getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])

# Archive and Export DNB Sets
df_dnb <- data.frame(matrix(NA, nrow = max(lengths(CTS.Lib)), ncol = length(CTS.Lib)))
for (i in 1:length(CTS.Lib)) { df_dnb[, i] <- c(CTS.Lib[[i]], rep(NA, nrow(df_dnb) - length(CTS.Lib[[i]]))) }
colnames(df_dnb) <- names(CTS.Lib)
write.csv(df_dnb, file = "TCell_DNB_Targets.csv", row.names = TRUE)

# 3.3 Significance Simulation & Index Calculation (C=500 passes)
M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
simuMCI <- lapply(CTS.Lib, function(cts) simulationMCI(length(cts), samplesL, logmat, B = 500, fun = "BioTIP", M = M))

BioTIP_scores <- SimResults_g <- list()
for (i in 1:length(CTS.Lib)) {
  BioTIP_scores[[i]] <- getIc(logmat, samplesL, CTS.Lib[[i]], fun = "BioTIP", shrink = TRUE)
  SimResults_g[[i]]  <- simulation_Ic(length(CTS.Lib[[i]]), samplesL, logmat, B = 500, fun = "BioTIP", shrink = TRUE)
}
names(BioTIP_scores) <- names(SimResults_g) <- names(CTS.Lib)

# 3.4 Trajectory Scoring Rendering (Axis 6-0-2)
x <- factor(c(6, 0, 2), levels = c(6, 0, 2))
y <- c(0, 4.1663762, 0.1273221)
plot(x, y, type = "b", col = "black", lwd = 2, xlab = "T-Cell Cluster", ylab = "DNB Score")
text(x, y, labels = c(0, 84, 71), pos = 3, col = "red", font = 2)

# --------------------------------------------------------------------------------------------------
# #### [SECTION 4] BIVARIATE CO-EXPRESSION & TOPOLOGY (networkD3) ####
# --------------------------------------------------------------------------------------------------
library(networkD3); library(htmlwidgets); library(ggrepel)

# 4.1 Log2FC Mapping for 2D Evaluation
expr_mat_use  <- get_cluster_expr_matrix(rownames(seurat_obj), seurat_obj, c("6", "0", "2"))
log2fc_vector <- log2((expr_mat_use[, "2"] + 1e-6) / (expr_mat_use[, "0"] + 1e-6))
names(log2fc_vector) <- rownames(expr_mat_use)

# 4.2 Main DNB Loop: Quadrant Assignment & Permutation Test
for (g in DNB_vec) {
  result <- get_top_correlated(g, seurat_obj, cluster_id = "0", top_n = 50, allowed_genes = deg_genes)
  if (is.null(result) || length(result$genes) < 2) next
  
  # Execute direction-specific permutation tests
  pos_res <- calculate_mean_and_pvalue(result$genes[result$signs == "Positive"], log2fc_vector, universe_genes)
  neg_res <- calculate_mean_and_pvalue(result$genes[result$signs == "Negative"], log2fc_vector, universe_genes)
  
  # Archive structured results for current target gene
  # ... (Archival logic consistent with Figure 6 Macro script)
}

# 4.3 Multi-layered Sankey Render: Prolif vs Exhaustion
for (dnb in c("Cd6", "Foxp3", "Pkp3", "Ccr5", "Tnfrsf1b", "Ssh2")[4]) {
  dnb_data <- result_table[result_table$DNB == dnb, ]
  nodes    <- data.frame(name = c(dnb, paste0("Cluster_", sort(unique(dnb_data$TYPE))), dnb_data$TARGET))
  
  layer1 <- dnb_data %>% group_by(DNB, TYPE) %>% summarise(value = n(), .groups = 'drop') %>%
            mutate(source = match(DNB, nodes$name)-1, target = match(paste0("Cluster_", TYPE), nodes$name)-1, group = as.factor(TYPE))
  layer2 <- dnb_data %>% mutate(source = match(paste0("Cluster_", TYPE), nodes$name)-1, target = match(TARGET, nodes$name)-1, 
                                value = 1, group = as.factor(COR))
  
  sankey <- sankeyNetwork(Links = bind_rows(layer1, layer2), Nodes = nodes, Source = "source", Target = "target", 
                          Value = "value", NodeID = "name", LinkGroup = "group", fontSize = 14)
  saveWidget(sankey, paste0("Figure6_Tcell_Sankey_", dnb, ".html"))
}

# 4.4 Synchronized Expression Heatmap
pheatmap(heatmap_data, cluster_rows = FALSE, scale = "row",
         color = colorRampPalette(c("#2E75B6", "white", "#C00000"))(100),
         annotation_row = data.frame(Correlation = dnb_data$COR, Cluster = as.factor(dnb_data$TYPE)) %>% set_rownames(make.unique(dnb_data$TARGET, sep = "-")))

# --------------------------------------------------------------------------------------------------
# #### [SECTION 5] QUANTITATIVE PCC & SD MATRICES EXPORT ####
# --------------------------------------------------------------------------------------------------
# Construction of Node Attributes for Cytoscape/Network mapping
pcc_results <- list(); sd_results <- list()
for (cluster in c(6, 0, 2)) {
  cluster_exprs <- exprs_matrix_subset[, which(cluster_info == cluster), drop = FALSE]
  pcc_results[[as.character(cluster)]] <- as.data.frame(as.table(cor(t(cluster_exprs), method = "pearson")))
  sd_results[[as.character(cluster)]]  <- data.frame(gene = rownames(cluster_exprs), sd = apply(cluster_exprs, 1, sd))
}

write.csv(do.call(rbind, pcc_results), "./Figure6/T/pcc_sd_results/TCell_PCC_Master.csv")
cat("Master PCC and SD matrices for T-Cell clusters successfully exported!\n")
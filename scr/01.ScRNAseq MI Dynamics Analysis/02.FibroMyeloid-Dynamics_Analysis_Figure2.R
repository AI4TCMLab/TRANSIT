# ==============================================================================
# Script: Figure2_Global_Preprocessing_and_Doublets.R
# Purpose: Global data preprocessing, doublet detection via DoubletFinder, 
#          and global cell type proportion visualization (Figure 2.A-E candidate).
# Dataset: Global Seurat Object (seurat_obj.rds)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment Setup & Dependencies
# ------------------------------------------------------------------------------
setwd("/home/tsh/E-MTAB-7895/")

# Active Dependencies (Strictly required for this script)
library(Seurat)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(ggalluvial)
library(scales) # Required for percent_format()

# Historically loaded dependencies (Commented out to prevent namespace masking)
# pkgs_legacy <- c("batchelor","cluster","corrplot","dynamicTreeCut","edgeR",
#                  "gplots","gridExtra","igraph","limma","org.Mm.eg.db",
#                  "pheatmap","psych","scater","scran","SingleCellExperiment",
#                  "stringr","BioTIP","progress","tidyverse","DropletUtils",
#                  "patchwork","slingshot","RColorBrewer","clustree")

# ------------------------------------------------------------------------------
# 2. Data Loading & Global Subsetting
# ------------------------------------------------------------------------------
obj <- readRDS("./seurat_obj.rds")
Idents(obj) <- "orig.ident"

# Ensure pre-annotated labels exist
stopifnot("anno_human" %in% colnames(obj@meta.data))

# Merge Monocytes I & II into a unified 'Monocytes' consensus group
obj$sub_anno2 <- as.character(obj$anno_human)
obj$sub_anno2[obj$sub_anno2 %in% c("Monocytes I", "Monocytes II")] <- "Monocytes"
obj$sub_anno2 <- factor(obj$sub_anno2)

# Define target lineages for downstream analysis
keep_types <- c(
  "Monocytes", "Macrophages",
  "Fibroblast I", "Fibroblast II", "Fibroblast III", "Myofibroblasts",
  "Endothelium", "Granulocytes", "B cells", "NK cell", "DCs"
)

# Subset to specific temporal stages and target cell identities
obj_sub <- subset(obj, subset = orig.ident %in% c("0day", "1day", "3day", "5day", "7day"))
obj_sub <- subset(obj_sub, subset = sub_anno2 %in% keep_types)
obj_sub$sub_anno2 <- droplevels(obj_sub$sub_anno2)

# ------------------------------------------------------------------------------
# 3. Global Dimensionality Reduction & Clustering
# ------------------------------------------------------------------------------
DefaultAssay(obj_sub) <- "RNA"
obj_sub <- NormalizeData(obj_sub) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunUMAP(dims = 1:30, reduction.name = "umap.sub", verbose = FALSE) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.6, verbose = FALSE)

# ==============================================================================
# 4. Doublet Detection (DoubletFinder applied per sample)
# ==============================================================================
# Split object by sample to account for batch-specific doublet rates
list_sce <- SplitObject(obj_sub, split.by = "orig.ident")

processed_list <- lapply(list_sce, function(seu_temp) {
  s_name <- unique(seu_temp$orig.ident)
  n_cells <- ncol(seu_temp)
  
  # Skip samples with insufficient cells for reliable doublet modeling
  if (n_cells < 50) { 
    message("Skipping sample ", s_name, ": only ", n_cells, " cells")
    return(NULL) 
  }
  
  # Re-process individual sample
  seu_temp <- NormalizeData(seu_temp, verbose = FALSE) %>% 
    FindVariableFeatures(verbose = FALSE) %>% 
    ScaleData(verbose = FALSE)
  
  current_pcs <- min(15, n_cells - 1)
  seu_temp <- RunPCA(seu_temp, npcs = current_pcs, verbose = FALSE)
  
  # Hyperparameter tuning (pK identification)
  sweep.res <- paramSweep(seu_temp, PCs = 1:current_pcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  opt_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # Calculate homotypic doublet proportion to adjust expected doublet rate (Assuming 5% generic rate)
  homo_prop <- if ("sub_anno2" %in% colnames(seu_temp@meta.data)) modelHomotypic(seu_temp$sub_anno2) else 0
  nExp <- round(0.050 * n_cells)
  nExp_adj <- round(nExp * (1 - homo_prop))
  
  # Run DoubletFinder
  seu_temp <- doubletFinder(seu_temp, PCs = 1:current_pcs, pN = 0.25, pK = opt_pK, nExp = nExp_adj, sct = FALSE)
  
  # Standardize metadata column names dynamically
  meta <- seu_temp@meta.data
  meta$pANN_final <- meta[[grep("^pANN_", colnames(meta), value = TRUE)]]
  meta$DF_final <- meta[[grep("^DF.classifications_", colnames(meta), value = TRUE)]]
  seu_temp@meta.data <- meta
  
  return(seu_temp)
})

# Filter out skipped samples and merge doublet metadata back to the global object
processed_list <- Filter(Negate(is.null), processed_list)
all_doublet_meta <- bind_rows(lapply(processed_list, function(x) x@meta.data))
obj_sub <- AddMetaData(obj_sub, all_doublet_meta[, c("pANN_final", "DF_final")])

# ==============================================================================
# #### Figure 2.D & 2.E: Doublet Statistics Visualization ####
# ==============================================================================
plot_df <- obj_sub@meta.data %>% dplyr::select(sub_anno2, pANN = pANN_final, DF = DF_final)

# Summarize pANN (Proportion of Artificial Nearest Neighbors)
summary_d <- plot_df %>%
  group_by(sub_anno2) %>%
  summarise(
    mean_pANN = mean(as.numeric(pANN), na.rm = TRUE),
    sd_pANN   = sd(as.numeric(pANN), na.rm = TRUE), 
    n = n()
  ) %>%
  mutate(se_pANN = sd_pANN / sqrt(n))

# Summarize Doublet/Singlet Frequencies
summary_e <- plot_df %>%
  filter(!is.na(DF)) %>%
  group_by(sub_anno2, DF) %>% 
  tally() %>% 
  mutate(freq = n / sum(n))

# -- Figure 2.D: Average pANN Bar Plot --
p1 <- ggplot(summary_d, aes(x = sub_anno2, y = mean_pANN)) +
  geom_bar(stat = "identity", fill = "grey40", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_pANN - se_pANN, ymax = mean_pANN + se_pANN), width = 0.2) +
  theme_classic() +
  labs(y = "Average pANN", x = "", title = "Average pANN by Cell Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

# -- Figure 2.E: Doublet/Singlet Frequency Stacked Bar --
summary_e$DF <- factor(summary_e$DF, levels = c("Doublet", "Singlet"))
p2 <- ggplot(summary_e, aes(x = sub_anno2, y = freq, fill = DF)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = c("Doublet" = "#E41A1C", "Singlet" = "#D9D9D9")) +
  theme_classic() +
  labs(y = "Frequency (%)", x = "", fill = "Classification") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold")
  )

# Note: Requires patchwork to render side-by-side (library added at top)
print(p1 + p2)

# ==============================================================================
# #### Figure 2.C: Alluvial Plot of Major Cell Lineages ####
# ==============================================================================
# Classify detailed annotations into major functional lineages
obj_sub$major_type <- dplyr::case_when(
  obj_sub$sub_anno2 %in% c("Monocytes", "Macrophages", "DCs", "Granulocytes") ~ "Myeloid",
  obj_sub$sub_anno2 %in% c("Fibroblast I", "Fibroblast II", "Fibroblast III", "Myofibroblasts") ~ "Fibroblast",
  obj_sub$sub_anno2 %in% c("Endothelium", "Lymphatic endothelium") ~ "Endothelium",
  obj_sub$sub_anno2 %in% c("B cells", "NK cell", "T cells", "T/NK", "Lymphocytes") ~ "Lymphoid",
  TRUE ~ "Other"
) |> factor(levels = c("Myeloid", "Fibroblast", "Endothelium", "Lymphoid", "Other"))

meta_df2 <- obj_sub@meta.data %>%
  dplyr::select(orig.ident, major_type) %>%
  dplyr::group_by(orig.ident, major_type) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::mutate(major_type = factor(major_type, levels = levels(obj_sub$major_type))) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::arrange(orig.ident, major_type, .by_group = TRUE) %>%
  dplyr::mutate(prop = count / sum(count), cumprop = cumsum(prop)) %>% 
  ungroup()

# FIXME [Environment Dependency]: `pal_sub` is called here but is NEVER defined in this script.
# You need to define `pal_sub <- c("Myeloid" = "#...", "Fibroblast" = "#...", ...)` before running this plot.
p_alluv <- ggplot(meta_df2, aes(x = orig.ident, y = prop, fill = major_type, stratum = major_type, alluvium = major_type)) +
  geom_col(position = "stack", width = 0.8) +
  geom_alluvium(aes(stratum = major_type), width = 0.5, alpha = 0.5, color = "white", linewidth = 1, curve_type = "sigmoid") +
  geom_stratum(width = 0.8, alpha = 0.75, color = "white") +
  scale_fill_manual(values = pal_sub, drop = FALSE) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Percentage of cells", title = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  )

print(p_alluv)

# ==============================================================================
# Script: Figure2_Part2_Myeloid_Heterogeneity_and_Trajectory.R
# Purpose: Downstream visualization and pseudotime trajectory inference for 
#          myeloid/inflammation subpopulations (Figure 2.D, E, F, G).
# Dataset: E-MTAB-7895 (Myeloid subset across timepoints: 0d to 7d)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Dependencies
# ------------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(monocle3)
library(patchwork)
# library(ggrastr) # Loaded dynamically via :: in the code

# ------------------------------------------------------------------------------
# 2. Data Loading & Annotations (Figure 2.D Setup)
# ------------------------------------------------------------------------------
# Load the pre-processed Seurat object containing myeloid/inflammation cells
temp_sce <- readRDS("/home/tsh/E-MTAB-7895/inflammation_temp_sce.rds")
temp_sce <- subset(temp_sce, subset = orig.ident %in% c("0day", "1day", "3day", "5day", "7day"))

# Define standardized biological order and color palette for myeloid subtypes
sub_order <- c("Monocytes", "Pro_inflammation", "Proliferation", "Reparative", "Antigen-presenting", "Collagens")

cell_colors <- c(
  "Monocytes"          = "#CC79A7",
  "Pro_inflammation"   = "#D55E00",
  "Proliferation"      = "#E69F00",
  "Reparative"         = "#009E73",
  "Antigen-presenting" = "#0072B2",
  "Collagens"          = "#56B4E9"
)

temp_sce$inflammation_human <- factor(temp_sce$inflammation_human, levels = sub_order)

# Extract UMAP coordinates and merge with metadata for ggplot2 rendering
emb <- as.data.frame(Embeddings(temp_sce, "umap"))
colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
meta <- temp_sce@meta.data %>% mutate(cell = rownames(.))
emb$cell <- rownames(emb)
df <- left_join(emb, meta, by = "cell")

# ------------------------------------------------------------------------------
# 3. Custom Plotting Utility Functions
# Note: These functions are duplicated from the Fibroblast pipeline.
# ------------------------------------------------------------------------------
theme_umap_cns <- function(base_size = 9, base_family = "Helvetica") {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.key.height = unit(4.0, "mm"),
      legend.key.width  = unit(4.0, "mm"),
      plot.margin = margin(10, 10, 30, 30)
    )
}

add_umap_axes_arrows_outside <- function(p, df, dx_frac = 0.12, dy_frac = 0.12, x_pad_frac = 0.10, y_pad_frac = 0.10) {
  xr <- range(df$UMAP_1); yr <- range(df$UMAP_2)
  dx <- dx_frac * diff(xr); dy <- dy_frac * diff(yr)
  x0 <- xr[1] - x_pad_frac * diff(xr); y0 <- yr[1] - y_pad_frac * diff(yr)
  
  p +
    coord_equal(clip = "off") +
    annotate("segment", x = x0, y = y0, xend = x0 + dx, yend = y0, linewidth = 0.35, arrow = arrow(length = unit(7, "pt"), type = "closed")) +
    annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + dy, linewidth = 0.35, arrow = arrow(length = unit(7, "pt"), type = "closed")) +
    annotate("text", x = x0 + dx, y = y0, label = "UMAP1", hjust = -0.10, vjust = 0.50, size = 2.7) +
    annotate("text", x = x0, y = y0 + dy, label = "UMAP2", hjust = 0.50, vjust = -0.60, size = 2.7)
}

# ==============================================================================
# #### Figure 2.D: Global UMAP of Myeloid Subtypes ####
# ==============================================================================
p_umap <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(color = "grey65", size = 0.10, alpha = 0.08, shape = 16, stroke = 0) +
  geom_point(aes(color = inflammation_human), size = 0.10, alpha = 0.55, shape = 16, stroke = 0) +
  scale_color_manual(values = cell_colors, breaks = sub_order, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
  labs(title = "Myeloid subtypes UMAP") +
  theme_umap_cns()

p_umap <- add_umap_axes_arrows_outside(p_umap, df)
print(p_umap)

# ==============================================================================
# #### Figure 2.E: Faceted UMAP by Timepoint ####
# ==============================================================================
df$timepoint <- factor(df$orig.ident, levels = sort(unique(df$orig.ident)))

label_df <- df %>% 
  group_by(timepoint, inflammation_human) %>% 
  summarise(x = median(UMAP_1), y = median(UMAP_2), n = dplyr::n(), .groups = "drop")

p_split <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  ggrastr::rasterise(geom_point(color = "grey65", size = 0.10, alpha = 0.08, shape = 16, stroke = 0), dpi = 600) +
  ggrastr::rasterise(geom_point(aes(color = inflammation_human), size = 0.08, alpha = 0.80, shape = 16, stroke = 0), dpi = 600) +
  scale_color_manual(values = cell_colors, drop = FALSE) +
  coord_equal() +
  facet_wrap(~ timepoint, ncol = 7) +
  theme_umap_cns() +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
print(p_split)

# ==============================================================================
# #### Figure 2.F: Alluvial Plot of Subpopulation Dynamics ####
# Note: In Part 1, Figure 2.F was a faceted UMAP. Proceeding with Alluvial plot.
# ==============================================================================
meta_stack <- temp_sce@meta.data %>%
  dplyr::select(orig.ident, cell_type = inflammation_human) %>%
  mutate(cell_type = factor(cell_type, levels = sub_order)) %>%
  group_by(orig.ident, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  arrange(orig.ident, cell_type, .by_group = TRUE) %>%
  mutate(cumsum = cumsum(count), prop = count / sum(count), cumprop = cumsum(prop)) %>% 
  ungroup()

p_alluv <- ggplot(meta_stack, aes(x = orig.ident, y = prop, fill = cell_type, stratum = cell_type, alluvium = cell_type)) +
  geom_col(position = "stack", width = 0.8) +
  geom_alluvium(aes(stratum = cell_type), width = 0.5, alpha = 0.2, color = "white", linewidth = 1, curve_type = "sigmoid") +
  geom_stratum(width = 0.8, alpha = 0.75, color = "white") +
  scale_fill_manual(values = cell_colors, drop = FALSE) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Percentage of cells", title = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  )
print(p_alluv)

# ==============================================================================
# #### Figure 2.G: Monocle3 Trajectory Inference (Myeloid Compartment) ####
# ==============================================================================
DefaultAssay(temp_sce) <- "RNA"
expression_matrix <- temp_sce[["RNA"]]$counts
cell_metadata <- data.frame(temp_sce@meta.data)
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))

# Construct Monocle3 CellDataSet (CDS)
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# Preprocess and dimensionality reduction
cds <- preprocess_cds(cds, num_dim = 50, norm_method = c("log"))
cds <- reduce_dimension(cds, reduction_method = "UMAP", cores = 5)
plot_cells(cds, color_cells_by = "seurat_clusters")

# Clustering
cds <- cluster_cells(cds, cluster_method = "louvain", reduction_method = "UMAP", k = 25)

# Integrate native Seurat UMAP coordinates
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(temp_sce, reduction = "umap")
cis.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- cis.embed

# Trajectory graph learning
cds <- learn_graph(cds, use_partition = TRUE, close_loop = TRUE)

#' Helper function to specify root nodes for pseudotime calculation
myselect <- function(cds, select.classify, my_select) {
  cell_ids <- which(colData(cds)[, select.classify] == my_select)
  closest_vertex <- as.matrix(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))]
  return(root_pr_nodes)
}

# Order cells using Monocytes as the developmental root
cds <- order_cells(cds, root_pr_nodes = myselect(cds, select.classify = 'inflammation_human', my_select = "Monocytes"))
save(cds, file = "/home/tsh/E-MTAB-7895/大论文图和数据/Figure 2/Macro_monocle3.Rdata", compress = TRUE)

# ------------------------------------------------------------------------------
# Evaluate Gene Dynamics Along Pseudotime
# ------------------------------------------------------------------------------
selected_cells <- temp_sce[, temp_sce$seurat_clusters %in% c(13, 6, 12, 3, 11, 8)]

# FIXME [Environment Dependency]: `selected_cells$pseudotime` is called below, but was 
# not explicitly assigned from the `cds` object in this code block. Ensure `temp_sce$pseudotime` 
# is populated prior to this step, e.g., temp_sce$pseudotime <- pseudotime(cds)

# Define target gene modules (Note: Genes are MHC-II related, but title says "Proliferation")
target_genes <- c("H2-DMa", "H2-Aa", "H2-Ab1", "H2-DMb1", "Ighm")
existing_genes <- intersect(target_genes, rownames(selected_cells))
expr_mat <- as.matrix(GetAssayData(selected_cells, slot = "data")[existing_genes, , drop = FALSE])

cell_stats <- data.frame(
  cell = colnames(expr_mat), 
  pseudotime = selected_cells$pseudotime, 
  idx_mean = colMeans(expr_mat), 
  idx_min = apply(expr_mat, 2, min), 
  idx_max = apply(expr_mat, 2, max)
) %>% filter(!is.na(pseudotime))

plot_grid <- data.frame(pseudotime = seq(min(cell_stats$pseudotime), max(cell_stats$pseudotime), length.out = 500))

# Perform external Loess smoothing for index prediction limits
plot_grid$y_mean <- predict(loess(idx_mean ~ pseudotime, data = cell_stats, span = 0.75), plot_grid)
plot_grid$y_min  <- predict(loess(idx_min ~ pseudotime, data = cell_stats, span = 0.75), plot_grid)
plot_grid$y_max  <- predict(loess(idx_max ~ pseudotime, data = cell_stats, span = 0.75), plot_grid)

# Plot Gene Set Index Trend
p_trend <- ggplot(plot_grid, aes(x = pseudotime)) +
  geom_ribbon(aes(ymin = y_min, ymax = y_max), fill = "grey80", alpha = 0.4) +
  geom_line(aes(y = y_mean), color = "#D55E00", linewidth = 1.2) +
  theme_classic(base_size = 14) +
  labs(x = "Pseudotime", y = "Gene Set Index", title = "Proliferation Program") +
  theme(
    axis.line = element_line(color = "black"), 
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

# Combine Density and Trend Plots using patchwork
# FIXME [Environment Dependency]: `p_density` is not defined in this script. 
# It likely persists in your global environment from the Fibroblast analysis (Figure 2.K). 
# If run in a fresh R session, this will throw an "object 'p_density' not found" error.
final_plot <- p_density / p_trend + plot_layout(heights = c(1, 4))
print(final_plot)


# ==============================================================================
# Script: Figure2_Fibroblast_Heterogeneity_and_Trajectory.R
# Purpose: Downstream visualization and pseudotime trajectory inference for 
#          fibroblast subpopulations (Figure 2.H, I, F, J, K).
# Dataset: E-MTAB-7895 (Fibroblast subset across timepoints: 0d to 7d)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Dependencies
# ------------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(colorspace)
library(scales)
library(ggalluvial)
library(devtools)
library(monocle3)
library(leidenbase)
library(tidyr)
library(patchwork)
library(ggpubr)
library(ggsci)
# library(ggrastr) # Loaded dynamically via :: in the code

# ------------------------------------------------------------------------------
# 2. Data Loading & Preprocessing
# ------------------------------------------------------------------------------
# Load the pre-processed Seurat object containing fibroblast cells
temp_sce <- readRDS("/home/tsh/E-MTAB-7895/temp_sce.rds")

# Subset data to specific experimental timepoints
temp_sce <- subset(temp_sce, subset = orig.ident %in% c("0day", "1day", "3day", "5day", "7day"))

# Initial visualization for quality control
DimPlot(temp_sce, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), label = TRUE)

# ------------------------------------------------------------------------------
# 3. Define Cell Annotations & Palettes (Figure 2.H Setup)
# ------------------------------------------------------------------------------
# Map generic Seurat clusters to functional major cell types
temp_sce$major_type <- dplyr::case_when(
  temp_sce$seurat_clusters %in% c("1", "4", "14", "12", "11") ~ "IR",
  temp_sce$seurat_clusters %in% c("9", "5", "8", "10", "13") ~ "Proliferation",
  temp_sce$seurat_clusters %in% c("3", "6", "2", "7", "0") ~ "Myofib",
  TRUE ~ "Other"
) |> factor(levels = c("IR", "Proliferation", "Myofib"))

# Define standardized color palette for fibroblast subtypes across all plots
cell_colors <- c(
  "IR"            = "#CC79A7",  # Sky blue (Note: Hex is historically pink/purple, keeping as is)
  "Proliferation" = "#E69F00",  # Orange
  "Myofib"        = "#009E73"   # Bluish green
)

# Extract UMAP coordinates and merge with metadata for ggplot2 custom plotting
emb <- as.data.frame(Embeddings(temp_sce, "umap"))
colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
meta <- temp_sce@meta.data %>% mutate(cell = rownames(.))
emb$cell <- rownames(emb)
df <- left_join(emb, meta, by = "cell")

# ------------------------------------------------------------------------------
# 4. Custom Plotting Utility Functions
# ------------------------------------------------------------------------------
#' Custom ggplot theme for clean UMAP visualizations
theme_umap_cns <- function(base_size = 9, base_family = "Helvetica") {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.key.height = unit(4.0, "mm"),
      legend.key.width  = unit(4.0, "mm"),
      plot.margin = margin(10, 10, 30, 30)
    )
}

#' Add aesthetic UMAP directional arrows outside the plot area
add_umap_axes_arrows_outside <- function(p, df,
                                         dx_frac = 0.12, dy_frac = 0.12,
                                         x_pad_frac = 0.10, y_pad_frac = 0.10) {
  xr <- range(df$UMAP_1); yr <- range(df$UMAP_2)
  dx <- dx_frac * diff(xr)
  dy <- dy_frac * diff(yr)
  x0 <- xr[1] - x_pad_frac * diff(xr)
  y0 <- yr[1] - y_pad_frac * diff(yr)
  
  p +
    coord_equal(clip = "off") +
    annotate("segment", x = x0, y = y0, xend = x0 + dx, yend = y0,
             linewidth = 0.35, arrow = arrow(length = unit(7, "pt"), type = "closed")) +
    annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + dy,
             linewidth = 0.35, arrow = arrow(length = unit(7, "pt"), type = "closed")) +
    annotate("text", x = x0 + dx, y = y0, label = "UMAP1",
             hjust = -0.10, vjust = 0.50, size = 2.7) +
    annotate("text", x = x0, y = y0 + dy, label = "UMAP2",
             hjust = 0.50, vjust = -0.60, size = 2.7)
}

# ==============================================================================
# #### Figure 2.H: Global UMAP of Fibroblast Subtypes ####
# ==============================================================================
p_umap <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(color = "grey65", size = 0.50, alpha = 0.8, shape = 16, stroke = 0) +
  geom_point(aes(color = major_type),
             size = 0.50, alpha = 0.8, shape = 16, stroke = 0) +
  scale_color_manual(values = cell_colors,  drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
  labs(title = "Fibro subtypes UMAP") +
  theme_umap_cns()

p_umap <- add_umap_axes_arrows_outside(p_umap, df)
print(p_umap)

# ==============================================================================
# #### Figure 2.I: Faceted UMAP by Timepoint (Rasterized, Fine points) ####
# ==============================================================================
df$timepoint <- factor(df$orig.ident, levels = sort(unique(df$orig.ident)))

label_df <- df %>%
  group_by(timepoint, major_type) %>%
  summarise(
    x = median(UMAP_1),
    y = median(UMAP_2),
    n = dplyr::n(),
    .groups = "drop"
  )

p_split_I <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  ggrastr::rasterise(
    geom_point(color = "grey65", size = 0.10, alpha = 0.08, shape = 16, stroke = 0),
    dpi = 600
  ) +
  ggrastr::rasterise(
    geom_point(aes(color = major_type), size = 0.3, alpha = 0.80, shape = 16, stroke = 0),
    dpi = 600
  ) +
  scale_color_manual(values = cell_colors, drop = FALSE) +
  coord_equal() +
  facet_wrap(~ timepoint, ncol = 7) +
  theme_umap_cns() +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
print(p_split_I)

# ==============================================================================
# #### Figure 2.F: Faceted UMAP by Timepoint (Rasterized, Coarse points) ####
# Note: Conceptually identical to 2.I, but uses different size/alpha parameters.
# ==============================================================================
p_split_F <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  ggrastr::rasterise(
    geom_point(color = "grey65", size = 0.50, alpha = 0.8, shape = 16, stroke = 0),
    dpi = 600
  ) +
  ggrastr::rasterise(
    geom_point(aes(color = major_type), size = 0.5, alpha = 0.80, shape = 16, stroke = 0),
    dpi = 600
  ) +
  scale_color_manual(values = cell_colors,  drop = FALSE) +
  coord_equal() +
  facet_wrap(~ timepoint, ncol = 7) +
  theme_umap_cns() +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
print(p_split_F)

# ==============================================================================
# #### Figure 2.J: Alluvial Plot of Cell Type Proportions Across Time ####
# ==============================================================================
meta_stack <- temp_sce@meta.data %>%
  dplyr::select(orig.ident, cell_type = major_type) %>%
  mutate(cell_type = factor(cell_type)) %>%
  group_by(orig.ident, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  arrange(orig.ident, cell_type, .by_group = TRUE) %>%
  mutate(
    cumsum   = cumsum(count),
    prop     = count / sum(count),
    cumprop  = cumsum(prop)
  ) %>%
  ungroup()

p_alluv <- ggplot(
  meta_stack,
  aes(
    x = orig.ident,
    y = prop,
    fill = cell_type,
    stratum = cell_type,
    alluvium = cell_type
  )
) +
  geom_col(position = "stack", width = 0.8) +
  geom_alluvium(aes(stratum = cell_type), width = 0.5, alpha = 0.2, color = "white", linewidth = 1, curve_type = "sigmoid") +
  geom_stratum(width = 0.8, alpha = 0.75, color = "white") +
  scale_fill_manual(values = cell_colors, drop = FALSE) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Percentage of cells", title = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.title = element_blank()
  )
print(p_alluv)

# ==============================================================================
# #### Figure 2.K: Monocle3 Trajectory Inference & Gene Dynamics ####
# ==============================================================================

DefaultAssay(temp_sce) <- "RNA"

# 1. Construct Monocle3 CellDataSet (CDS)
expression_matrix <- temp_sce[["RNA"]]$counts
cell_metadata <- data.frame(temp_sce@meta.data)

gene_annotation <- data.frame(expression_matrix[, 1])
gene_annotation[, 1] <- row.names(gene_annotation)
colnames(gene_annotation) <- c("gene_short_name")

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# 2. Preprocessing & Dimensionality Reduction
cds <- preprocess_cds(cds, num_dim = 50, norm_method = c("log"))
cds <- reduce_dimension(cds, reduction_method = "UMAP", cores = 5)
plot_cells(cds, color_cells_by = "seurat_clusters")

# 3. Clustering & Trajectory Graph Learning
# Note: Resolution is set low to enforce specific clustering representation
cds <- cluster_cells(cds, cluster_method = "louvain", reduction_method = "UMAP", k = 23)
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = FALSE, 
           label_leaves = FALSE, label_branch_points = FALSE)

# 4. Integrate Seurat's Original UMAP Embeddings into Monocle CDS
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(temp_sce, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- learn_graph(
  cds,
  use_partition = TRUE,
  close_loop = TRUE
)

# 5. Define Root Nodes & Order Cells
#' Helper function to select the root node based on a specific cell annotation
myselect <- function(cds, select.classify, my_select) {
  cell_ids <- which(colData(cds)[, select.classify] == my_select)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))]
  return(root_pr_nodes)
}

# Order cells starting from the 'IR' major type / Cluster '1'
cds <- order_cells(cds, root_pr_nodes = myselect(cds, select.classify = 'major_type', my_select = "IR"))
cds <- order_cells(cds, root_pr_nodes = myselect(cds, select.classify = 'seurat_clusters', my_select = "1"))

# Save the intermediate Monocle3 object
save(cds, file = "/home/tsh/E-MTAB-7895/大论文图和数据/Figure 2/Fibro_monocle3.Rdata", compress = TRUE)

# Append pseudotime values back to the Seurat object
temp_sce$pseudotime <- pseudotime(cds, reduction_method = "UMAP")

# 6. Evaluate Gene Dynamics Along Pseudotime
# Note: `selected_cells` is reassigned immediately. Preserving history for context.
selected_cells <- temp_sce[, temp_sce$seurat_clusters %in% c(1, 4, 9, 5, 8)]
selected_cells <- temp_sce[, temp_sce$seurat_clusters %in% c(1, 4, 9, 6, 7)] # Active override

# Define target gene modules
inflammation_genes <- c("Mt2", "Cxcl5", "Ccl2", "Angptl4", "Cxcl3")
ECM_genes <- c("Comp", "Col11a1", "Sfrp2", "Eln", "Cthrc1") 
genes_use <- c("Comp", "Col11a1", "Sfrp2", "Eln", "Cthrc1") 

# Extract expression matrix and align pseudotime
expr_mat <- GetAssayData(selected_cells, slot = "data")[genes_use, , drop = FALSE]
pt <- selected_cells$pseudotime
pt <- pt[colnames(expr_mat)]   

# Reshape data for ggplot2 visualization
plot_df <- as.data.frame(t(as.matrix(expr_mat)))
plot_df$cell <- rownames(plot_df)
plot_df$pseudotime <- pt[plot_df$cell]

plot_long <- plot_df %>%
  pivot_longer(
    cols = all_of(genes_use),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  filter(!is.na(pseudotime)) %>%
  mutate(log_expr = log1p(expression))

# Plot Individual Gene Trajectories
p_genes_traj <- ggplot(plot_long, aes(x = pseudotime, y = log_expr, color = gene)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.75, linewidth = 1.4) +
  theme_classic(base_size = 14) +
  labs(x = "Pseudotime", y = "Log(expression + 1)", color = NULL) +
  theme(legend.position = "right", axis.line = element_line(color = "black"))
print(p_genes_traj)

# 7. Pseudotime Density & ECM Program Index Analysis
df_cds <- pData(cds) 
df_cds$Pseudotime <- pseudotime(cds, reduction_method = "UMAP")

p_density <- ggplot(df_cds, aes(Pseudotime, colour = major_type, fill = major_type)) +     
  geom_density(bw = 0.5, size = 0.8, alpha = 0.5) + 
  theme_classic2()

# Calculate Extracellular Matrix (ECM) Program Index Statistics
existing_genes <- intersect(ECM_genes, rownames(selected_cells))
expr_mat_ecm <- GetAssayData(selected_cells, slot = "data")[existing_genes, , drop = FALSE]
expr_mat_ecm <- as.matrix(expr_mat_ecm)

cell_stats <- data.frame(
  cell = colnames(expr_mat_ecm),
  pseudotime = selected_cells$pseudotime,
  idx_mean = colMeans(expr_mat_ecm),                     # Mean expression (Trend line)
  idx_min = apply(expr_mat_ecm, 2, min),                 # Lower bound
  idx_max = apply(expr_mat_ecm, 2, max)                  # Upper bound
) %>% filter(!is.na(pseudotime))

# Perform external Loess smoothing for index prediction limits
plot_grid <- data.frame(pseudotime = seq(min(cell_stats$pseudotime), max(cell_stats$pseudotime), length.out = 500))

fit_mean <- loess(idx_mean ~ pseudotime, data = cell_stats, span = 0.75)
fit_min  <- loess(idx_min ~ pseudotime, data = cell_stats, span = 0.75)
fit_max  <- loess(idx_max ~ pseudotime, data = cell_stats, span = 0.75)

plot_grid$y_mean <- predict(fit_mean, plot_grid)
plot_grid$y_min  <- predict(fit_min, plot_grid)
plot_grid$y_max  <- predict(fit_max, plot_grid)

# Plot ECM Gene Set Index Dynamics
p_trend <- ggplot(plot_grid, aes(x = pseudotime)) +
  geom_ribbon(aes(ymin = y_min, ymax = y_max), fill = "grey80", alpha = 0.4) +
  geom_line(aes(y = y_mean), color = "#D55E00", linewidth = 1.2) +
  theme_classic(base_size = 14) +
  labs(x = "Pseudotime", y = "Gene Set Index", title = "ECM Program") +
  theme(
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

# Combine Density and Trend Plots using patchwork
final_plot <- p_density / p_trend + plot_layout(heights = c(1, 4))
print(final_plot)

# 8. Visualization of Trajectory on UMAP
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE) + 
  plot_cells(cds,
             color_cells_by = "major_type",
             label_cell_groups = FALSE,
             label_leaves = FALSE,
             label_branch_points = FALSE,
             graph_label_size = 1) + 
  scale_color_manual(values = pal_jco("default", alpha = 0.6)(19))

# 9. Additional Marker Identifications
DimPlot(temp_sce, reduction = "umap", split.by = c("major_type"), label = TRUE)
markers2 <- FindAllMarkers(temp_sce, group.by = "major_type")

dotplot_markers <- c(
  "Mt2", "Cxcl5", "Ccl2", "Angptl4", "Cxcl3",
  "Timp1", "Spp1", "Acta2", "Tpm2", "Serpinb2",
  "Cdc20", "Cenpe", "Stmn1", "Cks1b", "Hmgb2",
  "Comp", "Col11a1", "Sfrp2", "Eln", "Cthrc1"
) 


#### Figure2.L ####
temp_sce <- readRDS("../E-MTAB-7895/Granulocytes_temp_sce.rds")
temp_sce <- subset(temp_sce, subset = orig.ident %in% c("0day","1day","3day","5day","7day" ))
DimPlot(temp_sce, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"),label = T)
factor(temp_sce@meta.data$seurat_clusters)

temp_sce$major_type <- dplyr::case_when(
  temp_sce$seurat_clusters %in% c("0","3","4") ~ "Pro_inflammation",
  temp_sce$seurat_clusters %in% c("1","2") ~ "Reparative",
  temp_sce$seurat_clusters %in% c("5","6") ~ "Other",
  TRUE ~ "Other"
) |> factor(levels = c("Pro_inflammation","Reparative","Other"))

cell_colors <- c(
  "Pro_inflammation"          = "#CC79A7",  # sky blue
  
  "Reparative"      = "#009E73",  # orange
  "Other"         = "grey"  # bluish green
  
)
# 取坐标
emb <- as.data.frame(Embeddings(temp_sce, "umap"))
colnames(emb)[1:2] <- c("UMAP_1","UMAP_2")
meta <- temp_sce@meta.data %>% mutate(cell = rownames(.))
emb$cell <- rownames(emb)
df <- left_join(emb, meta, by = "cell")

theme_umap_cns <- function(base_size = 9, base_family = "Helvetica") {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.key.height = unit(4.0, "mm"),
      legend.key.width  = unit(4.0, "mm"),
      plot.margin = margin(10, 10, 30, 30)
    )
}

add_umap_axes_arrows_outside <- function(p, df,
                                         dx_frac = 0.12, dy_frac = 0.12,
                                         x_pad_frac = 0.10, y_pad_frac = 0.10){
  xr <- range(df$UMAP_1); yr <- range(df$UMAP_2)
  dx <- dx_frac * diff(xr)
  dy <- dy_frac * diff(yr)
  x0 <- xr[1] - x_pad_frac * diff(xr)
  y0 <- yr[1] - y_pad_frac * diff(yr)
  
  p +
    coord_equal(clip = "off") +
    annotate("segment", x = x0, y = y0, xend = x0 + dx, yend = y0,
             linewidth = 0.35, arrow = arrow(length = unit(7, "pt"), type = "closed")) +
    annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + dy,
             linewidth = 0.35, arrow = arrow(length = unit(7, "pt"), type = "closed")) +
    annotate("text", x = x0 + dx, y = y0, label = "UMAP1",
             hjust = -0.10, vjust = 0.50, size = 2.7) +
    annotate("text", x = x0, y = y0 + dy, label = "UMAP2",
             hjust = 0.50, vjust = -0.60, size = 2.7)
}

p_umap <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(color = "grey65", size = 0.50, alpha = 0.8, shape = 16, stroke = 0) +
  geom_point(aes(color = major_type),
             size = 0.50, alpha = 0.8, shape = 16, stroke = 0)+
  scale_color_manual(values = cell_colors,  drop = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
  labs(title = "Gran subtypes UMAP") +
  theme_umap_cns()

p_umap <- add_umap_axes_arrows_outside(p_umap, df)
p_umap

#### Figure2.M ####
df$timepoint <- factor(df$orig.ident, levels = sort(unique(df$orig.ident)))  # 可选：固定顺序
library(dplyr)

label_df <- df %>%
  group_by(timepoint, major_type) %>%
  summarise(
    x = median(UMAP_1),
    y = median(UMAP_2),
    n = dplyr::n(),
    .groups = "drop"
  )
p_split <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  ggrastr::rasterise(
    geom_point(color = "grey65", size = 0.50, alpha = 0.8, shape = 16, stroke = 0),
    dpi = 600
  ) +
  ggrastr::rasterise(
    geom_point(aes(color = major_type), size = 0.5, alpha = 0.80, shape = 16, stroke = 0),
    dpi = 600
  ) +
  # geom_text(
  #   data = subset(label_df, n >= 100),
  #   aes(x = x, y = y, label = major_type),
  #   inherit.aes = FALSE,
  #   size = 2.4,
  #   fontface = "bold"
  # ) +
  scale_color_manual(values = cell_colors, breaks = sub_order, drop = FALSE) +
  coord_equal() +
  facet_wrap(~ timepoint, ncol = 7) +
  theme_umap_cns() +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
p_split
#### Figure2.N ####
library(dplyr)

meta_stack <- temp_sce@meta.data %>%
  dplyr::select(orig.ident, cell_type = major_type) %>%
  mutate(cell_type = factor(cell_type)) %>%
  group_by(orig.ident, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  arrange(orig.ident, cell_type, .by_group = TRUE) %>%
  mutate(
    cumsum   = cumsum(count),
    prop     = count / sum(count),
    cumprop  = cumsum(prop)
  ) %>%
  ungroup()
library(ggplot2)
library(scales)

p_alluv <- ggplot(
  meta_stack,
  aes(
    x = orig.ident,
    y = prop,
    fill = cell_type,
    stratum = cell_type,
    alluvium = cell_type
  )
) +
  geom_col(position = "stack", width = 0.8) +
  geom_alluvium(aes(stratum = cell_type), width = 0.5, alpha = 0.6, color = "white", linewidth = 1, curve_type = "sigmoid") +
  geom_stratum(width = 0.8, alpha = 0.75, color = "white") +
  scale_fill_manual(values = cell_colors, drop = FALSE) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Percentage of cells", title = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.title = element_blank()
  )

p_alluv
#### Figure2.O ####
library(devtools)
library(monocle3)
library(Seurat)
library(leidenbase)
DefaultAssay(temp_sce) <- "RNA"
#构建Monocle3所需数据
#expression_matrix：行是gene，列是cell的表达矩阵
expression_matrix = temp_sce[["RNA"]]$counts

#cell_metadata:对应Seurat对象的meta.data信息
cell_metadata = data.frame(temp_sce@meta.data)

#gene_metadata:基因名称信息，可用Gene symbol
gene_annotation = data.frame(expression_matrix[,1])
gene_annotation[,1] = row.names(gene_annotation)
colnames(gene_annotation)=c("gene_short_name")
##构建Monocle3 cds对象
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#数据预处理
cds <- preprocess_cds(cds, num_dim = 50,norm_method = c("log"))

## 降维，默认是"Umap"方式
cds <- reduce_dimension(cds,reduction_method="UMAP",cores=5)
plot_cells(cds,color_cells_by = "seurat_clusters")
## 聚类分群，分辨率调小，是为了让细胞是一群可以更好展示
#cds <- cluster_cells(cds,resolution = 0.0000001)
cds <- cluster_cells(cds,cluster_method = "louvain", reduction_method ="UMAP")
## 拟时序
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster=FALSE, label_leaves=FALSE,label_branch_points=FALSE)

##选择特定细胞作为起点
##这里我们假定以HEpiD细胞为起点
myselect <- function(cds,select.classify,my_select){
  cell_ids <- which(colData(cds)[,select.classify] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes}


cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'major_type',my_select = "IR"))
##使用Seurat的UMAP信息，这样可以与Seurat对象的细胞分布保持一致
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(temp_sce , reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

#root_group = colnames(cds)[clusters(cds) == 1]
#cds = order_cells(cds, root_cells = root_group)

##不同细胞类型拟时序数值
##拟时序值越高表示细胞分化程度越高，这里仅为演示，并非真实分化情况
library(ggsci)
plot_cells(cds,color_cells_by = "pseudotime",
           #group_cells_by = "new.labels",
           #color_cells_by = "orig.ident",
           show_trajectory_graph=F) + plot_cells(cds,
                                                 #color_cells_by = "seurat_clusters",
                                                 color_cells_by = "major_type",
                                                 #color_cells_by = "label",
                                                 #color_cells_by = "orig.ident",
                                                 label_cell_groups=FALSE,
                                                 label_leaves=FALSE,
                                                 label_branch_points=FALSE,
                                                 graph_label_size=1)+ scale_color_manual(values = pal_jco("default", alpha = 0.6)(19))


# ==============================================================================
# Script: Figure2_Orthogonal_Validation_and_Heatmap.R
# Purpose: Complex Heatmap visualization (Fig 2.P) and orthogonal trajectory 
#          validation using Monocle 2 / CytoTRACE2 (Fig 2.Q, 2.S).
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Dependencies
# ------------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(monocle)      # Monocle 2
library(igraph)
library(CytoTRACE2)
library(magrittr)     # For set_rownames

# ==============================================================================
# #### Figure 2.P: Pathway Gene Set Heatmap ####
# ==============================================================================
# Define target functional gene sets
genesets_list <- list(
  inflammation     = c("Serpina3n", "Angptl4", "Cxcl2", "Il1b", "S100a8", "Mif", "Nfkb1"),
  oxidative_stress = c("Mt2", "Timp1", "Mt1", "Hmox1", "Sod1", "Hmox1", "Atf4"),
  apoptosis        = c("Slc16a3", "Cdkn1a", "Lgals3", "Xbp1", "Ctsc", "Ier3"),
  proliferation    = c("Aprt", "Sdcbp", "Stmn1", "Ranbp1", "Ccnd1", "Nme1", "Mki67"),
  angiogenesis     = c("Postn", "Igf1", "Fstl1", "Cdh2", "Anxa2"),
  phagocytosis     = c("Ms4a7", "Trem2", "C1qb", "Msr1", "Fabp5"),
  fibrosis         = c("Col4a1", "Eln", "Ltbp1", "Tgfb3", "Ccn5", "Acta2", "Col1a1")
)

geneset_colors <- c(
  "inflammation"     = "#66C2A5",
  "oxidative_stress" = "#FC8D62",
  "apoptosis"        = "#8DA0CB",
  "proliferation"    = "#E78AC3",
  "angiogenesis"     = "#A6D854",
  "phagocytosis"     = "#FF6347",
  "fibrosis"         = "#FFD92F"
)

# --- Configuration ---
# FIXME: Ensure `obj_sub` is loaded in the environment from the upstream script.
seu <- obj_sub              
assay_name <- "RNA"          
slot_to_use <- "data"        
group.by <- "orig.ident"     
group_order <- c("0day", "1day", "3day", "5day", "7day") 
per_gene_z <- TRUE           
clip_val <- 3                

# --- Gene Name Mapping & Validation ---
assay_rows <- rownames(seu[[assay_name]])
map_gene <- function(g) {
  # Attempt common capitalization variants to match assay features
  candidates <- unique(c(g, toupper(g), tolower(g), paste0(toupper(substr(g,1,1)), tolower(substr(g,2,nchar(g))))))
  ok <- candidates[candidates %in% assay_rows]
  if (length(ok) > 0) return(ok[1]) else return(NA_character_)
}

genesets_mapped <- lapply(genesets_list, function(gs) {
  gmap <- sapply(gs, map_gene, USE.NAMES = FALSE)
  gmap <- gmap[!is.na(gmap)]
  unique(gmap)
})

# Report missing/unmapped genes
missing_report <- unlist(lapply(names(genesets_list), function(nm) {
  miss <- setdiff(genesets_list[[nm]], genesets_mapped[[nm]])
  if (length(miss) > 0) paste0(nm, ": ", paste(miss, collapse = ", ")) else NULL
}))
if (length(missing_report) > 0) message("Missing genes (ignored):\n", paste(missing_report, collapse = "\n"))

ordered_genes <- unlist(genesets_mapped[names(genesets_list)])
if (length(ordered_genes) == 0) stop("No mapped genes found in assay. Check gene names and assay.")

# --- Aggregate Expression by Groups ---
avg <- Seurat::AverageExpression(seu, assays = assay_name, features = ordered_genes, group.by = group.by, slot = slot_to_use)
mat_genes_x_groups <- as.matrix(avg[[assay_name]]) 

groups_present <- colnames(mat_genes_x_groups)
groups_use <- intersect(group_order, groups_present)
if (length(groups_use) == 0) groups_use <- groups_present

mat_sub <- mat_genes_x_groups[ordered_genes, groups_use, drop = FALSE] 

# --- Row-wise Z-score Normalization ---
if (isTRUE(per_gene_z)) {
  z_mat <- t(apply(mat_sub, 1, function(x) {
    if (all(is.na(x))) return(rep(NA_real_, length(x)))
    s <- sd(x, na.rm = TRUE)
    if (s == 0 || is.na(s)) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE))/s
  }))
  if (!is.null(clip_val)) {
    z_mat[z_mat > clip_val] <- clip_val
    z_mat[z_mat < -clip_val] <- -clip_val
  }
  plot_mat <- z_mat
  legend_title <- "z-score"
} else {
  plot_mat <- mat_sub
  legend_title <- "avg expression"
}

# --- Aesthetic Setup & Row Splitting ---
row2geneset <- sapply(rownames(plot_mat), function(g) {
  nm <- names(which(sapply(genesets_mapped, function(x) g %in% x)))
  if (length(nm) > 0) nm[1] else NA_character_
})
row_split <- factor(row2geneset, levels = names(genesets_list))

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b")) 
if (!isTRUE(per_gene_z)) {
  col_fun <- colorRamp2(quantile(plot_mat, probs = c(0, 0.5, 1), na.rm=TRUE), c("#f7fbff", "#91bfdb", "#08306b"))
}

height_per_gene_cm <- 0.35
min_height_cm <- 8   
plot_width_cm <- 12  
genes_n <- nrow(plot_mat)
total_height_cm <- max(min_height_cm, genes_n * height_per_gene_cm)
row_fontsize <- ifelse(genes_n <= 40, 10, ifelse(genes_n <= 80, 8, 6))

ha_left <- rowAnnotation(
  Geneset = row_split,
  col = list(Geneset = geneset_colors),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(Geneset = list(title = "Geneset", at = names(geneset_colors), labels = names(geneset_colors)))
)

# dev.off() # Historical artifact: likely used to clear plotting device during interactive sessions.

ht <- Heatmap(
  plot_mat,
  name = legend_title,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_row_slices = TRUE,         # Independent clustering within each geneset slice
  clustering_distance_rows = "pearson", 
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_split = row_split,
  row_title = NULL,
  column_title = "Timepoint",
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = row_fontsize),
  column_order = groups_use,
  row_gap = unit(3, "mm"),
  left_annotation = ha_left,         
  heatmap_legend_param = list(title = legend_title, legend_height = unit(4, "cm")),
  use_raster = TRUE,
  raster_device = "png",
  width = unit(plot_width_cm, "cm"),
  height = unit(total_height_cm, "cm")
)
print(ht)



# ==============================================================================
# #### Figure 2.Q: Myeloid Trajectory via Monocle 2 & CytoTRACE2 ####
# ==============================================================================
# Note: Immediate overwrite is preserved to maintain historical context.
temp_sce <- readRDS("/home/tsh/E-MTAB-7895/inflammation_temp_sce.rds")
temp_sce <- readRDS("/home/tsh/E-MTAB-7895/filtered_inflammation_temp_sce.rds")

temp_sce <- subset(temp_sce, subset = orig.ident %in% c("0day", "1day", "3day", "5day", "7day"))
temp_sce <- subset(temp_sce, subset = seurat_clusters %in% c("10", "13", "4", "9", "6", "12", "3", "5", "2", "0", "1", "15", "14"))
DotPlot(temp_sce, features = c("Cdh2", "Pdgfra"), group.by = "orig.ident")

stopifnot("inflammation_human" %in% colnames(temp_sce@meta.data))
Idents(temp_sce) <- "inflammation_human"

# --- Monocle 2 CDS Construction ---
sparse_data <- as(as.matrix(temp_sce[["RNA"]]$counts), "sparseMatrix")
mdata <- new("AnnotatedDataFrame", data = temp_sce@meta.data)
fData <- data.frame(gene_short_name = row.names(sparse_data), row.names = row.names(sparse_data))
fd <- new("AnnotatedDataFrame", data = fData)

monocle_cds <- newCellDataSet(
  cellData = sparse_data,
  phenoData = mdata,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 10)

# Identify highly variable genes using Seurat's selection
expressed_genes <- VariableFeatures(temp_sce)
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes, ], fullModelFormulaStr = "~inflammation_human")

# Filter genes by strict q-value (0.01 instead of default 0.1)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01)) 
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)

# DDRTree Dimensionality Reduction & Ordering
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, reduction_method = "DDRTree")
monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds, root_state = "2")
monocle_cds <- orderCells(monocle_cds, root_state = "1") # Active override for root state

save(monocle_cds, file = "./大论文图和数据/Figure 2/marco_monocle2.Rdata", compress = TRUE)

# --- Visualization ---
color_myeloid <- c(
  "Monocytes"          = "#CC79A7", 
  "Pro_inflammation"   = "#D55E00", 
  "Proliferation"      = "#E69F00", 
  "Reparative"         = "#009E73", 
  "Antigen-presenting" = "#0072B2", 
  "Collagens"          = "#56B4E9"  
)

p <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 0.6, show_backbone = TRUE) + facet_wrap("~orig.ident", nrow = 1)
p0 <- plot_cell_trajectory(monocle_cds, color_by = "orig.ident", cell_size = 0.5, show_backbone = TRUE) + facet_wrap("~orig.ident", nrow = 1)
p1 <- plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters", cell_size = 0.6, show_backbone = TRUE) +
  theme(legend.position = 'none', panel.border = element_blank()) +
  scale_color_manual(values = color_myeloid) 

p2 <- plot_complex_cell_trajectory(monocle_cds, x = 1, y = 2, color_by = "label") + theme(legend.title = element_blank())

# --- CytoTRACE2 (Myeloid) ---
cytotrace2_sce <- cytotrace2(
  temp_sce, parallelize_models = FALSE, parallelize_smoothing = FALSE, ncores = 30, 
  is_seurat = TRUE, slot_type = "data", species = 'mouse', seed = 1234
)

annotation <- data.frame(phenotype = temp_sce@meta.data$seurat_clusters) %>% set_rownames(., colnames(temp_sce))
plots <- plotData(cytotrace2_result = cytotrace2_sce, annotation = annotation, is_seurat = TRUE)
save(cytotrace2_sce, file = "./大论文图和数据/Figure 2/marco_cytotrace2_sce.Rdata", compress = TRUE)


# ==============================================================================
# #### Figure 2.S: Fibroblast Trajectory via Monocle 2 & CytoTRACE2 ####
# ==============================================================================
temp_sce <- readRDS("/home/tsh/E-MTAB-7895/temp_sce.rds") 
# Notice the addition of late timepoints (14day, 28day)
temp_sce <- subset(temp_sce, subset = orig.ident %in% c("0day", "1day", "3day", "5day", "7day", "14day", "28day"))

temp_sce$major_type <- dplyr::case_when(
  temp_sce$seurat_clusters %in% c("1", "4", "14", "12", "11") ~ "IR",
  temp_sce$seurat_clusters %in% c("9", "5", "8", "10", "13") ~ "Proliferation",
  temp_sce$seurat_clusters %in% c("3", "6", "7", "0", "2") ~ "Myofib",
  TRUE ~ "Other"
) |> factor(levels = c("IR", "Proliferation", "Myofib"))

cell_colors_fibro <- c(
  "IR"            = "#CC79A7", 
  "Proliferation" = "#E69F00", 
  "Myofib"        = "#009E73"  
)

stopifnot("seurat_clusters" %in% colnames(temp_sce@meta.data))
Idents(temp_sce) <- "seurat_clusters"

# --- Monocle 2 CDS Construction ---
sparse_data <- as(as.matrix(temp_sce[["RNA"]]$counts), "sparseMatrix")
mdata <- new("AnnotatedDataFrame", data = temp_sce@meta.data)
fData <- data.frame(gene_short_name = row.names(sparse_data), row.names = row.names(sparse_data))
fd <- new("AnnotatedDataFrame", data = fData)

monocle_cds <- newCellDataSet(
  cellData = sparse_data, phenoData = mdata, featureData = fd,
  lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size()
)

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 10)

expressed_genes <- VariableFeatures(temp_sce)
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes, ], fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01)) 
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)

# DDRTree Dimensionality Reduction & Ordering
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, reduction_method = "DDRTree")
monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds, root_state = "2")
monocle_cds <- orderCells(monocle_cds, root_state = "1")

save(monocle_cds, file = "./大论文图和数据/Figure 2/Fibro_monocle2.Rdata", compress = TRUE)

# --- Visualization ---
p <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 0.8, show_backbone = TRUE) + facet_wrap("~orig.ident", nrow = 1)
p0 <- plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters", cell_size = 0.5, show_backbone = TRUE) + facet_wrap("~orig.ident", nrow = 1)

# FIXME: Original code used `scale_color_manual(values= color)` which was defined for Myeloid cells. 
# Corrected to use `cell_colors_fibro` defined locally for Fibroblasts.
p1 <- plot_cell_trajectory(monocle_cds, color_by = "major_type", cell_size = 0.8, show_backbone = TRUE, show_tree = TRUE) +
  theme(legend.position = 'none', panel.border = element_blank()) +
  scale_color_manual(values = cell_colors_fibro) 

p2 <- plot_complex_cell_trajectory(monocle_cds, x = 1, y = 2, color_by = "label") + theme(legend.title = element_blank())

# --- CytoTRACE2 (Fibroblasts) ---
cytotrace2_sce <- cytotrace2(
  temp_sce, parallelize_models = FALSE, parallelize_smoothing = FALSE, ncores = 30, 
  is_seurat = TRUE, slot_type = "count", species = 'mouse', seed = 1234
)

# Note: Overwritten immediately. Preserved for context.
annotation <- data.frame(phenotype = temp_sce@meta.data$major_type) %>% set_rownames(., colnames(temp_sce))
annotation <- data.frame(phenotype = temp_sce@meta.data$seurat_clusters) %>% set_rownames(., colnames(temp_sce))

plots <- plotData(cytotrace2_result = cytotrace2_sce, annotation = annotation, is_seurat = TRUE)
save(cytotrace2_sce, file = "./大论文图和数据/Figure 2/Fibro_cytotrace2_sce.Rdata", compress = TRUE)
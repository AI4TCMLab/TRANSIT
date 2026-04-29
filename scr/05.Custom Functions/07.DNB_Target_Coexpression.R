# Mandatory libraries
if (!requireNamespace("ggtern")) install.packages("ggtern") # Ternary plotting
library(ggplot2); library(dplyr); library(tidyr); library(ggtern)

#### Macrophage subset UMAP marker annotation and abundance shifts ###
temp_sce <- read_rds("./filtered_inflammation_temp_sce.rds")
#### Fibroblast subset UMAP marker annotation and abundance shifts ###
temp_sce <- readRDS("/home/tsh/E-MTAB-7895/temp_sce.rds")
#### Neutrophil subset UMAP marker annotation and abundance shifts ###
temp_sce <- readRDS("/home/tsh/E-MTAB-7895/Granulocytes_temp_sce.rds")

# ---------- 1. Configuration: Gene Sets and Parameters ----------
# Merging from existing geneset tables
geneset <- readxl::read_xlsx("./大论文图和数据/TableS4 Geneset score.xlsx", sheet = "三象限图")
gs1 <- unique(na.omit(c(geneset$Inflammation))) # Inflammation + Oxidative Stress
gs2 <- unique(na.omit(c(geneset$Repair))) # Proliferation / Angiogenesis / Antigen Presentation
topN <- 50
trajectory <- c("12","5","8","3","4","1","0")
trajectory <- c("1","4","9","6","7","0")
seurat_obj <- temp_sce 
table(seurat_obj$seurat_clusters)

# Function: get_cluster_expr_matrix (Summary expression per cluster)
get_cluster_expr_matrix <- function(genes, seurat_obj, trajectory) {
  expr_mat <- matrix(NA, nrow = length(genes), ncol = length(trajectory))
  rownames(expr_mat) <- genes
  colnames(expr_mat) <- as.character(trajectory)
  for (i in seq_along(trajectory)) {
    cells <- colnames(seurat_obj)[seurat_obj$seurat_clusters == trajectory[i]]
    if (length(cells) == 0) {
      expr_mat[, i] <- NA
    } else {
      avg_expr <- rowMeans(as.matrix(seurat_obj[["RNA"]]$data[genes, cells, drop = FALSE]), na.rm = TRUE)
      expr_mat[, i] <- avg_expr
    }
  }
  return(expr_mat)
}

# === Extract Expression Matrix (Stable Interface) ===
expr_all <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "data"))

# === Simplified: Calculate correlation across all cells and select topN (Spearman) ===
# Identify top N co-expressed genes within a specific cluster
get_top_correlated <- function(
    gene,
    seurat_obj,
    cluster_id = "5",
    top_n = 50,
    pval_cutoff = 0.05,
    allowed_genes = NULL,
    balance_sign = T,          # Ensure balance between positive and negative correlations
    method = c("pearson", "spearman"),
    adjust_p = FALSE,              # Apply multiple testing correction
    p_adjust_method = "BH"         # p.adjust method (Benjamini-Hochberg)
) {
  method <- match.arg(method)
  
  # 1) Extract cells and matrix for the target cluster
  cells <- colnames(seurat_obj)[seurat_obj$seurat_clusters == cluster_id]
  if (length(cells) == 0L) {
    warning("No cells found for cluster_id = ", cluster_id)
    return(list(genes = gene, correlations = NA, signs = "Center"))
  }
  data_mat <- as.matrix(seurat_obj[["RNA"]]$data[, cells, drop = FALSE])
  
  # 2) Target gene detection
  if (!gene %in% rownames(data_mat)) {
    warning("gene not found in RNA assay: ", gene)
    return(list(genes = gene, correlations = NA, signs = "Center"))
  }
  
  # 3) Calculate correlation and significance for the target gene
  gene_expr <- data_mat[gene, ]
  # Row-wise correlation (Suppressed warnings for stability)
  cors  <- apply(data_mat, 1, function(x) suppressWarnings(cor(x, gene_expr, method = method)))
  pvals <- apply(data_mat, 1, function(x) suppressWarnings(cor.test(x, gene_expr, method = method)$p.value))
  
  # 4) Optional p-value adjustment
  p_use <- if (adjust_p) p.adjust(pvals, method = p_adjust_method) else pvals
  
  # 5) Initial filtering: Significance + Self-exclusion + Allowed Gene List
  valid_genes <- names(p_use)[p_use < pval_cutoff]
  valid_genes <- setdiff(valid_genes, gene)
  if (!is.null(allowed_genes)) valid_genes <- intersect(valid_genes, allowed_genes)
  
  # Remove genes with NA correlations
  valid_genes <- valid_genes[!is.na(cors[valid_genes])]
  if (length(valid_genes) == 0L) {
    warning("No valid genes passed the filters.")
    return(list(genes = gene, correlations = NA, signs = "Center"))
  }
  
  # 6) Selection of top_n
  take_n <- min(top_n, length(valid_genes))
  
  if (!balance_sign) {
    # Rank purely by absolute correlation magnitude
    sorted <- valid_genes[order(abs(cors[valid_genes]), decreasing = TRUE)]
    top_genes <- head(sorted, take_n)
  } else {
    # Balanced Sign: Selection of ceil(n/2) positive and floor(n/2) negative correlations
    k_pos <- take_n %/% 2 + (take_n %% 2)  # ceil
    k_neg <- take_n - k_pos                # floor
    
    pos <- valid_genes[cors[valid_genes] > 0]
    neg <- valid_genes[cors[valid_genes] < 0]
    zer <- valid_genes[cors[valid_genes] == 0]
    
    # Sort positive descending; negative by absolute magnitude (most negative first)
    pos_sorted <- pos[order(cors[pos], decreasing = TRUE)]
    neg_sorted <- neg[order(abs(cors[neg]), decreasing = TRUE)]
    
    pick_pos <- head(pos_sorted, k_pos)
    pick_neg <- head(neg_sorted, k_neg)
    picked   <- c(pick_pos, pick_neg)
    top_genes <- picked
  }
  # 6) Assemble Output (Center gene / Correlation value / Directional sign)
  correlation_signs <- ifelse(cors[top_genes] > 0, "Positive",
                              ifelse(cors[top_genes] < 0, "Negative", "Zero"))
  signs <- c("Center", correlation_signs)
  
  return(list(
    genes = c(gene, top_genes),
    correlations = c(NA, cors[top_genes]),
    pvalues = c(NA, pvals[top_genes]),  # Return p-values
    signs = signs
  ))
}

universe_genes <- rownames(expr_all)
table(dnb_class$regulator_type)
dnb_class <- read.csv("./大论文图和数据/TableS7 Macro_DNB_annotation_result.csv")
dnb_class <- read.csv("./大论文图和数据/TableS9 Fibro_DNB_annotation_result.csv")
DNB_vec <- subset(dnb_class,regulator_type == "dual_regulator") %>% .[["gene"]]
DNB_vec <- dnb_class$gene
DNB_vec <- CTS.Lib[["4"]]
DNB_vec <- CTS.Lib[["5"]]
DNB_vec <- DNB

# Retrieve log2FC values
expr_mat_use <- get_cluster_expr_matrix(rownames(seurat_obj), seurat_obj, trajectory)
expr_mat_use <- expr_mat_use[rowMeans(expr_mat_use[, c("3", "12")]) != 0, ]
log2fc_vector <- log2((expr_mat_use[, "3"] + 1e-6) / 
                        (expr_mat_use[, "12"] + 1e-6))
names(log2fc_vector) <- rownames(expr_mat_use)
Log2fc_vector <- log2fc_vector
Bg_genes <- universe_genes
hist(log2fc_vector)

degs <- FindMarkers(seurat_obj, ident.1 = "3", ident.2 = "12", group.by = "seurat_clusters",only.pos = FALSE, min.pct =0, logfc.threshold = 0,test.use	= "wilcox")
deg_genes <- rownames(degs)[abs(degs$avg_log2FC) > 0.5  & degs$p_val_adj < 0.05]

g <- "Lyar"
# Traverse DNB list
res_list <- vector("list", length(DNB_vec))
names(res_list) <- DNB_vec

# Function: calculate_mean_and_pvalue (Non-parametric permutation test)
calculate_mean_and_pvalue <- function(coexpr_genes, 
                                      Log2fc_vector, 
                                      Bg_genes, 
                                      n_perm = 1000) {
  # Parameter verification before execution
  print(head(coexpr_genes))  # Verify positive correlation genes
  print(head(Log2fc_vector))   # Verify log2FC vector
  print(head(Bg_genes))        # Verify background pool
  
  # Calculate observed mean log2FC
  observed_mean <- mean(Log2fc_vector[coexpr_genes], na.rm = TRUE)
  
  # Determine background pool directionality based on observed mean
  if (observed_mean >= 0) {
    bg_pool <- Bg_genes[Log2fc_vector[Bg_genes] >= 0]
  } else {
    bg_pool <- Bg_genes[Log2fc_vector[Bg_genes] < 0]
  }
  
  # Exclude outliers and ensure valid vectors
  bg_pool <- bg_pool[abs(Log2fc_vector[bg_pool]) < 5 & !is.na(Log2fc_vector[bg_pool])]
  
  # Ensure background pool size is sufficient for sampling
  if (length(bg_pool) < length(coexpr_genes) || length(bg_pool) == 0) {
    warning(paste("Insufficient background pool: required", length(coexpr_genes), "but only", length(bg_pool), "available"))
    return(list(observed_mean = observed_mean, p_value = NA))
  }
  
  # Execute permutation test
  permuted_means <- replicate(n_perm, {
    sampled_genes <- sample(bg_pool, length(coexpr_genes))
    mean(Log2fc_vector[sampled_genes], na.rm = TRUE)
  })
  
  # Calculate empirical p-value
  if (observed_mean >= 0) {
    p_val <- mean(permuted_means >= observed_mean)
  } else {
    p_val <- mean(permuted_means <= observed_mean)
  }

  return(list(observed_mean = observed_mean, p_value = p_val))
}

g <- "Mdk"
# Initialize output dataframes
res_list <- vector("list", length(DNB_vec))
names(res_list) <- DNB_vec
result_table <- data.frame(DNB = character(), TARGET = character(), TYPE = character(), stringsAsFactors = FALSE)
dnb_coexpr_genes_list <- list()

# Main Loop: Iterate through each DNB molecule
for (g in DNB_vec) {
  cat("Processing:", g, "...\n")
  bg_genes <- setdiff(names(log2fc_vector), g)  # Exclude self from background
  # Identify co-expressed genes for the target DNB
  result <- get_top_correlated(g, seurat_obj, cluster_id = "5", top_n = 50, pval_cutoff = 0.05, 
                               allowed_genes = deg_genes, method = "spearman", balance_sign = TRUE)
  
  if (is.null(result) || length(result$genes) < 2) {
    # Skip if no significant correlations found
    res_list[[g]] <- data.frame(gene = g, n_top = 0, positive_mean_log2FC = NA, negative_mean_log2FC = NA,
                                positive_p_value = NA, negative_p_value = NA,
                                stringsAsFactors = FALSE)
    next
  }
  # Save co-expressed genes if identification was successful
  if (!is.null(result) && length(result$genes) > 1) {
    top_genes <- result$genes[-1] # Exclude target gene itself
    dnb_coexpr_genes_list[[g]] <- top_genes
  } else {
    dnb_coexpr_genes_list[[g]] <- NULL
  }
  
  top_genes <- result$genes[-1]
  
  # Separate positive and negative correlation modules
  positive_genes <- top_genes[result$signs[top_genes] == "Positive"]
  negative_genes <- top_genes[result$signs[top_genes] == "Negative"]
  
  # Significance testing for Positive modules
  positive_results <- calculate_mean_and_pvalue(positive_genes, log2fc_vector, bg_genes)
  positive_mean_log2FC <- positive_results$observed_mean
  positive_p_value <- positive_results$p_value
  
  # Significance testing for Negative modules
  negative_results <- calculate_mean_and_pvalue(negative_genes, log2fc_vector, bg_genes)
  negative_mean_log2FC <- negative_results$observed_mean
  negative_p_value <- negative_results$p_value
  
  # Construct result dataframe for the current DNB
  df <- data.frame(
    gene = g,
    n_top = length(top_genes), 
    positive_mean_log2FC = positive_mean_log2FC, 
    negative_mean_log2FC = negative_mean_log2FC, 
    positive_p_value = positive_p_value, 
    negative_p_value = negative_p_value, 
    stringsAsFactors = FALSE
  )
  
  res_list[[g]] <- df
  
  # Mapping DNB to target gene correlation attributes
  cor_values <- result$correlations[-1]
  cor_direction <- result$signs[-1]
  cor_p_values <- result$pvalues[-1]
  cor_log2fc <- log2fc_vector[result$genes[-1]]
  
  # Create mapping for current DNB and targets
  current_data <- data.frame(DNB = rep(g, length(top_genes)), 
                             TARGET = top_genes, 
                             COR = cor_direction, 
                             Cor_values = cor_values, 
                             Cor_p_values = cor_p_values, 
                             Cor_log2fc = cor_log2fc,  
                             stringsAsFactors = FALSE)
  
  result_table <- rbind(result_table, current_data)
}

# Consolidate results into a single dataframe
res_df <- bind_rows(res_list)

head(res_df)
result_table <- read_csv("./大论文图和数据/TableS9.Macro_dnb_coexp.csv")
res_df <- read_csv("./大论文图和数据/TableS10.Macro_dnb_coexp_plot.csv")
res_df <- read_csv("./大论文图和数据/TableS12.Fibro_dnb_coexp_plot.csv")

### 2D Quadrant Plot Generation
library(ggplot2)
res_df <- res_df %>%
  left_join(dnb_class, by = c("gene" = "gene"))

library(dplyr)
library(ggrepel)
# Log transformation of p-values for node scaling
res_df <- res_df %>%
  mutate(
    log_positive_p_value = -log10(positive_p_value+1e-6),
    log_negative_p_value = -log10(negative_p_value+1e-6 ),
    size_factor = pmin(log_positive_p_value, log_negative_p_value),
    
    # Categorization of significance domains
    color_factor = case_when(
      positive_p_value < 0.2 & negative_p_value < 0.2 ~ "red",  # Significant Dual
      negative_p_value < 0.2 ~ "green",                      # Significant Negative
      positive_p_value < 0.2 ~ "blue",                       # Significant Positive
      TRUE ~ "grey"                                            # Non-significant
    )
  )

table(res_df$color_factor)
selected_DNB  <- res_df$gene[res_df$color_factor != "grey"]

# Scatter plot: 2D Quadrant Evaluation using ggrepel for annotation stability
ggplot(res_df, aes(x = negative_mean_log2FC, y = positive_mean_log2FC)) +
  geom_point(aes(color = color_factor, size = size_factor), alpha = 0.6) + 
  geom_text_repel(aes(label = gene), size = 3.6, box.padding = 0.35, point.padding = 0.5, 
                  max.overlaps = 10, color = "black") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + 
  scale_size_continuous(range = c(3, 10), name = "-log10(P-value)") + 
  scale_color_manual(values = c("red" = "red", "green" = "green", "blue" = "blue", "grey" = "grey"), 
                     name = "Significance", 
                     labels = c("Pos Significant","Neg Significant", "Non-Significant","Dual Significant")) + 
  labs(title = "2D Quadrant Evaluation of Target Gene Regulation",
       x = "Inflammation Mean Log2FC",
       y = "Repair Mean Log2FC",
       color = "Gene Type", 
       size = "Node Size (-log10(p-value))") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.title = element_text(size = 12), legend.text = element_text(size = 10))

# Consolidate unique co-expressed targets across all DNBs
all_coexpr_genes <- unique(unlist(dnb_coexpr_genes_list))
all_coexpr_genes <- unique(result_table$TARGET)
library(ClusterGVis)
cat("Total number of unique coexpressed genes across all DNBs:", length(all_coexpr_genes), "\n")
expr_mat <- expr_mat_use[all_coexpr_genes,]
expr_mat <- expr_mat[apply(expr_mat, 1, function(x) sd(x, na.rm = TRUE) > 0), , drop = FALSE]
all_genes <- rownames(expr_mat)

# Mfuzz clustering for temporal trend discovery
ck <- clusterData(exp = expr_mat,
                  cluster.method = "mfuzz",
                  cluster.num = 5)
visCluster(object = ck,
           plot.type = "both",
           add.mline = F,
           add.line = F,
           add.bar = F,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           show_row_names = F,
           column_names_rot = 45)
write.csv(ck$long.res,"./大论文图和数据/TableS_Fibroblasts_Coexp_Clusters.csv")

# Identify functional gene types from clusters
inflammation_genes <- ck$wide.res %>% filter(cluster == 1) %>% .[["gene"]]
repair_genes <- ck$wide.res %>% filter(cluster == 2) %>% .[["gene"]]
inflammation_genes <- ck$wide.res %>% filter(cluster == 5) %>% .[["gene"]]
ECM_genes <- ck$wide.res %>% filter(cluster == 2) %>% .[["gene"]]
Prolif_genes <- ck$wide.res %>% filter(cluster %in% c(1,3)) %>% .[["gene"]]
Vascular_genes <- ck$wide.res %>% filter(cluster == 4) %>% .[["gene"]]

### Late-stage Fibroblast Classification ###
inflammation_genes <- ck$wide.res %>% filter(cluster == c(4)) %>% .[["gene"]]
repair_genes <- ck$wide.res %>% filter(cluster %in% c(2,5)) %>% .[["gene"]]
ECM_genes <- ck$wide.res %>% filter(cluster %in% c(1,3)) %>% .[["gene"]]

# Define functional types in result table
result_table <- result_table %>%
  mutate(TYPE = case_when(
    TARGET %in% inflammation_genes ~ "Inflam", 
    TARGET %in% ECM_genes ~ "Fibrosis",     
    TARGET %in% repair_genes ~ "Repair",     
    TRUE ~ "Other"                           
  ))

expr_mat_use["Atf4",]

# Exploratory visualizations for feature validation
FeaturePlot(temp_sce,features = "Atf4",split.by = "orig.ident",reduction = "umap")
VlnPlot(temp_sce,features = "Cdh2",group.by =  "orig.ident")
DotPlot(temp_sce,features = "Irg1",group.by =  "orig.ident")
DotPlot(temp_sce,features = "Nfkb1",group.by =  "orig.ident")
DotPlot(temp_sce,features = c("Nos2"),group.by =  "seurat_clusters")
FeaturePlot(temp_sce,features = "Ly6c2",split.by = "orig.ident",reduction = "umap")

### Late-stage Macrophage Classification ###
inflammation_genes <- ck$wide.res %>% filter(cluster == c(1)) %>% .[["gene"]]
repair_genes <- ck$wide.res %>% filter(cluster == c(4,5)) %>% .[["gene"]]
ECM_genes <- ck$wide.res %>% filter(cluster == 2) %>% .[["gene"]]

result_table <- result_table %>%
  mutate(TYPE = case_when(
    TARGET %in% inflammation_genes ~ "Inflam",
    TARGET %in% ECM_genes ~ "Fibrosis",     
    TARGET %in% repair_genes ~ "Repair",     
    TRUE ~ "Other"                           
  ))

result_table <- result_table %>%
  mutate(TYPE = case_when(
    TARGET %in% inflammation_genes ~ "Inflam",
    TARGET %in% ECM_genes ~ "Fibrosis",     
    TARGET %in% Prolif_genes ~ "Prolif",     
    TARGET %in% Vascular_genes ~ "Vascular",     
    TRUE ~ "Other"                           
  ))

# GO Enrichment: Biological Process analysis for Inflammatory modules
inflammation_enrich <- enrichGO(
  gene = inflammation_genes, 
  OrgDb = org.Mm.eg.db,       # Species: Mouse
  keyType = "SYMBOL",         
  ont = "BP",                 
  pAdjustMethod = "BH",       
  qvalueCutoff = 0.05         
)

head(inflammation_enrich)
dotplot(inflammation_enrich, showCategory = 10)
barplot(inflammation_enrich, showCategory = 5) 

repair_enrich <- enrichGO(
  gene = Prolif_genes, 
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL", 
  ont = "BP", 
  pAdjustMethod = "BH", 
  qvalueCutoff = 0.05
)

head(repair_enrich)
barplot(repair_enrich, showCategory = 5) 

result_table <- result_table %>%
  mutate(TYPE = case_when(
    TARGET %in% inflammation_genes ~ "Inflam", 
    TARGET %in% repair_genes ~ "Repair",      
    TRUE ~ "Other"                           
  ))

# Save archival DNB execution data
save(DNB_vec,result_table, expr_mat,ck,res_df,dnb_coexpr_genes_list,expr_mat_use,file = "./Figure3/Fibor_early_analysis_DNB_data.RData",compress = T)

head(result_table)
result_table <- result_table %>%
  mutate(DNB_log2fc = log2fc_vector[match(DNB, names(log2fc_vector))])
head(result_table)
getwd()
write.csv(result_table, file = "Macro_analysis_DNB_data.csv", row.names = FALSE)
selected_result_table <- subset(result_table,DNB %in% selected_DNB)
write.csv(selected_result_table, file = "Fibro_analysis_DNB_p0.2_early_data.csv", row.names = FALSE)

result_table <- read.csv("./Fibro_analysis_DNB_data.csv")
result_table <- read.csv("./Macro_analysis_DNB_data.csv")
unique_dnbs <- DNB_vec
dnb <- unique_dnbs[c(1,4)]
library(dplyr)

# Specific DNB analysis loop for topology generation
dnb <- c("Fabp4")
dnb <- c("Tnfaip6")
dnb <- c("Ubash3b")
dnb <- c("Sc5d","F10","Tes","Prdx6","Atf4","Lyar")[6]

for (dnb in unique_dnbs) {
  # Load Sankey mapping data
  dnb_data <- read_xlsx("./大论文图和数据/ASII_Rd_Sankeyplot.xlsx")
  colnames(dnb_data)[2] <- "TARGET" 
  colnames(dnb_data)[5] <- "COR" 
  dnb <- unique(dnb_data$DNB)
  dnb_data <- result_table[result_table$DNB %in% dnb, ]
  
  # Determine node order for current DNB module
  dnb_gene_order <- dnb_data %>%
    count(TYPE, TARGET) %>%
    arrange(TYPE, desc(n)) %>%
    pull(TARGET)
  
  cat(paste0("\n=== DNB: ", dnb, " - Gene Order and Correlation Summary ===\n"))
  correlation_summary <- dnb_data %>%
    dplyr::select(TARGET, TYPE, COR) %>%
    arrange(TYPE, TARGET)
  print(correlation_summary)
  
  # Construct hierarchical nodes
  dnb_nodes <- dnb
  cluster_nodes <- paste0("Cluster_", sort(unique(dnb_data$TYPE)))
  coexgene_nodes <- dnb_gene_order
  
  nodes <- data.frame(name = c(dnb_nodes, cluster_nodes, coexgene_nodes))
  
  # Layer 1: DNB -> Cluster linkage
  layer1_links_dnb <- dnb_data %>%
    group_by(DNB, TYPE) %>%
    summarise(Value = n(), .groups = 'drop') %>%
    mutate(
      source = match(DNB, nodes$name) - 1,
      target = match(paste0("Cluster_", TYPE), nodes$name) - 1,
      group = as.factor(paste0("Cluster_", TYPE))
    ) %>%
    dplyr::select(source, target, value = Value, group)
  
  # Layer 2: Cluster -> Co-expressed Target linkage
  layer2_links_dnb <- dnb_data %>%
    mutate(
      source = match(paste0("Cluster_", TYPE), nodes$name) - 1,
      target = match(TARGET, nodes$name) - 1,
      value = 1,
      group = as.factor(COR) 
    ) %>%
    dplyr::select(source, target, value, group)
  
  links_dnb <- bind_rows(layer1_links_dnb, layer2_links_dnb)
  
  # Configuration of D3 color scales for Sankey topology
  dnb_all_groups <- unique(c(as.character(layer1_links_dnb$group), 
                             as.character(layer2_links_dnb$group)))
  dnb_cluster_groups <- grep("^Cluster_", dnb_all_groups, value = TRUE)
  dnb_cluster_colors <- rainbow(length(dnb_cluster_groups))
  
  dnb_direction_colors <- c()
  if ("Positive" %in% dnb_all_groups) dnb_direction_colors["Positive"] <- "red"
  if ("positive" %in% dnb_all_groups) dnb_direction_colors["positive"] <- "red"
  if ("Negative" %in% dnb_all_groups) dnb_direction_colors["Negative"] <- "green"
  if ("negative" %in% dnb_all_groups) dnb_direction_colors["negative"] <- "green"
  if ("unknown" %in% dnb_all_groups) dnb_direction_colors["unknown"] <- "gray"
  
  dnb_all_colors <- c(dnb_cluster_colors, dnb_direction_colors)
  dnb_all_domains <- c(dnb_cluster_groups, names(dnb_direction_colors))
  
  dnb_color_scale_js <- paste0('d3.scaleOrdinal() .domain(["', 
                               paste(dnb_all_domains, collapse = '", "'), 
                               '"]) .range(["', 
                               paste(dnb_all_colors, collapse = '", "'), 
                               '"])')
  
  # Render Sankey Diagram
  sankey_plot <- sankeyNetwork(
    Links = links_dnb,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    LinkGroup = "group",
    units = "",
    fontSize = 14,
    nodeWidth = 20,
    height = 800,
    width = 900,
    nodePadding = 8,
    sinksRight = FALSE,
    iterations = 0, # Fixed node positioning
    colourScale = JS(dnb_color_scale_js)
  )
  
  output_dir <- "/home/tsh/E-MTAB-7895/大论文图和数据/"
  html_file <- file.path(output_dir, paste0("Figure2-3A Lyar_Sankey", gsub("[^A-Za-z0-9]", "_", dnb), ".html"))
  pdf_file <- file.path(output_dir, paste0("Figure2-3A F10_Sankey", gsub("[^A-Za-z0-9]", "_", dnb), ".pdf"))
  
  library(plotly)
  library(htmlwidgets)
  saveWidget(sankey_plot, html_file)
  
  # Generate synchronized heatmap for current DNB module
  existing_genes <- dnb_gene_order[dnb_gene_order %in% rownames(expr_mat)]
  
  if (length(existing_genes) > 0) {
    heatmap_data <- expr_mat[existing_genes, , drop = FALSE]
    heatmap_data <- heatmap_data[dnb_gene_order[dnb_gene_order %in% existing_genes], , drop = FALSE]
    heatmap_data <- dnb_data[,c("TARGET","logFC","Cor_values")] %>% column_to_rownames("TARGET") %>% .[dnb_gene_order,,drop = FALSE]
    
    annotation_row <- data.frame(
      Correlation = dnb_data$COR[match(rownames(heatmap_data), dnb_data$TARGET)],
      Cluster = as.factor(dnb_data$TYPE[match(rownames(heatmap_data), dnb_data$TARGET)])
    )
    rownames(annotation_row) <- rownames(heatmap_data)
    
    annotation_colors <- list(
      Correlation = c(Positive = "red", Negative = "green", unknown = "gray"),
      Cluster = setNames(rainbow(length(unique(annotation_row$Cluster))), 
                         unique(annotation_row$Cluster))
    )
    
    pheatmap(heatmap_data,
             main = paste("Expression Heatmap for DNB:", dnb, "\n(Red: Positive, Green: Negative)"),
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             scale = "row", border_color = "white",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             show_rownames = TRUE,
             show_colnames = TRUE,
             annotation_row = annotation_row,
             annotation_colors = annotation_colors)
    dev.off()
  }
  
  # Archive gene order and correlation metrics for publication tables
  dnb_order_df <- dnb_data %>%
    select(DNB, CoexGene, Cluster, Coexp_Direction, correlation) %>%
    arrange(Cluster, CoexGene)
  
  write.table(dnb_order_df,
              file.path(output_dir, paste0("DNB_", gsub("[^A-Za-z0-9]", "_", dnb), "_gene_order.txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

write.csv(heatmap_data,"adam12_coexp_heatmap.csv")
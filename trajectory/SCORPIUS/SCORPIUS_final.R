#!/usr/bin/env Rscript

############################################################
## EW2 SCORPIUS trajectory reconstruction
## Expression-only exploratory trajectory analysis
##
## This script harmonises SCORPIUS input preparation with the
## Monocle3, Slingshot and PAGA expression-only TI analyses:
##   - EW2 cells are extracted from the same annotated Seurat object.
##   - The RNA assay is used.
##   - BCR/TCR-like metadata are removed before trajectory inference.
##   - Immunoglobulin and TCR receptor genes are excluded from features.
##   - PCA/UMAP used for visualisation are recomputed after receptor-gene exclusion.
##   - SCORPIUS pseudotime is oriented post hoc using early B-cell states.
############################################################

############################
## 0. Setup
############################

setwd("/localdisk/home/s2845297/dissertation/trajectory/SCORPIUS_expression_only")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SCORPIUS)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  library(patchwork)
  library(grid)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

############################
## 1. Paths and parameters
############################

seurat_rds <- "/localdisk/home/s2845297/dissertation/trajectory/EW2_annotated.rds"

outdir  <- "SCORPIUS_expression_only_final"
fig_dir <- file.path(outdir, "figures")
tab_dir <- file.path(outdir, "tables")
rds_dir <- file.path(outdir, "rds")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

## Dataset / assay settings
sample_col <- "orig.ident"
target_sample <- "EW2"
assay_use <- "RNA"
celltype_col <- "celltype"

## Optional cell-type restriction. Keep NULL to use all EW2 cells.
bcell_types <- NULL

## Expression-only feature filtering, harmonised with other TI scripts
exclude_receptor_genes <- TRUE
min_cells_expressed <- 10
n_hvg <- 2000

## Recompute expression-only PCA/UMAP for plotting and cross-method consistency
rerun_expression_only_umap <- TRUE
num_dim <- 100

## SCORPIUS pseudotime orientation
## SCORPIUS does not use a root during trajectory inference; direction is arbitrary.
## The final pseudotime direction is oriented so that these early B-cell states
## have lower pseudotime values where possible.
root_celltypes <- c("Transitional", "NAIVE 1", "NAIVE 2")

## Plotting parameters
n_pt_bins <- 60
cell_point_size <- 0.9
cell_alpha <- 0.85
label_size <- 3.5
plot_width <- 8.5
plot_height <- 6.5
plot_dpi <- 300

############################
## 2. Save package versions
############################

packages_used <- c(
  "Seurat", "SeuratObject", "SCORPIUS", "Matrix", "dplyr", "tibble",
  "readr", "ggplot2", "ggrepel", "viridis", "patchwork"
)

package_versions <- data.frame(
  package = packages_used,
  version = sapply(packages_used, function(pkg) as.character(packageVersion(pkg))),
  row.names = NULL
)

write.csv(
  package_versions,
  file.path(outdir, "SCORPIUS_package_versions.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "sessionInfo_SCORPIUS.txt")
)

############################
## 3. Colour-blind friendly palette
############################

celltype_order <- c(
  "Transitional",
  "NAIVE 1",
  "NAIVE 2",
  "Early-Activation",
  "IgD+ve Memory",
  "IgD-ve Memory",
  "HSP+",
  "DN2",
  "ASC 1",
  "ASC 2",
  "ASC 3",
  "Cell-Cycling"
)

celltype_cols <- c(
  "Transitional"     = "#332288",
  "NAIVE 1"          = "#88CCEE",
  "NAIVE 2"          = "#44AA99",
  "Early-Activation" = "#117733",
  "IgD+ve Memory"    = "#999933",
  "IgD-ve Memory"    = "#DDCC77",
  "HSP+"             = "#EE7733",
  "DN2"              = "#CC6677",
  "ASC 1"            = "#882255",
  "ASC 2"            = "#AA4499",
  "ASC 3"            = "#CC3311",
  "Cell-Cycling"     = "#661100"
)

############################
## 4. Helper functions
############################

save_plot <- function(p, filename, width = plot_width, height = plot_height) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = plot_dpi,
    bg = "white"
  )
}

get_assay_counts <- function(seu, assay = "RNA") {
  mat <- tryCatch(
    {
      GetAssayData(seu, assay = assay, layer = "counts")
    },
    error = function(e) {
      GetAssayData(seu, assay = assay, slot = "counts")
    }
  )
  as(mat, "dgCMatrix")
}

get_assay_data_safe <- function(seu, assay = "RNA", layer_or_slot = "data") {
  mat <- tryCatch(
    {
      GetAssayData(seu, assay = assay, layer = layer_or_slot)
    },
    error = function(e) {
      GetAssayData(seu, assay = assay, slot = layer_or_slot)
    }
  )
  mat
}

is_receptor_gene_symbol <- function(gene_symbol) {
  gene_symbol <- as.character(gene_symbol)

  receptor_patterns <- c(
    "^IGHV", "^IGHD[0-9]", "^IGHJ", "^IGHC",
    "^IGHA[0-9]*$", "^IGHG[0-9]*$", "^IGHM$", "^IGHD$", "^IGHE$",
    "^IGKV", "^IGKJ", "^IGKC$",
    "^IGLV", "^IGLJ", "^IGLC[0-9]*$",
    "^TRAV", "^TRAJ", "^TRAC$",
    "^TRBV", "^TRBD", "^TRBJ", "^TRBC[0-9]*$",
    "^TRGV", "^TRGJ", "^TRGC[0-9]*$",
    "^TRDV", "^TRDD", "^TRDJ", "^TRDC$"
  )

  grepl(
    paste(receptor_patterns, collapse = "|"),
    gene_symbol
  )
}

select_top_variable_features <- function(seu, assay, candidate_genes, n_features) {
  DefaultAssay(seu) <- assay

  ## Request more variable features than ultimately needed, because receptor genes
  ## and genes outside the expression-only candidate pool are removed afterwards.
  n_request <- min(
    max(n_features * 3, n_features + 1000),
    length(candidate_genes)
  )

  seu <- FindVariableFeatures(
    seu,
    assay = assay,
    selection.method = "vst",
    nfeatures = n_request,
    verbose = FALSE
  )

  hvg <- VariableFeatures(seu)
  hvg <- hvg[hvg %in% candidate_genes]

  if (length(hvg) < min(100, n_features)) {
    message("Too few candidate genes were recovered from VariableFeatures(); ranking candidate genes by variance in log-normalised data.")
    expr_data <- get_assay_data_safe(seu, assay = assay, layer_or_slot = "data")
    candidate_genes <- intersect(candidate_genes, rownames(expr_data))
    expr_sub <- as.matrix(expr_data[candidate_genes, , drop = FALSE])
    gene_var <- apply(expr_sub, 1, var, na.rm = TRUE)
    gene_var <- sort(gene_var, decreasing = TRUE)
    hvg <- names(gene_var)[seq_len(min(length(gene_var), n_features))]
  }

  hvg <- hvg[seq_len(min(length(hvg), n_features))]
  list(seu = seu, hvg = hvg)
}

############################
## 5. Load Seurat object and extract EW2
############################

if (!file.exists(seurat_rds)) {
  stop("Seurat RDS file not found: ", seurat_rds)
}

seu <- readRDS(seurat_rds)

cat("Total cells in Seurat object:", ncol(seu), "\n")
cat("Available assays:", paste(SeuratObject::Assays(seu), collapse = ", "), "\n")
cat("Available reductions:", paste(Reductions(seu), collapse = ", "), "\n")

DefaultAssay(seu) <- assay_use

if (sample_col %in% colnames(seu@meta.data)) {
  sample_values <- as.character(seu@meta.data[[sample_col]])

  if (target_sample %in% sample_values) {
    seu <- subset(
      seu,
      cells = colnames(seu)[sample_values == target_sample]
    )
  } else {
    message("target_sample not found in orig.ident. Using the full Seurat object.")
  }

} else {
  message("sample_col not found in metadata. Using the full Seurat object.")
}

cat("Cells after EW2 subsetting:", ncol(seu), "\n")

DefaultAssay(seu) <- assay_use

if (inherits(seu[[assay_use]], "Assay5")) {
  seu[[assay_use]] <- JoinLayers(seu[[assay_use]])
}

if (!celltype_col %in% colnames(seu@meta.data)) {
  stop("celltype_col not found in metadata: ", celltype_col)
}

if (!is.null(bcell_types)) {
  keep_cells <- rownames(seu@meta.data)[
    as.character(seu@meta.data[[celltype_col]]) %in% bcell_types
  ]
  if (length(keep_cells) == 0) {
    stop("No cells matched bcell_types.")
  }
  bcell <- subset(seu, cells = keep_cells)
} else {
  bcell <- seu
}

present_celltypes <- unique(as.character(bcell@meta.data[[celltype_col]]))
ordered_present <- c(
  celltype_order[celltype_order %in% present_celltypes],
  setdiff(sort(present_celltypes), celltype_order)
)

bcell@meta.data[[celltype_col]] <- factor(
  as.character(bcell@meta.data[[celltype_col]]),
  levels = ordered_present
)

celltype_cols_present <- celltype_cols[ordered_present]
missing_cols <- setdiff(ordered_present, names(celltype_cols_present))
if (length(missing_cols) > 0) {
  extra_cols <- viridis::viridis(length(missing_cols), option = "D")
  names(extra_cols) <- missing_cols
  celltype_cols_present <- c(celltype_cols_present, extra_cols)
}

cat("Cells used for SCORPIUS:", ncol(bcell), "\n")
print(table(bcell@meta.data[[celltype_col]], useNA = "ifany"))

write_csv(
  tibble(celltype = names(table(bcell@meta.data[[celltype_col]])),
         n_cells = as.integer(table(bcell@meta.data[[celltype_col]]))),
  file.path(tab_dir, "EW2_SCORPIUS_celltype_counts.csv")
)

saveRDS(
  bcell,
  file.path(rds_dir, "EW2_seurat_for_SCORPIUS_before_expression_filtering.rds")
)

############################
## 6. Remove receptor-derived metadata
############################

bcr_like_cols <- c(
  "clone_id", "cloneId", "clone", "clonotype", "CloneID",
  "clone_subgroup_id", "sub_clone_id",
  "isotype", "c_call", "constant_call",
  "v_call", "j_call", "d_call",
  "cdr3", "cdr3_aa", "junction", "junction_aa",
  "shm", "SHM", "mu_freq", "v_mutation_rate", "v_identity",
  "clone_size", "clone_expanded",
  "productive", "locus", "sequence_id"
)

bcr_like_cols_present <- intersect(
  bcr_like_cols,
  colnames(bcell@meta.data)
)

message("BCR/TCR-like metadata removed before SCORPIUS inference:")
print(bcr_like_cols_present)

write_csv(
  tibble(removed_metadata_column = bcr_like_cols_present),
  file.path(tab_dir, "EW2_SCORPIUS_removed_BCR_TCR_metadata_columns.csv")
)

bcell@meta.data <- bcell@meta.data[
  ,
  setdiff(colnames(bcell@meta.data), bcr_like_cols_present),
  drop = FALSE
]

############################
## 7. Select expression-only feature pool
############################

counts_mat <- get_assay_counts(bcell, assay = assay_use)

if (nrow(counts_mat) == 0 || ncol(counts_mat) == 0) {
  stop("Counts matrix is empty.")
}

gene_df <- tibble(
  gene_id = rownames(counts_mat),
  gene_short_name = rownames(counts_mat),
  num_cells_expressed = Matrix::rowSums(counts_mat > 0)
)

expressed_genes <- gene_df %>%
  filter(num_cells_expressed >= min_cells_expressed)

if (exclude_receptor_genes) {
  expressed_genes <- expressed_genes %>%
    mutate(is_receptor_gene = is_receptor_gene_symbol(gene_short_name))

  receptor_genes_removed <- expressed_genes %>%
    filter(is_receptor_gene)

  trajectory_gene_pool <- expressed_genes %>%
    filter(!is_receptor_gene) %>%
    pull(gene_id)

  write_csv(
    receptor_genes_removed,
    file.path(tab_dir, "EW2_SCORPIUS_receptor_genes_removed_from_feature_pool.csv")
  )

} else {
  expressed_genes <- expressed_genes %>%
    mutate(is_receptor_gene = FALSE)
  trajectory_gene_pool <- expressed_genes$gene_id
}

cat("Number of expressed genes:", nrow(expressed_genes), "\n")
cat("Number of expressed non-receptor genes in SCORPIUS feature pool:", length(trajectory_gene_pool), "\n")

if (length(trajectory_gene_pool) < 500) {
  stop("Too few expression-only trajectory genes. Please check gene symbols or filtering threshold.")
}

write_csv(
  tibble(gene = trajectory_gene_pool),
  file.path(tab_dir, "EW2_SCORPIUS_expression_only_gene_pool.csv")
)

############################
## 8. Recompute expression-only normalisation, HVGs, PCA and UMAP
############################

DefaultAssay(bcell) <- assay_use

bcell <- NormalizeData(
  bcell,
  assay = assay_use,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

hvg_result <- select_top_variable_features(
  seu = bcell,
  assay = assay_use,
  candidate_genes = trajectory_gene_pool,
  n_features = n_hvg
)

bcell <- hvg_result$seu
hvg <- hvg_result$hvg
hvg <- intersect(hvg, rownames(bcell))

if (length(hvg) < 100) {
  stop("Too few HVGs available for SCORPIUS after filtering: ", length(hvg))
}

hvg_table <- tibble(
  gene = hvg,
  is_receptor_gene = is_receptor_gene_symbol(hvg)
)

if (any(hvg_table$is_receptor_gene)) {
  stop("Receptor genes are still present in the final SCORPIUS HVG set.")
}

write_csv(
  hvg_table,
  file.path(tab_dir, "EW2_SCORPIUS_expression_only_HVGs.csv")
)

cat("Genes used for SCORPIUS:", length(hvg), "\n")

if (rerun_expression_only_umap) {
  message("Recomputing PCA and UMAP for SCORPIUS visualisation using receptor-gene-excluded features.")

  pca_features <- intersect(trajectory_gene_pool, rownames(bcell))
  if (length(pca_features) < 500) {
    stop("Too few genes available for expression-only PCA/UMAP.")
  }

  bcell <- ScaleData(
    bcell,
    assay = assay_use,
    features = pca_features,
    verbose = FALSE
  )

  npcs_use <- min(num_dim, ncol(bcell) - 1, length(pca_features) - 1)
  if (npcs_use < 2) {
    stop("Too few cells/features to compute PCA.")
  }

  bcell <- RunPCA(
    bcell,
    assay = assay_use,
    features = pca_features,
    npcs = npcs_use,
    verbose = FALSE
  )

  bcell <- RunUMAP(
    bcell,
    reduction = "pca",
    dims = seq_len(npcs_use),
    reduction.name = "umap",
    reduction.key = "UMAP_",
    verbose = FALSE
  )

  saveRDS(
    bcell,
    file.path(rds_dir, "EW2_seurat_expression_only_for_SCORPIUS_visualisation.rds")
  )
}

if (!"umap" %in% Reductions(bcell)) {
  stop("UMAP reduction not found. Please enable rerun_expression_only_umap or provide an existing UMAP reduction.")
}

############################
## 9. Build SCORPIUS input matrix
############################

expr_mat <- get_assay_data_safe(
  bcell,
  assay = assay_use,
  layer_or_slot = "data"
)

expr_use <- expr_mat[hvg, , drop = FALSE]
expr_scorpius <- t(as.matrix(expr_use))

gene_sd <- apply(expr_scorpius, 2, sd, na.rm = TRUE)
expr_scorpius <- expr_scorpius[, gene_sd > 0, drop = FALSE]
expr_scorpius <- scale(expr_scorpius)
expr_scorpius[is.na(expr_scorpius)] <- 0
expr_scorpius[is.infinite(expr_scorpius)] <- 0

cat("SCORPIUS input matrix dimensions:\n")
print(dim(expr_scorpius))

saveRDS(
  expr_scorpius,
  file.path(rds_dir, "EW2_SCORPIUS_input_expression_matrix.rds")
)

write_csv(
  tibble(gene = colnames(expr_scorpius)),
  file.path(tab_dir, "EW2_SCORPIUS_final_input_genes_after_zero_sd_filter.csv")
)

############################
## 10. Run SCORPIUS
############################

space <- SCORPIUS::reduce_dimensionality(
  expr_scorpius,
  "spearman"
)

traj <- SCORPIUS::infer_trajectory(space)

if (is.null(traj$time)) {
  stop("SCORPIUS did not return pseudotime in traj$time.")
}

pseudotime_raw <- as.numeric(traj$time)
names(pseudotime_raw) <- rownames(expr_scorpius)

bcell$SCORPIUS_pseudotime_raw <- pseudotime_raw[colnames(bcell)]

############################
## 11. Orient and normalise pseudotime
############################

bcell$SCORPIUS_pseudotime <- bcell$SCORPIUS_pseudotime_raw

root_celltypes_present <- intersect(
  root_celltypes,
  as.character(unique(bcell@meta.data[[celltype_col]]))
)

write_csv(
  tibble(root_celltype_used_for_orientation = root_celltypes_present),
  file.path(tab_dir, "EW2_SCORPIUS_root_celltypes_used_for_orientation.csv")
)

if (length(root_celltypes_present) > 0) {
  root_cells <- rownames(bcell@meta.data)[
    as.character(bcell@meta.data[[celltype_col]]) %in% root_celltypes_present
  ]

  root_median <- median(
    bcell$SCORPIUS_pseudotime_raw[root_cells],
    na.rm = TRUE
  )

  global_median <- median(
    bcell$SCORPIUS_pseudotime_raw,
    na.rm = TRUE
  )

  if (is.finite(root_median) && is.finite(global_median) && root_median > global_median) {
    max_time <- max(bcell$SCORPIUS_pseudotime_raw, na.rm = TRUE)
    min_time <- min(bcell$SCORPIUS_pseudotime_raw, na.rm = TRUE)
    bcell$SCORPIUS_pseudotime <- max_time - bcell$SCORPIUS_pseudotime_raw + min_time
    message(
      "SCORPIUS pseudotime direction was reversed so that early B-cell root states are earlier: ",
      paste(root_celltypes_present, collapse = ", ")
    )
  } else {
    message("SCORPIUS pseudotime direction was kept unchanged.")
  }

} else {
  warning("No specified root_celltypes were found. Pseudotime direction was not manually oriented.")
}

pt <- bcell$SCORPIUS_pseudotime

if (max(pt, na.rm = TRUE) == min(pt, na.rm = TRUE)) {
  stop("SCORPIUS pseudotime has no variation.")
}

bcell$SCORPIUS_pseudotime_01 <- (
  pt - min(pt, na.rm = TRUE)
) / (
  max(pt, na.rm = TRUE) - min(pt, na.rm = TRUE)
)

saveRDS(
  bcell,
  file.path(rds_dir, "EW2_seurat_with_SCORPIUS_pseudotime.rds")
)

saveRDS(
  list(space = space, trajectory = traj),
  file.path(rds_dir, "EW2_SCORPIUS_space_and_trajectory.rds")
)

############################
## 12. Prepare plotting and summary data
############################

umap_mat <- Embeddings(bcell, reduction = "umap")[, 1:2, drop = FALSE]
colnames(umap_mat) <- c("UMAP_1", "UMAP_2")

meta <- bcell@meta.data %>%
  rownames_to_column("cell_id_seurat") %>%
  mutate(
    celltype_plot = factor(
      as.character(.data[[celltype_col]]),
      levels = ordered_present
    )
  )

umap_df <- as.data.frame(umap_mat) %>%
  rownames_to_column("cell_id_seurat") %>%
  left_join(meta, by = "cell_id_seurat") %>%
  filter(!is.na(SCORPIUS_pseudotime_01))

write_csv(
  umap_df,
  file.path(tab_dir, "EW2_SCORPIUS_cell_metadata_with_pseudotime.csv")
)

pt_summary <- umap_df %>%
  group_by(celltype_plot) %>%
  summarise(
    n = n(),
    min_pt = min(SCORPIUS_pseudotime_01, na.rm = TRUE),
    q25_pt = quantile(SCORPIUS_pseudotime_01, 0.25, na.rm = TRUE),
    median_pt = median(SCORPIUS_pseudotime_01, na.rm = TRUE),
    mean_pt = mean(SCORPIUS_pseudotime_01, na.rm = TRUE),
    q75_pt = quantile(SCORPIUS_pseudotime_01, 0.75, na.rm = TRUE),
    max_pt = max(SCORPIUS_pseudotime_01, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_pt)

write_csv(
  pt_summary,
  file.path(tab_dir, "EW2_SCORPIUS_pseudotime_summary_by_celltype.csv")
)

label_df <- umap_df %>%
  group_by(celltype_plot) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(celltype_plot))

write_csv(
  label_df,
  file.path(tab_dir, "EW2_SCORPIUS_celltype_label_positions.csv")
)

## This is only a visual summary of SCORPIUS pseudotime over UMAP.
## The actual SCORPIUS trajectory is inferred in SCORPIUS reduced space.
umap_path_df <- umap_df %>%
  mutate(
    pt_bin = cut(
      SCORPIUS_pseudotime_01,
      breaks = seq(0, 1, length.out = n_pt_bins + 1),
      include.lowest = TRUE,
      labels = FALSE
    )
  ) %>%
  group_by(pt_bin) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    mean_pt = mean(SCORPIUS_pseudotime_01, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(pt_bin), n >= 3) %>%
  arrange(mean_pt)

write_csv(
  umap_path_df,
  file.path(tab_dir, "EW2_SCORPIUS_umap_pseudotime_trend_line.csv")
)

############################
## 13. Plot: UMAP coloured by SCORPIUS pseudotime
############################

p_umap_pt <- ggplot(
  umap_df,
  aes(x = UMAP_1, y = UMAP_2, color = SCORPIUS_pseudotime_01)
) +
  geom_point(size = cell_point_size, alpha = cell_alpha) +
  geom_path(
    data = umap_path_df,
    aes(x = UMAP_1, y = UMAP_2),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.1,
    alpha = 0.85,
    arrow = arrow(length = unit(0.16, "cm"), type = "closed")
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = UMAP_1, y = UMAP_2, label = celltype_plot),
    inherit.aes = FALSE,
    color = "black",
    size = label_size,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_color_viridis_c(
    option = "cividis",
    na.value = "grey85",
    name = "SCORPIUS\npseudotime"
  ) +
  coord_equal() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    title = "SCORPIUS pseudotime on expression-only UMAP",
    subtitle = "Black line shows binned UMAP trend for visualisation only",
    x = "UMAP 1",
    y = "UMAP 2"
  )

save_plot(
  p_umap_pt,
  "01_SCORPIUS_pseudotime_on_expression_only_UMAP.png"
)

ggsave(
  file.path(fig_dir, "01_SCORPIUS_pseudotime_on_expression_only_UMAP.pdf"),
  p_umap_pt,
  width = plot_width,
  height = plot_height,
  bg = "white"
)

############################
## 14. Plot: UMAP coloured by cell type
############################

p_umap_celltype <- ggplot(
  umap_df,
  aes(x = UMAP_1, y = UMAP_2, color = celltype_plot)
) +
  geom_point(size = cell_point_size, alpha = cell_alpha) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = UMAP_1, y = UMAP_2, label = celltype_plot),
    inherit.aes = FALSE,
    color = "black",
    size = label_size,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_color_manual(
    values = celltype_cols_present,
    name = "Cell type",
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 4, alpha = 1)
    )
  ) +
  coord_equal() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    title = "Annotated B-cell states on expression-only UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
  )

save_plot(
  p_umap_celltype,
  "02_SCORPIUS_celltype_on_expression_only_UMAP.png"
)

ggsave(
  file.path(fig_dir, "02_SCORPIUS_celltype_on_expression_only_UMAP.pdf"),
  p_umap_celltype,
  width = plot_width,
  height = plot_height,
  bg = "white"
)

############################
## 15. Plot: UMAP cell type with SCORPIUS trend line
############################

p_umap_celltype_path <- p_umap_celltype +
  geom_path(
    data = umap_path_df,
    aes(x = UMAP_1, y = UMAP_2),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.1,
    alpha = 0.9,
    arrow = arrow(length = unit(0.16, "cm"), type = "closed")
  ) +
  labs(
    title = "SCORPIUS trajectory trend over annotated B-cell states",
    subtitle = "Black line shows binned UMAP trend for visualisation only"
  )

save_plot(
  p_umap_celltype_path,
  "03_SCORPIUS_celltype_with_pseudotime_trend_on_expression_only_UMAP.png"
)

ggsave(
  file.path(fig_dir, "03_SCORPIUS_celltype_with_pseudotime_trend_on_expression_only_UMAP.pdf"),
  p_umap_celltype_path,
  width = plot_width,
  height = plot_height,
  bg = "white"
)

############################
## 16. Plot: SCORPIUS reduced space with inferred path
############################

space_df <- as.data.frame(space)

if (ncol(space_df) >= 2) {
  colnames(space_df)[1:2] <- c("SCORPIUS_dim1", "SCORPIUS_dim2")
  space_df <- space_df %>%
    rownames_to_column("cell_id_seurat") %>%
    left_join(
      meta %>% select(cell_id_seurat, celltype_plot, SCORPIUS_pseudotime_01),
      by = "cell_id_seurat"
    )

  path_df <- as.data.frame(traj$path)

  if (!is.null(path_df) && ncol(path_df) >= 2) {
    colnames(path_df)[1:2] <- c("SCORPIUS_dim1", "SCORPIUS_dim2")

    p_scorpius_space <- ggplot(
      space_df,
      aes(x = SCORPIUS_dim1, y = SCORPIUS_dim2, color = celltype_plot)
    ) +
      geom_point(size = cell_point_size, alpha = cell_alpha) +
      geom_path(
        data = path_df,
        aes(x = SCORPIUS_dim1, y = SCORPIUS_dim2),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 1.2,
        alpha = 0.9,
        arrow = arrow(length = unit(0.16, "cm"), type = "closed")
      ) +
      scale_color_manual(
        values = celltype_cols_present,
        name = "Cell type",
        drop = FALSE
      ) +
      guides(
        color = guide_legend(
          override.aes = list(size = 4, alpha = 1)
        )
      ) +
      coord_equal() +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
      ) +
      labs(
        title = "SCORPIUS inferred trajectory in reduced space",
        x = "SCORPIUS dimension 1",
        y = "SCORPIUS dimension 2"
      )

    save_plot(
      p_scorpius_space,
      "04_SCORPIUS_reduced_space_trajectory_celltype.png"
    )

    ggsave(
      file.path(fig_dir, "04_SCORPIUS_reduced_space_trajectory_celltype.pdf"),
      p_scorpius_space,
      width = plot_width,
      height = plot_height,
      bg = "white"
    )

    write_csv(
      space_df,
      file.path(tab_dir, "EW2_SCORPIUS_reduced_space_coordinates.csv")
    )

    write_csv(
      path_df,
      file.path(tab_dir, "EW2_SCORPIUS_path_coordinates.csv")
    )
  }
}

############################
## 17. Plot: cell type pseudotime distribution
############################

p_celltype_pt <- ggplot(
  umap_df,
  aes(x = celltype_plot, y = SCORPIUS_pseudotime_01, fill = celltype_plot)
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.8,
    color = "grey30",
    linewidth = 0.3
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.size = 0.2,
    alpha = 0.75,
    color = "grey20"
  ) +
  scale_fill_manual(
    values = celltype_cols_present,
    name = "Cell type",
    drop = FALSE
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  labs(
    x = "Annotated B-cell state",
    y = "SCORPIUS pseudotime",
    title = "SCORPIUS pseudotime distribution across annotated B-cell states"
  )

save_plot(
  p_celltype_pt,
  "05_celltype_vs_SCORPIUS_pseudotime.png",
  width = 9,
  height = 5.5
)

ggsave(
  file.path(fig_dir, "05_celltype_vs_SCORPIUS_pseudotime.pdf"),
  p_celltype_pt,
  width = 9,
  height = 5.5,
  bg = "white"
)

############################
## 18. Combined figure
############################

combined_plot <- p_umap_celltype_path | p_umap_pt

ggsave(
  file.path(fig_dir, "06_combined_SCORPIUS_celltype_pseudotime.pdf"),
  combined_plot,
  width = 16,
  height = 6.5,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "06_combined_SCORPIUS_celltype_pseudotime.png"),
  combined_plot,
  width = 16,
  height = 6.5,
  dpi = plot_dpi,
  bg = "white"
)

############################
## 19. Save final metadata and objects
############################

write_csv(
  bcell@meta.data %>% rownames_to_column("cell_id_seurat"),
  file.path(tab_dir, "EW2_final_metadata_with_SCORPIUS_pseudotime.csv")
)

saveRDS(
  bcell,
  file.path(rds_dir, "EW2_SCORPIUS_FINAL.rds")
)

cat("SCORPIUS expression-only exploratory trajectory reconstruction completed.\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", tab_dir, "\n")
cat("RDS files saved to:", rds_dir, "\n")

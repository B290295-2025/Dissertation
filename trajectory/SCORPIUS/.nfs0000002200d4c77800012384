############################################################
## EW2 SCORPIUS trajectory reconstruction
## Expression-only exploratory trajectory analysis
############################################################

############################
## 0. Setup
############################

setwd("/localdisk/home/s2845297/dissertation/trajectory/SCORPIUS")

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

outdir  <- "SCORPIUS_final"
fig_dir <- file.path(outdir, "figures")
tab_dir <- file.path(outdir, "tables")
rds_dir <- file.path(outdir, "rds")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

celltype_col <- "celltype"
assay_use <- "RNA"
bcell_types <- NULL
n_hvg <- 2000
exclude_receptor_genes <- TRUE
root_celltype <- "Transitional"
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

write.csv(package_versions, file.path(outdir, "SCORPIUS_package_versions.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_SCORPIUS.txt"))

############################
## 3. Colour-blind friendly palette
############################

celltype_order <- c(
  "Transitional", "NAIVE 1", "NAIVE 2", "Early-Activation",
  "IgD+ve Memory", "IgD-ve Memory", "HSP+", "DN2",
  "ASC 1", "ASC 2", "ASC 3", "Cell-Cycling"
)

celltype_cols <- c(
  "Transitional"     = "#332",  
  "NAIVE 1"          = "#28CCEE", 
  "NAIVE 2"          = "#149", 
  "Early-Activation" = "#117733", 
  "IgD+ve Memory"    = "#999933", 
  "IgD-ve Memory"    = "#D06677", 
  "HSP+"             = "#E71133",  
  "DN2"              = "blue",  
  "ASC 1"            = "#882255",  
  "ASC 2"            = "#AAAAA9", 
  "ASC 3"            = "#CC3311",  
  "Cell-Cycling"     = "#661176"
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

get_assay_data_safe <- function(obj, assay, layer_or_slot = "data") {
  mat <- tryCatch(
    GetAssayData(obj, assay = assay, layer = layer_or_slot),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer_or_slot)
  )
  mat
}

is_receptor_gene_symbol <- function(gene_symbol) {
  receptor_patterns <- c(
    "^IGHV", "^IGHD[0-9]", "^IGHJ", "^IGHC",
    "^IGHA[0-9]*$", "^IGHG[0-9]*$", "^IGHM$", "^IGHD$", "^IGHE$",
    "^IGKV", "^IGKJ", "^IGKC$", "^IGLV", "^IGLJ", "^IGLC[0-9]*$",
    "^TRAV", "^TRAJ", "^TRAC$", "^TRBV", "^TRBD", "^TRBJ", "^TRBC[0-9]*$",
    "^TRGV", "^TRGJ", "^TRGC[0-9]*$", "^TRDV", "^TRDD", "^TRDJ", "^TRDC$"
  )
  grepl(paste(receptor_patterns, collapse = "|"), as.character(gene_symbol))
}

############################
## 5. Load Seurat object
############################

if (!file.exists(seurat_rds)) stop("Seurat RDS file not found: ", seurat_rds)

seu <- readRDS(seurat_rds)
cat("Loaded Seurat object with", ncol(seu), "cells.\n")
cat("Available assays:", paste(SeuratObject::Assays(seu), collapse = ", "), "\n")
cat("Available reductions:", paste(Reductions(seu), collapse = ", "), "\n")

if (!celltype_col %in% colnames(seu@meta.data)) stop("celltype_col not found: ", celltype_col)
if (!"umap" %in% Reductions(seu)) stop("The Seurat object does not contain a UMAP reduction.")

DefaultAssay(seu) <- assay_use
if (inherits(seu[[assay_use]], "Assay5")) seu[[assay_use]] <- JoinLayers(seu[[assay_use]])

if (!is.null(bcell_types)) {
  keep_cells <- rownames(seu@meta.data)[as.character(seu@meta.data[[celltype_col]]) %in% bcell_types]
  if (length(keep_cells) == 0) stop("No cells matched bcell_types.")
  bcell <- subset(seu, cells = keep_cells)
} else {
  bcell <- seu
}

present_celltypes <- unique(as.character(bcell@meta.data[[celltype_col]]))
ordered_present <- c(celltype_order[celltype_order %in% present_celltypes], setdiff(sort(present_celltypes), celltype_order))
bcell@meta.data[[celltype_col]] <- factor(as.character(bcell@meta.data[[celltype_col]]), levels = ordered_present)

celltype_cols_present <- celltype_cols[ordered_present]
missing_cols <- setdiff(ordered_present, names(celltype_cols_present))
if (length(missing_cols) > 0) {
  extra_cols <- viridis::viridis(length(missing_cols), option = "D")
  names(extra_cols) <- missing_cols
  celltype_cols_present <- c(celltype_cols_present, extra_cols)
}

cat("Cells used for SCORPIUS:", ncol(bcell), "\n")
print(table(bcell@meta.data[[celltype_col]], useNA = "ifany"))
saveRDS(bcell, file.path(rds_dir, "EW2_seurat_for_SCORPIUS.rds"))

############################
## 6. Build expression-only SCORPIUS input matrix
############################

expr_mat <- get_assay_data_safe(bcell, assay = assay_use, layer_or_slot = "data")
if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) stop("Expression matrix is empty.")

hvg <- VariableFeatures(bcell)
if (length(hvg) == 0) {
  message("No VariableFeatures detected. Running FindVariableFeatures().")
  bcell <- FindVariableFeatures(bcell, assay = assay_use, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
  hvg <- VariableFeatures(bcell)
}

hvg <- hvg[seq_len(min(length(hvg), n_hvg))]
hvg <- intersect(hvg, rownames(expr_mat))

hvg_table <- tibble(gene = hvg, is_receptor_gene = is_receptor_gene_symbol(hvg))
if (exclude_receptor_genes) {
  receptor_removed <- hvg_table %>% filter(is_receptor_gene)
  hvg <- hvg_table %>% filter(!is_receptor_gene) %>% pull(gene)
  write.csv(receptor_removed, file.path(tab_dir, "EW2_SCORPIUS_receptor_genes_removed.csv"), row.names = FALSE)
}
if (length(hvg) < 100) stop("Too few genes available for SCORPIUS after filtering: ", length(hvg))

write.csv(tibble(gene = hvg), file.path(tab_dir, "EW2_SCORPIUS_expression_only_genes.csv"), row.names = FALSE)
cat("Genes used for SCORPIUS:", length(hvg), "\n")

expr_use <- expr_mat[hvg, , drop = FALSE]
expr_scorpius <- t(as.matrix(expr_use))
gene_sd <- apply(expr_scorpius, 2, sd, na.rm = TRUE)
expr_scorpius <- expr_scorpius[, gene_sd > 0, drop = FALSE]
expr_scorpius <- scale(expr_scorpius)
expr_scorpius[is.na(expr_scorpius)] <- 0
expr_scorpius[is.infinite(expr_scorpius)] <- 0

cat("SCORPIUS input matrix dimensions:\n")
print(dim(expr_scorpius))
saveRDS(expr_scorpius, file.path(rds_dir, "EW2_SCORPIUS_input_expression_matrix.rds"))

############################
## 7. Run SCORPIUS
############################

space <- SCORPIUS::reduce_dimensionality(expr_scorpius, "spearman")
traj <- SCORPIUS::infer_trajectory(space)
if (is.null(traj$time)) stop("SCORPIUS did not return pseudotime in traj$time.")

pseudotime_raw <- as.numeric(traj$time)
names(pseudotime_raw) <- rownames(expr_scorpius)
bcell$SCORPIUS_pseudotime_raw <- pseudotime_raw[colnames(bcell)]

############################
## 8. Orient and normalise pseudotime
############################

bcell$SCORPIUS_pseudotime <- bcell$SCORPIUS_pseudotime_raw
if (!is.null(root_celltype) && root_celltype %in% as.character(bcell@meta.data[[celltype_col]])) {
  root_cells <- rownames(bcell@meta.data)[as.character(bcell@meta.data[[celltype_col]]) == root_celltype]
  root_median <- median(bcell$SCORPIUS_pseudotime_raw[root_cells], na.rm = TRUE)
  global_median <- median(bcell$SCORPIUS_pseudotime_raw, na.rm = TRUE)
  if (is.finite(root_median) && is.finite(global_median) && root_median > global_median) {
    max_time <- max(bcell$SCORPIUS_pseudotime_raw, na.rm = TRUE)
    min_time <- min(bcell$SCORPIUS_pseudotime_raw, na.rm = TRUE)
    bcell$SCORPIUS_pseudotime <- max_time - bcell$SCORPIUS_pseudotime_raw + min_time
    message("SCORPIUS pseudotime direction was reversed so that root_celltype is early: ", root_celltype)
  } else {
    message("SCORPIUS pseudotime direction was kept unchanged.")
  }
} else {
  warning("root_celltype not found. Pseudotime direction was not manually oriented.")
}

pt <- bcell$SCORPIUS_pseudotime
if (max(pt, na.rm = TRUE) == min(pt, na.rm = TRUE)) stop("SCORPIUS pseudotime has no variation.")
bcell$SCORPIUS_pseudotime_01 <- (pt - min(pt, na.rm = TRUE)) / (max(pt, na.rm = TRUE) - min(pt, na.rm = TRUE))

saveRDS(bcell, file.path(rds_dir, "EW2_seurat_with_SCORPIUS_pseudotime.rds"))
saveRDS(list(space = space, trajectory = traj), file.path(rds_dir, "EW2_SCORPIUS_space_and_trajectory.rds"))

############################
## 9. Prepare plotting data
############################

umap_mat <- Embeddings(bcell, reduction = "umap")[, 1:2, drop = FALSE]
colnames(umap_mat) <- c("UMAP_1", "UMAP_2")

meta <- bcell@meta.data %>%
  rownames_to_column("cell_id_seurat") %>%
  mutate(celltype_plot = factor(as.character(.data[[celltype_col]]), levels = ordered_present))

umap_df <- as.data.frame(umap_mat) %>%
  rownames_to_column("cell_id_seurat") %>%
  left_join(meta, by = "cell_id_seurat") %>%
  filter(!is.na(SCORPIUS_pseudotime_01))

write.csv(umap_df, file.path(tab_dir, "EW2_SCORPIUS_cell_metadata_with_pseudotime.csv"), row.names = FALSE)

label_df <- umap_df %>%
  group_by(celltype_plot) %>%
  summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE), UMAP_2 = median(UMAP_2, na.rm = TRUE), n = n(), .groups = "drop") %>%
  filter(!is.na(celltype_plot))

write.csv(label_df, file.path(tab_dir, "EW2_SCORPIUS_celltype_label_positions.csv"), row.names = FALSE)

umap_path_df <- umap_df %>%
  mutate(pt_bin = cut(SCORPIUS_pseudotime_01, breaks = seq(0, 1, length.out = n_pt_bins + 1), include.lowest = TRUE, labels = FALSE)) %>%
  group_by(pt_bin) %>%
  summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE), UMAP_2 = median(UMAP_2, na.rm = TRUE), mean_pt = mean(SCORPIUS_pseudotime_01, na.rm = TRUE), n = n(), .groups = "drop") %>%
  filter(!is.na(pt_bin), n >= 3) %>%
  arrange(mean_pt)

write.csv(umap_path_df, file.path(tab_dir, "EW2_SCORPIUS_umap_pseudotime_trend_line.csv"), row.names = FALSE)

############################
## 10. Plot: UMAP coloured by SCORPIUS pseudotime
############################

p_umap_pt <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = SCORPIUS_pseudotime_01)) +
  geom_point(size = cell_point_size, alpha = cell_alpha) +
  geom_path(data = umap_path_df, aes(x = UMAP_1, y = UMAP_2), inherit.aes = FALSE, color = "black", linewidth = 1.1, alpha = 0.85, arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
  ggrepel::geom_text_repel(data = label_df, aes(x = UMAP_1, y = UMAP_2, label = celltype_plot), inherit.aes = FALSE, color = "black", size = label_size, fontface = "bold", box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_viridis_c(option = "cividis", na.value = "grey85", name = "SCORPIUS\npseudotime") +
  coord_equal() +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right") +
  labs(title = "SCORPIUS pseudotime on Seurat UMAP", subtitle = "Black line shows binned UMAP trend for visualisation only", x = "UMAP 1", y = "UMAP 2")

save_plot(p_umap_pt, "01_SCORPIUS_pseudotime_on_UMAP.png")
ggsave(file.path(fig_dir, "01_SCORPIUS_pseudotime_on_UMAP.pdf"), p_umap_pt, width = plot_width, height = plot_height, bg = "white")

############################
## 11. Plot: UMAP coloured by cell type
############################

p_umap_celltype <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltype_plot)) +
  geom_point(size = cell_point_size, alpha = cell_alpha) +
  ggrepel::geom_text_repel(data = label_df, aes(x = UMAP_1, y = UMAP_2, label = celltype_plot), inherit.aes = FALSE, color = "black", size = label_size, fontface = "bold", box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_manual(values = celltype_cols_present, name = "Cell type", drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  coord_equal() +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right") +
  labs(title = "Annotated B-cell states on Seurat UMAP", x = "UMAP 1", y = "UMAP 2")

save_plot(p_umap_celltype, "02_SCORPIUS_celltype_on_UMAP.png")
ggsave(file.path(fig_dir, "02_SCORPIUS_celltype_on_UMAP.pdf"), p_umap_celltype, width = plot_width, height = plot_height, bg = "white")

############################
## 12. Plot: UMAP cell type with SCORPIUS trend line
############################

p_umap_celltype_path <- p_umap_celltype +
  geom_path(data = umap_path_df, aes(x = UMAP_1, y = UMAP_2), inherit.aes = FALSE, color = "black", linewidth = 1.1, alpha = 0.9, arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
  labs(title = "SCORPIUS trajectory trend over annotated B-cell states", subtitle = "Black line shows binned UMAP trend for visualisation only")

save_plot(p_umap_celltype_path, "03_SCORPIUS_celltype_with_pseudotime_trend_on_UMAP.png")
ggsave(file.path(fig_dir, "03_SCORPIUS_celltype_with_pseudotime_trend_on_UMAP.pdf"), p_umap_celltype_path, width = plot_width, height = plot_height, bg = "white")

############################
## 13. Plot: SCORPIUS reduced space with inferred path
############################

space_df <- as.data.frame(space)
if (ncol(space_df) >= 2) {
  colnames(space_df)[1:2] <- c("SCORPIUS_dim1", "SCORPIUS_dim2")
  space_df <- space_df %>%
    rownames_to_column("cell_id_seurat") %>%
    left_join(meta %>% select(cell_id_seurat, celltype_plot, SCORPIUS_pseudotime_01), by = "cell_id_seurat")

  path_df <- as.data.frame(traj$path)
  if (!is.null(path_df) && ncol(path_df) >= 2) {
    colnames(path_df)[1:2] <- c("SCORPIUS_dim1", "SCORPIUS_dim2")

    p_scorpius_space <- ggplot(space_df, aes(x = SCORPIUS_dim1, y = SCORPIUS_dim2, color = celltype_plot)) +
      geom_point(size = cell_point_size, alpha = cell_alpha) +
      geom_path(data = path_df, aes(x = SCORPIUS_dim1, y = SCORPIUS_dim2), inherit.aes = FALSE, color = "black", linewidth = 1.2, alpha = 0.9, arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
      scale_color_manual(values = celltype_cols_present, name = "Cell type", drop = FALSE) +
      guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
      coord_equal() +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right") +
      labs(title = "SCORPIUS inferred trajectory in reduced space", x = "SCORPIUS dimension 1", y = "SCORPIUS dimension 2")

    save_plot(p_scorpius_space, "04_SCORPIUS_reduced_space_trajectory_celltype.png")
    ggsave(file.path(fig_dir, "04_SCORPIUS_reduced_space_trajectory_celltype.pdf"), p_scorpius_space, width = plot_width, height = plot_height, bg = "white")

    write.csv(space_df, file.path(tab_dir, "EW2_SCORPIUS_reduced_space_coordinates.csv"), row.names = FALSE)
    write.csv(path_df, file.path(tab_dir, "EW2_SCORPIUS_path_coordinates.csv"), row.names = FALSE)
  }
}

############################
## 14. Plot: cell type pseudotime distribution
############################

p_celltype_pt <- ggplot(umap_df, aes(x = celltype_plot, y = SCORPIUS_pseudotime_01, fill = celltype_plot)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8, color = "grey30", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.75, color = "grey20") +
  scale_fill_manual(values = celltype_cols_present, name = "Cell type", drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(color = "black"), axis.title = element_text(color = "black"), plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none") +
  labs(x = "Annotated B-cell state", y = "SCORPIUS pseudotime", title = "SCORPIUS pseudotime distribution across annotated B-cell states")

save_plot(p_celltype_pt, "05_celltype_vs_SCORPIUS_pseudotime.png", width = 9, height = 5.5)
ggsave(file.path(fig_dir, "05_celltype_vs_SCORPIUS_pseudotime.pdf"), p_celltype_pt, width = 9, height = 5.5, bg = "white")

############################
## 15. Combined figure
############################

combined_plot <- p_umap_celltype_path | p_umap_pt

ggsave(file.path(fig_dir, "06_combined_SCORPIUS_celltype_pseudotime.pdf"), combined_plot, width = 16, height = 6.5, bg = "white")
ggsave(file.path(fig_dir, "06_combined_SCORPIUS_celltype_pseudotime.png"), combined_plot, width = 16, height = 6.5, dpi = plot_dpi, bg = "white")

############################
## 16. Save final metadata and objects
############################

write.csv(bcell@meta.data, file.path(tab_dir, "EW2_final_metadata_with_SCORPIUS_pseudotime.csv"), row.names = TRUE)
saveRDS(bcell, file.path(rds_dir, "EW2_SCORPIUS_FINAL.rds"))

cat("SCORPIUS expression-only exploratory trajectory reconstruction completed.\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", tab_dir, "\n")
cat("RDS files saved to:", rds_dir, "\n")

############################################################
## EW2 Seurat analysis
## Project: Evaluating single-cell trajectories using BCR lineages
############################################################

############################
## 0. Setup
############################

setwd("/localdisk/home/s2845297/dissertation/seurat")

library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clustree)
library(tidyr)
library(scales)

## Record package versions for reproducibility
writeLines(capture.output(sessionInfo()), "sessionInfo_Seurat_analysis.txt")


############################
## 1. Read Cell Ranger outputs
############################

data_dir_ew1 <- "/localdisk/home/s2845297/dissertation/cellranger_EW1/EW1_output/outs/per_sample_outs/EW1_output/count/sample_filtered_feature_bc_matrix.h5"
data_dir_ew2 <- "/localdisk/home/s2845297/dissertation/cellranger_EW2/EW2_output/outs/per_sample_outs/EW2_output/count/sample_filtered_feature_bc_matrix.h5"
data_dir_ew3 <- "/localdisk/home/s2845297/dissertation/cellranger_EW3/EW3_output/outs/per_sample_outs/EW3_output/count/sample_filtered_feature_bc_matrix.h5"

counts_ew1 <- Read10X_h5(data_dir_ew1)
counts_ew2 <- Read10X_h5(data_dir_ew2)
counts_ew3 <- Read10X_h5(data_dir_ew3)

ew1 <- CreateSeuratObject(
  counts = counts_ew1,
  project = "EW-1",
  min.cells = 3,
  min.features = 200
)

ew2_all <- CreateSeuratObject(
  counts = counts_ew2,
  project = "EW-2",
  min.cells = 3,
  min.features = 200
)

ew3 <- CreateSeuratObject(
  counts = counts_ew3,
  project = "EW-3",
  min.cells = 3,
  min.features = 200
)


############################
## 2. Merge three samples for pre-QC overview
############################

merged_obj <- merge(
  ew1,
  y = c(ew2_all, ew3),
  add.cell.ids = c("EW-1", "EW-2", "EW-3"),
  project = "UoE_BCR_Trajectory"
)

merged_obj[["percent.mt"]] <- PercentageFeatureSet(
  merged_obj,
  pattern = "^MT-"
)

table(merged_obj$orig.ident)

sample_cols <- c(
  "EW-1" = "#E69F00",
  "EW-2" = "#56B4E9",
  "EW-3" = "#009E73"
)

p_gene <- VlnPlot(
  merged_obj,
  features = "nFeature_RNA",
  group.by = "orig.ident",
  pt.size = 0.1,
  cols = sample_cols
) +
  theme(legend.position = "none") +
  labs(title = "Number of Features", x = "", y = "nFeature_RNA")

p_count <- VlnPlot(
  merged_obj,
  features = "nCount_RNA",
  group.by = "orig.ident",
  pt.size = 0.1,
  cols = sample_cols
) +
  theme(legend.position = "none") +
  labs(title = "Total UMI Counts per Cell", x = "", y = "nCount_RNA")

p_mt <- VlnPlot(
  merged_obj,
  features = "percent.mt",
  group.by = "orig.ident",
  pt.size = 0.1,
  cols = sample_cols
) +
  theme(legend.position = "none") +
  labs(title = "Mitochondrial Transcript %", x = "", y = "percent.mt (%)")

cell_counts <- as.data.frame(table(merged_obj$orig.ident))
colnames(cell_counts) <- c("Sample", "CellCount")

p_cells <- ggplot(cell_counts, aes(x = Sample, y = CellCount, fill = Sample)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = CellCount), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = sample_cols) +
  theme_classic() +
  labs(title = "Total Cell Count", x = "", y = "Number of Cells") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

qc_combined_plot <- (p_cells | p_gene) / (p_count | p_mt) +
  plot_annotation(
    title = "Pre-QC Sample Summary and Metrics",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

ggsave(
  filename = "PreQC_Sample_Summary.png",
  plot = qc_combined_plot,
  width = 10,
  height = 8,
  dpi = 300
)


############################
## 3. Select EW2 for downstream analysis
############################

ew2 <- ew2_all

ew2[["percent.mt"]] <- PercentageFeatureSet(
  ew2,
  pattern = "^MT-"
)

VlnPlot(
  ew2,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

ew2 <- subset(
  ew2,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 5000 &
    percent.mt < 10
)


############################
## 4. Normalisation, feature selection and scaling
############################

ew2 <- NormalizeData(ew2)

ew2 <- FindVariableFeatures(
  ew2,
  selection.method = "vst",
  nfeatures = 3000
)

VariableFeaturePlot(ew2)

ew2 <- ScaleData(ew2)


############################
## 5. PCA
############################

ew2 <- RunPCA(
  ew2,
  features = VariableFeatures(ew2)
)

ElbowPlot(ew2)


############################
## 6. Resolution selection using clustree
############################

ew2 <- FindNeighbors(
  ew2,
  dims = 1:20
)

ew2 <- FindClusters(
  ew2,
  resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2)
)

clustree(
  ew2@meta.data,
  prefix = "RNA_snn_res."
)

ggsave(
  filename = "EW2_clustree_resolution_selection.png",
  width = 10,
  height = 10,
  dpi = 300
)


############################
## 7. Final clustering
############################

ew2 <- FindNeighbors(
  ew2,
  dims = 1:20
)

ew2 <- FindClusters(
  ew2,
  resolution = 0.6
)

table(Idents(ew2))


############################
## 8. UMAP visualisation
############################

ew2 <- RunUMAP(
  ew2,
  dims = 1:20
)

DimPlot(
  ew2,
  reduction = "umap",
  label = TRUE
)

ggsave(
  filename = "EW2_UMAP_clusters.png",
  width = 8,
  height = 6,
  dpi = 300
)


############################
## 9. Marker gene visualisation
############################

PanB <- FeaturePlot(
  ew2,
  features = c("CD79A", "CD79B", "MS4A1", "CD19", "CD74", "HLA-DRA"),
  reduction = "umap",
  ncol = 3,
  label = TRUE
)

Diff_B <- FeaturePlot(
  ew2,
  features = c("FCER2", "IL4R", "AICDA", "MKI67", "MYC", "XBP1", "SDC1", "PRDM1"),
  reduction = "umap",
  ncol = 3,
  label = TRUE
)

IG_region <- FeaturePlot(
  ew2,
  features = c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHA1", "IGKC", "IGLC2"),
  reduction = "umap",
  ncol = 3,
  label = TRUE
)

T_cell <- FeaturePlot(
  ew2,
  features = c("CD3D", "CD3E", "TRBC1"),
  reduction = "umap",
  label = TRUE
)

NK <- FeaturePlot(
  ew2,
  features = c("NKG7", "GNLY", "KLRD1"),
  reduction = "umap",
  label = TRUE
)

Mono <- FeaturePlot(
  ew2,
  features = c("LYZ", "S100A8", "FCN1"),
  reduction = "umap",
  label = TRUE
)

plasma <- FeaturePlot(
  ew2,
  features = c("XBP1", "SDC1", "MZB1", "JCHAIN"),
  reduction = "umap",
  label = TRUE
)


############################
## 10. Differential marker analysis
############################

markers <- FindAllMarkers(
  ew2,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers,
  "EW2_markers.csv",
  row.names = FALSE
)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

write.csv(
  top10,
  "EW2_top10_markers.csv",
  row.names = FALSE
)


############################
## 11. Manual annotation
############################

new.cluster.ids <- c(
  "IgD-ve Memory",
  "IgD+ve Memory",
  "NAIVE 1",
  "NAIVE 2",
  "ASC 2",
  "DN2",
  "HSP+",
  "Cell-Cycling",
  "Early-Activation",
  "ASC 3",
  "Transitional",
  "ASC 1"
)

names(new.cluster.ids) <- levels(ew2)
ew2 <- RenameIdents(ew2, new.cluster.ids)

ew2[["celltype"]] <- Idents(ew2)

celltype_order <- c(
  "Transitional",
  "NAIVE 1",
  "NAIVE 2",
  "Early-Activation",
  "IgD-ve Memory",
  "IgD+ve Memory",
  "HSP+",
  "DN2",
  "ASC 1",
  "ASC 2",
  "ASC 3",
  "Cell-Cycling"
)

ew2$celltype <- factor(
  ew2$celltype,
  levels = celltype_order
)

Idents(ew2) <- "celltype"

table(ew2$celltype)


############################
## 12. Dot plot for manual annotation markers
############################

manually_markers <- c(
  "CD19", "MS4A1", "MME", "CD83", "BACH2", "TCL1A", "FCER2",
  "IGHD", "CD27", "CD24", "LY9", "CD69", "CCR7", "IRF4",
  "EGR1", "EGR3", "CR2", "HSPA1A", "HSPA6", "HSPA1B",
  "SOX5", "ITGAX", "FCRL5", "ZEB2", "DUSP4", "CD86",
  "CD80", "CD38", "PRDM1", "IGHG1", "IGHG2", "IGHG3",
  "IGHG4", "IGHM", "IGHA1", "JCHAIN", "PPIB", "NEAT1",
  "MKI67", "NUSAP1"
)

dotplot_manual_markers <- DotPlot(
  ew2,
  features = manually_markers
) +
  scale_color_gradient2(
    low = "blue",
    mid = "#grey",
    high = "#darkred",
    midpoint = 0,
    limits = c(-1, 2),
    oob = scales::squish
  ) +
  RotatedAxis()

ggsave(
  filename = "EW2_manual_marker_dotplot.png",
  plot = dotplot_manual_markers,
  width = 16,
  height = 8,
  dpi = 300
)


############################
## 13. Annotated UMAP
############################

celltype_cols <- c(
  "Transitional"     = "#E69F00",
  "NAIVE 1"          = "#56B4E9",
  "NAIVE 2"          = "#009E73",
  "IgD-ve Memory"    = "#F0E442",
  "IgD+ve Memory"    = "#0072B2",
  "HSP+"             = "#D55E00",
  "Early-Activation" = "#999933",
  "DN2"              = "#999999",
  "ASC 1"            = "#882255",
  "ASC 2"            = "#AA4499",
  "ASC 3"            = "#6699CC",
  "Cell-Cycling"     = "#117733"
)

annotated_umap <- DimPlot(
  ew2,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  cols = celltype_cols
)

ggsave(
  filename = "EW2_annotated_UMAP.png",
  plot = annotated_umap,
  width = 10,
  height = 7,
  dpi = 300
)


############################
## 14. Violin plot for selected marker genes
############################

plot_df <- FetchData(
  ew2,
  vars = c("celltype", "MS4A1", "JCHAIN", "ITGAX", "CD86")
)

plot_long <- plot_df %>%
  pivot_longer(
    cols = c("MS4A1", "JCHAIN", "ITGAX", "CD86"),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(
    celltype = factor(celltype, levels = rev(celltype_order)),
    Gene = factor(
      Gene,
      levels = c("MS4A1", "JCHAIN", "ITGAX", "CD86"),
      labels = c("MS4A1 (CD20)", "JCHAIN", "ITGAX (CD11c)", "CD86")
    )
  )

marker_violin <- ggplot(
  plot_long,
  aes(x = celltype, y = Expression, fill = celltype)
) +
  geom_violin(scale = "width", trim = TRUE, color = "gray50", linewidth = 0.35) +
  coord_flip() +
  facet_wrap(~ Gene, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = celltype_cols) +
  labs(x = "Identity", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "right"
  )

ggsave(
  filename = "EW2_selected_marker_violin.png",
  plot = marker_violin,
  width = 12,
  height = 6,
  dpi = 300
)


############################
## 15. Save annotated Seurat object
############################

saveRDS(
  ew2,
  file = "EW2_annotated.rds"
)


############################
## 16. Export data for downstream PAGA analysis
############################

seu <- readRDS("EW2_annotated.rds")
DefaultAssay(seu) <- "RNA"

dir.create("paga_export", showWarnings = FALSE)

meta <- seu@meta.data
meta$cell_id <- rownames(meta)

write.csv(
  meta,
  file = "paga_export/EW2_metadata.csv",
  row.names = FALSE,
  quote = TRUE
)

if (!"pca" %in% names(seu@reductions)) {
  stop("Seurat object does not contain PCA reduction. Please run RunPCA first.")
}

pca <- Embeddings(seu, reduction = "pca")
pca <- as.data.frame(pca)
pca$cell_id <- rownames(pca)

write.csv(
  pca,
  file = "paga_export/EW2_pca.csv",
  row.names = FALSE,
  quote = TRUE
)

if ("umap" %in% names(seu@reductions)) {
  umap <- Embeddings(seu, reduction = "umap")
  umap <- as.data.frame(umap)
  umap$cell_id <- rownames(umap)

  write.csv(
    umap,
    file = "paga_export/EW2_umap.csv",
    row.names = FALSE,
    quote = TRUE
  )
}

print("PAGA export files completed.")

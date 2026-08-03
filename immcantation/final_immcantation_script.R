############################################################
## EW2 BCR analysis using Immcantation, SHM and Dowser
## Project: Evaluating single-cell trajectories using BCR lineages
############################################################

############################
## 0. Setup
############################

setwd("/data")

suppressPackageStartupMessages({
  library(airr)
  library(alakazam)
  library(scoper)
  library(shazam)
  library(dowser)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggtree)
  library(ape)
  library(patchwork)
  library(Seurat)
})

options(stringsAsFactors = FALSE)

## Save session information for reproducibility
writeLines(
  capture.output(sessionInfo()),
  "sessionInfo_immcantation_SHM_dowser.txt"
)

############################
## 1. Paths and parameters
############################

bcr_file <- "/data/cellranger_results/BCR_data_sequences_igblast_db-pass.tsv"
gex_rds  <- "/data/EW2_annotated.rds"

outdir      <- "/data/heavy_light"
tree_outdir <- file.path(outdir, "tree")
shm_outdir  <- file.path(outdir, "SHM_plots")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tree_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(shm_outdir, recursive = TRUE, showWarnings = FALSE)

## PHYLIP dnapars executable inside the Immcantation container
dnapars_exec <- "/usr/local/bin/dnapars"

## Clone definition threshold
clone_threshold <- 0.15

## Minimum number of IGH sequences required for Dowser tree construction
minseq_use <- 5

## IMGT human germline reference
references <- readIMGT(dir = "/usr/local/share/germlines/imgt/human/vdj")


############################
## 2. Read AIRR/Change-O BCR table
############################

bcr_data <- airr::read_rearrangement(
  bcr_file,
  aux_types = c(
    "v_germline_length" = "i",
    "d_germline_length" = "i",
    "j_germline_length" = "i",
    "day" = "i"
  )
)

cat("Initial BCR rows:", nrow(bcr_data), "\n")

write.csv(
  as.data.frame(head(bcr_data, 30)),
  file.path(outdir, "01_bcr_data_head.csv"),
  row.names = FALSE
)


############################
## 3. Basic filtering: productive rearrangements and valid chain calls
############################

required_bcr_cols <- c(
  "sequence_id",
  "productive",
  "v_call",
  "j_call",
  "c_call",
  "locus",
  "junction",
  "junction_aa",
  "sequence_alignment",
  "germline_alignment"
)

missing_bcr_cols <- setdiff(required_bcr_cols, colnames(bcr_data))

if (length(missing_bcr_cols) > 0) {
  warning(
    "The BCR table is missing the following expected columns: ",
    paste(missing_bcr_cols, collapse = ", ")
  )
}

## Productive filter compatible with logical or character input
bcr_data <- bcr_data %>%
  mutate(
    productive_bool = productive %in% c(TRUE, "T", "TRUE", "True", "true", "1")
  ) %>%
  filter(productive_bool)

cat("Rows after productive filtering:", nrow(bcr_data), "\n")

## Retain rows with V/J/C calls consistent with immunoglobulin locus
bcr_data <- bcr_data %>%
  filter(
    (grepl("^IGHV", v_call) & grepl("^IGHJ", j_call) & grepl("^IGH[MGADE]", c_call)) |
      (grepl("^IGKV", v_call) & grepl("^IGKJ", j_call) & grepl("^IGKC", c_call)) |
      (grepl("^IGLV", v_call) & grepl("^IGLJ", j_call) & grepl("^IGLC", c_call))
  )

cat("Rows after V/J/C-locus consistency filtering:", nrow(bcr_data), "\n")
cat("Chain distribution after filtering:\n")
print(table(bcr_data$locus))


############################
## 4. Extract cell identifiers
############################

## Remove only the contig suffix, preserving the full cell barcode
bcr_data <- bcr_data %>%
  mutate(
    cell_id = sub("_contig.*", "", sequence_id),
    cell_id_unique = cell_id
  )

cat("Unique cell_id_unique:", n_distinct(bcr_data$cell_id_unique), "\n")

write.csv(
  bcr_data %>%
    select(sequence_id, cell_id, cell_id_unique) %>%
    head(100),
  file.path(outdir, "02_cell_id_check.csv"),
  row.names = FALSE
)


############################
## 5. Retain strictly paired heavy-light cells
############################

multi_heavy <- table(filter(bcr_data, locus == "IGH")$cell_id_unique)
multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]

cat("Cells with multiple IGH chains:", length(multi_heavy_cells), "\n")

bcr_data <- bcr_data %>%
  filter(!cell_id_unique %in% multi_heavy_cells)

cat("Rows after removing cells with multiple IGH chains:", nrow(bcr_data), "\n")

paired_cells <- bcr_data %>%
  group_by(cell_id_unique) %>%
  filter(any(locus == "IGH")) %>%
  filter(any(locus == "IGK") | any(locus == "IGL")) %>%
  ungroup()

cat("Rows from cells with both IGH and light chain:", nrow(paired_cells), "\n")

chain_number <- paired_cells %>%
  group_by(cell_id_unique) %>%
  summarise(
    n_total = n(),
    IGH = sum(locus == "IGH"),
    IGK = sum(locus == "IGK"),
    IGL = sum(locus == "IGL"),
    .groups = "drop"
  )

write.csv(
  chain_number,
  file.path(outdir, "03_paired_chain_number.csv"),
  row.names = FALSE
)

good_cells <- chain_number %>%
  filter(IGH == 1, IGK + IGL == 1)

cat("Cells with exactly 1 IGH and 1 light chain:", nrow(good_cells), "\n")

bcr_data <- paired_cells %>%
  filter(cell_id_unique %in% good_cells$cell_id_unique)

cat("Rows after strict paired heavy-light filtering:", nrow(bcr_data), "\n")
cat("Chain distribution after strict paired filtering:\n")
print(table(bcr_data$locus))


############################
## 6. Add Seurat gene-expression annotation
############################

gex_db <- readRDS(gex_rds)

anno <- data.frame(
  cell_id_unique = Cells(gex_db),
  gex_annotation = as.character(Idents(gex_db))
)

bcr_data <- bcr_data %>%
  select(-any_of("gex_annotation"))

bcr_data <- left_join(
  bcr_data,
  anno,
  by = "cell_id_unique"
)

match_rate <- mean(!is.na(bcr_data$gex_annotation))
cat("BCR-Seurat annotation match rate:", match_rate, "\n")

if (match_rate < 0.5) {
  warning("Low BCR-Seurat cell ID matching rate. Check cell_id_unique and Seurat cell names.")
}

bcr_data <- bcr_data %>%
  filter(!is.na(gex_annotation))

cat("Rows after adding GEX annotation and removing NA:", nrow(bcr_data), "\n")

write.csv(
  as.data.frame(head(bcr_data, 50)),
  file.path(outdir, "04_bcr_data_annotated_head.csv"),
  row.names = FALSE
)


############################
## 7. Clone definition using SCOPer and light-chain resolution
############################

clone_results <- hierarchicalClones(
  bcr_data,
  cell_id = "cell_id_unique",
  threshold = clone_threshold,
  summarize_clones = FALSE
)

cat("Rows after hierarchicalClones:", nrow(clone_results), "\n")

clone_results <- resolveLightChains(clone_results)

cat("Rows after resolveLightChains:", nrow(clone_results), "\n")

write.csv(
  as.data.frame(head(clone_results, 100)),
  file.path(outdir, "05_results_after_resolveLightChains_head.csv"),
  row.names = FALSE
)

saveRDS(
  clone_results,
  file.path(outdir, "results_after_resolveLightChains.rds")
)


############################
## 8. Clone distribution across transcriptomic annotations
############################

cluster_df <- data.frame(
  cell_id_unique = Cells(gex_db),
  cluster = as.character(Idents(gex_db))
)

clone_df <- clone_results %>%
  distinct(cell_id_unique, clone_id)

clone_cluster <- left_join(
  clone_df,
  cluster_df,
  by = "cell_id_unique"
)

clone_spread <- clone_cluster %>%
  group_by(clone_id) %>%
  summarise(
    n_clusters = n_distinct(cluster),
    n_cells = n(),
    .groups = "drop"
  )

write.csv(
  clone_spread,
  file.path(outdir, "06_clone_spread.csv"),
  row.names = FALSE
)

cat("Distribution of number of transcriptomic clusters per clone:\n")
print(table(clone_spread$n_clusters))


############################
## 9. Germline reconstruction
############################

clone_results <- createGermlines(
  clone_results,
  references,
  nproc = 1
)

cat("Rows after createGermlines:", nrow(clone_results), "\n")

saveRDS(
  clone_results,
  file.path(outdir, "results_after_createGermlines.rds")
)


############################
## 10. Heavy-chain SHM calculation
############################

results_heavy <- clone_results %>%
  filter(locus == "IGH")

cat("IGH rows used for SHM calculation:", nrow(results_heavy), "\n")

data_mut <- observedMutations(
  results_heavy,
  sequenceColumn = "sequence_alignment",
  germlineColumn = "germline_alignment_d_mask",
  frequency = TRUE,
  combine = TRUE
)

drop_cols <- intersect(
  c("gex_annotation", "clone_subgroup_id"),
  colnames(data_mut)
)

if (length(drop_cols) > 0) {
  data_mut <- data_mut %>%
    select(-all_of(drop_cols))
}

data_mut2 <- left_join(
  data_mut,
  clone_results %>%
    select(
      sequence_id,
      cell_id_unique,
      clone_id,
      clone_subgroup_id,
      gex_annotation
    ),
  by = "sequence_id"
)

data_mut2 <- data_mut2 %>%
  mutate(
    clone_subgroup_id = as.character(clone_subgroup_id),
    gex_annotation = as.character(gex_annotation),
    gex_annotation = ifelse(
      is.na(gex_annotation) | gex_annotation == "",
      "Unknown",
      gex_annotation
    ),
    mu_freq = as.numeric(mu_freq)
  )

saveRDS(
  data_mut2,
  file.path(outdir, "data_mut2.rds")
)

write.csv(
  data_mut2,
  file.path(outdir, "07_data_mut2.csv"),
  row.names = FALSE
)

shm_summary <- data_mut2 %>%
  group_by(clone_subgroup_id, gex_annotation) %>%
  summarise(
    N = n(),
    mean_SHM = mean(mu_freq, na.rm = TRUE),
    median_SHM = median(mu_freq, na.rm = TRUE),
    sd_SHM = sd(mu_freq, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  shm_summary,
  file.path(outdir, "08_shm_summary.csv"),
  row.names = FALSE
)


############################
## 11. Annotation order and colours
############################

ann_order <- c(
  "NAIVE 2",    
  "NAIVE 1",
  "Transitional", 
  "Cell-Cycling",
  "Early-Activation",
  "HSP+",     
  "IgD-ve Memory",    
  "IgD+ve Memory",  
  "DN2",
  "ASC 1",
  "ASC 2",
  "ASC 3",           
  "Germline"
)

cb_palette <- c(
  "Transitional"     = "#E69F00",
  "NAIVE 1"          = "#56B4E9",
  "NAIVE 2"          = "#009E73",
  "IgD-ve Memory"    = "#F0E442",
  "IgD+ve Memory"    = "#0072B2",
  "HSP+"             = "#D55E00",
  "Early-Activation" = "#999933",
  "DN2"              = "#999999",
  "ASC 1"            = "#117733",
  "ASC 2"            = "#AA4499",
  "ASC 3"            = "#6699CC",
  "Cell-Cycling"     = "#882255",
  "Germline"         = "black"
)

all_ann_in_data <- sort(unique(c(
  as.character(data_mut2$gex_annotation),
  as.character(clone_results$gex_annotation)
)))

extra_ann <- setdiff(all_ann_in_data, ann_order)

if (length(extra_ann) > 0) {
  ann_order <- c(
    setdiff(ann_order, "Germline"),
    extra_ann,
    "Germline"
  )
}

missing_cols <- setdiff(ann_order, names(cb_palette))

if (length(missing_cols) > 0) {
  extra_cols <- setNames(
    grDevices::rainbow(length(missing_cols)),
    missing_cols
  )
  cb_palette <- c(cb_palette, extra_cols)
}

cb_palette <- cb_palette[ann_order]

write.csv(
  data.frame(
    gex_annotation = names(cb_palette),
    colour = as.character(cb_palette)
  ),
  file.path(outdir, "09_annotation_colour_map.csv"),
  row.names = FALSE
)


############################
## 12. Global SHM boxplot
############################

plot_df <- data_mut2 %>%
  filter(!is.na(mu_freq)) %>%
  mutate(
    clone_subgroup_id = as.character(clone_subgroup_id),
    gex_annotation = factor(
      as.character(gex_annotation),
      levels = ann_order
    )
  )

y_max <- as.numeric(quantile(plot_df$mu_freq, 0.99, na.rm = TRUE))

if (!is.finite(y_max) || y_max <= 0) {
  y_max <- max(plot_df$mu_freq, na.rm = TRUE)
}

if (!is.finite(y_max) || y_max <= 0) {
  y_max <- 1
}

plot_df <- plot_df %>%
  mutate(
    mu_freq_plot = pmin(mu_freq, y_max)
  )

shm_ann_order <- setdiff(ann_order, "Germline")
shm_palette <- cb_palette[shm_ann_order]

p_global_shm <- ggplot(
  plot_df,
  aes(
    x = gex_annotation,
    y = mu_freq_plot,
    fill = gex_annotation
  )
) +
  geom_boxplot(
    width = 0.75,
    outlier.shape = 16,
    outlier.size = 0.8,
    alpha = 0.9
  ) +
  scale_fill_manual(
    values = shm_palette,
    drop = FALSE
  ) +
  scale_x_discrete(
    limits = shm_ann_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, y_max),
    expand = c(0.02, 0.02)
  ) +
  labs(
    title = "EW2 global SHM frequency",
    x = "B cell state",
    y = "SHM frequency"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(shm_outdir, "EW2_global_SHM_boxplot_ordered.png"),
  plot = p_global_shm,
  width = 8,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(shm_outdir, "EW2_global_SHM_boxplot_ordered.pdf"),
  plot = p_global_shm,
  width = 8,
  height = 5,
  bg = "white"
)


############################
## 13. Clone-level SHM panel for clones spanning >=3 B-cell states
############################

clone_cluster_summary <- plot_df %>%
  group_by(clone_subgroup_id) %>%
  summarise(
    size = n(),
    n_cluster = n_distinct(gex_annotation),
    clusters = paste(sort(unique(as.character(gex_annotation))), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(size))

write.csv(
  clone_cluster_summary,
  file.path(outdir, "10_clone_cluster_summary_from_SHM.csv"),
  row.names = FALSE
)

candidate_clones <- clone_cluster_summary %>%
  filter(n_cluster >= 3) %>%
  pull(clone_subgroup_id)

cat("Number of clone_subgroup_id spanning >=3 B-cell states:", length(candidate_clones), "\n")

panel_df <- plot_df %>%
  filter(clone_subgroup_id %in% candidate_clones)

clone_order <- panel_df %>%
  count(clone_subgroup_id, name = "size") %>%
  arrange(desc(size)) %>%
  pull(clone_subgroup_id)

panel_df <- panel_df %>%
  mutate(
    clone_subgroup_id = factor(
      clone_subgroup_id,
      levels = clone_order
    ),
    gex_annotation = factor(
      as.character(gex_annotation),
      levels = shm_ann_order
    )
  )

p_clone_panel_shm <- ggplot(
  panel_df,
  aes(
    x = gex_annotation,
    y = mu_freq_plot,
    fill = gex_annotation
  )
) +
  geom_boxplot(
    width = 0.75,
    outlier.shape = 16,
    outlier.size = 0.5,
    alpha = 0.9
  ) +
  facet_wrap(
    ~ clone_subgroup_id,
    scales = "fixed",
    ncol = 4
  ) +
  scale_fill_manual(
    values = shm_palette,
    drop = FALSE
  ) +
  scale_x_discrete(
    limits = shm_ann_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, y_max),
    expand = c(0.02, 0.02)
  ) +
  labs(
    title = "SHM frequency of clone_subgroup_id with n_cluster >= 3",
    x = "B cell state",
    y = "SHM frequency"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(shm_outdir, "EW2_clone_panel_SHM_boxplot_ordered.png"),
  plot = p_clone_panel_shm,
  width = 16,
  height = 12,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(shm_outdir, "EW2_clone_panel_SHM_boxplot_ordered.pdf"),
  plot = p_clone_panel_shm,
  width = 16,
  height = 12,
  bg = "white"
)

p_shm_combined <- p_global_shm / p_clone_panel_shm +
  plot_layout(heights = c(1, 3))

ggsave(
  filename = file.path(shm_outdir, "EW2_SHM_boxplot_combined_ordered.png"),
  plot = p_shm_combined,
  width = 16,
  height = 16,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(shm_outdir, "EW2_SHM_boxplot_combined_ordered.pdf"),
  plot = p_shm_combined,
  width = 16,
  height = 16,
  bg = "white"
)


############################
## 14. Dowser lineage-tree construction using IGH sequences
############################

results_tree <- clone_results %>%
  filter(locus == "IGH") %>%
  mutate(
    sequence_id = as.character(sequence_id),
    clone_subgroup_id = as.character(clone_subgroup_id),
    gex_annotation = as.character(gex_annotation),
    gex_annotation = ifelse(
      is.na(gex_annotation) | gex_annotation == "",
      "Unknown",
      gex_annotation
    )
  )

cat("IGH rows used for Dowser tree construction:", nrow(results_tree), "\n")
cat("Number of clone_subgroup_id used for Dowser:", n_distinct(results_tree$clone_subgroup_id), "\n")

clones <- formatClones(
  results_tree,
  clone = "clone_subgroup_id",
  traits = c("gex_annotation", "clone_subgroup_id"),
  minseq = minseq_use
)

cat("Number of clones after formatClones:", length(clones$clone_id), "\n")

saveRDS(
  clones,
  file.path(outdir, "clones.rds")
)

trees <- getTrees(
  clones,
  build = "dnapars",
  exec = dnapars_exec
)

cat("Number of trees after getTrees:", length(trees$clone_id), "\n")
print(head(trees$clone_id))

saveRDS(
  trees,
  file.path(outdir, "trees.rds")
)

tree_clone_cluster_summary <- results_tree %>%
  group_by(clone_subgroup_id) %>%
  summarise(
    size = n(),
    n_cluster = n_distinct(gex_annotation),
    clusters = paste(sort(unique(gex_annotation)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(size))

write.csv(
  tree_clone_cluster_summary,
  file.path(tree_outdir, "tree_clone_cluster_summary.csv"),
  row.names = FALSE
)

candidate_tree <- tree_clone_cluster_summary %>%
  filter(n_cluster >= 3) %>%
  mutate(clone_subgroup_id = as.character(clone_subgroup_id))

tree_ids <- as.character(trees$clone_id)
clone_ids <- as.character(clones$clone_id)

candidate_check <- candidate_tree %>%
  mutate(
    in_formatClones = clone_subgroup_id %in% clone_ids,
    in_getTrees = clone_subgroup_id %in% tree_ids
  ) %>%
  arrange(in_getTrees, in_formatClones, desc(size))

write.csv(
  candidate_check,
  file.path(tree_outdir, "candidate_tree_check.csv"),
  row.names = FALSE
)

cat("Candidate clones spanning >=3 B-cell states:", nrow(candidate_tree), "\n")
cat("Candidates included by formatClones:", sum(candidate_check$in_formatClones), "\n")
cat("Candidates successfully reconstructed by getTrees:", sum(candidate_check$in_getTrees), "\n")


############################
## 15. Plot Dowser trees for clone_subgroup_id spanning >=3 B-cell states
############################

tree_to_plot <- candidate_check %>%
  filter(in_getTrees) %>%
  pull(clone_subgroup_id) %>%
  unique()

tree_index <- match(tree_to_plot, tree_ids)
tree_index <- tree_index[!is.na(tree_index)]

if (length(tree_index) == 0) {
  warning("No clone with n_cluster >= 3 successfully reconstructed by Dowser. Skipping tree plotting.")
} else {

  cat("Number of Dowser trees to plot:", length(tree_index), "\n")

  mu_tbl <- data_mut2 %>%
    transmute(
      sequence_id = as.character(sequence_id),
      mu_freq = as.numeric(mu_freq)
    ) %>%
    distinct(sequence_id, .keep_all = TRUE)

  mu_source <- as.numeric(mu_tbl$mu_freq)
  mu_source <- mu_source[is.finite(mu_source)]

  if (length(mu_source) == 0) {
    mu_limits <- c(0, 1)
  } else {
    mu_max <- max(mu_source, na.rm = TRUE)
    if (!is.finite(mu_max) || mu_max <= 0) {
      mu_limits <- c(0, 1)
    } else {
      mu_limits <- c(0, mu_max)
    }
  }

  mu_breaks <- pretty(mu_limits, n = 4)
  mu_breaks <- mu_breaks[
    mu_breaks >= mu_limits[1] &
      mu_breaks <= mu_limits[2]
  ]

  if (length(mu_breaks) == 0) {
    mu_breaks <- mu_limits
  }

  get_tree_xmax <- function(tr) {
    if (is.null(tr) || !inherits(tr, "phylo")) {
      return(NA_real_)
    }

    x <- tryCatch({
      if (!is.null(tr$edge.length)) {
        ape::node.depth.edgelength(tr)
      } else {
        ape::node.depth(tr)
      }
    }, error = function(e) {
      NA_real_
    })

    x <- x[is.finite(x)]

    if (length(x) == 0) {
      return(NA_real_)
    }

    max(x, na.rm = TRUE)
  }

  tree_xmax <- vapply(
    trees$trees[tree_index],
    get_tree_xmax,
    numeric(1)
  )

  if (!any(is.finite(tree_xmax))) {
    global_xmax <- 1
  } else {
    global_xmax <- max(tree_xmax, na.rm = TRUE)
  }

  if (!is.finite(global_xmax) || global_xmax <= 0) {
    global_xmax <- 1
  }

  global_xmax <- global_xmax * 1.08
  x_left_pad <- global_xmax * 0.12
  x_right_pad <- global_xmax * 0.03
  x_limits <- c(-x_left_pad, global_xmax + x_right_pad)

  x_breaks <- pretty(c(0, global_xmax), n = 5)
  x_breaks <- x_breaks[
    x_breaks >= 0 &
      x_breaks <= global_xmax
  ]

  if (length(x_breaks) < 2) {
    x_breaks <- c(0, global_xmax)
  }

  x_break_diffs <- diff(x_breaks)
  x_break_diffs <- x_break_diffs[x_break_diffs > 0]

  if (length(x_break_diffs) == 0) {
    scale_bar_width <- global_xmax / 5
  } else {
    scale_bar_width <- x_break_diffs[1]
  }

  scale_bar_x <- x_limits[1] + x_left_pad * 0.15
  scale_bar_y <- 0

  plot_log <- vector("list", length(tree_index))

  for (k in seq_along(tree_index)) {

    i <- tree_index[k]

    tree_id <- as.character(trees$clone_id[i])
    tree <- trees$trees[[i]]

    message("Plotting Dowser tree ", k, "/", length(tree_index), ": ", tree_id)

    res <- tryCatch({

      clone_i <- match(tree_id, clone_ids)

      if (is.na(clone_i)) {
        stop("Clone not found in formatClones output: ", tree_id)
      }

      clone <- clones$data[[clone_i]]
      clone_dat <- as_tibble(clone@data)

      if (!"sequence_id" %in% colnames(clone_dat)) {
        stop("clone@data does not contain sequence_id: ", tree_id)
      }

      if (!"gex_annotation" %in% colnames(clone_dat)) {
        clone_dat$gex_annotation <- "Unknown"
      }

      if (!"clone_subgroup_id" %in% colnames(clone_dat)) {
        clone_dat$clone_subgroup_id <- tree_id
      }

      clone_dat <- clone_dat %>%
        mutate(
          sequence_id = as.character(sequence_id),
          clone_subgroup_id = as.character(clone_subgroup_id),
          gex_annotation = as.character(gex_annotation),
          gex_annotation = ifelse(
            is.na(gex_annotation) | gex_annotation == "",
            "Unknown",
            gex_annotation
          )
        )

      if ("mu_freq" %in% colnames(clone_dat)) {
        ann <- clone_dat %>%
          select(sequence_id, gex_annotation, clone_subgroup_id, mu_freq)
      } else {
        ann <- clone_dat %>%
          left_join(mu_tbl, by = "sequence_id") %>%
          select(sequence_id, gex_annotation, clone_subgroup_id, mu_freq)
      }

      ann <- ann %>%
        mutate(
          label = as.character(sequence_id),
          gex_annotation = as.character(gex_annotation),
          gex_annotation = ifelse(
            is.na(gex_annotation) | gex_annotation == "",
            "Unknown",
            gex_annotation
          ),
          mu_freq = as.numeric(mu_freq),
          mu_freq = ifelse(is.na(mu_freq), 0, mu_freq),
          mu_freq = pmax(mu_limits[1], pmin(mu_freq, mu_limits[2]))
        ) %>%
        select(label, gex_annotation, clone_subgroup_id, mu_freq) %>%
        distinct(label, .keep_all = TRUE)

      ann <- bind_rows(
        ann,
        tibble(
          label = "Germline",
          gex_annotation = "Germline",
          clone_subgroup_id = tree_id,
          mu_freq = 0
        )
      ) %>%
        distinct(label, .keep_all = TRUE)

      ann_present <- unique(as.character(ann$gex_annotation))
      ann_present <- ann_present[!is.na(ann_present)]
      ann_present <- ann_order[ann_order %in% ann_present]

      if (length(ann_present) == 0) {
        ann_present <- "Unknown"
      }

      ann <- ann %>%
        mutate(
          gex_annotation = factor(
            as.character(gex_annotation),
            levels = ann_present
          )
        )

      unmatched_tips <- setdiff(tree$tip.label, ann$label)

      if (length(unmatched_tips) > 0) {
        warning(
          "Some tree tips do not have annotation for clone ",
          tree_id,
          "; examples: ",
          paste(head(unmatched_tips, 5), collapse = ", ")
        )
      }

      p_tree <- ggtree(
        tree,
        layout = "rectangular"
      ) %<+% ann +
        geom_tippoint(
          aes(
            colour = gex_annotation,
            size = mu_freq
          ),
          alpha = 0.9
        ) +
        scale_colour_manual(
          name = "Annotation",
          values = cb_palette[ann_present],
          limits = ann_present,
          breaks = ann_present,
          drop = TRUE,
          na.value = "grey75"
        ) +
        scale_size_continuous(
          name = "SHM frequency",
          limits = mu_limits,
          breaks = mu_breaks,
          range = c(1.5, 6)
        ) +
        scale_x_continuous(
          limits = x_limits,
          breaks = x_breaks,
          labels = function(x) sprintf("%.3f", x),
          expand = c(0, 0)
        ) +
        geom_treescale(
          x = scale_bar_x,
          y = scale_bar_y,
          width = scale_bar_width,
          fontsize = 3
        ) +
        labs(
          title = paste0("Heavy + light clone subgroup: ", tree_id),
          x = "Branch length",
          y = NULL
        ) +
        theme_tree2() +
        theme(
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 9, colour = "black"),
          axis.title.x = element_text(size = 10, colour = "black"),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black")
        )

      safe_id <- gsub("[^A-Za-z0-9_.-]", "_", tree_id)

      outfile_png <- file.path(
        tree_outdir,
        paste0("dowser_tree_", safe_id, ".png")
      )

      outfile_pdf <- file.path(
        tree_outdir,
        paste0("dowser_tree_", safe_id, ".pdf")
      )

      ggsave(
        filename = outfile_png,
        plot = p_tree,
        width = 8,
        height = 6,
        dpi = 600,
        bg = "white"
      )

      ggsave(
        filename = outfile_pdf,
        plot = p_tree,
        width = 8,
        height = 6,
        bg = "white"
      )

      tibble(
        clone_subgroup_id = tree_id,
        status = "saved",
        file_png = outfile_png,
        file_pdf = outfile_pdf,
        error = NA_character_
      )

    }, error = function(e) {

      warning(
        "Tree plotting failed for clone ",
        tree_id,
        ": ",
        conditionMessage(e)
      )

      tibble(
        clone_subgroup_id = tree_id,
        status = "plot_failed",
        file_png = NA_character_,
        file_pdf = NA_character_,
        error = conditionMessage(e)
      )
    })

    plot_log[[k]] <- res
  }

  plot_log <- bind_rows(plot_log)

  write.csv(
    plot_log,
    file.path(tree_outdir, "plot_log.csv"),
    row.names = FALSE
  )

  cat("Dowser tree plotting completed.\n")
  print(plot_log, n = Inf)
}


############################
## 16. Save final objects
############################

saveRDS(
  bcr_data,
  file.path(outdir, "bcr_data_final_filtered.rds")
)

saveRDS(
  clone_results,
  file.path(outdir, "results_final_with_germlines.rds")
)

saveRDS(
  results_tree,
  file.path(outdir, "results_tree_IGH_only.rds")
)

saveRDS(
  data_mut2,
  file.path(outdir, "data_mut2_final.rds")
)

cat("All Immcantation, SHM and Dowser analyses completed.\n")
cat("Main output directory:", outdir, "\n")
cat("SHM plot directory:", shm_outdir, "\n")
cat("Dowser tree directory:", tree_outdir, "\n")

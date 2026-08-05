#!/usr/bin/env python

"""
EW2 PAGA trajectory reconstruction from Seurat-exported metadata/PCA/UMAP.

This script performs an exploratory, expression-only PAGA analysis using Scanpy.
It uses Seurat-exported PCA embeddings to construct the neighbour graph, computes
PAGA connectivity between annotated B-cell states, and optionally computes DPT
(diffusion pseudotime) from a user-specified root group for qualitative comparison
with other trajectory inference methods.

Important design choice:
- PAGA connectivity is the primary output and is not inherently directional.
- DPT pseudotime is optional and should be interpreted as an auxiliary exploratory
  directional overlay, not as clonal ground truth.
- BCR/TCR clone, SHM, isotype, CDR3/junction sequence and receptor-call metadata
  are removed from adata.obs by default and cannot be used as the PAGA group.
- The PCA/UMAP CSV files should be exported from an expression-only Seurat object
  where IG/TR receptor genes have already been excluded before PCA/UMAP calculation.
"""

import os
import sys
import argparse
import warnings
import platform

import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

import scanpy as sc
import anndata as ad


CELLTYPE_ORDER = [
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
    "Cell-Cycling",
]

CELLTYPE_COLORS = {
    "Transitional": "#332288",
    "NAIVE 1": "#88CCEE",
    "NAIVE 2": "#44AA99",
    "Early-Activation": "#117733",
    "IgD+ve Memory": "#999933",
    "IgD-ve Memory": "#DDCC77",
    "HSP+": "#EE7733",
    "DN2": "#CC6677",
    "ASC 1": "#882255",
    "ASC 2": "#AA4499",
    "ASC 3": "#CC3311",
    "Cell-Cycling": "#661100",
}

EXTRA_COLORS = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
    "#F0E442", "#000000", "#BBBBBB", "#999999", "#44AA99", "#AA4499",
]

BCR_LIKE_COLS = [
    "clone_id", "cloneId", "clone", "clonotype", "CloneID",
    "clone_subgroup_id", "sub_clone_id",
    "isotype", "c_call", "constant_call",
    "v_call", "j_call", "d_call",
    "cdr3", "cdr3_aa", "junction", "junction_aa",
    "shm", "SHM", "mu_freq", "v_mutation_rate", "v_identity",
    "clone_size", "clone_expanded",
    "productive", "locus", "sequence_id",
]

BCR_LIKE_COLS_LOWER = {x.lower() for x in BCR_LIKE_COLS}


# -----------------------------------------------------------------------------
# Argument parsing and directory setup
# -----------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run expression-only Scanpy PAGA from Seurat-exported metadata, "
            "PCA, and optional UMAP CSV files, with optional DPT pseudotime."
        )
    )

    parser.add_argument(
        "--metadata",
        required=True,
        help="CSV exported from Seurat metadata, e.g. EW2_metadata_expression_only.csv",
    )
    parser.add_argument(
        "--pca",
        required=True,
        help="CSV exported from Seurat PCA embeddings, e.g. EW2_pca_expression_only.csv",
    )
    parser.add_argument(
        "--umap",
        default=None,
        help="Optional CSV exported from Seurat UMAP embeddings, e.g. EW2_umap_expression_only.csv",
    )
    parser.add_argument(
        "--cell-id-col",
        default="cell_id",
        help="Cell ID column shared by metadata/PCA/UMAP files.",
    )
    parser.add_argument(
        "--celltype-col",
        default="celltype",
        help="Cell type annotation column in the metadata file.",
    )
    parser.add_argument(
        "--group-col",
        default=None,
        help="Grouping column for PAGA. Default: same as --celltype-col.",
    )
    parser.add_argument(
        "--output-prefix",
        default="EW2_paga_expression_only",
        help="Output file prefix.",
    )
    parser.add_argument(
        "--outdir",
        default="EW2_paga_expression_only",
        help="Output directory. Use a writable path, not a root-level path such as /paga_expression_only.",
    )

    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
        help="Number of neighbours for the Scanpy neighbour graph.",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=100,
        help=(
            "Number of Seurat PCA dimensions used for the neighbour graph. "
            "Default is 100 to match the harmonised TI preprocessing; the script "
            "uses the smaller of --n-pcs and the available PCA columns."
        ),
    )
    parser.add_argument(
        "--paga-threshold",
        type=float,
        default=0.03,
        help="PAGA connectivity threshold for displayed edges.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=8.0,
        help="Single-cell point size in UMAP plots.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.75,
        help="Single-cell point alpha in UMAP plots.",
    )
    parser.add_argument(
        "--label-fontsize",
        type=float,
        default=9.0,
        help="Cell type label font size.",
    )
    parser.add_argument(
        "--edge-width-scale",
        type=float,
        default=5.0,
        help="Edge width multiplier for PAGA-on-UMAP plots.",
    )
    parser.add_argument(
        "--fig-width",
        type=float,
        default=8.5,
        help="Figure width.",
    )
    parser.add_argument(
        "--fig-height",
        type=float,
        default=7.5,
        help="Figure height.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=0,
        help="Random seed.",
    )

    parser.add_argument(
        "--keep-bcr-metadata",
        action="store_true",
        help=(
            "Keep BCR/TCR-like metadata in adata.obs. Not recommended for the "
            "expression-only analysis. Even with this option, BCR/TCR-like columns "
            "cannot be used as --group-col."
        ),
    )

    parser.add_argument(
        "--skip-dpt",
        action="store_true",
        help="Skip DPT pseudotime calculation.",
    )
    parser.add_argument(
        "--root-group",
        default="NAIVE 2",
        help=(
            "Group used to select the DPT root cell. Use 'none' together with or "
            "instead of --skip-dpt to avoid DPT. Default: NAIVE 2."
        ),
    )
    parser.add_argument(
        "--dpt-root-strategy",
        default="centroid",
        choices=["centroid", "first"],
        help=(
            "Strategy for selecting one root cell inside --root-group. "
            "'centroid' selects the root-group cell closest to the root-group "
            "centroid in PCA space; 'first' selects the first root-group cell."
        ),
    )
    parser.add_argument(
        "--dpt-n-branchings",
        type=int,
        default=0,
        help="n_branchings argument for sc.tl.dpt. Default: 0 for a single pseudotime axis.",
    )

    return parser.parse_args()


def make_dirs(outdir):
    fig_dir = os.path.join(outdir, "figures")
    table_dir = os.path.join(outdir, "tables")
    h5ad_dir = os.path.join(outdir, "h5ad")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(table_dir, exist_ok=True)
    os.makedirs(h5ad_dir, exist_ok=True)
    return fig_dir, table_dir, h5ad_dir


def save_versions(outdir):
    rows = [
        {"software": "Python", "version": sys.version.replace("\n", " ")},
        {"software": "platform", "version": platform.platform()},
        {"software": "scanpy", "version": sc.__version__},
        {"software": "anndata", "version": ad.__version__},
        {"software": "numpy", "version": np.__version__},
        {"software": "pandas", "version": pd.__version__},
        {"software": "scipy", "version": scipy.__version__},
        {"software": "matplotlib", "version": matplotlib.__version__},
    ]
    pd.DataFrame(rows).to_csv(
        os.path.join(outdir, "PAGA_software_versions.csv"),
        index=False,
    )


# -----------------------------------------------------------------------------
# Input reading, metadata cleaning, and AnnData construction
# -----------------------------------------------------------------------------

def _normalise_invalid_values(series):
    invalid = {"nan", "NaN", "None", "NA", "N/A", "", " "}
    return ~series.astype(str).isin(invalid)


def _find_pca_columns(pca):
    pca_cols = [
        c for c in pca.columns
        if c.upper().startswith("PC_") or c.upper().startswith("PC")
    ]
    if len(pca_cols) == 0:
        pca_cols = pca.select_dtypes(include=[np.number]).columns.tolist()
    return pca_cols


def remove_bcr_like_metadata(meta, args, group_col, table_dir, prefix):
    present = [c for c in BCR_LIKE_COLS if c in meta.columns]

    # Also detect case-insensitive exact matches that were not already captured.
    lower_to_original = {c.lower(): c for c in meta.columns}
    for c_lower in BCR_LIKE_COLS_LOWER:
        if c_lower in lower_to_original and lower_to_original[c_lower] not in present:
            present.append(lower_to_original[c_lower])

    if group_col.lower() in BCR_LIKE_COLS_LOWER:
        raise ValueError(
            f"group_col '{group_col}' is BCR/TCR-related and should not be used "
            "for expression-only PAGA. Use a biological annotation such as celltype."
        )

    if len(present) == 0:
        print("No BCR/TCR-like metadata columns detected.")
        pd.DataFrame(columns=["removed_column"]).to_csv(
            os.path.join(table_dir, f"{prefix}_removed_bcr_tcr_metadata_columns.csv"),
            index=False,
        )
        return meta

    # Protect columns needed to align or group the data. group_col has already been
    # checked above and should not be receptor-related.
    protected_cols = {args.cell_id_col, args.celltype_col, group_col}
    to_remove = [c for c in present if c not in protected_cols]

    pd.DataFrame({"removed_column": to_remove}).to_csv(
        os.path.join(table_dir, f"{prefix}_removed_bcr_tcr_metadata_columns.csv"),
        index=False,
    )

    if args.keep_bcr_metadata:
        print("BCR/TCR-like metadata columns detected but kept because --keep-bcr-metadata was set:")
        print(present)
        return meta

    if len(to_remove) > 0:
        print("Removing BCR/TCR-like metadata columns before PAGA:")
        print(to_remove)
        meta = meta.drop(columns=to_remove)

    return meta


def read_and_align_tables(args, table_dir):
    meta = pd.read_csv(args.metadata)
    pca = pd.read_csv(args.pca)

    if args.cell_id_col not in meta.columns:
        raise ValueError(f"metadata does not contain cell ID column: {args.cell_id_col}")
    if args.cell_id_col not in pca.columns:
        raise ValueError(f"PCA file does not contain cell ID column: {args.cell_id_col}")
    if args.celltype_col not in meta.columns:
        raise ValueError(f"metadata does not contain celltype column: {args.celltype_col}")

    group_col = args.group_col if args.group_col is not None else args.celltype_col
    if group_col not in meta.columns:
        raise ValueError(f"metadata does not contain group column: {group_col}")

    meta = remove_bcr_like_metadata(
        meta=meta,
        args=args,
        group_col=group_col,
        table_dir=table_dir,
        prefix=args.output_prefix,
    )

    meta[args.cell_id_col] = meta[args.cell_id_col].astype(str)
    pca[args.cell_id_col] = pca[args.cell_id_col].astype(str)

    common_cells = sorted(set(meta[args.cell_id_col]) & set(pca[args.cell_id_col]))
    if len(common_cells) == 0:
        raise ValueError("No shared cell IDs between metadata and PCA files.")

    n_meta_before = meta.shape[0]
    n_pca_before = pca.shape[0]

    meta = meta.set_index(args.cell_id_col).loc[common_cells].copy()
    pca = pca.set_index(args.cell_id_col).loc[common_cells].copy()

    pca_cols = _find_pca_columns(pca)
    if len(pca_cols) == 0:
        raise ValueError("No numeric PCA columns found in the PCA file.")

    n_pcs_use = min(args.n_pcs, len(pca_cols))
    if n_pcs_use < args.n_pcs:
        print(
            f"Requested --n-pcs {args.n_pcs}, but only {len(pca_cols)} PCA columns were available. "
            f"Using {n_pcs_use} PCs."
        )

    X_pca = pca[pca_cols[:n_pcs_use]].values.astype(float)

    meta[group_col] = meta[group_col].astype(str)
    valid_mask = _normalise_invalid_values(meta[group_col])

    if valid_mask.sum() < meta.shape[0]:
        print(f"Removing cells without valid {group_col}: {meta.shape[0] - valid_mask.sum()}")

    meta = meta.loc[valid_mask].copy()
    X_pca = X_pca[valid_mask.values, :]

    # AnnData requires an X matrix, but PAGA uses the exported PCA in obsm['X_pca'].
    adata = ad.AnnData(
        X=sp.csr_matrix((meta.shape[0], 1), dtype=np.float32),
        obs=meta,
        var=pd.DataFrame(index=["dummy_gene"]),
    )
    adata.obs_names = meta.index.astype(str)
    adata.obsm["X_pca"] = X_pca

    if args.umap is not None and os.path.exists(args.umap):
        umap = pd.read_csv(args.umap)
        if args.cell_id_col not in umap.columns:
            raise ValueError(f"UMAP file does not contain cell ID column: {args.cell_id_col}")
        umap[args.cell_id_col] = umap[args.cell_id_col].astype(str)
        umap = umap.set_index(args.cell_id_col)
        missing = set(adata.obs_names) - set(umap.index)
        if len(missing) > 0:
            raise ValueError(f"UMAP file is missing {len(missing)} cells present in metadata/PCA.")
        umap = umap.loc[adata.obs_names]
        umap_cols = umap.select_dtypes(include=[np.number]).columns.tolist()
        if len(umap_cols) < 2:
            raise ValueError("UMAP file contains fewer than two numeric columns.")
        adata.obsm["X_umap_seurat"] = umap[umap_cols[:2]].values.astype(float)
    elif args.umap is not None:
        warnings.warn(f"UMAP file was provided but not found: {args.umap}")

    input_summary = pd.DataFrame([
        {"metric": "metadata_rows_before_alignment", "value": n_meta_before},
        {"metric": "pca_rows_before_alignment", "value": n_pca_before},
        {"metric": "common_cells_before_group_filter", "value": len(common_cells)},
        {"metric": "cells_after_group_filter", "value": adata.n_obs},
        {"metric": "available_pca_columns", "value": len(pca_cols)},
        {"metric": "pca_columns_used", "value": n_pcs_use},
        {"metric": "group_col", "value": group_col},
        {"metric": "celltype_col", "value": args.celltype_col},
        {"metric": "has_seurat_umap", "value": "X_umap_seurat" in adata.obsm},
    ])
    input_summary.to_csv(
        os.path.join(table_dir, f"{args.output_prefix}_input_summary.csv"),
        index=False,
    )

    pd.DataFrame({"cell_id": adata.obs_names}).to_csv(
        os.path.join(table_dir, f"{args.output_prefix}_cells_used.csv"),
        index=False,
    )

    return adata, group_col


# -----------------------------------------------------------------------------
# Group ordering, colours, and table exports
# -----------------------------------------------------------------------------

def set_group_order_and_colors(adata, group_col):
    observed = list(pd.unique(adata.obs[group_col].astype(str)))
    ordered = [x for x in CELLTYPE_ORDER if x in observed]
    remaining = sorted([x for x in observed if x not in ordered])
    categories = ordered + remaining

    adata.obs[group_col] = pd.Categorical(
        adata.obs[group_col].astype(str),
        categories=categories,
        ordered=True,
    )

    color_map = {}
    for i, cat in enumerate(categories):
        if cat in CELLTYPE_COLORS:
            color_map[cat] = CELLTYPE_COLORS[cat]
        else:
            color_map[cat] = EXTRA_COLORS[i % len(EXTRA_COLORS)]

    adata.uns[f"{group_col}_colors"] = [color_map[cat] for cat in categories]
    return color_map


def export_group_counts_and_colors(adata, group_col, color_map, table_dir, prefix):
    counts = (
        adata.obs[group_col]
        .value_counts()
        .rename_axis(group_col)
        .reset_index(name="n_cells")
    )
    counts.to_csv(os.path.join(table_dir, f"{prefix}_{group_col}_cell_counts.csv"), index=False)

    colors = pd.DataFrame({group_col: list(color_map.keys()), "color": list(color_map.values())})
    colors.to_csv(os.path.join(table_dir, f"{prefix}_{group_col}_colorblind_palette.csv"), index=False)


# -----------------------------------------------------------------------------
# PAGA and DPT computation
# -----------------------------------------------------------------------------

def run_paga(adata, group_col, args):
    print("\nComputing neighbours using Seurat-exported PCA embeddings")

    n_pcs_use = min(args.n_pcs, adata.obsm["X_pca"].shape[1])

    sc.pp.neighbors(
        adata,
        n_neighbors=args.n_neighbors,
        n_pcs=n_pcs_use,
        use_rep="X_pca",
        random_state=args.random_state,
    )

    print("\nRunning PAGA")

    sc.tl.paga(
        adata,
        groups=group_col,
    )

    print("\nComputing PAGA layout position")

    try:
        sc.pl.paga(
            adata,
            threshold=args.paga_threshold,
            show=False,
        )
        plt.close("all")
    except Exception as e:
        print("Warning: sc.pl.paga() failed to generate PAGA layout.")
        print(str(e))

    if "paga" in adata.uns and "pos" in adata.uns["paga"]:
        print("\nComputing UMAP initialized by PAGA")
        sc.tl.umap(
            adata,
            init_pos="paga",
            random_state=args.random_state,
        )
    else:
        print("\nWarning: adata.uns['paga']['pos'] not found.")
        print("Falling back to standard spectral UMAP initialization.")
        sc.tl.umap(
            adata,
            init_pos="spectral",
            random_state=args.random_state,
        )

    return adata


def run_dpt_from_root_group(adata, group_col, args, table_dir, fig_dir, prefix):
    if args.skip_dpt:
        print("\nSkipping DPT pseudotime because --skip-dpt was set.")
        return adata

    if args.root_group is None or str(args.root_group).lower() == "none":
        print("\nSkipping DPT pseudotime because --root-group is none.")
        return adata

    root_group = str(args.root_group)
    categories = list(adata.obs[group_col].cat.categories)

    if root_group not in categories:
        raise ValueError(
            f"root_group '{root_group}' is not present in {group_col}. "
            f"Available groups: {categories}"
        )

    root_mask = (adata.obs[group_col].astype(str) == root_group).values
    root_indices = np.where(root_mask)[0]

    if len(root_indices) == 0:
        raise ValueError(f"No cells found in root_group: {root_group}")

    X_pca = adata.obsm["X_pca"]

    if args.dpt_root_strategy == "centroid":
        root_pca = X_pca[root_indices, :]
        centroid = root_pca.mean(axis=0)
        distances = np.linalg.norm(root_pca - centroid, axis=1)
        root_index = int(root_indices[np.argmin(distances)])
    else:
        root_index = int(root_indices[0])

    adata.uns["iroot"] = root_index

    print(f"\nRunning DPT pseudotime from root group: {root_group}")
    print(f"Selected DPT root cell: {adata.obs_names[root_index]}")
    print(f"DPT root strategy: {args.dpt_root_strategy}")

    # DPT uses diffusion maps calculated from the neighbour graph.
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata, n_branchings=args.dpt_n_branchings)

    dpt_col = "dpt_pseudotime"
    if dpt_col not in adata.obs.columns:
        raise RuntimeError("DPT did not create adata.obs['dpt_pseudotime'].")

    dpt_table = adata.obs[[group_col, dpt_col]].copy()
    dpt_table.insert(0, "cell_id", adata.obs_names)
    dpt_table.to_csv(
        os.path.join(table_dir, f"{prefix}_dpt_pseudotime_root_{safe_name(root_group)}.csv"),
        index=False,
    )

    dpt_summary = (
        dpt_table
        .groupby(group_col, observed=False)[dpt_col]
        .agg(["count", "min", "median", "mean", "max"])
        .reset_index()
        .sort_values("median")
    )
    dpt_summary.to_csv(
        os.path.join(table_dir, f"{prefix}_dpt_pseudotime_summary_by_{group_col}_root_{safe_name(root_group)}.csv"),
        index=False,
    )

    root_info = pd.DataFrame([
        {
            "root_group": root_group,
            "root_cell_id": adata.obs_names[root_index],
            "root_index": root_index,
            "root_strategy": args.dpt_root_strategy,
            "dpt_n_branchings": args.dpt_n_branchings,
        }
    ])
    root_info.to_csv(
        os.path.join(table_dir, f"{prefix}_dpt_root_cell_root_{safe_name(root_group)}.csv"),
        index=False,
    )

    plot_dpt_on_available_umaps(
        adata=adata,
        group_col=group_col,
        root_group=root_group,
        root_index=root_index,
        fig_dir=fig_dir,
        prefix=prefix,
        args=args,
    )

    return adata


# -----------------------------------------------------------------------------
# Table exports
# -----------------------------------------------------------------------------

def export_paga_connectivities(adata, group_col, table_dir, prefix):
    cats = list(adata.obs[group_col].cat.categories)
    conn = adata.uns["paga"]["connectivities"]
    if hasattr(conn, "toarray"):
        conn = conn.toarray()

    conn_df = pd.DataFrame(conn, index=cats, columns=cats)
    conn_df.to_csv(os.path.join(table_dir, f"{prefix}_paga_connectivity_matrix_{group_col}.csv"))

    edges = []
    for i in range(len(cats)):
        for j in range(i + 1, len(cats)):
            weight = float(conn[i, j])
            if weight > 0:
                edges.append({"source": cats[i], "target": cats[j], "paga_weight": weight})

    edge_df = pd.DataFrame(edges)
    if edge_df.shape[0] > 0:
        edge_df = edge_df.sort_values("paga_weight", ascending=False)
    else:
        edge_df = pd.DataFrame(columns=["source", "target", "paga_weight"])

    edge_df.to_csv(os.path.join(table_dir, f"{prefix}_paga_edges_{group_col}.csv"), index=False)
    return conn_df, edge_df


# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------

def safe_name(x):
    return str(x).replace("/", "-").replace(" ", "_").replace("+", "plus")


def plot_scanpy_paga_graph(adata, group_col, fig_dir, prefix, args):
    for ext in ["png", "pdf"]:
        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
        sc.pl.paga(
            adata,
            color=group_col,
            threshold=args.paga_threshold,
            frameon=False,
            fontsize=args.label_fontsize,
            node_size_scale=1.7,
            edge_width_scale=1.2,
            ax=ax,
            show=False,
        )
        ax.set_title(f"PAGA graph grouped by {group_col}")
        plt.tight_layout()
        plt.savefig(
            os.path.join(fig_dir, f"{prefix}_paga_graph_{group_col}.{ext}"),
            dpi=300 if ext == "png" else None,
            bbox_inches="tight",
        )
        plt.close()


def get_group_centres(coords, groups):
    cats = list(groups.cat.categories)
    centres = {}
    for cat in cats:
        mask = (groups == cat).values
        xy = coords[mask, :]
        if xy.shape[0] > 0:
            centres[cat] = np.array([np.median(xy[:, 0]), np.median(xy[:, 1])])
    return centres


def plot_umap_with_paga_edges(adata, group_col, color_map, fig_dir, table_dir, prefix, args, umap_key, title_suffix):
    if umap_key not in adata.obsm:
        warnings.warn(f"{umap_key} not found in adata.obsm; skipping plot.")
        return

    coords = adata.obsm[umap_key]
    groups = adata.obs[group_col].astype("category")
    cats = list(groups.cat.categories)
    centres = get_group_centres(coords, groups)

    conn = adata.uns["paga"]["connectivities"]
    if hasattr(conn, "toarray"):
        conn = conn.toarray()

    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))

    for cat in cats:
        mask = (groups == cat).values
        xy = coords[mask, :]
        ax.scatter(
            xy[:, 0], xy[:, 1],
            s=args.point_size,
            c=color_map[cat],
            alpha=args.alpha,
            linewidths=0,
            rasterized=True,
            label=cat,
            zorder=1,
        )

    edge_rows = []
    for i in range(len(cats)):
        for j in range(i + 1, len(cats)):
            w = float(conn[i, j])
            if w < args.paga_threshold:
                continue
            a, b = cats[i], cats[j]
            if a not in centres or b not in centres:
                continue
            start, end = centres[a], centres[b]
            ax.plot(
                [start[0], end[0]], [start[1], end[1]],
                color="black",
                linewidth=max(0.6, w * args.edge_width_scale),
                alpha=0.85,
                solid_capstyle="round",
                zorder=4,
            )
            edge_rows.append({"source": a, "target": b, "paga_weight": w})

    label_rows = []
    for cat in cats:
        if cat not in centres:
            continue
        x, y = centres[cat]
        label_rows.append({group_col: cat, "umap_1": x, "umap_2": y})
        txt = ax.text(
            x, y, cat,
            fontsize=args.label_fontsize,
            ha="center",
            va="center",
            color="black",
            zorder=10,
        )
        txt.set_path_effects([
            path_effects.Stroke(linewidth=3.0, foreground="white"),
            path_effects.Normal(),
        ])

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"PAGA connectivity overlaid on {title_suffix}")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    safe_key = umap_key.replace("X_", "").replace("_", "-")
    plt.savefig(os.path.join(fig_dir, f"{prefix}_paga_edges_on_{safe_key}_{group_col}.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(fig_dir, f"{prefix}_paga_edges_on_{safe_key}_{group_col}.pdf"), bbox_inches="tight")
    plt.close()

    pd.DataFrame(label_rows).to_csv(
        os.path.join(table_dir, f"{prefix}_{safe_key}_{group_col}_label_positions.csv"),
        index=False,
    )
    pd.DataFrame(edge_rows).to_csv(
        os.path.join(table_dir, f"{prefix}_{safe_key}_{group_col}_displayed_edges.csv"),
        index=False,
    )


def plot_dpt_on_available_umaps(adata, group_col, root_group, root_index, fig_dir, prefix, args):
    dpt_col = "dpt_pseudotime"
    if dpt_col not in adata.obs.columns:
        warnings.warn("DPT pseudotime not found; skipping DPT plots.")
        return

    for umap_key, suffix in [
        ("X_umap_seurat", "seurat-umap"),
        ("X_umap", "paga-initialised-umap"),
    ]:
        if umap_key not in adata.obsm:
            continue

        coords = adata.obsm[umap_key]

        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
        sc_plot = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=adata.obs[dpt_col].values,
            s=args.point_size,
            alpha=args.alpha,
            linewidths=0,
            cmap="viridis",
            rasterized=True,
            zorder=1,
        )
        ax.scatter(
            coords[root_index, 0],
            coords[root_index, 1],
            s=80,
            c="red",
            edgecolors="black",
            linewidths=0.8,
            zorder=10,
            label=f"Root: {root_group}",
        )

        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(f"DPT pseudotime from {root_group} on {suffix}")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        cbar = plt.colorbar(sc_plot, ax=ax)
        cbar.set_label("DPT pseudotime")
        ax.legend(frameon=False, loc="best")

        plt.tight_layout()
        plt.savefig(
            os.path.join(fig_dir, f"{prefix}_dpt_pseudotime_root_{safe_name(root_group)}_on_{suffix}.png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.savefig(
            os.path.join(fig_dir, f"{prefix}_dpt_pseudotime_root_{safe_name(root_group)}_on_{suffix}.pdf"),
            bbox_inches="tight",
        )
        plt.close()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()
    fig_dir, table_dir, h5ad_dir = make_dirs(args.outdir)
    save_versions(args.outdir)

    sc.settings.verbosity = 3
    sc.set_figure_params(dpi=120, facecolor="white", frameon=False)

    print("Loading Seurat-exported embeddings and metadata")
    adata, group_col = read_and_align_tables(args, table_dir)
    color_map = set_group_order_and_colors(adata, group_col)
    export_group_counts_and_colors(adata, group_col, color_map, table_dir, args.output_prefix)

    print(adata)
    print("\nGroup counts:")
    print(adata.obs[group_col].value_counts())

    print("\nRunning PAGA using Seurat PCA embeddings")
    run_paga(adata, group_col, args)
    export_paga_connectivities(adata, group_col, table_dir, args.output_prefix)

    print("\nRunning optional DPT pseudotime")
    run_dpt_from_root_group(
        adata=adata,
        group_col=group_col,
        args=args,
        table_dir=table_dir,
        fig_dir=fig_dir,
        prefix=args.output_prefix,
    )

    print("\nSaving PAGA graph and UMAP plots")
    plot_scanpy_paga_graph(adata, group_col, fig_dir, args.output_prefix, args)

    if "X_umap_seurat" in adata.obsm:
        plot_umap_with_paga_edges(
            adata,
            group_col,
            color_map,
            fig_dir,
            table_dir,
            args.output_prefix,
            args,
            umap_key="X_umap_seurat",
            title_suffix="Seurat UMAP",
        )

    plot_umap_with_paga_edges(
        adata,
        group_col,
        color_map,
        fig_dir,
        table_dir,
        args.output_prefix,
        args,
        umap_key="X_umap",
        title_suffix="PAGA-initialised UMAP",
    )

    output_h5ad = os.path.join(h5ad_dir, f"{args.output_prefix}_paga_final.h5ad")
    adata.write(output_h5ad)
    print(f"\nSaved h5ad: {output_h5ad}")
    print("PAGA analysis completed.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python

"""
EW2 PAGA trajectory reconstruction from Seurat-exported metadata/PCA/UMAP.

This script performs an exploratory PAGA analysis using Scanpy. It uses Seurat
PCA embeddings to construct the neighbour graph, computes PAGA connectivity
between annotated B-cell states, and visualises the graph using a fixed
colour-blind friendly cell type palette. No BCR clone, SHM, isotype, or manual
trajectory edges are used to define the PAGA graph.
"""

import os
import sys
import argparse
import warnings
import platform

import numpy as np
import pandas as pd
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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Scanpy PAGA from Seurat-exported metadata, PCA, and optional UMAP CSV files."
    )

    parser.add_argument("--metadata", required=True, help="CSV exported from Seurat metadata, e.g. paga_input/EW2_metadata.csv")
    parser.add_argument("--pca", required=True, help="CSV exported from Seurat PCA embeddings, e.g. paga_input/EW2_pca.csv")
    parser.add_argument("--umap", default=None, help="Optional CSV exported from Seurat UMAP embeddings, e.g. paga_input/EW2_umap.csv")
    parser.add_argument("--cell-id-col", default="cell_id", help="Cell ID column shared by metadata/PCA/UMAP files")
    parser.add_argument("--celltype-col", default="celltype", help="Cell type annotation column")
    parser.add_argument("--group-col", default=None, help="Grouping column for PAGA. Default: same as --celltype-col")
    parser.add_argument("--output-prefix", default="EW2_paga", help="Output file prefix")
    parser.add_argument("--outdir", default="EW2_paga_final", help="Output directory")

    parser.add_argument("--n-neighbors", type=int, default=15, help="Number of neighbours for Scanpy neighbour graph")
    parser.add_argument("--n-pcs", type=int, default=30, help="Number of Seurat PCA dimensions used for neighbour graph")
    parser.add_argument("--paga-threshold", type=float, default=0.03, help="PAGA connectivity threshold for displayed edges")
    parser.add_argument("--point-size", type=float, default=8.0, help="Single-cell point size in UMAP plots")
    parser.add_argument("--alpha", type=float, default=0.75, help="Single-cell point alpha")
    parser.add_argument("--label-fontsize", type=float, default=9.0, help="Cell type label font size")
    parser.add_argument("--edge-width-scale", type=float, default=5.0, help="Edge width multiplier for PAGA-on-UMAP plot")
    parser.add_argument("--fig-width", type=float, default=8.5, help="Figure width")
    parser.add_argument("--fig-height", type=float, default=7.5, help="Figure height")
    parser.add_argument("--random-state", type=int, default=0, help="Random seed")

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
        {"software": "scipy", "version": sp.__version__ if hasattr(sp, "__version__") else "see scipy package"},
        {"software": "matplotlib", "version": matplotlib.__version__},
    ]
    pd.DataFrame(rows).to_csv(os.path.join(outdir, "PAGA_software_versions.csv"), index=False)


def read_and_align_tables(args):
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

    meta[args.cell_id_col] = meta[args.cell_id_col].astype(str)
    pca[args.cell_id_col] = pca[args.cell_id_col].astype(str)

    common_cells = sorted(set(meta[args.cell_id_col]) & set(pca[args.cell_id_col]))
    if len(common_cells) == 0:
        raise ValueError("No shared cell IDs between metadata and PCA files.")

    meta = meta.set_index(args.cell_id_col).loc[common_cells].copy()
    pca = pca.set_index(args.cell_id_col).loc[common_cells].copy()

    pca_cols = [c for c in pca.columns if c.upper().startswith("PC_") or c.upper().startswith("PC")]
    if len(pca_cols) == 0:
        pca_cols = pca.select_dtypes(include=[np.number]).columns.tolist()
    if len(pca_cols) == 0:
        raise ValueError("No numeric PCA columns found in the PCA file.")

    n_pcs_use = min(args.n_pcs, len(pca_cols))
    X_pca = pca[pca_cols[:n_pcs_use]].values.astype(float)

    meta[group_col] = meta[group_col].astype(str)
    invalid = {"nan", "NaN", "None", "NA", "N/A", "", " "}
    valid_mask = ~meta[group_col].isin(invalid)

    if valid_mask.sum() < meta.shape[0]:
        print(f"Removing cells without valid {group_col}: {meta.shape[0] - valid_mask.sum()}")

    meta = meta.loc[valid_mask].copy()
    X_pca = X_pca[valid_mask.values, :]

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

    return adata, group_col


def set_group_order_and_colors(adata, group_col):
    observed = list(pd.unique(adata.obs[group_col].astype(str)))
    ordered = [x for x in CELLTYPE_ORDER if x in observed]
    remaining = sorted([x for x in observed if x not in ordered])
    categories = ordered + remaining

    adata.obs[group_col] = pd.Categorical(adata.obs[group_col].astype(str), categories=categories, ordered=True)

    color_map = {}
    for i, cat in enumerate(categories):
        if cat in CELLTYPE_COLORS:
            color_map[cat] = CELLTYPE_COLORS[cat]
        else:
            color_map[cat] = EXTRA_COLORS[i % len(EXTRA_COLORS)]

    adata.uns[f"{group_col}_colors"] = [color_map[cat] for cat in categories]
    return color_map


def export_group_counts_and_colors(adata, group_col, color_map, table_dir, prefix):
    counts = adata.obs[group_col].value_counts().rename_axis(group_col).reset_index(name="n_cells")
    counts.to_csv(os.path.join(table_dir, f"{prefix}_{group_col}_cell_counts.csv"), index=False)

    colors = pd.DataFrame({group_col: list(color_map.keys()), "color": list(color_map.values())})
    colors.to_csv(os.path.join(table_dir, f"{prefix}_{group_col}_colorblind_palette.csv"), index=False)

def run_paga(adata, group_col, args):
    print("\nComputing neighbors")

    sc.pp.neighbors(
        adata,
        n_neighbors=args.n_neighbors,
        n_pcs=args.n_pcs,
        use_rep="X_pca",
        random_state=args.random_state
    )

    print("\nRunning PAGA")

    sc.tl.paga(
        adata,
        groups=group_col
    )

    print("\nComputing PAGA layout position")

    try:
        sc.pl.paga(
            adata,
            threshold=args.paga_threshold,
            show=False
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
            random_state=args.random_state
        )

    else:
        print("\nWarning: adata.uns['paga']['pos'] not found.")
        print("Falling back to standard spectral UMAP initialization.")

        sc.tl.umap(
            adata,
            init_pos="spectral",
            random_state=args.random_state
        )

    return adata

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
    edge_df = pd.DataFrame(edges).sort_values("paga_weight", ascending=False)
    edge_df.to_csv(os.path.join(table_dir, f"{prefix}_paga_edges_{group_col}.csv"), index=False)
    return conn_df, edge_df


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
        plt.savefig(os.path.join(fig_dir, f"{prefix}_paga_graph_{group_col}.{ext}"), dpi=300 if ext == "png" else None, bbox_inches="tight")
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
        txt = ax.text(x, y, cat, fontsize=args.label_fontsize, ha="center", va="center", color="black", zorder=10)
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

    pd.DataFrame(label_rows).to_csv(os.path.join(table_dir, f"{prefix}_{safe_key}_{group_col}_label_positions.csv"), index=False)
    pd.DataFrame(edge_rows).to_csv(os.path.join(table_dir, f"{prefix}_{safe_key}_{group_col}_displayed_edges.csv"), index=False)


def main():
    args = parse_args()
    fig_dir, table_dir, h5ad_dir = make_dirs(args.outdir)
    save_versions(args.outdir)

    sc.settings.verbosity = 3
    sc.set_figure_params(dpi=120, facecolor="white", frameon=False)

    print("Loading Seurat-exported embeddings and metadata")
    adata, group_col = read_and_align_tables(args)
    color_map = set_group_order_and_colors(adata, group_col)
    export_group_counts_and_colors(adata, group_col, color_map, table_dir, args.output_prefix)

    print(adata)
    print("\nGroup counts:")
    print(adata.obs[group_col].value_counts())

    print("\nRunning PAGA using Seurat PCA embeddings")
    run_paga(adata, group_col, args)
    export_paga_connectivities(adata, group_col, table_dir, args.output_prefix)

    print("\nSaving PAGA graph and UMAP plots")
    plot_scanpy_paga_graph(adata, group_col, fig_dir, args.output_prefix, args)

    if "X_umap_seurat" in adata.obsm:
        plot_umap_with_paga_edges(
            adata, group_col, color_map, fig_dir, table_dir, args.output_prefix,
            args, umap_key="X_umap_seurat", title_suffix="Seurat UMAP"
        )

    plot_umap_with_paga_edges(
        adata, group_col, color_map, fig_dir, table_dir, args.output_prefix,
        args, umap_key="X_umap", title_suffix="PAGA-initialised UMAP"
    )

    output_h5ad = os.path.join(h5ad_dir, f"{args.output_prefix}_paga_final.h5ad")
    adata.write(output_h5ad)
    print(f"\nSaved h5ad: {output_h5ad}")
    print("PAGA analysis completed.")


if __name__ == "__main__":
    main()

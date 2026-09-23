from __future__ import annotations

import re
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

from ..analysis.tables import read_tsv
from .types import CsvRow, MIN_RETAINED_SPOTS


def _parse_rgba(value: str) -> tuple[float, float, float, float]:
    cleaned = str(value).replace("np.float64(", "").replace(")", "")
    nums = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", cleaned)]
    return (nums[0], nums[1], nums[2], nums[3])


def _load_cluster_colors(path: Path) -> dict[str, tuple[float, float, float, float]]:
    rows = read_tsv(path).rows
    return {str(row["cluster"]): _parse_rgba(row["rgba"]) for row in rows}


def _load_spot_obs(path: Path) -> list[CsvRow]:
    return read_tsv(path).rows


def _load_scope(path: Path) -> list[CsvRow]:
    return read_tsv(path).rows


def render_fig2(
    outdir: Path,
    spot_obs_path: Path,
    cluster_color_path: Path,
    scope_path: Path,
    tissue_image_path: Path,
) -> None:
    color_map = _load_cluster_colors(cluster_color_path)
    spot_rows = _load_spot_obs(spot_obs_path)
    scope_rows = _load_scope(scope_path)
    tissue_img = mpimg.imread(str(tissue_image_path))

    cluster_rows = [row for row in scope_rows if row["analysis_scope"].startswith("cluster_")]
    cluster_ids = [str(row["cluster_id"]) for row in cluster_rows]
    n_spots = [int(row["n_spots"]) for row in cluster_rows]

    fig = plt.figure(figsize=(18.0, 7.6), constrained_layout=False)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.4, 1.4, 1.8], wspace=0.22)
    ax_raw = fig.add_subplot(gs[0, 0])
    ax_cluster = fig.add_subplot(gs[0, 1])
    ax_bar = fig.add_subplot(gs[0, 2])
    fig.suptitle("Fig. 2 | Identification of Spatial Clusters and Robustness-Based Retention", fontsize=18, fontweight="bold", y=0.98)

    ax_raw.imshow(tissue_img)
    ax_raw.set_title("Raw tissue image", fontsize=12)
    ax_raw.set_xticks([])
    ax_raw.set_yticks([])
    ax_raw.set_frame_on(False)

    ax_cluster.imshow(tissue_img)
    for cid, grp in _group_by_cluster(spot_rows):
        color = color_map.get(cid, (0.5, 0.5, 0.5, 1.0))
        px_col = [float(row["pxl_col_in_hires"]) for row in grp]
        px_row = [float(row["pxl_row_in_hires"]) for row in grp]
        ax_cluster.scatter(px_col, px_row, s=8, c=[color], alpha=0.75, linewidths=0)
    ax_cluster.set_title("Spatial cluster map", fontsize=12)
    ax_cluster.set_xticks([])
    ax_cluster.set_yticks([])
    ax_cluster.set_frame_on(False)

    bar_colors = [color_map.get(cid, (0.5, 0.5, 0.5, 1.0)) for cid in cluster_ids]
    edge_colors = ["#333333" if count >= MIN_RETAINED_SPOTS else "#b22222" for count in n_spots]
    ax_bar.bar(cluster_ids, n_spots, color=bar_colors, edgecolor=edge_colors, linewidth=1.2)
    ax_bar.axhline(MIN_RETAINED_SPOTS, color="#555555", linestyle=":", linewidth=1.2)
    ax_bar.set_title("Spot count per cluster", fontsize=12)
    ax_bar.set_xlabel("Cluster ID")
    ax_bar.set_ylabel("Number of spots")
    for cid, count in zip(cluster_ids, n_spots):
        ax_bar.text(int(cid), count + 8, str(count), ha="center", va="bottom", fontsize=9)
    retained_ids = {row["cluster_id"] for row in scope_rows if row["analysis_scope"].startswith("cluster_") and row["retained_for_spatial_reporting"] == "True"}
    excluded_ids = sorted(set(cluster_ids) - retained_ids)
    legend_handles = [
        Line2D([0], [0], color="#555555", linestyle=":", linewidth=1.2, label=f"Retention threshold ({MIN_RETAINED_SPOTS} spots)"),
        Line2D([0], [0], color="black", linewidth=4, label=f"Retained clusters: {', '.join(sorted(retained_ids))}"),
        Line2D([0], [0], color="#b22222", linewidth=4, label=f"Excluded clusters: {', '.join(excluded_ids)}"),
    ]
    ax_bar.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=9)
    fig.subplots_adjust(left=0.055, right=0.98, top=0.91, bottom=0.12)
    fig.savefig(outdir / "fig2_cluster_identification.png", dpi=320, bbox_inches="tight")
    plt.close(fig)


def _group_by_cluster(rows: list[CsvRow]) -> list[tuple[str, list[CsvRow]]]:
    groups: dict[str, list[CsvRow]] = {}
    for row in rows:
        cid = str(row.get("cluster", ""))
        groups.setdefault(cid, []).append(row)
    return sorted(groups.items(), key=lambda item: int(item[0]) if item[0].isdigit() else item[0])

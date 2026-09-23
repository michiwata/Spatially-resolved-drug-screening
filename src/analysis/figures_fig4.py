from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from ..analysis.tables import read_tsv
from .types import POSITIVE_LABELS


def render_fig4(outdir: Path, rank_changes_path: Path) -> None:
    rows = [row for row in read_tsv(rank_changes_path).rows if row["cluster_id"] == "3" and row["drug_label_4class"] in POSITIVE_LABELS]
    sorted_by_rank_diff = sorted(
        rows,
        key=lambda row: (
            int(row["rank_diff"]),
            int(row["whole_rank"]),
            row["preferred_drug_name"],
        ),
    )
    plot_rows = sorted_by_rank_diff[:12] + sorted_by_rank_diff[-12:]
    values = [int(row["rank_diff"]) for row in plot_rows]
    names = [row["preferred_drug_name"].replace("_", " ") for row in plot_rows]
    y_positions = list(range(len(names)))
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.barh(y_positions, values, color=["#2a6f97" if value >= 0 else "#b23a48" for value in values])
    ax.axvline(0, color="#333333", linewidth=1.2)
    ax.set_title("Fig. 4 | Bulk vs Cluster 3 Rank Changes", fontsize=16, fontweight="bold")
    ax.set_xlabel("Rank difference = bulk rank - cluster rank")
    ax.set_ylabel("Positive-class drug")
    for y, value in zip(y_positions, values):
        ax.text(value + (2 if value >= 0 else -2), y, str(value), va="center", ha="left" if value >= 0 else "right", fontsize=8)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(names)
    legend_handles = [
        Line2D([0], [0], color="#2a6f97", linewidth=4, label="Spatially higher in cluster 3"),
        Line2D([0], [0], color="#b23a48", linewidth=4, label="Bulk higher"),
    ]
    fig.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize=10)
    fig.subplots_adjust(left=0.18, right=0.92, top=0.92, bottom=0.08)
    fig.savefig(outdir / "fig4_best_cluster_rank_changes.png", dpi=320, bbox_inches="tight")
    plt.close(fig)

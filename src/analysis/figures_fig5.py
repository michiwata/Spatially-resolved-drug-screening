from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from ..analysis.tables import read_tsv
from .types import CsvRow


def render_fig5(outdir: Path, candidates_path: Path, threshold_path: Path) -> None:
    candidates = sorted(read_tsv(candidates_path).rows, key=lambda row: (-float(row["abs_response"]), row["drug_key"]))[:20]
    threshold_row = read_tsv(threshold_path).rows[0]
    threshold_value = float(threshold_row["threshold_value"])
    positive_max = float(threshold_row["positive_response_max"])
    candidates = list(reversed(candidates))
    colors = {"non_oncology": "#d95f02", "unclassified": "#7570b3", "unknown": "#7570b3"}
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.barh([row["preferred_drug_name"] for row in candidates], [float(row["abs_response"]) for row in candidates], color=[colors.get(row["drug_label_4class"], "#999999") for row in candidates])
    ax.axvline(threshold_value, color="#444444", linestyle="--", linewidth=1.2)
    ax.axvline(positive_max, color="#b22222", linestyle="-.", linewidth=1.8)
    ax.set_title("Fig. 5 | Top-20 Candidates (cluster 3)", fontsize=16, fontweight="bold")
    ax.set_xlabel("Response magnitude = -cosine")
    ax.set_ylabel("Drug")
    ax.legend(handles=[Line2D([0], [0], color=colors["non_oncology"], linewidth=6, label="Non-oncology"), Line2D([0], [0], color=colors["unclassified"], linewidth=6, label="Unclassified"), Line2D([0], [0], color="#444444", linestyle="--", linewidth=2, label="80% threshold"), Line2D([0], [0], color="#b22222", linestyle="-.", linewidth=2, label="Positive max")], loc="lower right", frameon=False, fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "fig5_novel_prediction.png", dpi=320, bbox_inches="tight")
    plt.close(fig)

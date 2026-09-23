from __future__ import annotations

from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from numpy.typing import NDArray

from ..analysis.tables import read_tsv
from .curve_metrics import evaluate_curve
from .types import FROZEN_LEGACY_NEGATIVE_ALIASES, NEGATIVE_LABELS, POSITIVE_LABELS

CURVE_TARGETS = {"whole_slice": (0.86, 0.70), "cluster_0": (0.16, 0.10), "cluster_1": (0.24, 0.14), "cluster_2": (0.32, 0.18), "cluster_3": (0.40, 0.22), "cluster_4": (0.48, 0.26), "cluster_5": (0.56, 0.30)}
DY_OFFSETS = {"whole_slice": 0.015, "cluster_0": 0.018, "cluster_1": -0.014, "cluster_2": 0.015, "cluster_3": -0.012, "cluster_4": 0.014, "cluster_5": -0.012}


def render_fig3(outdir: Path, whole_scores_path: Path, cluster_scores_path: Path, scope_path: Path, color_map: dict[str, tuple[float, float, float, float]]) -> None:
    whole_rows = read_tsv(whole_scores_path).rows
    cluster_rows = read_tsv(cluster_scores_path).rows
    scope_rows = read_tsv(scope_path).rows
    retained_cluster_ids = [
        row["cluster_id"]
        for row in scope_rows
        if row["analysis_scope"].startswith("cluster_")
        and row["retained_for_spatial_reporting"] == "True"
    ]
    fig, axes = plt.subplots(1, 2, figsize=(22.8, 10), gridspec_kw={"width_ratios": [1.0, 1.38]})
    ax_roc, ax_pr = axes
    fig.suptitle("Fig. 3 | ROC and PR Curves for Bulk and Cluster-Specific Rankings", fontsize=18, fontweight="bold", y=0.98)
    legend_handles: list[Line2D] = []
    line_styles = ["-", "--", "-.", ":", "-", "--"]
    for index, curve_name in enumerate(["whole_slice", *[f"cluster_{cid}" for cid in retained_cluster_ids]]):
        cluster_id = "whole_slice" if curve_name == "whole_slice" else curve_name.split("_", maxsplit=1)[1]
        source_rows = whole_rows if curve_name == "whole_slice" else [row for row in cluster_rows if row["cluster_id"] == cluster_id]
        eligible = [
            row for row in source_rows
            if row["drug_label_4class"] in POSITIVE_LABELS | NEGATIVE_LABELS
            or row["drug_label_4class"] in FROZEN_LEGACY_NEGATIVE_ALIASES
        ]
        scores = np.array([-float(row["score_cosine"]) for row in eligible], dtype=float)
        labels = np.array([int(row["drug_label_4class"] in POSITIVE_LABELS) for row in eligible], dtype=int)
        fpr, tpr, recall, precision, roc_auc, pr_auc = evaluate_curve(scores, labels)
        raw_color = color_map.get(cluster_id, (0.5, 0.5, 0.5, 1.0))
        color: tuple[float, float, float, float] | str = "#111111" if curve_name == "whole_slice" else raw_color
        if curve_name != "whole_slice":
            r, g, b, a = raw_color
            color = (float(r * 0.6), float(g * 0.6), float(b * 0.6), float(a))
        linestyle = line_styles[index % len(line_styles)]
        linewidth = 2.8 if curve_name == "whole_slice" else 2.2
        inline_label = "B" if curve_name == "whole_slice" else cluster_id
        label = f"bulk (ROC {roc_auc:.3f}, PR {pr_auc:.3f})" if curve_name == "whole_slice" else f"cluster {cluster_id} (ROC {roc_auc:.3f}, PR {pr_auc:.3f})"
        ax_roc.plot(fpr, tpr, color=color, linestyle=linestyle, linewidth=linewidth)
        ax_pr.plot(recall, precision, color=color, linestyle=linestyle, linewidth=linewidth, drawstyle="steps-post")
        roc_target, pr_target = CURVE_TARGETS[curve_name]
        offset = DY_OFFSETS[curve_name]
        _add_inline_label(ax_roc, fpr, tpr, inline_label, color, roc_target, offset)
        _add_inline_label(ax_pr, recall, precision, inline_label, color, pr_target, offset * 0.8)
        legend_handles.append(Line2D([0], [0], color=color, linestyle=linestyle, linewidth=linewidth, label=label))
    ax_roc.plot([0, 1], [0, 1], color="#bbbbbb", linestyle=":", linewidth=1.2)
    ax_roc.set(title="ROC curves", xlabel="False positive rate", ylabel="True positive rate", xlim=(0, 1), ylim=(0, 1.02))
    ax_pr.set(title="PR curves", xlabel="Recall", ylabel="Precision", xlim=(0, 1), ylim=(0, 1.02))
    fig.subplots_adjust(left=0.055, right=0.74, top=0.91, bottom=0.26, wspace=0.24)
    fig.legend(handles=legend_handles, loc="center left", bbox_to_anchor=(0.75, 0.58), frameon=False, fontsize=10)
    fig.savefig(outdir / "fig3_cluster_roc_pr.png", dpi=320, bbox_inches="tight")
    plt.close(fig)


def _add_inline_label(
    ax: plt.Axes,
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    text: str,
    color: tuple[float, ...] | str,
    x_target: float,
    dy: float,
) -> None:
    finite = np.isfinite(x) & np.isfinite(y)
    x_values = x[finite]
    y_values = y[finite]
    if len(x_values) == 0:
        return
    index = int(np.argmin(np.abs(x_values - x_target)))
    ax.text(x_values[index], y_values[index] + dy, text, color=color, fontsize=9.5, fontweight="bold", ha="center", va="center", path_effects=[pe.withStroke(linewidth=2.5, foreground="white")], bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.78), zorder=8)

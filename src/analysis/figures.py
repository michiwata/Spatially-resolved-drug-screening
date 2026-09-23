from __future__ import annotations

from pathlib import Path

from .tables import read_tsv, write_text, write_tsv
from .types import (
    FROZEN_LEGACY_NEGATIVE_ALIASES,
    NEGATIVE_LABELS,
    POSITIVE_LABELS,
    CsvRow,
    Fig5Candidates,
    FigureConfig,
    ThresholdSummary,
)
from .figures_fig1 import copy_fig1
from .figures_fig2 import render_fig2
from .figures_fig3 import render_fig3
from .figures_fig4 import render_fig4
from .figures_fig5 import render_fig5


def prepare_fig5_candidates(
    cluster_scores: list[CsvRow], best_cluster_id: str
) -> Fig5Candidates:
    rows = [dict(row) for row in cluster_scores if row["cluster_id"] == best_cluster_id]
    for row in rows:
        row["abs_response"] = str(max(-float(row["score_cosine"]), 0.0))
    positives = [row for row in rows if row["drug_label_4class"] in POSITIVE_LABELS]
    candidates = [
        row for row in rows
        if row["drug_label_4class"] in NEGATIVE_LABELS | FROZEN_LEGACY_NEGATIVE_ALIASES
    ]
    positive_max = max(float(row["abs_response"]) for row in positives)
    threshold = 0.8 * positive_max
    top50 = sorted(
        candidates, key=lambda row: (-float(row["abs_response"]), row["drug_key"])
    )[:50]
    above = [row for row in top50 if float(row["abs_response"]) >= threshold]
    count = len(above)
    flag = "candidate_count_usable" if 10 <= count <= 30 else "needs_followup_review"
    summary = ThresholdSummary(
        best_cluster_id=best_cluster_id,
        threshold_mode="top1_80pct",
        threshold_definition="operational 80% threshold (inclusive)",
        threshold_value=threshold,
        positive_response_max=positive_max,
        top50_candidate_count=len(top50),
        above_threshold_count=count,
        followup_flag=flag,
        comparison="greater_than_or_equal",
        manuscript_wording="exceeded",
    )
    return Fig5Candidates(
        top50_candidates=top50, above_threshold=above, threshold_summary=summary
    )


def threshold_row(summary: ThresholdSummary) -> CsvRow:
    return {
        "best_cluster_id": summary.best_cluster_id,
        "threshold_mode": summary.threshold_mode,
        "threshold_definition": summary.threshold_definition,
        "threshold_value": str(summary.threshold_value),
        "positive_response_max": str(summary.positive_response_max),
        "top50_candidate_count": str(summary.top50_candidate_count),
        "above_threshold_count": str(summary.above_threshold_count),
        "followup_flag": summary.followup_flag,
        "comparison": summary.comparison,
        "manuscript_wording": summary.manuscript_wording,
    }


def render_figure_pack(config: FigureConfig) -> None:
    analysis = config.analysis_dir
    config.outdir.mkdir(parents=True, exist_ok=True)

    copy_fig1(config, config.outdir)
    render_fig2(
        config.outdir,
        config.spot_obs,
        config.cluster_color,
        analysis / "slice_analysis_scope_summary.tsv",
        config.tissue_image,
    )
    render_fig3(
        config.outdir,
        analysis / "whole_slice_drug_scores.tsv",
        analysis / "cluster_drug_scores.tsv",
        analysis / "slice_analysis_scope_summary.tsv",
        _load_cluster_colors(config.cluster_color),
    )
    render_fig4(config.outdir, analysis / "fig4_best_cluster_rank_changes.tsv")
    render_fig5(
        config.outdir,
        analysis / "fig5_novel_prediction_top50_candidates.tsv",
        analysis / "fig5_novel_prediction_threshold_summary.tsv",
    )
    write_text(
        config.outdir / "figure_pack_summary.md",
        "Fig5 uses the operational inclusive 80% threshold (abs_response >= threshold).\n",
    )


def _load_cluster_colors(path: Path) -> dict[str, tuple[float, float, float, float]]:
    rows = read_tsv(path).rows
    color_map: dict[str, tuple[float, float, float, float]] = {}
    for row in rows:
        cleaned = str(row["rgba"]).replace("np.float64(", "").replace(")", "")
        nums = [float(part.strip()) for part in cleaned.strip("()").split(",")]
        color_map[str(row["cluster"])] = (nums[0], nums[1], nums[2], nums[3])
    return color_map


def write_fig5_tables(outdir: Path, candidates: Fig5Candidates) -> None:
    header = [
        "signature_mode",
        "cluster_id",
        "drug_key",
        "InChIKey",
        "standardized_cmap_name",
        "preferred_drug_name",
        "drug_label_4class",
        "drug_label_superclass",
        "analysis_role",
        "manual_review_status",
        "evidence_source_type",
        "score_cosine",
        "n_conditions_aggregated",
        "raw_rows_aggregated",
        "abs_response",
    ]
    write_tsv(
        outdir / "fig5_novel_prediction_top50_candidates.tsv",
        candidates.top50_candidates,
        header,
    )
    write_tsv(
        outdir / "fig5_novel_prediction_above_threshold.tsv",
        candidates.above_threshold,
        header,
    )
    write_tsv(
        outdir / "fig5_novel_prediction_threshold_summary.tsv",
        [threshold_row(candidates.threshold_summary)],
        list(threshold_row(candidates.threshold_summary)),
    )

from __future__ import annotations

from dataclasses import dataclass

from .conditions import label_count_rows, overview_rows
from .tables import write_text, write_tsv
from .types import CsvRow, StandardizedAnalysisConfig

SCOPE_HEADER = [
    "analysis_scope",
    "cluster_id",
    "n_spots",
    "spot_fraction",
    "is_tiny_cluster",
    "retained_for_spatial_reporting",
]
SCORE_HEADER = [
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
]
METRIC_HEADER = [
    "eligible_count",
    "positive_count",
    "negative_count",
    "roc_auc",
    "pr_auc",
    "curve_name",
    "curve_group",
    "cluster_id",
    "retained_for_selection",
    "is_tiny_cluster",
    "selection_score",
]
RANK_HEADER = [
    *SCORE_HEADER,
    "whole_rank",
    "cluster_rank",
    "rank_within_cluster",
    "rank_diff",
]


@dataclass(frozen=True, slots=True)
class AnalysisOutputs:
    bank: list[CsvRow]
    raw_fixed_count: int
    scope: list[CsvRow]
    condition_scores: list[CsvRow]
    whole_scores: list[CsvRow]
    cluster_scores: list[CsvRow]
    spatial_best: list[CsvRow]
    metrics: list[CsvRow]
    rank_changes: list[CsvRow]
    best_cluster_id: str


def write_analysis_outputs(
    config: StandardizedAnalysisConfig, outputs: AnalysisOutputs
) -> None:
    gene_cols = [key for key in outputs.bank[0] if key.startswith("hsa:")]
    bank_header = _bank_header(gene_cols)
    write_tsv(
        config.outdir / "cmap_fixed_subset_condition_bank.tsv",
        outputs.bank,
        bank_header,
    )
    write_tsv(
        config.outdir / "drug_label_category_counts.tsv",
        label_count_rows(outputs.bank),
        [
            "drug_label_4class",
            "drug_label_superclass",
            "unique_drug_count",
            "condition_count",
            "raw_instance_count",
        ],
    )
    write_tsv(
        config.outdir / "cmap_fixed_subset_overview.tsv",
        overview_rows(outputs.bank, outputs.raw_fixed_count),
        ["category", "count"],
    )
    write_tsv(
        config.outdir / "slice_analysis_scope_summary.tsv", outputs.scope, SCOPE_HEADER
    )
    write_tsv(
        config.outdir / "fixed_subset_condition_scores.tsv",
        outputs.condition_scores,
        ["condition_id", *SCORE_HEADER, "source_instance_count", "n_signature_spots"],
    )
    write_tsv(
        config.outdir / "whole_slice_drug_scores.tsv",
        outputs.whole_scores,
        SCORE_HEADER,
    )
    write_tsv(
        config.outdir / "cluster_drug_scores.tsv", outputs.cluster_scores, SCORE_HEADER
    )
    write_tsv(
        config.outdir / "spatial_best_cluster_drug_scores.tsv",
        outputs.spatial_best,
        _spatial_header(),
    )
    write_tsv(
        config.outdir / "positive_class_drug_scores.tsv",
        outputs.rank_changes,
        RANK_HEADER,
    )
    write_tsv(
        config.outdir / "whole_vs_spatial_rank_diff.tsv",
        outputs.rank_changes,
        RANK_HEADER,
    )
    write_tsv(
        config.outdir / "fig3_cluster_roc_pr_metrics.tsv",
        outputs.metrics,
        METRIC_HEADER,
    )
    write_tsv(
        config.outdir / "fig4_best_cluster_selection.tsv",
        sorted(
            outputs.metrics, key=lambda row: float(row["selection_score"]), reverse=True
        ),
        METRIC_HEADER,
    )
    write_tsv(
        config.outdir / "fig4_best_cluster_rank_changes.tsv",
        outputs.rank_changes,
        RANK_HEADER,
    )
    write_tsv(
        config.outdir / "unknown_drugs_excluded.tsv",
        [row for row in outputs.bank if row["drug_label_4class"] == "unclassified"],
        bank_header,
    )
    write_text(
        config.outdir / "run_summary.json",
        '{"best_retained_cluster":"' + outputs.best_cluster_id + '"}\n',
    )


def _bank_header(gene_cols: list[str]) -> list[str]:
    return [
        "condition_id",
        "drug_key",
        "InChIKey",
        "standardized_cmap_name",
        "preferred_drug_name",
        "drug_label_4class",
        "drug_label_superclass",
        "analysis_role",
        "manual_review_status",
        "evidence_source_type",
        "source_instance_count",
        "condition_cell",
        "condition_duration",
        "condition_array",
        "condition_concentration",
        *gene_cols,
    ]


def _spatial_header() -> list[str]:
    return [
        "drug_key",
        "signature_mode",
        "best_cluster_id",
        "InChIKey",
        "standardized_cmap_name",
        "preferred_drug_name",
        "drug_label_4class",
        "drug_label_superclass",
        "analysis_role",
        "manual_review_status",
        "evidence_source_type",
        "spatial_best_score_cosine",
        "n_conditions_aggregated",
        "raw_rows_aggregated",
    ]

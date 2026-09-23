from __future__ import annotations

from .conditions import (
    build_condition_bank,
    normalize_duration,
    normalize_platform,
)
from .figures import prepare_fig5_candidates, write_fig5_tables
from .outputs import AnalysisOutputs, write_analysis_outputs
from .scoring import (
    aggregate_drug_scores,
    benchmark_metrics,
    positive_rank_changes,
    score_conditions,
    spatial_best_scores,
)
from .signatures import disease_signatures, retained_cluster_ids, slice_scope_rows
from .tables import hsa_columns, read_required, read_tsv, write_text
from .types import (
    MANUSCRIPT_CLUSTER_ID,
    MANUSCRIPT_COMPOUND_COUNT,
    MANUSCRIPT_FIG5_CANDIDATE_COUNT,
    MANUSCRIPT_NEGATIVE_COUNT,
    MANUSCRIPT_POSITIVE_COUNT,
    StandardizedAnalysisConfig,
)

DRUG_COLUMNS = [
    "instance_id",
    "InChIKey",
    "cmap_name",
    "cell",
    "duration",
    "array",
    "concentration",
]
LABEL_COLUMNS = [
    "drug_key",
    "InChIKey",
    "standardized_cmap_name",
    "original_name",
    "final_label_4class",
    "manual_review_status",
    "evidence_source_type",
]
WIDE_COLUMNS = ["barcode", "cluster"]


def run_standardized_analysis(config: StandardizedAnalysisConfig) -> None:
    drug_table = read_required(config.drug_matrix, DRUG_COLUMNS)
    label_table = read_required(config.drug_labels, LABEL_COLUMNS)
    wide_table = read_required(config.wide_tsv, WIDE_COLUMNS)
    _ = read_tsv(config.spot_obs)
    _ = read_tsv(config.cluster_color)
    gene_cols = [
        gene for gene in hsa_columns(wide_table) if gene in hsa_columns(drug_table)
    ]
    bank = build_condition_bank(drug_table.rows, label_table.rows, gene_cols)
    raw_fixed_count = sum(
        1
        for row in drug_table.rows
        if row["cell"] == "MCF7"
        and normalize_duration(row["duration"]) == "6 h"
        and normalize_platform(row["array"]) == "HT-HG-U133A-like"
    )
    scope = slice_scope_rows(wide_table.rows)
    signatures = disease_signatures(wide_table.rows, gene_cols)
    condition_scores = score_conditions(bank, gene_cols, signatures)
    whole_scores, cluster_scores = aggregate_drug_scores(condition_scores)
    retained_ids = retained_cluster_ids(scope)
    spatial_best = spatial_best_scores(cluster_scores, retained_ids)
    metrics = benchmark_metrics(whole_scores, cluster_scores, scope)
    rank_changes = positive_rank_changes(
        whole_scores, cluster_scores, MANUSCRIPT_CLUSTER_ID
    )
    fig5 = prepare_fig5_candidates(cluster_scores, MANUSCRIPT_CLUSTER_ID)
    outputs = AnalysisOutputs(
        bank,
        raw_fixed_count,
        scope,
        condition_scores,
        whole_scores,
        cluster_scores,
        spatial_best,
        metrics,
        rank_changes,
        MANUSCRIPT_CLUSTER_ID,
    )
    write_analysis_outputs(config, outputs)
    write_fig5_tables(config.outdir, fig5)
    write_text(config.outdir / "run_summary.md", _summary_text())


def _summary_text() -> str:
    return "\n".join(
        [
            "# Standardized 092A Drug Perturbation Summary",
            "Fixed CMap subset: MCF7 + HT-HG-U133A-like + 6 h.",
            "Cosine scores are ranked ascending; more negative means stronger reversal tendency.",
            f"Fig3 manuscript benchmark: {MANUSCRIPT_COMPOUND_COUNT} drugs, {MANUSCRIPT_POSITIVE_COUNT} positives, and {MANUSCRIPT_NEGATIVE_COUNT} non_oncology + unclassified negatives.",
            f"Fig4/Fig5 manuscript cluster: {MANUSCRIPT_CLUSTER_ID}.",
            f"Fig5 operational 80% threshold is inclusive (>=); frozen reference: {MANUSCRIPT_FIG5_CANDIDATE_COUNT} candidates.",
        ]
    )

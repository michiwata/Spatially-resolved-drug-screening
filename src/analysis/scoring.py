from __future__ import annotations

import math
from collections import defaultdict

import numpy as np
from numpy.typing import NDArray

from .tables import sort_numeric_text
from .types import FROZEN_LEGACY_NEGATIVE_ALIASES, NEGATIVE_LABELS, POSITIVE_LABELS, CsvRow, Signature

META_COLUMNS = [
    "drug_key",
    "InChIKey",
    "standardized_cmap_name",
    "preferred_drug_name",
    "drug_label_4class",
    "drug_label_superclass",
    "analysis_role",
    "manual_review_status",
    "evidence_source_type",
]


def cosine(vector: NDArray[np.float64], row: NDArray[np.float64]) -> float:
    denom = max(float(np.linalg.norm(vector) * np.linalg.norm(row)), 1e-12)
    return float(row @ vector / denom)


def score_conditions(
    bank: list[CsvRow], gene_cols: list[str], signatures: list[Signature]
) -> list[CsvRow]:
    rows: list[CsvRow] = []
    for signature in signatures:
        for condition in bank:
            values = np.array(
                [float(condition[gene]) for gene in gene_cols], dtype=np.float64
            )
            output = {
                key: condition[key]
                for key in ["condition_id", *META_COLUMNS, "source_instance_count"]
            }
            output.update(
                {
                    "signature_mode": signature.mode,
                    "cluster_id": signature.cluster_id,
                    "score_cosine": str(cosine(signature.vector, values)),
                    "n_signature_spots": str(signature.n_spots),
                }
            )
            rows.append(output)
    return rows


def aggregate_drug_scores(
    condition_scores: list[CsvRow],
) -> tuple[list[CsvRow], list[CsvRow]]:
    grouped: dict[tuple[str, ...], list[CsvRow]] = defaultdict(list)
    for row in condition_scores:
        grouped[
            tuple(row[key] for key in ["signature_mode", "cluster_id", *META_COLUMNS])
        ].append(row)
    whole: list[CsvRow] = []
    clusters: list[CsvRow] = []
    for key, rows in sorted(grouped.items()):
        output = {
            name: value
            for name, value in zip(
                ["signature_mode", "cluster_id", *META_COLUMNS], key, strict=True
            )
        }
        output["score_cosine"] = str(
            sum(float(row["score_cosine"]) for row in rows) / len(rows)
        )
        output["n_conditions_aggregated"] = str(len(rows))
        output["raw_rows_aggregated"] = str(
            sum(int(row["source_instance_count"]) for row in rows)
        )
        if output["signature_mode"] == "whole_slice_mean":
            whole.append(output)
        else:
            clusters.append(output)
    return whole, clusters


def spatial_best_scores(
    cluster_scores: list[CsvRow], retained_ids: set[str]
) -> list[CsvRow]:
    grouped: dict[str, list[CsvRow]] = defaultdict(list)
    for row in cluster_scores:
        if row["cluster_id"] in retained_ids:
            grouped[row["drug_key"]].append(row)
    output_rows: list[CsvRow] = []
    for drug_key, rows in sorted(grouped.items()):
        best = min(
            rows, key=lambda row: (float(row["score_cosine"]), row["cluster_id"])
        )
        output = {key: best[key] for key in META_COLUMNS}
        output.update(
            {
                "drug_key": drug_key,
                "signature_mode": best["signature_mode"],
                "best_cluster_id": best["cluster_id"],
                "spatial_best_score_cosine": best["score_cosine"],
                "n_conditions_aggregated": best["n_conditions_aggregated"],
                "raw_rows_aggregated": best["raw_rows_aggregated"],
            }
        )
        output_rows.append(output)
    return output_rows


def positive_rank_changes(
    whole_scores: list[CsvRow], cluster_scores: list[CsvRow], best_cluster_id: str
) -> list[CsvRow]:
    whole_pos = [
        row for row in whole_scores if row["drug_label_4class"] in POSITIVE_LABELS
    ]
    cluster_pos = [
        row
        for row in cluster_scores
        if row["cluster_id"] == best_cluster_id
        and row["drug_label_4class"] in POSITIVE_LABELS
    ]
    whole_rank = {
        row["drug_key"]: index
        for index, row in enumerate(
            sorted(
                whole_pos, key=lambda row: (float(row["score_cosine"]), row["drug_key"])
            ),
            start=1,
        )
    }
    cluster_ranked = sorted(
        cluster_pos, key=lambda row: (float(row["score_cosine"]), row["drug_key"])
    )
    rows: list[CsvRow] = []
    for cluster_rank, row in enumerate(cluster_ranked, start=1):
        whole = whole_rank[row["drug_key"]]
        output = dict(row)
        output.update(
            {
                "whole_rank": str(whole),
                "cluster_rank": str(cluster_rank),
                "rank_within_cluster": str(cluster_rank),
                "rank_diff": str(whole - cluster_rank),
            }
        )
        rows.append(output)
    return rows


def roc_auc(y_true: list[int], scores: list[float]) -> float:
    pos = sum(y_true)
    neg = len(y_true) - pos
    if pos == 0 or neg == 0:
        return math.nan
    ordered = sorted(zip(scores, y_true, strict=True), key=lambda item: item[0])
    rank_sum = sum(index for index, item in enumerate(ordered, start=1) if item[1] == 1)
    return (rank_sum - pos * (pos + 1) / 2) / (pos * neg)


def pr_auc(y_true: list[int], scores: list[float]) -> float:
    pos = sum(y_true)
    if pos == 0:
        return math.nan
    ordered = sorted(
        zip(scores, y_true, strict=True), key=lambda item: item[0], reverse=True
    )
    hits = 0
    prev_recall = 0.0
    area = 0.0
    for index, item in enumerate(ordered, start=1):
        hits += item[1]
        recall = hits / pos
        precision = hits / index
        area += (recall - prev_recall) * precision
        prev_recall = recall
    return area


def benchmark_metrics(
    whole_scores: list[CsvRow], cluster_scores: list[CsvRow], scope_rows: list[CsvRow]
) -> list[CsvRow]:
    tiny = {row["cluster_id"]: row["is_tiny_cluster"] == "True" for row in scope_rows}
    retained = {
        row["cluster_id"]: row["retained_for_spatial_reporting"] == "True"
        for row in scope_rows
    }
    groups = [("whole_slice", "whole_slice", whole_scores)] + [
        (
            f"cluster_{cid}",
            cid,
            [row for row in cluster_scores if row["cluster_id"] == cid],
        )
        for cid in sort_numeric_text({row["cluster_id"] for row in cluster_scores})
    ]
    metrics: list[CsvRow] = []
    for curve_name, cluster_id, rows in groups:
        eligible = [
            row
            for row in rows
            if row["drug_label_4class"] in POSITIVE_LABELS | NEGATIVE_LABELS
            or row["drug_label_4class"] in FROZEN_LEGACY_NEGATIVE_ALIASES
        ]
        y_true = [
            1 if row["drug_label_4class"] in POSITIVE_LABELS else 0 for row in eligible
        ]
        eval_scores = [-float(row["score_cosine"]) for row in eligible]
        roc = roc_auc(y_true, eval_scores)
        pr = pr_auc(y_true, eval_scores)
        curve_group = "whole_slice" if cluster_id == "whole_slice" else "cluster"
        metrics.append(
            {
                "eligible_count": str(len(eligible)),
                "positive_count": str(sum(y_true)),
                "negative_count": str(len(y_true) - sum(y_true)),
                "roc_auc": str(roc),
                "pr_auc": str(pr),
                "curve_name": curve_name,
                "curve_group": curve_group,
                "cluster_id": cluster_id,
                "retained_for_selection": str(retained.get(cluster_id, True)),
                "is_tiny_cluster": str(tiny.get(cluster_id, False)),
                "selection_score": str((roc + pr) / 2),
            }
        )
    return metrics


def best_retained_cluster(metrics: list[CsvRow]) -> str:
    eligible = [
        row
        for row in metrics
        if row["curve_group"] == "cluster" and row["retained_for_selection"] == "True"
    ]
    return max(
        eligible,
        key=lambda row: (
            float(row["selection_score"]),
            float(row["roc_auc"]),
            float(row["pr_auc"]),
            -int(row["cluster_id"]),
        ),
    )["cluster_id"]

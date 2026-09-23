from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

from .tables import sort_numeric_text
from .types import MIN_RETAINED_SPOTS, CsvRow, Signature

SIGNATURE_STD_FLOOR = 1e-12


def _gene_standardize(vector: NDArray[np.float64]) -> NDArray[np.float64]:
    values: list[np.float64] = list(vector)
    mean = sum(values) / len(values)
    centered = [value - mean for value in values]
    scale = max(
        math.sqrt(sum(value * value for value in centered) / len(centered)),
        SIGNATURE_STD_FLOOR,
    )
    return np.array([value / scale for value in centered], dtype=np.float64)


def slice_scope_rows(wide_rows: list[CsvRow]) -> list[CsvRow]:
    cluster_ids = sort_numeric_text({row["cluster"] for row in wide_rows})
    total = len(wide_rows)
    rows = [
        {
            "analysis_scope": "whole_slice",
            "cluster_id": "whole_slice",
            "n_spots": str(total),
            "spot_fraction": "1.0",
            "is_tiny_cluster": "False",
            "retained_for_spatial_reporting": "True",
        }
    ]
    for cluster_id in cluster_ids:
        count = sum(1 for row in wide_rows if row["cluster"] == cluster_id)
        fraction = count / total
        tiny = count < MIN_RETAINED_SPOTS
        rows.append(
            {
                "analysis_scope": f"cluster_{cluster_id}",
                "cluster_id": cluster_id,
                "n_spots": str(count),
                "spot_fraction": str(fraction),
                "is_tiny_cluster": str(tiny),
                "retained_for_spatial_reporting": str(not tiny),
            }
        )
    return rows


def retained_cluster_ids(scope_rows: list[CsvRow]) -> set[str]:
    return {
        row["cluster_id"]
        for row in scope_rows
        if row["analysis_scope"].startswith("cluster_")
        and row["retained_for_spatial_reporting"] == "True"
    }


def disease_signatures(
    wide_rows: list[CsvRow], gene_cols: list[str]
) -> list[Signature]:
    raw_counts = np.array(
        [[float(row[gene]) for gene in gene_cols] for row in wide_rows],
        dtype=np.float64,
    )
    matrix = np.log1p(raw_counts)
    whole_vector = _gene_standardize(
        np.asarray(matrix.mean(axis=0), dtype=np.float64)
    )
    signatures = [
        Signature(
            mode="whole_slice_mean",
            cluster_id="whole_slice",
            n_spots=len(wide_rows),
            vector=whole_vector,
        )
    ]
    for cluster_id in sort_numeric_text({row["cluster"] for row in wide_rows}):
        indexes = [
            index for index, row in enumerate(wide_rows) if row["cluster"] == cluster_id
        ]
        rest_indexes = [
            index for index, row in enumerate(wide_rows) if row["cluster"] != cluster_id
        ]
        cluster_matrix = matrix[indexes]
        rest_matrix = matrix[rest_indexes]
        vector = _gene_standardize(
            np.asarray(
                cluster_matrix.mean(axis=0) - rest_matrix.mean(axis=0),
                dtype=np.float64,
            )
        )
        signatures.append(
            Signature(
                mode="cluster_mean_vs_rest",
                cluster_id=cluster_id,
                n_spots=len(indexes),
                vector=vector,
            )
        )
    return signatures

from __future__ import annotations

import re
from collections import defaultdict

from .types import NEGATIVE_LABELS, POSITIVE_LABELS, CsvRow, EmptySelectionError


def clean_text(value: str) -> str:
    text = value.strip()
    return (
        "<missing>"
        if text.lower() in {"", "<missing>", "-", "na", "n/a", "nan", "none", "null"}
        else text
    )


def normalize_platform(value: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "", clean_text(value).lower())
    if "hthgu133a" in normalized or "hgu133a" in normalized:
        return "HT-HG-U133A-like"
    if "u133plus2" in normalized:
        return "HG-U133_Plus_2-like"
    return clean_text(value)


def normalize_duration(value: str) -> str:
    return (
        "6 h"
        if clean_text(value).lower().replace(" ", "") == "6h"
        else clean_text(value)
    )


def drug_key(row: CsvRow) -> str:
    inchikey = clean_text(row["InChIKey"])
    return (
        inchikey
        if inchikey != "<missing>"
        else f"missing::{clean_text(row['cmap_name'])}"
    )


def standardize_label(value: str) -> str:
    labels = {
        "breast_cancer": "breast_cancer_related",
        "oncology_other": "other_oncology_related",
        "non_oncology": "non_oncology",
        "unknown": "unclassified",
        "unclassified": "unclassified",
    }
    return labels.get(clean_text(value), "unclassified")


def analysis_role(label: str) -> str:
    if label in POSITIVE_LABELS:
        return "positive"
    if label in NEGATIVE_LABELS:
        return "negative"
    return "excluded"


def superclass(label: str) -> str:
    return "oncology_related" if label in POSITIVE_LABELS else "non_oncology_related"


def build_condition_bank(
    drug_rows: list[CsvRow], label_rows: list[CsvRow], gene_cols: list[str]
) -> list[CsvRow]:
    fixed_rows = [
        row
        for row in drug_rows
        if row["cell"] == "MCF7"
        and normalize_duration(row["duration"]) == "6 h"
        and normalize_platform(row["array"]) == "HT-HG-U133A-like"
    ]
    if not fixed_rows:
        raise EmptySelectionError(name="fixed CMap subset")
    label_by_key = {row["drug_key"]: row for row in label_rows}
    grouped: dict[tuple[str, ...], list[CsvRow]] = defaultdict(list)
    for row in fixed_rows:
        key = (
            drug_key(row),
            clean_text(row["InChIKey"]),
            clean_text(row["cmap_name"]),
            clean_text(row["cell"]),
            normalize_duration(row["duration"]),
            normalize_platform(row["array"]),
            clean_text(row["concentration"]),
        )
        grouped[key].append(row)
    bank: list[CsvRow] = []
    for index, key in enumerate(sorted(grouped), start=1):
        rows = grouped[key]
        label_row = label_by_key.get(key[0], {})
        label = standardize_label(label_row.get("final_label_4class", "unclassified"))
        output = {
            "condition_id": f"fixed_condition_{index:05d}",
            "drug_key": key[0],
            "InChIKey": label_row.get("InChIKey", key[1]),
            "standardized_cmap_name": label_row.get("standardized_cmap_name", key[2]),
            "preferred_drug_name": label_row.get("original_name", key[2]),
            "drug_label_4class": label,
            "drug_label_superclass": superclass(label),
            "analysis_role": analysis_role(label),
            "manual_review_status": label_row.get("manual_review_status", "unresolved"),
            "evidence_source_type": label_row.get("evidence_source_type", "unresolved"),
            "source_instance_count": str(len(rows)),
            "condition_cell": key[3],
            "condition_duration": key[4],
            "condition_array": key[5],
            "condition_concentration": key[6],
        }
        for gene in gene_cols:
            output[gene] = str(sum(float(row[gene]) for row in rows) / len(rows))
        bank.append(output)
    return bank


def label_count_rows(bank: list[CsvRow]) -> list[CsvRow]:
    rows: list[CsvRow] = []
    for label in [
        "breast_cancer_related",
        "other_oncology_related",
        "non_oncology",
        "unclassified",
    ]:
        keys = {row["drug_key"] for row in bank if row["drug_label_4class"] == label}
        conditions = [row for row in bank if row["drug_label_4class"] == label]
        rows.append(
            {
                "drug_label_4class": label,
                "drug_label_superclass": superclass(label),
                "unique_drug_count": str(len(keys)),
                "condition_count": str(len(conditions)),
                "raw_instance_count": str(
                    sum(int(row["source_instance_count"]) for row in conditions)
                ),
            }
        )
    return [row for row in rows if row["unique_drug_count"] != "0"]


def overview_rows(bank: list[CsvRow], raw_count: int) -> list[CsvRow]:
    labels = {
        label: {row["drug_key"] for row in bank if row["drug_label_4class"] == label}
        for label in POSITIVE_LABELS | NEGATIVE_LABELS
    }
    return [
        {"category": "raw_fixed_subset_instances", "count": str(raw_count)},
        {"category": "deduplicated_fixed_conditions", "count": str(len(bank))},
        {
            "category": "unique_drugs_all",
            "count": str(len({row["drug_key"] for row in bank})),
        },
        {
            "category": "positive_drugs",
            "count": str(sum(len(labels[label]) for label in POSITIVE_LABELS)),
        },
        {"category": "negative_drugs", "count": str(len(labels["non_oncology"]))},
        {"category": "unclassified_drugs", "count": str(len(labels["unclassified"]))},
    ]

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Final, TypeAlias, override

import numpy as np
from numpy.typing import NDArray

CsvRow: TypeAlias = dict[str, str]
CsvCell: TypeAlias = str | int | float | bool
OutputRow: TypeAlias = Mapping[str, CsvCell]

POSITIVE_LABELS: Final = frozenset({"breast_cancer_related", "other_oncology_related"})
NEGATIVE_LABELS: Final = frozenset({"non_oncology", "unclassified"})
FROZEN_LEGACY_NEGATIVE_ALIASES: Final = frozenset({"unknown"})
MIN_RETAINED_SPOTS: Final = 64
MANUSCRIPT_CLUSTER_ID: Final = "3"
MANUSCRIPT_COMPOUND_COUNT: Final = 1291
MANUSCRIPT_POSITIVE_COUNT: Final = 142
MANUSCRIPT_NEGATIVE_COUNT: Final = 1149
MANUSCRIPT_FIG5_CANDIDATE_COUNT: Final = 15


class AnalysisError(Exception):
    pass


@dataclass(frozen=True, slots=True)
class MissingColumnError(AnalysisError):
    path: Path
    columns: tuple[str, ...]

    @override
    def __str__(self) -> str:
        return f"{self.path}: missing required columns: {', '.join(self.columns)}"


@dataclass(frozen=True, slots=True)
class EmptySelectionError(AnalysisError):
    name: str

    @override
    def __str__(self) -> str:
        return f"{self.name}: selection produced no rows"


@dataclass(frozen=True, slots=True)
class Table:
    path: Path
    header: tuple[str, ...]
    rows: list[CsvRow]


@dataclass(frozen=True, slots=True)
class StandardizedAnalysisConfig:
    drug_matrix: Path
    drug_labels: Path
    wide_tsv: Path
    spot_obs: Path
    cluster_color: Path
    outdir: Path

    @classmethod
    def from_toml(cls, path: Path) -> StandardizedAnalysisConfig:
        values = read_toml_section(path, "analysis")
        base = path.parent
        return cls(
            drug_matrix=base / values["drug_matrix"],
            drug_labels=base / values["drug_labels"],
            wide_tsv=base / values["wide_tsv"],
            spot_obs=base / values["spot_obs"],
            cluster_color=base / values["cluster_color"],
            outdir=base / values["outdir"],
        )


@dataclass(frozen=True, slots=True)
class FigureConfig:
    analysis_dir: Path
    spot_obs: Path
    cluster_color: Path
    outdir: Path
    figure1_static: Path
    tissue_image: Path

    @classmethod
    def from_toml(cls, path: Path) -> FigureConfig:
        values = read_toml_section(path, "figures")
        base = path.parent
        return cls(
            analysis_dir=base / values["analysis_dir"],
            spot_obs=base / values["spot_obs"],
            cluster_color=base / values["cluster_color"],
            outdir=base / values["outdir"],
            figure1_static=base / values["figure1_static"],
            tissue_image=base / values["tissue_image"],
        )


@dataclass(frozen=True, slots=True)
class Signature:
    mode: str
    cluster_id: str
    n_spots: int
    vector: NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class ThresholdSummary:
    best_cluster_id: str
    threshold_mode: str
    threshold_definition: str
    threshold_value: float
    positive_response_max: float
    top50_candidate_count: int
    above_threshold_count: int
    followup_flag: str
    comparison: str
    manuscript_wording: str


@dataclass(frozen=True, slots=True)
class Fig5Candidates:
    top50_candidates: list[CsvRow]
    above_threshold: list[CsvRow]
    threshold_summary: ThresholdSummary


def read_toml_section(path: Path, section_name: str) -> dict[str, str]:
    section = ""
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped.removeprefix("[").removesuffix("]")
            continue
        if section == section_name and "=" in stripped:
            key, raw_value = stripped.split("=", maxsplit=1)
            values[key.strip()] = raw_value.strip().strip('"')
    return values

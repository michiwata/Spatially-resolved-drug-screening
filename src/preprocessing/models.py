from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import NewType, Protocol, TypeAlias

import numpy as np
from numpy.typing import NDArray

Barcode = NewType("Barcode", str)
ClusterId = NewType("ClusterId", str)
EntrezId = NewType("EntrezId", int)
JsonScalar: TypeAlias = str | int | float | bool | None
JsonMap: TypeAlias = dict[str, JsonScalar]


@dataclass(frozen=True, slots=True)
class SamplePaths:
    sample_id: str
    h5: Path
    positions: Path
    scalefactors: Path
    hires_image: Path | None


@dataclass(frozen=True, slots=True)
class GeneFeature:
    feature_id: str
    symbol: str


@dataclass(frozen=True, slots=True)
class CountMatrix:
    barcodes: tuple[Barcode, ...]
    features: tuple[GeneFeature, ...]
    counts: NDArray[np.int64]


@dataclass(frozen=True, slots=True)
class PositionRecord:
    barcode: Barcode
    in_tissue: int
    array_row: int
    array_col: int
    pxl_row_in_fullres: float
    pxl_col_in_fullres: float


@dataclass(frozen=True, slots=True)
class ScaleFactors:
    tissue_hires_scalef: float | None


@dataclass(frozen=True, slots=True)
class SpatialRecord:
    sample_id: str
    barcode: Barcode
    in_tissue: int
    array_row: int
    array_col: int
    pxl_row_in_fullres: float
    pxl_col_in_fullres: float
    pxl_row_in_hires: float
    pxl_col_in_hires: float
    cluster: ClusterId
    total_counts: float
    n_genes_by_counts: int
    pct_counts_mt: float


@dataclass(frozen=True, slots=True)
class ClusterRequest:
    matrix: CountMatrix
    barcodes: tuple[Barcode, ...]
    resolution: float
    random_state: int


class Clusterer(Protocol):
    def assign(self, request: ClusterRequest) -> dict[Barcode, ClusterId]: ...


@dataclass(frozen=True, slots=True)
class SpatialBuildInputs:
    matrix: CountMatrix
    positions: tuple[PositionRecord, ...]
    scalefactors: ScaleFactors
    sample_id: str
    clusterer: Clusterer
    in_tissue_only: bool
    resolution: float
    random_state: int


@dataclass(frozen=True, slots=True)
class GeneMaps:
    ensembl_to_entrez: dict[str, EntrezId]
    symbol_to_entrez: dict[str, EntrezId]
    discontinued_to_current: dict[EntrezId, EntrezId]


@dataclass(frozen=True, slots=True)
class CMapGeneColumn:
    name: str
    entrez_id: EntrezId


@dataclass(frozen=True, slots=True)
class AggregateStats:
    n_genes_raw: int
    n_genes_with_entrez: int
    n_genes_in_cmap_after_map: int
    n_unique_entrez_in_cmap: int
    n_duplicate_features_collapsed: int
    n_nonzero_cmap_counts: int


@dataclass(frozen=True, slots=True)
class AggregatedCounts:
    gene_columns: list[str]
    counts: NDArray[np.int64]
    stats: AggregateStats


@dataclass(frozen=True, slots=True)
class ExportConfig:
    sample: SamplePaths
    matrix: CountMatrix
    cmap_signature_path: Path
    mapping_dir: Path
    out_root: Path
    clusterer: Clusterer
    in_tissue_only: bool = True
    leiden_resolution: float = 0.6
    random_state: int = 0


@dataclass(frozen=True, slots=True)
class ExportResult:
    wide_tsv: Path
    spot_obs_tsv: Path
    report_json: Path
    stats: AggregateStats

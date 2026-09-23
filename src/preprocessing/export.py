from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Final

import numpy as np
from numpy.typing import NDArray

from .errors import InputSchemaError
from .genes import aggregate_counts_to_cmap, map_features_to_entrez
from .io import load_gene_maps, load_positions, load_scalefactors
from .models import (
    AggregateStats,
    CMapGeneColumn,
    EntrezId,
    ExportConfig,
    ExportResult,
    SpatialBuildInputs,
    SpatialRecord,
)
from .spatial import build_spatial_table

META_COLUMNS: Final = [
    "sample_id",
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
    "pxl_row_in_hires",
    "pxl_col_in_hires",
    "cluster",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
]


def cmap_gene_columns_from_header(path: Path) -> list[CMapGeneColumn]:
    with path.open(encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    columns = [_cmap_column(path, name) for name in header if name.startswith("hsa:")]
    if not columns:
        raise InputSchemaError(path=path, reason="header has no hsa:<EntrezID> CMap gene columns")
    return columns


def export_cmapwide(config: ExportConfig) -> ExportResult:
    cmap_columns = cmap_gene_columns_from_header(config.cmap_signature_path)
    gene_maps = load_gene_maps(config.mapping_dir)
    spatial_rows = build_spatial_table(SpatialBuildInputs(
        matrix=config.matrix,
        positions=load_positions(config.sample.positions),
        scalefactors=load_scalefactors(config.sample.scalefactors),
        sample_id=config.sample.sample_id,
        clusterer=config.clusterer,
        in_tissue_only=config.in_tissue_only,
        resolution=config.leiden_resolution,
        random_state=config.random_state,
    ))
    mapped = map_features_to_entrez(config.matrix.features, gene_maps)
    aggregated = aggregate_counts_to_cmap(config.matrix, mapped, cmap_columns)
    sample_root = config.out_root / config.sample.sample_id
    wide_tsv = sample_root / "wide" / f"{config.sample.sample_id}__spot_by_gene__cmapwide.tsv"
    spot_obs_tsv = sample_root / "spot_obs.tsv"
    report_json = sample_root / "wide_export_report.json"
    wide_tsv.parent.mkdir(parents=True, exist_ok=True)
    _write_wide_tsv(wide_tsv, spatial_rows, aggregated.gene_columns, aggregated.counts)
    _write_spot_obs(spot_obs_tsv, spatial_rows)
    _write_report(report_json, config, spatial_rows, aggregated.stats)
    return ExportResult(wide_tsv=wide_tsv, spot_obs_tsv=spot_obs_tsv, report_json=report_json, stats=aggregated.stats)


def _cmap_column(path: Path, name: str) -> CMapGeneColumn:
    try:
        entrez = EntrezId(int(name.split(":", maxsplit=1)[1]))
    except ValueError as error:
        raise InputSchemaError(path=path, reason=f"invalid CMap gene column {name!r}") from error
    return CMapGeneColumn(name=name, entrez_id=entrez)


def _write_wide_tsv(path: Path, rows: list[SpatialRecord], gene_columns: list[str], counts: NDArray[np.int64]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([*META_COLUMNS, *gene_columns])
        for index, record in enumerate(rows):
            writer.writerow([*_record_values(record), *[str(int(value)) for value in counts[index, :]]])


def _write_spot_obs(path: Path, rows: list[SpatialRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(META_COLUMNS)
        for record in rows:
            writer.writerow(_record_values(record))


def _write_report(path: Path, config: ExportConfig, rows: list[SpatialRecord], stats: AggregateStats) -> None:
    path.write_text(
        json.dumps({
            **asdict(stats),
            "sample_id": config.sample.sample_id,
            "n_spots_exported": len(rows),
            "in_tissue_only": config.in_tissue_only,
            "leiden_resolution": config.leiden_resolution,
            "random_state": config.random_state,
        }, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _record_values(record: SpatialRecord) -> list[str | int | float]:
    return [
        record.sample_id,
        str(record.barcode),
        record.in_tissue,
        record.array_row,
        record.array_col,
        record.pxl_row_in_fullres,
        record.pxl_col_in_fullres,
        record.pxl_row_in_hires,
        record.pxl_col_in_hires,
        str(record.cluster),
        record.total_counts,
        record.n_genes_by_counts,
        record.pct_counts_mt,
    ]

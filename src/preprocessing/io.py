from __future__ import annotations

import csv
import gzip
import json
from collections.abc import Iterable
from pathlib import Path
from typing import TextIO, TypeAlias, TypeGuard

import numpy as np

from .errors import InputSchemaError, OptionalDependencyError
from .models import (
    Barcode,
    CountMatrix,
    EntrezId,
    GeneFeature,
    GeneMaps,
    JsonScalar,
    PositionRecord,
    SamplePaths,
    ScaleFactors,
)

JsonMap: TypeAlias = dict[str, JsonScalar]

H5_SUFFIX = "_filtered_feature_bc_matrix.h5"
POSITION_COLUMNS = ("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")


def find_samples(raw_root: Path) -> list[SamplePaths]:
    root = raw_root.expanduser().resolve()
    samples: list[SamplePaths] = []
    for h5_path in sorted(root.glob(f"*{H5_SUFFIX}")):
        sample_id = h5_path.name.removesuffix(H5_SUFFIX)
        samples.append(SamplePaths(
            sample_id=sample_id,
            h5=h5_path,
            positions=_first_existing(root, sample_id, "_tissue_positions_list.csv", zipped=True),
            scalefactors=_first_existing(root, sample_id, "_scalefactors_json.json", zipped=True),
            hires_image=_optional_existing(root, sample_id, "_tissue_hires_image.png", zipped=True),
        ))
    return samples


def load_positions(path: Path) -> tuple[PositionRecord, ...]:
    with _open_text(path) as handle:
        first = handle.readline().rstrip("\n\r")
        handle.seek(0)
        reader = csv.reader(handle)
        rows = list(reader)
    data_rows = rows[1:] if _has_header(first) else rows
    return tuple(_position_from_cells(path, row) for row in data_rows)


def load_scalefactors(path: Path) -> ScaleFactors:
    with _open_text(path) as handle:
        payload: JsonScalar | JsonMap = json.load(handle)
    if not _is_json_map(payload):
        raise InputSchemaError(path=path, reason="scalefactors JSON must be an object")
    value = _json_number(payload.get("tissue_hires_scalef"))
    fallback = _json_number(payload.get("tissue_hires_scalefactor"))
    return ScaleFactors(tissue_hires_scalef=value if value is not None else fallback)


def load_counts_tsv(path: Path) -> CountMatrix:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) < 2 or header[0] != "barcode":
            raise InputSchemaError(path=path, reason="counts TSV header must start with barcode and feature columns")
        barcodes: list[Barcode] = []
        values: list[list[int]] = []
        for row in reader:
            if len(row) != len(header):
                raise InputSchemaError(path=path, reason=f"row for {row[0] if row else '<empty>'} has wrong column count")
            barcodes.append(Barcode(row[0]))
            values.append([int(cell) for cell in row[1:]])
    features = tuple(GeneFeature(feature_id=name, symbol=name) for name in header[1:])
    return CountMatrix(barcodes=tuple(barcodes), features=features, counts=np.array(values, dtype=np.int64))


def load_visium_h5(path: Path) -> CountMatrix:
    try:
        scanpy = __import__("scanpy")
    except ModuleNotFoundError as error:
        raise OptionalDependencyError(package="scanpy", purpose="load 10x Visium H5 matrices") from error
    adata = scanpy.read_10x_h5(str(path))
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    feature_ids = tuple(str(name) for name in adata.var_names)
    symbols = feature_ids
    raw_counts = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    return CountMatrix(
        barcodes=tuple(Barcode(str(name)) for name in adata.obs_names),
        features=tuple(GeneFeature(feature_id=feature_id, symbol=symbol) for feature_id, symbol in zip(feature_ids, symbols, strict=True)),
        counts=np.asarray(raw_counts, dtype=np.int64),
    )


def load_gene_maps(mapping_dir: Path) -> GeneMaps:
    root = mapping_dir.expanduser().resolve()
    ensembl_to_entrez: dict[str, EntrezId] = {}
    symbol_to_entrez: dict[str, EntrezId] = {}
    discontinued_to_current: dict[EntrezId, EntrezId] = {}
    _read_hgnc(root / "hgnc_complete_set.txt", ensembl_to_entrez, symbol_to_entrez)
    _read_gene2ensembl(root / "gene2ensembl.gz", ensembl_to_entrez)
    _read_gene_history(root / "gene_history.gz", discontinued_to_current)
    return GeneMaps(ensembl_to_entrez=ensembl_to_entrez, symbol_to_entrez=symbol_to_entrez, discontinued_to_current=discontinued_to_current)


def _first_existing(root: Path, sample_id: str, suffix: str, *, zipped: bool) -> Path:
    gz_path = root / f"{sample_id}{suffix}.gz"
    plain_path = root / f"{sample_id}{suffix}"
    if zipped and gz_path.exists():
        return gz_path
    return plain_path


def _optional_existing(root: Path, sample_id: str, suffix: str, *, zipped: bool) -> Path | None:
    path = _first_existing(root, sample_id, suffix, zipped=zipped)
    return path if path.exists() else None


def _open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def _has_header(first_line: str) -> bool:
    cells = first_line.split(",")
    return "barcode" in cells or "in_tissue" in cells


def _position_from_cells(path: Path, row: list[str]) -> PositionRecord:
    if len(row) != len(POSITION_COLUMNS):
        raise InputSchemaError(path=path, reason="positions row must have six 092A Visium columns")
    return PositionRecord(
        barcode=Barcode(row[0]),
        in_tissue=int(row[1]),
        array_row=int(row[2]),
        array_col=int(row[3]),
        pxl_row_in_fullres=float(row[4]),
        pxl_col_in_fullres=float(row[5]),
    )


def _json_number(value: JsonScalar | None) -> float | None:
    return float(value) if isinstance(value, int | float) and not isinstance(value, bool) else None


def _is_json_map(value: JsonScalar | JsonMap) -> TypeGuard[JsonMap]:
    return isinstance(value, dict)


def _read_hgnc(path: Path, ensembl: dict[str, EntrezId], symbols: dict[str, EntrezId]) -> None:
    if not path.exists():
        return
    for row in _dict_rows(path):
        entrez = _parse_entrez(row.get("entrez_id", ""))
        if entrez is None:
            continue
        if row.get("ensembl_gene_id", "").strip():
            ensembl[row["ensembl_gene_id"].strip()] = entrez
        if row.get("symbol", "").strip():
            symbols[row["symbol"].strip()] = entrez


def _read_gene2ensembl(path: Path, ensembl: dict[str, EntrezId]) -> None:
    if not path.exists():
        return
    for row in _dict_rows(path):
        entrez = _parse_entrez(row.get("GeneID", ""))
        gene_id = row.get("Ensembl_gene_identifier", "").strip()
        if row.get("tax_id", "") == "9606" and entrez is not None and gene_id:
            ensembl.setdefault(gene_id, entrez)


def _read_gene_history(path: Path, history: dict[EntrezId, EntrezId]) -> None:
    if not path.exists():
        return
    for row in _dict_rows(path):
        current = _parse_entrez(row.get("GeneID", ""))
        old = _parse_entrez(row.get("Discontinued_GeneID", ""))
        if current is not None and old is not None:
            history[old] = current


def _dict_rows(path: Path) -> Iterable[dict[str, str]]:
    with _open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _parse_entrez(value: str) -> EntrezId | None:
    text = value.strip()
    if text == "" or text.lower() == "nan":
        return None
    try:
        return EntrezId(int(float(text)))
    except ValueError:
        return None

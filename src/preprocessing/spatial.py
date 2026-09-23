from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import NDArray

from .errors import BarcodeMismatchError, OptionalDependencyError
from .models import (
    Barcode,
    ClusterId,
    ClusterRequest,
    CountMatrix,
    PositionRecord,
    SpatialBuildInputs,
    SpatialRecord,
)

MT_PREFIX: Final = "MT-"


@dataclass(frozen=True, slots=True)
class StaticClusterer:
    assignments: dict[str, str]

    def assign(self, request: ClusterRequest) -> dict[Barcode, ClusterId]:
        return {barcode: ClusterId(self.assignments[str(barcode)]) for barcode in request.barcodes}


class OrderedFixtureClusterer:
    def assign(self, request: ClusterRequest) -> dict[Barcode, ClusterId]:
        return {barcode: ClusterId(str(index)) for index, barcode in enumerate(request.barcodes)}


class ScanpyLeidenClusterer:
    def assign(self, request: ClusterRequest) -> dict[Barcode, ClusterId]:
        try:
            anndata = __import__("anndata")
            scanpy = __import__("scanpy")
        except ModuleNotFoundError as error:
            raise OptionalDependencyError(package="scanpy anndata", purpose="compute Leiden clusters") from error
        adata = anndata.AnnData(X=request.matrix.counts)
        adata.obs_names = [str(barcode) for barcode in request.matrix.barcodes]
        adata.var_names = [feature.symbol for feature in request.matrix.features]
        scanpy.pp.highly_variable_genes(adata, n_top_genes=min(3000, len(request.matrix.features)), flavor="seurat_v3")
        scanpy.pp.normalize_total(adata, target_sum=10000.0)
        scanpy.pp.log1p(adata)
        adata = adata[:, adata.var["highly_variable"]].copy()
        scanpy.pp.scale(adata, max_value=10)
        scanpy.tl.pca(adata, svd_solver="arpack", random_state=request.random_state)
        scanpy.pp.neighbors(adata, n_neighbors=15, n_pcs=min(30, len(request.matrix.features) - 1), random_state=request.random_state)
        scanpy.tl.leiden(adata, resolution=request.resolution, random_state=request.random_state, flavor="leidenalg", directed=False)
        return {Barcode(str(barcode)): ClusterId(str(cluster)) for barcode, cluster in zip(adata.obs_names, adata.obs["leiden"], strict=True)}


def build_spatial_table(inputs: SpatialBuildInputs) -> list[SpatialRecord]:
    position_by_barcode = {record.barcode: record for record in inputs.positions}
    selected = tuple(
        barcode for barcode in inputs.matrix.barcodes
        if barcode in position_by_barcode and (not inputs.in_tissue_only or position_by_barcode[barcode].in_tissue == 1)
    )
    if not selected:
        raise BarcodeMismatchError(sample_id=inputs.sample_id, reason="no overlapping barcodes between matrix and positions")
    clusters = inputs.clusterer.assign(ClusterRequest(
        matrix=inputs.matrix,
        barcodes=selected,
        resolution=inputs.resolution,
        random_state=inputs.random_state,
    ))
    barcode_to_index = {barcode: index for index, barcode in enumerate(inputs.matrix.barcodes)}
    return [
        _spatial_record(
            matrix=inputs.matrix,
            row_index=barcode_to_index[barcode],
            position=position_by_barcode[barcode],
            sample_id=inputs.sample_id,
            cluster=clusters[barcode],
            hires_scale=inputs.scalefactors.tissue_hires_scalef,
        )
        for barcode in selected
    ]


def _spatial_record(
    matrix: CountMatrix,
    row_index: int,
    position: PositionRecord,
    sample_id: str,
    cluster: ClusterId,
    hires_scale: float | None,
) -> SpatialRecord:
    counts: NDArray[np.int64] = matrix.counts[row_index, :]
    total = float(np.sum(counts))
    mt_count = float(np.sum(counts[_mitochondrial_mask(matrix)]))
    scale = 0.0 if hires_scale is None else hires_scale
    return SpatialRecord(
        sample_id=sample_id,
        barcode=position.barcode,
        in_tissue=position.in_tissue,
        array_row=position.array_row,
        array_col=position.array_col,
        pxl_row_in_fullres=position.pxl_row_in_fullres,
        pxl_col_in_fullres=position.pxl_col_in_fullres,
        pxl_row_in_hires=position.pxl_row_in_fullres * scale,
        pxl_col_in_hires=position.pxl_col_in_fullres * scale,
        cluster=cluster,
        total_counts=total,
        n_genes_by_counts=int(np.count_nonzero(counts)),
        pct_counts_mt=0.0 if total == 0.0 else mt_count / total * 100.0,
    )


def _mitochondrial_mask(matrix: CountMatrix) -> NDArray[np.bool_]:
    return np.array([feature.symbol.upper().startswith(MT_PREFIX) for feature in matrix.features], dtype=bool)

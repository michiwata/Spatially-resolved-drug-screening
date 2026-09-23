from __future__ import annotations

from collections import Counter

import numpy as np

from .errors import GeneMappingError
from .models import (
    AggregatedCounts,
    AggregateStats,
    CMapGeneColumn,
    CountMatrix,
    EntrezId,
    GeneFeature,
    GeneMaps,
)


def map_features_to_entrez(features: tuple[GeneFeature, ...], gene_maps: GeneMaps) -> list[EntrezId | None]:
    return [_feature_entrez(feature, gene_maps) for feature in features]


def aggregate_counts_to_cmap(
    matrix: CountMatrix,
    mapped_entrez: list[EntrezId | None],
    cmap_columns: list[CMapGeneColumn],
) -> AggregatedCounts:
    cmap_ids = [column.entrez_id for column in cmap_columns]
    cmap_index = {entrez: index for index, entrez in enumerate(cmap_ids)}
    output = np.zeros((matrix.counts.shape[0], len(cmap_columns)), dtype=np.int64)
    kept: list[EntrezId] = []
    for feature_index, entrez in enumerate(mapped_entrez):
        if entrez is not None and entrez in cmap_index:
            output[:, cmap_index[entrez]] += matrix.counts[:, feature_index]
            kept.append(entrez)
    if not kept:
        raise GeneMappingError(reason="no features mapped to CMap Entrez columns")
    unique_count = len(set(kept))
    stats = AggregateStats(
        n_genes_raw=len(matrix.features),
        n_genes_with_entrez=sum(entrez is not None for entrez in mapped_entrez),
        n_genes_in_cmap_after_map=len(kept),
        n_unique_entrez_in_cmap=unique_count,
        n_duplicate_features_collapsed=sum(count - 1 for count in Counter(kept).values() if count > 1),
        n_nonzero_cmap_counts=int(np.count_nonzero(output)),
    )
    return AggregatedCounts(gene_columns=[column.name for column in cmap_columns], counts=output, stats=stats)


def _feature_entrez(feature: GeneFeature, gene_maps: GeneMaps) -> EntrezId | None:
    ensembl = _strip_ensembl_version(feature.feature_id)
    entrez = gene_maps.ensembl_to_entrez.get(ensembl)
    if entrez is None:
        entrez = gene_maps.symbol_to_entrez.get(feature.symbol)
    if entrez is None:
        return None
    return gene_maps.discontinued_to_current.get(entrez, entrez)


def _strip_ensembl_version(identifier: str) -> str:
    return identifier.split(".", maxsplit=1)[0]

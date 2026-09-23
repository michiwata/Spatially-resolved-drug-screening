# STASCAN: Spatially resolved transcriptomic correlation improves the prediction of anticancer drug candidates

This repository accompanies the final proof, `ci-2026-02416r_Proof_hi.pdf`, by Zetian Song and Michio Iwata. It pairs spatial transcriptomics with drug-induced gene expression signatures to predict anticancer drug candidates.

## Method

The framework compares a bulk disease signature and cluster-specific spatial signatures with Connectivity Map perturbation signatures. A negative cosine score indicates inverse correlation and a reversal candidate. The spatial signatures improve ranking of known anticancer drugs and surface compounds for manual review.

The study uses GEO `GSE210616`, sample `GSM6433585`, section `092A`, with 1,289 in-tissue spots and 12,399 shared `hsa:` genes. Clustering records eight clusters. Clusters with fewer than 64 spots are excluded, leaving clusters 0 to 5 for downstream spatial analysis.

## Logic and order

1. Preprocessing reads Visium counts, applies the spot and cluster rules, and exports the raw-count CMap-wide spot matrix and metadata.
2. The maintained public analysis path applies log1p and per-signature gene standardization, then aggregates the fixed MCF7, HT-HG-U133A-like, 6 h CMap context before cosine scoring.
3. Cosine scores are ranked ascending, with more negative values indicating stronger reversal.
4. Figures 2 to 5 are rendered from the resulting tables. Figure 1 is the preserved static overview PNG.

Two historical disease-signature lineages are documented separately. Lineage A is the archived raw-rebuild matching path represented by the maintained public `analysis/signatures.py` implementation above. Lineage B is the April 20/21 paper-figure rerun, which used raw disease-expression means with cosine scoring. These are distinct historical paths, not sequential stages of one unified run. The full mapping is in `docs/logic_chain.md` and `docs/methods_mapping.md`.

## Inputs and outputs

Inputs are the GEO Visium files, restricted CMap build 02 matrix, gene mappings, and drug labels. Outputs include the raw-count CMap-wide matrix, spot metadata, standardized signatures, drug scores, benchmark metrics, rank tables, and figure source tables. The fixed condition aggregation contains 1,291 compounds.

The proof benchmark contains 142 positive drugs and 1,149 negatives, including `unclassified`. The `unclassified` group is scientifically unresolved and is treated as negative only for this benchmark.

## Commands

Run from the repository root with `PYTHONPATH=src`.

```bash
# Preprocessing: export one 092A sample as a CMap-wide spot table
PYTHONPATH=src python -m st_virtual_screening.preprocess_cli export-wide \
  --raw-root data/external/geo/GSM6433585 \
  --sample-id GSM6433585_092A \
  --cmap-signatures data/restricted/cmap_build_02/BC_gene_expression_signatures_mas5.txt \
  --mapping-dir data/external/gene_mappings \
  --out-root data/generated \
  --leiden-resolution 0.6 \
  --random-state 0

# Standardized analysis
PYTHONPATH=src python -m st_virtual_screening.analyze_cli --config config/analysis.example.toml

# Figure pack
PYTHONPATH=src python -m st_virtual_screening.figures_cli --config config/figures.example.toml
```

## Reference results and figures

Cluster 3 is the fixed downstream cluster used for spatial reporting. Figure 3 reports bulk ROC AUC/AUPR of 0.416/0.096 and cluster 3 values of 0.575/0.146. Figure 4 shows the 12 smallest and 12 largest positive-class `rank_diff` values (24 bars); in the frozen table this is exactly the same drug set as `whole_rank` 1–12 and 131–142. Blue bars indicate drugs spatially higher in cluster 3 and red bars indicate drugs with higher bulk rank. Figure 5 displays the top 20 candidates, uses the inclusive `>=` threshold rule, and retains 15 candidates above threshold.

Nine figure source TSVs plus the static Figure 1 PNG make 10 figure-reference files. Canonical Figure 1 is preserved static artwork supplied for the final manuscript/proof (652×352, 83,701 bytes, SHA-256 `d48459af895967a5f3ad53060dfbaa758c85dce2d83895afd27d978b35b8bfa1`); no generating code is included or claimed for it. The release contains 22 mapped artifacts. See `docs/figure_mapping.md` and `docs/release_inventory.md`.

## Data access and rights

Raw GEO files and the restricted CMap matrix are not bundled. No official public direct-download URL or redistribution right is asserted for the exact CMap file; an institutional acquisition link supplied by the maintainer is recorded. Ownership and licensing for this release remain pending. See `data/README.md` and `LICENSE_PENDING.md`.

## Validation boundary

Lightweight tests and synthetic or frozen rendering checks were run. No full GEO or CMap workflow was rerun, and no environment was reconstructed. The code and reference results support reproduction, but the workflow was not independently rerun in a new environment.

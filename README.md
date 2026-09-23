# Spatially resolved transcriptomic correlation improves the prediction of anticancer drug candidates

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview
This repository contains the Python code for data analysis and figure generation used in the following paper:

> **"Spatially resolved transcriptomic correlation improves the prediction of anticancer drug candidates"**  
> Zetian Song, Keitaro Izumi, and Michio Iwata

## Requirements
This pipeline requires **Python 3.10+** (due to modern type hinting and dataclass features) and relies on standard spatial transcriptomics analysis libraries.

**Key Dependencies:**
- `scanpy` (for Leiden clustering and spatial data handling)
- `anndata` / `h5py` (for reading 10x Visium `.h5` matrices)
- `pandas` / `numpy` / `scipy` (for core data manipulation and scoring)
- `matplotlib` / `seaborn` (for rendering figures)

**Installation:**
It is recommended to use a virtual environment (e.g., `venv` or `conda`). Clone the repository and install the required packages using the provided `requirements.txt`:

```bash
git clone [https://github.com/yourusername/your-repo-name.git](https://github.com/yourusername/your-repo-name.git)
cd your-repo-name
pip install -r requirements.txt
```
*(If you are using Conda, you may prefer to install `scanpy` and its dependencies via `conda install -c conda-forge scanpy` before running pip.)*

## Repository Structure
```text
.
├── data/                       # Input data and reference metadata
│   ├── gene_mappings/          # Downloaded from NCBI and HGNC databases (i.e., hgnc_complete_set.txt, gene2ensembl.gz, and gene_history.gz)
│   ├── raw/
│       ├── CMap/               # Downloaded CMap signature matrices (e.g., mas5.txt)
│       ├── GEO/                # Downloaded raw Visium data (e.g., H5, tissue images)
│   ├── processed/
├── results/                    # Generated outputs and figure renderings
│   ├── analysis/               # Analysis outputs (scores, tables, run_summary.md)
│   └── figures/                # Final rendered manuscript figures
├── src/                        # CLI scripts and pipeline submodules
│   ├── 01_preprocess.py        # CLI for raw data discovery and CMap-wide formatting
│   ├── 02_analyze.py           # CLI for standardized spatial screening analysis
│   ├── 03_plot_results.py      # CLI for figure generation
│   └── analysis/               # Submodules supporting the main pipeline
│       ├── pipeline.py         # Core analysis execution function
│       ├── conditions.py, figures.py, outputs.py, scoring.py, signatures.py, tables.py, types.py
├── analysis.example.toml       # Example configuration for 02_analyze.py
├── figures.example.toml        # Example configuration for 03_plot_results.py
├── requirements.txt            # Python dependencies
└── README.md
```


## Usage / Pipeline
The analysis is divided into three steps. Please execute the scripts in numerical order.

### 1. Data Preprocessing (`01_preprocess.py`)
This script handles the initial preprocessing of spatial transcriptomics (10x Visium) raw data. It maps the spatial coordinates, performs Leiden clustering using Scanpy, and integrates Connectivity Map (CMap) signatures to evaluate each tissue spot. Finally, it exports the processed data as a comprehensive, wide-format spot table (CMap-wide spot table) ready for downstream virtual screening.

The script provides two main subcommands depending on your goal:

* **`audit-raw` (Data Discovery):** 
  Scans a specified local directory (e.g., downloaded GEO data) to automatically discover available Visium RAW samples. It outputs a list of sample IDs along with their essential file paths (H5 expression matrix, spatial positions, and scale factors).
  
* **`export-wide` (Main Processing):** 
  Executes the core preprocessing pipeline for a specific sample. It loads the H5 data and spatial metadata, runs Leiden clustering, applies the specified CMap gene signatures, and exports the fully integrated data as a wide-format text file.
  
**To list discoverable Visium RAW samples in a GEO directory:**
```bash
python -m src.01_preprocess audit-raw \
    --raw-root data/raw_geo_dir/
```
For the preprocessing script to successfully detect and process a Visium sample, the specified raw directory (e.g., `data/raw/GEO/GSM6433585/`) must contain the following standard Space Ranger output files:

*   `*_filtered_feature_bc_matrix.h5` (Expression matrix)
*   `*_tissue_positions_list.csv.gz` or `.csv` (Spot coordinates)
*   `*_scalefactors_json.json.gz` or `.json` (Scale factors)
*   `*_tissue_hires_image.png` (High-resolution tissue image)

**To export a specific sample as a CMap-wide spot table:**
```bash
python -m src.01_preprocess export-wide \
    --raw-root data/raw_geo_dir \
    --sample-id  \
    --cmap-signatures data/cmap_signatures.gmt \
    --mapping-dir data/mapping_dir \
    --out-root results \
    --leiden-resolution 0.6 \
    --random-state 0
```
*(Optional arguments: `--keep-out-of-tissue` to retain out-of-tissue spots, or `--counts-tsv` to use a small fixture dataset instead of the 10x H5 file.)*

**Example Command:**
```bash
python -m src.01_preprocess audit-raw \
    --raw-root data/raw/GEO/GSM6433585/

python -m src.01_preprocess export-wide \
    --raw-root data/raw/GEO/GSM6433585 \
    --sample-id GSM6433585_092A \
    --cmap-signatures data/raw/CMap/BC_gene_expression_signatures_mas5.txt \
    --mapping-dir data/gene_mappings \
    --out-root data/processed \
    --leiden-resolution 0.6 \
    --random-state 0
```

### 2. Main Analysis (`02_analyze.py`)
This script executes the standardized spatial virtual screening pipeline. It calculates disease signatures and drug reversal tendency scores (cosine similarity) for both the whole tissue and specific Leiden clusters. The pipeline also automatically generates the benchmarking metrics and candidate tables required for Manuscript Figures.

To ensure reproducibility, this script is entirely driven by a TOML configuration file rather than complex command-line arguments. This configuration links the CMap-wide spot table (generated in Step 1) with CMap drug reference matrices and spatial metadata.

**To run the analysis:**
```bash
python src/02_analyze.py --config path/to/analysis.example.toml
```

**Example Configuration (`analysis.example.toml`):**
```toml
[analysis]
drug_matrix = "../data/restricted/cmap_build_02/BC_gene_expression_signatures_mas5.txt"
drug_labels = "../data/reference/canonical_relabelled_drug_table.tsv"
wide_tsv = "../data/generated/GSM6433585_092A/wide/GSM6433585_092A__spot_by_gene__cmapwide.tsv"
spot_obs = "../data/reference/GSM6433585/spot_obs.tsv"
cluster_color = "../data/reference/GSM6433585/GSM6433585_092A__cluster_color_map.tsv"
outdir = "../results/output/analysis"
```

*(Note: This pipeline relies on submodules within the `analysis/` directory (e.g., `scoring.py`, `figures.py`). Upon successful execution, it outputs the analysis results, Figure 5 candidate tables, and a `run_summary.md` to the specified `outdir`.)*

**Submodules:**
The main analysis relies on several specialized submodules within the `src/analysis/` directory to modularize the workflow:

*   **`conditions.py`**: Normalizes experimental conditions (e.g., duration, platform) and builds the condition bank for drug perturbations.
*   **`figures.py`**: Prepares and exports the specific candidate tables required for Manuscript Figure 5.
*   **`outputs.py`**: Manages the data structures (`AnalysisOutputs`) and handles the exporting of final analysis results.
*   **`scoring.py`**: The core computational module. It calculates drug reversal scores, aggregates spatial best scores, and computes benchmark metrics and rank changes.
*   **`signatures.py`**: Computes disease signatures based on spatial gene expression and manages retained cluster IDs.
*   **`tables.py`**: Provides utility functions for data I/O, such as reading TSVs, extracting required columns, and writing text files.
*   **`types.py`**: Defines essential data classes (e.g., `StandardizedAnalysisConfig`) and hardcoded constants used in the manuscript (e.g., cluster IDs, compound counts).

### 3. Result Visualization (`03_plot_results.py`)
This script generates the final figures (and related data representations) for the manuscript based on the results from the main analysis pipeline. Similar to the analysis step, it is driven by a TOML configuration file that links the intermediate analysis outputs, spatial metadata, and high-resolution tissue images to the figure generation submodule.

**To generate the figures:**
```bash
python scripts/03_plot_results.py --config path/to/figures.example.toml
```

**Example Configuration (`figures.example.toml`):**
```toml
[figures]
analysis_dir = "../results/output/analysis"
spot_obs = "../data/reference/GSM6433585/spot_obs.tsv"
cluster_color = "../data/reference/GSM6433585/GSM6433585_092A__cluster_color_map.tsv"
outdir = "../results/output/figures"
figure1_static = "../results/reference/figures/fig1_overview.png"
tissue_image = "../data/external/geo/GSM6433585/GSM6433585_092A_tissue_hires_image.png"
```
*(Upon successful execution, the script will render the figure pack and print the path to the specified `outdir` in the console.)*

**Submodules:**
The  analysis also relies on several specialized submodules within the `src/analysis/` directory to modularize the workflow:

*   **`figures.py`**
*   **`types.py`**

## Expected Output
Upon successful execution of the entire pipeline (Steps 1-3), the following key directories and files will be generated in your specified output paths:

*   **1. Preprocessing Outputs (e.g., `data/generated/.../wide/`)**
    *   `*__spot_by_gene__cmapwide.tsv`: The comprehensive spatial expression matrix integrating CMap signatures and Leiden clusters.

*   **2. Analysis Outputs (`results/output/analysis/`)**
    *   `run_summary.md`: An automated summary report detailing the run parameters and manuscript constants.
    *   `whole_slice_drug_scores.tsv` & `cluster_drug_scores.tsv`: The calculated drug reversal tendency (cosine) scores for the whole tissue and per cluster.
    *   `spatial_best_cluster_drug_scores.tsv`: Aggregated best scores across spatial clusters.
    *   Various benchmarking and candidate selection tables (e.g., Figure 5 top-50 candidates).

*   **3. Figure Outputs (`results/output/figures/`)**
    *   Rendered image files (e.g., PNG/PDF formats) representing the manuscript figures (Figures 1-5), including spatial cluster maps, ROC/PR curves, and drug rank-change plots.

## Citation
To be updated.

## License
This project is licensed under the [MIT License](LICENSE).

## Contact
If you have any questions or bug reports, please contact Michio Iwata at Iwata.michio121@mail.kyutech.jp.

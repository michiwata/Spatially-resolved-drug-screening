# Spatially resolved transcriptomic correlation improves the prediction of anticancer drug candidates

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview
This repository contains the Python code for data analysis used in the following paper:

> **"Spatially resolved transcriptomic correlation improves the prediction of anticancer drug candidates"**  
> Zetian Song, Keitaro Izumi, and Michio Iwata

## Requirements
The following environment is required to run the analysis:

- Python 3.8+
- numpy >= 1.21.0
- pandas >= 1.3.0
- matplotlib >= 3.4.0

After cloning the repository, install the required packages using the following commands:

```bash
git clone https://github.com/yourusername/your-repo-name.git
cd your-repo-name
pip install -r requirements.txt
```

## Repository Structure
```text
.
├── data/               # Directory for raw and intermediate data
│   ├── row/
│   ├── processed/
├── results/            # Directory for output results (figures and tables)
├── src/            # Analysis scripts
│   ├── 01_preprocess.py
│   ├── 02_analyze.py
│   └── 03_plot_results.py
├── requirements.txt    # List of required Python packages
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
python scripts/01_preprocess.py audit-raw \
    --raw-root data/raw_geo_dir/
```

**To export a specific sample as a CMap-wide spot table:**
```bash
python scripts/01_preprocess.py export-wide \
    --raw-root data/raw_geo_dir/ \
    --sample-id  \
    --cmap-signatures data/cmap_signatures.gmt \
    --mapping-dir data/mapping_dir/ \
    --out-root results/preprocessed/ \
    --leiden-resolution 0.6 \
    --random-state 0
```
*(Optional arguments: `--keep-out-of-tissue` to retain out-of-tissue spots, or `--counts-tsv` to use a small fixture dataset instead of the 10x H5 file.)*

**Example Command:**
```bash
python scripts/01_preprocess.py export-wide \
    --raw-root data/raw/geo/GSM6433585 \
    --sample-id GSM6433585_092A \
    --cmap-signatures data/raw/BC_gene_expression_signatures_mas5.txt \
    --mapping-dir data/processed/gene_mappings/ \
    --out-root data/processed/ \
    --leiden-resolution 0.6 \
    --random-state 0
```

### 2. Main Analysis (`02_analyze.py`)
Performs analysis using the XX model on the preprocessed data.
```bash
python scripts/02_analyze.py --input data/processed_data.csv --output results/analysis_summary.json
```

### 3. Result Visualization (`03_plot_results.py`)
Generates Figure 1 and Figure 2 for the paper based on the analysis results.
```bash
python scripts/03_plot_results.py --input results/analysis_summary.json --outdir results/figures/
```

## Expected Output
Upon successful execution, the following files will be generated in the `results/figures/` directory:
- `figure_1_distribution.png`
- `figure_2_correlation.pdf`

## Citation
If you use this code or data in your research, please cite our paper:

```bibtex
@article{YourName202X,
  title={Your Paper Title},
  author={Lastname, Firstname and Lastname, Firstname},
  journal={Journal Name},
  volume={1},
  number={1},
  pages={1--10},
  year={202X},
  publisher={Publisher Name}
}
```

## License
This project is licensed under the [MIT License](LICENSE).

## Contact
If you have any questions or bug reports, please open an [Issue](https://github.com/yourusername/your-repo-name/issues) or contact [Your Name] at email@example.com.

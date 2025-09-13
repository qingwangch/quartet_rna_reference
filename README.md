# Quartet_ratio_ground_truth
> A ratio-based framework using Quartet reference materials for integrating long- and short-read RNA-seq

[![DOI](https://zenodo.org/badge/DOI/TODO.svg)](https://doi.org/TODO) <!-- 如暂未生成 DOI 可先注释 -->
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## Table of Contents
1. [Overview](#overview)
2. [Directory Structure](#directory-structure)
3. [Quick Start](#quick-start)
4. [Analysis Workflow](#analysis-workflow)
5. [Input Data](#input-data)
6. [Environment & Dependencies](#environment--dependencies)
7. [Results & Figures](#results--figures)
8. [Reproducing the Study](#reproducing-the-study)
9. [Troubleshooting](#troubleshooting)
10. [Contributing](#contributing)
11. [License](#license)
12. [Citation](#citation)
13. [Contact](#contact)

---

## Overview
Briefly introduce  
* **Research goal** – e.g. benchmarking isoform quantification across sequencing platforms using the Quartet reference.  
* **Key contribution** – ratio-based normalization, unified ground-truth set, reproducible pipeline.  
* Link to the corresponding manuscript (bioRxiv/Journal, DOI).

---

## Directory Structure
```text
quartet_rna_reference/
├── data_analysis/                 # Main data analysis
│   ├── 01_longvsshort/            # long- vs short-read comparison
│   ├── 02_reference_description/  # Quartet reference construction
│   └── 03_performance_eva/        # Accuracy, reproducibility, and other metrics
│
├── figures/                       # Publication-ready figures
│   ├── fig2_longvsshort/          # Fig. 2 — long- vs short-read comparison
│   ├── fig3_ratio/                # Fig. 3 — ratio-based normalization results
│   ├── fig4_refData/              # Fig. 4 — reference dataset characteristics
│   └── fig5_application/          # Fig. 5 — downstream application example
│
├── ref_construction/              # Scripts for building reference transcriptomes
│   ├── LO/                        # Long-only protocol–specific reference
│   ├── SO/                        # Short-only protocol–specific reference
│   └── src/                       # Helper code, configs, and logs
│
├── upstream/                      # Main upstream analysis pipeline (scripts, and workflows)
│   ├── 01_preprocessing/          # Raw FASTQ QC, adapter trimming, mapping
│   ├── 02_quantification/         # Isoform/AS quantification tool wrappers
│   ├── 03_statistics/             # Shared statistical and plotting functions
│   └── 03_others/                 # Miscellaneous utilities
│
└── README.md                      # Project overview (this file)
```

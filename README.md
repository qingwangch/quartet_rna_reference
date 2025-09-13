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
├── data_analysis/        # Scripts, notebooks, and workflow definitions
│   ├── 01_preprocessing/     # QC, trimming, mapping
│   ├── 02_quantification/    # Isoform/gene quant
│   └── 03_statistics/        # Ratio calculation, plots
├── figures/              # Final figures for the manuscript
├── upstream/             # External resources (e.g. annotation, helper scripts)
└── README.md
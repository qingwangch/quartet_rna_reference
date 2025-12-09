# Quartet\_ratio\_ground\_truth

> **A ratio-based framework using Quartet reference materials for integrating long- and short-read RNA-seq**



---

## Table of Contents

- [Quartet\_ratio\_ground\_truth](#quartet_ratio_ground_truth)
  - [Table of Contents](#table-of-contents)
  - [Overview](#overview)
  - [Directory Structure](#directory-structure)
  - [Quick Start](#quick-start)
  - [Analysis Workflow](#analysis-workflow)
  - [Input Data](#input-data)
  - [Environment \& Dependencies](#environment--dependencies)
    - [R environment (≥ 4.3)](#r-environment--43)
    - [Python (≥ 3.9)](#python--39)
    - [External tools](#external-tools)
  - [Results \& Figures](#results--figures)
  - [Reproducing the Study](#reproducing-the-study)
  - [Troubleshooting](#troubleshooting)
  - [Contributing](#contributing)
  - [License](#license)
  - [Citation](#citation)
  - [Contact](#contact)

---

## Overview
A manuscript link is available at: https://www.biorxiv.org/content/10.1101/2025.09.15.676287v1

This repository provides the complete analysis workflow for:

**A ratio-based, platform-independent framework for integrating long-read (lrRNA-seq) and short-read (srRNA-seq) RNA-seq at the isoform level, empowered by *****Quartet***** reference materials.**

Key contributions:

- **Unified isoform ground truth** derived from multi-center, multi-platform long-read sequencing.
- **Ratio-based quantification** leveraging the same Quartet samples in each batch.
- **A fully reproducible analytical pipeline**, from raw data → quantification → ratio-based quantification → reference construction → performance evaluation.


---

## Directory Structure

```text
quartet_rna_reference/
├── data_analysis/                 # Main downstream analyses
│   ├── 01_longvsshort/            # Long- vs short-read comparisons
│   ├── 02_reference_description/  # Quartet reference characterization
│   └── 03_performance_eva/        # Accuracy, reproducibility, and bias metrics
│
├── figures/                       # Manuscript figure scripts
│   ├── fig2_longvsshort/          # Fig. 2 — platform comparison
│   ├── fig3_ratioQuant/           # Fig. 3 — ratio normalization
│   ├── fig4_refData/              # Fig. 4 — reference dataset characteristics
│   └── fig5_application/          # Fig. 5 — downstream application
│
├── ref_construction/              # Building the isoform-level reference
│   ├── LO/                        # Long-read-only reference construction
│   ├── SO/                        # Short-read-only modules
│   └── src/                       # Shared functions, configs, logs
│
├── upstream/                      # Upstream processing pipeline
│   ├── 00_qc/                     # Initial quality control (FastQC, NanoPlot, etc.)
│   ├── 01_preprocessing/          # Adapter trimming, read filtering, and format conversion
│   ├── 02_mapping/                # Long- and short-read alignment to reference genome
│   ├── 03_quantification/         # Isoform and alternative splicing quantification
│   └── 04_others/                 # Auxiliary or miscellaneous helper scripts
│
└── README.md                      # This file
```

---

## Quick Start

Clone:

```bash
git clone https://github.com/qingwangch/quartet_rna_reference.git
cd quartet_rna_reference
```

---

## Analysis Workflow

1. **Raw data preprocessing**

   - Adapter trimming and QC
   - Long-read alignment: `minimap2`
   - Short-read alignment: `STAR`

2. **Isoform quantification**

   - Long-read mapping/quantification:  `StringTie2`, `Bambu`, `Oarfish`, `isoQuant`
   - Short-read quantification: `StringTie2`, `RSEM`, `Salmon`, and `kallisto`


3. **Ratio-based quantification**

   - Computes robust D5/F7/M8\:D6 ratios
   - Removes batch/platform effects
   - Generates platform-agnostic expression values

4. **Reference construction**

   - Selects high-confidence isoform datasets
   - Integrates the common isoform quantifications 
   - Produces reference datasets for benchmarking

5. **Performance evaluation**

   - Accuracy
   - Reproducibility
   - Bias quantification

---

## Input Data

Raw sequencing data are **not included** in this repository.

To run the full workflow, prepare:

- Long-read RNA-seq FASTQ/BAM (PacBio, ONT)
- Short-read RNA-seq FASTQ
- Quartet metadata
- GTF annotation (MANE / Ensembl / custom)
- Optional: spike-ins (SIRV-Set4)

---

## Environment & Dependencies

### R environment (≥ 4.3)

Required packages include:

- `dplyr`, `tidyr`, `ggplot2`, `data.table`
- `SummarizedExperiment`, `BiocParallel`
- `DRIMSeq`, `edgeR`

### Python (≥ 3.9)

- `pandas`, `numpy`, `scipy`, `matplotlib`

### External tools

- `minimap2`, `samtools`
- `SQANTI3`, `StringTie3`, `StringTie2`, `Bambu`, `Oarfish`, `isoQuant`, `MPACT`, `miniQuant`
- `STAR`, `Salmon`, `kallisto`

Environment files (`conda.yaml`, `renv.lock`) will be provided later.

---

## Results & Figures

All figure-generating scripts are in:

```
figures/
```

Running them will reproduce:

- Long vs short-read comparison plots
- Ratio normalization effects
- Reference data characteristics
- Downstream biological examples

Large figure files are **not stored** to keep the repo lightweight.

---

## Reproducing the Study

To reproduce all results:

1. Acquire raw Quartet RNA-seq (long-read + short-read).
2. Run `upstream/` for preprocessing and quantification.
3. Build the reference using `ref_construction/`.
4. Evaluate performance using `data_analysis/`.
5. Generate figures via `figures/`.

A full workflow diagram will be added.

---

## Troubleshooting

| Issue                      | Likely Cause                | Recommended Fix                       |
| -------------------------- | --------------------------- | ------------------------------------- |
| Missing isoform IDs        | Annotation mismatch         | Use the provided MANE/Ensembl version |
| Extreme ratio values       | Low expression              | Apply TPM/CPM filtering               |
| Inconsistent LR–SR results | SQANTI category differences | Restrict to FSM + ISM isoforms        |

---

## Contributing

Contributions are welcome.\
Please open an Issue or submit a Pull Request.

---

## License

This project is released under the **MIT License**.

---

## Citation

If you use this repository, please cite:

> Chen Q, et al. *A ratio-based framework using Quartet reference materials for integrating long- and short-read RNA-seq.* (Manuscript in revision)

DOI: https://www.biorxiv.org/content/10.1101/2025.09.15.676287v1
---

## Contact

Maintainer: **Qingwang Chen**\
Email: [qwchen20@fudan.edu.cn](mailto\:qwchen20@fudan.edu.cn)\
For issues and feature requests, please use GitHub Issues.


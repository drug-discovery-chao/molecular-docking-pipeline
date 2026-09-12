# Multi-Target Drug Screening Pipeline

> A unified, reusable computational pipeline for ChEMBL-based drug screening.
> Designed for computational biology / bioinformatics portfolio demonstration.

---

## Overview

This project implements a professional-grade drug screening pipeline that fetches bioactivity data from ChEMBL, computes molecular descriptors with RDKit, applies ADMET filters, and generates publication-ready visualizations.

**Target proteins screened:**
- **mTOR** (mechanistic Target of Rapamycin) — cancer / aging
- **SIRT1** (Sirtuin 1) — metabolism / longevity

**Key capability:** The pipeline is fully generic — add any ChEMBL target with a single line.

---

## Tech Stack

| Layer | Tool |
|-------|------|
| Molecular Cheminformatics | RDKit |
| Bioactivity Database | ChEMBL Web Resource API |
| Data Processing | pandas, numpy |
| Machine Learning (PCA) | scikit-learn |
| Visualization | matplotlib, seaborn |

---

## Core Features

### v2 (Current)

| Feature | Description |
|---------|-------------|
| **Generic Pipeline** | `DrugScreeningPipeline` class supports any ChEMBL target |
| **Smart Deduplication** | Groups by molecule ID, keeps **max** pChEMBL value |
| **PAINS Filter** | Removes Pan-Assay Interference Compounds (false positives) |
| **QED** | Quantitative Estimate of Drug-likeness (RDKit built-in) |
| **SA Score** | Synthetic Accessibility score (1=easy, 10=hard) |
| **Weighted Scoring** | 6-dimension composite (Lipinski/Veber/QED/SA/PAINS/MW sweet spot) |
| **Visualization Suite** | 5 chart types generated automatically |
| **Multi-Target Comparison** | Joint PCA, box plots, overlap detection across targets |
| **Production Logging** | Full exception handling + structured logging |

### v1 → v2 Evolution

See [`archive/`](archive/) for early standalone scripts. v2 consolidates duplicated code, adds scientific rigor (PAINS/QED/SA), and introduces software engineering practices (class architecture, error handling, logging).

---

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run full pipeline (mTOR + SIRT1)
python drug_screening_pipeline.py

# Single target
python drug_screening_pipeline.py --target MTOR

# Compare saved results (no re-fetch from ChEMBL)
python drug_screening_pipeline.py --compare
```

---

## Output Structure

```
results/
├── mtor_properties_full.csv
├── mtor_properties_filtered.csv
├── sirt1_properties_full.csv
├── sirt1_properties_filtered.csv
├── multi_target_comparison.csv
└── plots/
    ├── mtor_descriptor_distributions.png
    ├── mtor_activity_vs_druglikeness.png
    ├── mtor_pca_chemical_space.png
    ├── mtor_similarity_heatmap.png
    ├── mtor_screening_funnel.png
    ├── multi_target_pca.png
    └── multi_target_boxplots.png
```

---

## Scoring Model

The composite drug-likeness score (0–100) weights 6 dimensions:

| Dimension | Weight | Basis |
|-----------|--------|-------|
| Lipinski Rule of Five | 25 pts | Classical oral drug filter |
| Veber Rules | 15 pts | Oral bioavailability |
| QED | 20 pts | Quantitative drug-likeness (0–1 → 0–20) |
| SA Score | 15 pts | Synthetic accessibility (lower = easier) |
| PAINS Pass | 15 pts | No assay interference |
| MW Sweet Spot | 10 pts | 200–450 Da optimal range |

---

## Visualization Examples

| Chart | What It Shows |
|-------|---------------|
| Descriptor Distributions | MW/LogP/TPSA/RotatableBonds/QED/SA histograms |
| Activity vs Drug-likeness | Scatter of pChEMBL vs score (red = PAINS) |
| PCA Chemical Space | Morgan fingerprint PCA — target coverage |
| Tanimoto Heatmap | Structural similarity matrix |
| Screening Funnel | Compounds remaining after each filter stage |
| Multi-Target PCA | Overlap between mTOR and SIRT1 chemical spaces |

---

## File Structure

```
├── drug_screening_pipeline.py      # Main pipeline (v2)
├── requirements.txt                # Dependencies
├── molecular_properties.csv        # Reference data (5 common drugs)
├── archive/                        # Early versions for evolution demo
│   ├── README_ARCHIVE.md
│   ├── mtor_drug_screening_v1.py
│   └── sirt1_drug_screening_v1.py
└── README.md                         # This file
```

---

## Author

Built as a computational biology portfolio project demonstrating:
- Cheminformatics (RDKit, ChEMBL)
- Drug design & ADMET filtering
- Data science (pandas, visualization)
- Software engineering (class architecture, logging, CLI)

For academic / research use. Target audience: computational biology and bioinformatics MSc admissions.

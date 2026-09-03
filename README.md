# Molecular Docking Pipeline

An automated molecular docking pipeline built with Python, designed for virtual screening of small molecules against protein targets.

## 🎯 What This Project Does

- Batch prepares protein and ligand structures
- Runs molecular docking using DiffDock
- Analyzes binding affinity and poses
- Generates visualization reports

## 🛠️ Tech Stack

- Python 3.9+
- RDKit (cheminformatics)
- DiffDock (molecular docking)
- Pandas / NumPy (data processing)
- Matplotlib (visualization)

## 📂 Project Structure
├── data/ # Input molecules and protein structures
├── scripts/ # Python automation scripts
├── results/ # Docking outputs and reports
└── README.md
## 🧬 mTOR Screening Module

ChEMBL-based virtual screening pipeline:
- `mtor_drug_screening.py` — RDKit descriptor computation + Lipinski/Veber filtering
- `molecular_properties.csv` — Computed properties for mTOR inhibitors
- ## 🧬 SIRT1 Screening Module

ChEMBL-based virtual screening pipeline for SIRT1 (NAD+-dependent deacetylase sirtuin-1), a key longevity target linked to caloric restriction and metabolic health.

**Files:**
- `sirt1_drug_screening.py` — RDKit descriptor computation + Lipinski/Veber filtering
- `sirt1_properties_filtered.csv` — 48 lead compounds passing dual filters

**Filters Applied:**
- Lipinski's Rule of Five (MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10)
- Veber's Rules (Rotatable Bonds ≤ 10, TPSA ≤ 140)

**Scale:**
- Total compounds analyzed: 100
- Pass both filters + drug-likeness > 60: 48

**Target Significance:** SIRT1 activation mimics caloric restriction effects; inhibitors are investigated for metabolic disorders and potential anti-aging interventions.

Filters applied:
- Lipinski's Rule of Five (MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10)
- Veber's Rules (Rotatable Bonds ≤ 10, TPSA ≤ 140)
## 🚀 Quick Start

```bash
pip install -r requirements.txt
python run_docking.py --protein target.pdb --ligands ligands.sdf
📧 Contact
For collaboration or freelance inquiries: 508454132@qq.com


Note: This project is under active development. More features coming soon.

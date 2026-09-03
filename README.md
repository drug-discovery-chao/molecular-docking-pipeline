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

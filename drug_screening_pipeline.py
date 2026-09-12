#!/usr/bin/env python3
"""
Multi-Target Drug Screening Pipeline
====================================
A unified, reusable framework for ChEMBL-based drug screening.

Features:
  - Generic pipeline: any ChEMBL target by name or ID
  - Dedup by molecule ID, keeping the max pChEMBL value
  - PAINS filter (Pan-Assay Interference Compounds)
  - QED (Quantitative Estimate of Drug-likeness)
  - Synthetic Accessibility (SA) Score
  - Weighted composite drug-likeness score
  - Visualization: descriptor distributions, scatter, PCA, similarity heatmap
  - Multi-target comparison report

Usage:
  python drug_screening_pipeline.py                # run mTOR + SIRT1
  python drug_screening_pipeline.py --target MTOR  # single target
  python drug_screening_pipeline.py --compare       # compare saved results

Author: Computational Biology Portfolio Project
"""

import argparse
import logging
import os
import sys
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S'
)
log = logging.getLogger('DrugScreening')

# ---------------------------------------------------------------------------
# RDKit imports (deferred error handling)
# ---------------------------------------------------------------------------
try:
    from rdkit import Chem
    from rdkit.Chem import (
        AllChem,
        DataStructs,
        Descriptors,
        Lipinski,
        QED,
        FilterCatalog,
        rdMolDescriptors,
    )
    from rdkit.Chem.FilterCatalog import FilterCatalogParams
except ImportError:
    log.error("RDKit is required. Install with: pip install rdkit")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Matplotlib / Seaborn (deferred so script runs headless without them)
# ---------------------------------------------------------------------------
MPL_FONT_PATH = r'C:\Windows\Fonts\msyh.ttc'  # Microsoft YaHei
MPL_AVAILABLE = False
plt = None
sns = None
try:
    import matplotlib
    matplotlib.use('Agg')  # headless-safe
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    MPL_AVAILABLE = True
    if os.path.exists(MPL_FONT_PATH):
        fm.fontManager.addfont(MPL_FONT_PATH)
        _font_prop = fm.FontProperties(fname=MPL_FONT_PATH)
        plt.rcParams['font.family'] = _font_prop.get_name()
        plt.rcParams['axes.unicode_minus'] = False
    try:
        import seaborn as sns
        sns.set_style('whitegrid')
    except ImportError:
        sns = None
except ImportError:
    log.warning("matplotlib not available — skipping visualizations")


# ===================================================================
# 1.  ChEMBL Data Fetching
# ===================================================================

def fetch_chembl_data(target_name, target_chembl_id=None, limit=100):
    """
    Fetch activity data from ChEMBL for a given target.

    Args:
        target_name: e.g. 'MTOR', 'SIRT1' (used for logging / filenames)
        target_chembl_id: e.g. 'CHEMBL4506'. If None, look up by pref_name.
        limit: max unique molecules to keep.

    Returns:
        list of {chembl_id, smiles, pchembl_value}
    """
    try:
        from chembl_webresource_client.new_client import new_client
    except Exception as e:
        log.error(f"Cannot import ChEMBL client: {e}")
        log.error("Install with: pip install chembl_webresource_client")
        return []

    # Resolve target ID
    if target_chembl_id is None:
        log.info(f"Looking up ChEMBL target ID for '{target_name}' ...")
        try:
            target = new_client.target
            results = target.filter(
                pref_name__iexact=target_name
            ).only('target_chembl_id')
            targets = list(results)
            if not targets:
                log.error(f"No target found for name '{target_name}'")
                return []
            target_chembl_id = targets[0]['target_chembl_id']
        except Exception as e:
            log.error(f"Failed to look up target: {e}")
            return []

    log.info(f"Fetching activities for {target_name} (ID: {target_chembl_id}) ...")

    try:
        activity = new_client.activity
        activities = activity.filter(
            target_chembl_id=target_chembl_id,
            pchembl_value__isnull=False,
            molecule_chembl_id__isnull=False
        ).only('molecule_chembl_id', 'canonical_smiles', 'pchembl_value')
    except Exception as e:
        log.error(f"ChEMBL API request failed: {e}")
        return []

    # Group by molecule, keep max pChEMBL
    mol_data = defaultdict(lambda: {'smiles': None, 'pchembl_values': []})

    count = 0
    for act in activities:
        mol_id = act.get('molecule_chembl_id')
        smiles = act.get('canonical_smiles')
        pval = act.get('pchembl_value')

        if not mol_id or not smiles or len(smiles) <= 3 or pval is None:
            continue

        try:
            pval = float(pval)
        except (ValueError, TypeError):
            continue

        mol_data[mol_id]['smiles'] = smiles
        mol_data[mol_id]['pchembl_values'].append(pval)
        count += 1

    # Deduplicate: keep max pChEMBL per molecule
    data = []
    for mol_id, info in mol_data.items():
        if not info['smiles']:
            continue
        max_pval = max(info['pchembl_values'])
        data.append({
            'chembl_id': mol_id,
            'smiles': info['smiles'],
            'pchembl_value': max_pval
        })
        if len(data) >= limit:
            break

    log.info(f"Retrieved {count} raw activity records → {len(data)} unique compounds")
    return data


# ===================================================================
# 2.  Molecular Descriptors
# ===================================================================

def compute_properties(smiles):
    """Compute molecular descriptors + QED + SA score."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    props = {
        'MW': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'HBD': Lipinski.NumHDonors(mol),
        'HBA': Lipinski.NumHAcceptors(mol),
        'RotatableBonds': Lipinski.NumRotatableBonds(mol),
        'TPSA': rdMolDescriptors.CalcTPSA(mol),
        'NumRings': Lipinski.RingCount(mol),
        'NumAromRings': Lipinski.NumAromaticRings(mol),
        'QED': QED.qed(mol),
        'SA_Score': _compute_sa_score(mol),
    }
    return props


def _compute_sa_score(mol):
    """
    Synthetic Accessibility Score (1 = easy, 10 = hard).
    Simplified version based on ring count and complexity.
    """
    try:
        from rdkit.Chem import Crippen
        from rdkit.Chem.SA_Score import sascorer  # type: ignore
        return sascorer.calculateScore(mol)
    except Exception:
        # Fallback heuristic if SA module unavailable
        rings = Lipinski.RingCount(mol)
        arom = Lipinski.NumAromaticRings(mol)
        mw = Descriptors.MolWt(mol)
        score = 1.0 + rings * 0.5 + arom * 0.3 + (mw - 200) * 0.01
        return round(min(10.0, max(1.0, score)), 2)


# ===================================================================
# 3.  Filters
# ===================================================================

# PAINS catalog (lazy init)
_PAINS_CATALOG = None

def _get_pains_catalog():
    global _PAINS_CATALOG
    if _PAINS_CATALOG is None:
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(params.FilterCatalogs.PAINS)
        _PAINS_CATALOG = FilterCatalog.FilterCatalog(params)
    return _PAINS_CATALOG


def lipinski_filter(props):
    """Lipinski's Rule of Five."""
    return (
        props['MW'] <= 500 and
        props['LogP'] <= 5 and
        props['HBD'] <= 5 and
        props['HBA'] <= 10
    )


def veber_filter(props):
    """Veber's rules for oral bioavailability."""
    return (
        props['RotatableBonds'] <= 10 and
        props['TPSA'] <= 140
    )


def pains_filter(smiles):
    """Return True if molecule is NOT a PAINS compound (passes filter)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    catalog = _get_pains_catalog()
    return not catalog.HasMatch(mol)


# ===================================================================
# 4.  Composite Drug-likeness Score (weighted)
# ===================================================================

def druglikeness_score(props, smiles):
    """
    Weighted composite score (0-100).

    Components:
      - Lipinski compliance (25 pts)
      - Veber compliance  (15 pts)
      - QED score         (20 pts, scaled 0-1 → 0-20)
      - SA Score          (15 pts, lower = better)
      - PAINS pass        (15 pts)
      - MW in sweet spot  (10 pts, 200-450 Da)
    """
    score = 0.0

    # Lipinski (25)
    if lipinski_filter(props):
        score += 25
    else:
        score += max(0, 25 - 10 * sum([
            props['MW'] > 500,
            props['LogP'] > 5,
            props['HBD'] > 5,
            props['HBA'] > 10
        ]))

    # Veber (15)
    if veber_filter(props):
        score += 15
    else:
        score += max(0, 15 - 8 * sum([
            props['RotatableBonds'] > 10,
            props['TPSA'] > 140
        ]))

    # QED (20)
    score += props['QED'] * 20

    # SA Score (15) — 1=easy → full marks; 10=hard → 0
    sa = props['SA_Score']
    sa_pts = max(0, 15 * (10 - sa) / 9.0)
    score += sa_pts

    # PAINS (15)
    if pains_filter(smiles):
        score += 15

    # MW sweet spot (10) — 200-450 Da
    if 200 <= props['MW'] <= 450:
        score += 10
    elif 150 <= props['MW'] <= 500:
        score += 5

    return round(min(100, max(0, score)), 2)


# ===================================================================
# 5.  Pipeline (class-based)
# ===================================================================

class DrugScreeningPipeline:
    """Generic multi-target drug screening pipeline."""

    def __init__(self, target_name, target_chembl_id=None, limit=100,
                 output_dir='results'):
        self.target_name = target_name
        self.target_chembl_id = target_chembl_id
        self.limit = limit
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        self.df_full = None
        self.df_filtered = None

    def _out(self, suffix):
        return os.path.join(self.output_dir, f'{self.target_name.lower()}_{suffix}')

    def run(self):
        """Execute the full pipeline."""
        log.info(f"{'='*60}")
        log.info(f"Starting pipeline: {self.target_name}")
        log.info(f"{'='*60}")

        # 1. Fetch data
        raw_data = fetch_chembl_data(
            self.target_name,
            self.target_chembl_id,
            self.limit
        )
        if not raw_data:
            log.error(f"No data retrieved for {self.target_name}, skipping.")
            return None

        # 2. Compute properties + filters
        results = []
        n_pains_removed = 0
        for i, item in enumerate(raw_data):
            if (i + 1) % 10 == 0:
                log.info(f"  Processing {i+1}/{len(raw_data)} ...")

            props = compute_properties(item['smiles'])
            if props is None:
                log.warning(f"  Failed to parse SMILES: {item['chembl_id']}")
                continue

            props['chembl_id'] = item['chembl_id']
            props['smiles'] = item['smiles']
            props['pchembl_value'] = item['pchembl_value']
            props['pass_lipinski'] = lipinski_filter(props)
            props['pass_veber'] = veber_filter(props)
            props['pass_pains'] = pains_filter(item['smiles'])

            if not props['pass_pains']:
                n_pains_removed += 1

            props['druglikeness_score'] = druglikeness_score(
                props, item['smiles']
            )
            results.append(props)

        if not results:
            log.error(f"No valid compounds for {self.target_name}")
            return None

        # 3. Build DataFrame
        df = pd.DataFrame(results)
        col_order = [
            'chembl_id', 'smiles', 'pchembl_value',
            'MW', 'LogP', 'HBD', 'HBA', 'RotatableBonds', 'TPSA',
            'NumRings', 'NumAromRings', 'QED', 'SA_Score',
            'pass_lipinski', 'pass_veber', 'pass_pains',
            'druglikeness_score'
        ]
        df = df[col_order]

        # 4. Save full results
        full_path = self._out('properties_full.csv')
        df.to_csv(full_path, index=False, encoding='utf-8')
        log.info(f"Saved full results → {full_path} ({len(df)} compounds)")

        # 5. Filtered results
        filtered = df[
            df['pass_lipinski'] &
            df['pass_veber'] &
            df['pass_pains'] &
            (df['druglikeness_score'] > 60)
        ].sort_values('druglikeness_score', ascending=False)

        filt_path = self._out('properties_filtered.csv')
        filtered.to_csv(filt_path, index=False, encoding='utf-8')
        log.info(f"Saved filtered → {filt_path} ({len(filtered)} compounds)")

        # 6. Summary
        log.info(f"\n{'='*60}")
        log.info(f"Pipeline Summary: {self.target_name}")
        log.info(f"{'='*60}")
        log.info(f"  Total compounds:     {len(df)}")
        log.info(f"  Pass Lipinski:      {df['pass_lipinski'].sum()}")
        log.info(f"  Pass Veber:         {df['pass_veber'].sum()}")
        log.info(f"  Pass PAINS:          {df['pass_pains'].sum()}")
        log.info(f"  PAINS removed:       {n_pains_removed}")
        log.info(f"  Pass all + score>60: {len(filtered)}")
        log.info(f"\nTop 5 by drug-likeness score:")
        top5 = df.nlargest(5, 'druglikeness_score')[
            ['chembl_id', 'MW', 'LogP', 'QED', 'SA_Score', 'druglikeness_score']
        ]
        for _, row in top5.iterrows():
            log.info(
                f"  {row['chembl_id']}  MW={row['MW']:.1f}  "
                f"LogP={row['LogP']:.2f}  QED={row['QED']:.3f}  "
                f"SA={row['SA_Score']:.2f}  Score={row['druglikeness_score']}"
            )

        # 7. Visualization
        if MPL_AVAILABLE:
            self._generate_plots(df, filtered)

        self.df_full = df
        self.df_filtered = filtered
        return df


    # ---------------------------------------------------------------
    #  Visualization
    # ---------------------------------------------------------------

    def _generate_plots(self, df_full, df_filtered):
        """Generate 4 figures."""
        plots_dir = os.path.join(self.output_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        tname = self.target_name

        # --- Plot 1: Descriptor distributions ---
        fig, axes = plt.subplots(2, 3, figsize=(15, 9))
        fig.suptitle(f'{tname} - Molecular Descriptor Distributions',
                     fontsize=14, fontweight='bold')
        descriptors = ['MW', 'LogP', 'TPSA', 'RotatableBonds', 'QED', 'SA_Score']
        for ax, desc in zip(axes.flat, descriptors):
            if sns:
                sns.histplot(df_full[desc], bins=20, ax=ax, color='steelblue',
                             edgecolor='white', alpha=0.8)
            else:
                ax.hist(df_full[desc], bins=20, color='steelblue', alpha=0.8)
            ax.set_xlabel(desc)
            ax.set_ylabel('Count')
        plt.tight_layout()
        p = os.path.join(plots_dir, f'{tname.lower()}_descriptor_distributions.png')
        plt.savefig(p, dpi=150)
        plt.close()
        log.info(f"  Plot saved → {p}")

        # --- Plot 2: pChEMBL vs drug-likeness scatter ---
        fig, ax = plt.subplots(figsize=(8, 6))
        colors = np.where(df_full['pass_pains'], 'steelblue', 'red')
        ax.scatter(df_full['pchembl_value'], df_full['druglikeness_score'],
                   c=colors, alpha=0.7, edgecolors='navy', linewidths=0.5)
        ax.set_xlabel('pChEMBL Value (max)')
        ax.set_ylabel('Drug-likeness Score')
        ax.set_title(f'{tname} - Activity vs Drug-likeness\n(red = PAINS compound)')
        plt.tight_layout()
        p = os.path.join(plots_dir, f'{tname.lower()}_activity_vs_druglikeness.png')
        plt.savefig(p, dpi=150)
        plt.close()
        log.info(f"  Plot saved → {p}")

        # --- Plot 3: PCA of chemical space ---
        mols = []
        valid_idx = []
        for idx, row in df_full.iterrows():
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol is not None and mol.GetNumAtoms() > 0:
                mols.append(mol)
                valid_idx.append(idx)

        if len(mols) >= 2:
            from sklearn.decomposition import PCA  # optional
            from rdkit.Chem import rdFingerprintGenerator
            gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
            fps = np.zeros((len(mols), 1024), dtype=np.int8)
            for i, mol in enumerate(mols):
                fp = gen.GetFingerprint(mol)
                arr = np.zeros((1024,), dtype=np.int8)
                DataStructs.ConvertToNumpyArray(fp, arr)
                fps[i] = arr
            pca = PCA(n_components=2)
            pcs = pca.fit_transform(fps)
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(pcs[:, 0], pcs[:, 1], c='steelblue', alpha=0.7,
                       edgecolors='navy', linewidths=0.5)
            ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
            ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
            ax.set_title(f'{tname} - Chemical Space (PCA of Morgan FP)')
            plt.tight_layout()
            p = os.path.join(plots_dir, f'{tname.lower()}_pca_chemical_space.png')
            plt.savefig(p, dpi=150)
            plt.close()
            log.info(f"  Plot saved → {p}")
        else:
            log.warning("  Not enough valid molecules for PCA plot")

        # --- Plot 4: Tanimoto similarity heatmap ---
        # Generate fingerprints from molecules
        fps = []
        for mol in mols:
            fp = gen.GetFingerprint(mol)
            fps.append(fp)
        
        if 3 <= len(fps) <= 50:
            n = len(fps)
            sim_matrix = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    sim_matrix[i, j] = DataStructs.TanimotoSimilarity(
                        fps[i], fps[j]
                    )
            fig, ax = plt.subplots(figsize=(8, 7))
            if sns:
                sns.heatmap(sim_matrix, ax=ax, cmap='YlOrRd',
                            xticklabels=False, yticklabels=False,
                            cbar_kws={'label': 'Tanimoto Similarity'})
            else:
                im = ax.imshow(sim_matrix, cmap='YlOrRd')
                plt.colorbar(im, ax=ax, label='Tanimoto Similarity')
            ax.set_title(f'{tname} - Tanimoto Similarity Heatmap')
            plt.tight_layout()
            p = os.path.join(plots_dir, f'{tname.lower()}_similarity_heatmap.png')
            plt.savefig(p, dpi=150)
            plt.close()
            log.info(f"  Plot saved → {p}")

        # --- Plot 5 (bonus): Filter funnel ---
        stages = {
            'Total': len(df_full),
            'Lipinski': df_full['pass_lipinski'].sum(),
            'Veber': df_full['pass_veber'].sum(),
            'PAINS': df_full['pass_pains'].sum(),
            'Score>60': len(df_full[df_full['druglikeness_score'] > 60]),
            'Final': len(df_full[
                df_full['pass_lipinski'] &
                df_full['pass_veber'] &
                df_full['pass_pains'] &
                (df_full['druglikeness_score'] > 60)
            ])
        }
        fig, ax = plt.subplots(figsize=(8, 5))
        bars = ax.bar(stages.keys(), stages.values(),
                       color=['steelblue', 'cornflowerblue', 'skyblue',
                              'lightblue', 'lightskyblue', 'seagreen'],
                       edgecolor='navy')
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., h + 0.5,
                    f'{int(h)}', ha='center', fontsize=11)
        ax.set_ylabel('Number of Compounds')
        ax.set_title(f'{tname} - Screening Funnel')
        plt.tight_layout()
        p = os.path.join(plots_dir, f'{tname.lower()}_screening_funnel.png')
        plt.savefig(p, dpi=150)
        plt.close()
        log.info(f"  Plot saved → {p}")


# ===================================================================
# 6.  Multi-target Comparison
# ===================================================================

def compare_targets(targets_data, output_dir='results'):
    """
    Compare screening results across targets.

    Args:
        targets_data: dict of {target_name: DataFrame}
        output_dir: output directory
    """
    log.info(f"\n{'#'*60}")
    log.info(f"# Multi-Target Comparison Report")
    log.info(f"{'#'*60}")

    os.makedirs(os.path.join(output_dir, 'plots'), exist_ok=True)

    # Summary table
    rows = []
    for tname, df in targets_data.items():
        rows.append({
            'Target': tname,
            'Total_Compounds': len(df),
            'Pass_Lipinski': df['pass_lipinski'].sum() if 'pass_lipinski' in df else 0,
            'Pass_Veber': df['pass_veber'].sum() if 'pass_veber' in df else 0,
            'Pass_PAINS': df['pass_pains'].sum() if 'pass_pains' in df else 0,
            'Mean_MW': df['MW'].mean(),
            'Mean_LogP': df['LogP'].mean(),
            'Mean_TPSA': df['TPSA'].mean(),
            'Mean_QED': df['QED'].mean(),
            'Mean_SA': df['SA_Score'].mean(),
            'Mean_pChEMBL': df['pchembl_value'].mean(),
            'Mean_Score': df['druglikeness_score'].mean(),
        })
    summary = pd.DataFrame(rows)
    summary_path = os.path.join(output_dir, 'multi_target_comparison.csv')
    summary.to_csv(summary_path, index=False, encoding='utf-8')
    log.info(f"Comparison table saved → {summary_path}")
    log.info("\n" + summary.to_string(index=False))

    # Find overlapping compounds (by chembl_id)
    if len(targets_data) >= 2:
        id_sets = {t: set(df['chembl_id']) for t, df in targets_data.items()}
        all_ids = set()
        for s in id_sets.values():
            all_ids |= s
        overlap = all_ids.copy()
        for s in id_sets.values():
            overlap &= s
        if overlap:
            log.info(f"\nCompounds hitting multiple targets: {len(overlap)}")
            for cid in sorted(overlap):
                log.info(f"  {cid}")
        else:
            log.info(f"\nNo overlapping compounds across targets.")

    # Combined PCA plot
    if MPL_AVAILABLE and len(targets_data) >= 2:
        _plot_combined_pca(targets_data, output_dir)

    # Descriptor box plots
    if MPL_AVAILABLE and len(targets_data) >= 2:
        _plot_descriptor_boxplots(targets_data, output_dir)


def _plot_combined_pca(targets_data, output_dir):
    """PCA of all targets on one plot."""
    all_fps = []
    all_labels = []
    all_colors = []
    color_map = {}
    palette = ['steelblue', 'coral', 'seagreen', 'purple', 'darkorange']

    for i, (tname, df) in enumerate(targets_data.items()):
        color_map[tname] = palette[i % len(palette)]
        for _, row in df.iterrows():
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol and mol.GetNumAtoms() > 0:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
                all_fps.append(fp)
                all_labels.append(tname)
                all_colors.append(color_map[tname])

    if len(all_fps) < 3:
        log.warning("Not enough molecules for combined PCA")
        return

    try:
        from sklearn.decomposition import PCA
    except ImportError:
        log.warning("scikit-learn not available, skipping PCA")
        return

    fps = np.zeros((len(all_fps), 1024), dtype=np.int8)
    for i, fp in enumerate(all_fps):
        arr = np.zeros((1024,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fps[i] = arr

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(fps)

    fig, ax = plt.subplots(figsize=(9, 7))
    for tname in targets_data:
        mask = [l == tname for l in all_labels]
        xs = [pcs[j, 0] for j in range(len(mask)) if mask[j]]
        ys = [pcs[j, 1] for j in range(len(mask)) if mask[j]]
        ax.scatter(xs, ys, c=color_map[tname], alpha=0.7,
                   edgecolors='white', linewidths=0.5, label=tname)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
    ax.set_title('Multi-Target Chemical Space (PCA)')
    ax.legend()
    plt.tight_layout()
    p = os.path.join(output_dir, 'plots', 'multi_target_pca.png')
    plt.savefig(p, dpi=150)
    plt.close()
    log.info(f"  Combined PCA saved → {p}")


def _plot_descriptor_boxplots(targets_data, output_dir):
    """Box plots comparing descriptors across targets."""
    combined = pd.concat(
        [df.assign(Target=t) for t, df in targets_data.items()],
        ignore_index=True
    )
    descriptors = ['MW', 'LogP', 'TPSA', 'QED', 'SA_Score', 'druglikeness_score']
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle('Multi-Target Descriptor Comparison', fontsize=14,
                 fontweight='bold')
    for ax, desc in zip(axes.flat, descriptors):
        if sns:
            sns.boxplot(data=combined, x='Target', y=desc, ax=ax,
                        palette='Set2')
        else:
            targets = combined['Target'].unique()
            data = [combined[combined['Target'] == t][desc] for t in targets]
            ax.boxplot(data, labels=targets)
            ax.set_ylabel(desc)
    plt.tight_layout()
    p = os.path.join(output_dir, 'plots', 'multi_target_boxplots.png')
    plt.savefig(p, dpi=150)
    plt.close()
    log.info(f"  Box plots saved → {p}")


# ===================================================================
# 7.  Main
# ===================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Multi-Target Drug Screening Pipeline'
    )
    parser.add_argument('--target', type=str, default=None,
                        help='Single target to screen (e.g. MTOR, SIRT1)')
    parser.add_argument('--compare', action='store_true',
                        help='Compare previously saved results')
    parser.add_argument('--limit', type=int, default=100,
                        help='Max compounds per target')
    parser.add_argument('--output', type=str, default='results',
                        help='Output directory')
    args = parser.parse_args()

    # Target registry
    targets = [
        {'name': 'MTOR', 'chembl_id': None},      # dynamic lookup
        {'name': 'SIRT1', 'chembl_id': 'CHEMBL4506'},
    ]

    if args.target:
        targets = [{'name': args.target.upper(), 'chembl_id': None}]

    # If --compare, load saved CSVs instead of fetching
    if args.compare:
        targets_data = {}
        for t in targets:
            path = os.path.join(args.output, f"{t['name'].lower()}_properties_full.csv")
            if os.path.exists(path):
                targets_data[t['name']] = pd.read_csv(path)
                log.info(f"Loaded {t['name']}: {len(targets_data[t['name']])} compounds")
            else:
                log.warning(f"No saved data for {t['name']} at {path}")
        if len(targets_data) >= 1:
            compare_targets(targets_data, args.output)
        else:
            log.error("No saved data found. Run without --compare first.")
        return

    # Run pipeline for each target
    all_results = {}
    for t in targets:
        pipeline = DrugScreeningPipeline(
            target_name=t['name'],
            target_chembl_id=t['chembl_id'],
            limit=args.limit,
            output_dir=args.output
        )
        df = pipeline.run()
        if df is not None:
            all_results[t['name']] = df

    # Multi-target comparison
    if len(all_results) >= 2:
        compare_targets(all_results, args.output)

    log.info(f"\n{'='*60}")
    log.info(f"All pipelines complete. Results in: {args.output}/")
    log.info(f"{'='*60}")


if __name__ == '__main__':
    main()

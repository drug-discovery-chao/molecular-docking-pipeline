#!/usr/bin/env python3
"""
MTOR Inhibitor Drug Screening Pipeline
Uses RDKit + ChEMBL API to compute molecular properties and apply ADMET filters.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
import pandas as pd
from chembl_webresource_client.new_client import new_client
import warnings
warnings.filterwarnings('ignore')

def get_chembl_mtor_data(limit=50):
    """Fetch mTOR inhibitor data from ChEMBL."""
    target = new_client.target
    mtor = target.filter(pref_name__iexact='MTOR').only('target_chembl_id')
    mtor_id = list(mtor)[0]['target_chembl_id']

    activity = new_client.activity
    activities = activity.filter(
        target_chembl_id=mtor_id,
        pchembl_value__isnull=False,
        molecule_chembl_id__isnull=False
    ).only('molecule_chembl_id', 'canonical_smiles', 'pchembl_value')

    data = []
    seen = set()
    for act in activities:
        mol_id = act['molecule_chembl_id']
        if mol_id in seen:
            continue
        seen.add(mol_id)
        smiles = act['canonical_smiles']
        if smiles and len(smiles) > 3:
            data.append({
                'chembl_id': mol_id,
                'smiles': smiles,
                'pchembl_value': float(act['pchembl_value'])
            })
        if len(data) >= limit:
            break
    return data

def compute_properties(smiles):
    """Compute molecular descriptors using RDKit."""
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
        'NumAromRings': Lipinski.NumAromaticRings(mol)
    }
    return props

def lipinski_filter(props):
    """Lipinski's Rule of Five filter."""
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

def druglikeness_score(props):
    """Simple drug-likeness composite score (0-100)."""
    score = 100
    if props['MW'] > 500: score -= 15
    if props['LogP'] > 5: score -= 15
    if props['HBD'] > 5: score -= 10
    if props['HBA'] > 10: score -= 10
    if props['RotatableBonds'] > 10: score -= 10
    if props['TPSA'] > 140: score -= 10
    return max(0, score)

def main():
    print("Fetching mTOR inhibitor data from ChEMBL...")
    raw_data = get_chembl_mtor_data(limit=100)
    print(f"Retrieved {len(raw_data)} compounds.")

    results = []
    for item in raw_data:
        props = compute_properties(item['smiles'])
        if props is None:
            continue

        props['chembl_id'] = item['chembl_id']
        props['smiles'] = item['smiles']
        props['pchembl_value'] = item['pchembl_value']
        props['pass_lipinski'] = lipinski_filter(props)
        props['pass_veber'] = veber_filter(props)
        props['druglikeness_score'] = druglikeness_score(props)
        results.append(props)

    df = pd.DataFrame(results)
    cols = [
        'chembl_id', 'smiles', 'pchembl_value',
        'MW', 'LogP', 'HBD', 'HBA', 'RotatableBonds', 'TPSA',
        'NumRings', 'NumAromRings',
        'pass_lipinski', 'pass_veber', 'druglikeness_score'
    ]
    df = df[cols]

    # Save full results
    df.to_csv('molecular_properties_full.csv', index=False)
    print(f"Saved full results to molecular_properties_full.csv ({len(df)} compounds)")

    # Save filtered: Lipinski + Veber + drug-likeness > 60
    filtered = df[(df['pass_lipinski']) & (df['pass_veber']) & (df['druglikeness_score'] > 60)]
    filtered.to_csv('molecular_properties_filtered.csv', index=False)
    print(f"Saved filtered results to molecular_properties_filtered.csv ({len(filtered)} compounds)")

    # Summary
    print("\n=== Pipeline Summary ===")
    print(f"Total compounds analyzed: {len(df)}")
    print(f"Pass Lipinski: {df['pass_lipinski'].sum()}")
    print(f"Pass Veber: {df['pass_veber'].sum()}")
    print(f"Pass both + drug-likeness > 60: {len(filtered)}")
    print("\nTop 5 by drug-likeness score:")
    print(df.nlargest(5, 'druglikeness_score')[['chembl_id', 'MW', 'LogP', 'druglikeness_score']].to_string(index=False))

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Example usage of the atct library."""

import os
from atct import (
    healthcheck,
    get_species,
    search_species,
    get_species_by_smiles,
    get_species_by_casrn,
    get_species_by_formula,
    get_species_by_name,
    get_species_by_inchi,
    get_species_by_inchikey,
    get_species_covariance_by_atctid,
    get_species_covariance_matrix,
    as_dataframe
)

def main():
    """Demonstrate atct functionality."""
    print("=== atct Example Usage ===\n")
    
    # Set API base URL (optional - defaults to production)
    # os.environ["ATCT_API_BASE_URL"] = "https://atct.anl.gov/api/v1"
    
    # 1. Health check
    print("1. Health Check:")
    if healthcheck():
        print("   ✅ API is healthy")
    else:
        print("   ❌ API is not responding")
        return
    print()
    
    # 2. Get specific species by ATcT ID
    print("2. Get Species by ATcT ID:")
    species = get_species("67-56-1*0")  # Methanol (gas)
    print(f"   Name: {species.name}")
    print(f"   Formula: {species.formula}")
    print(f"   SMILES: {species.smiles}")
    print(f"   CASRN: {species.casrn}")
    print(f"   ΔHf(298K): {species.delta_h_298k} ± {species.delta_h_298k_uncertainty} kJ/mol")
    print()
    
    # 3. Search for species
    print("3. Search for Species:")
    results = search_species("methanol", limit=3)
    print(f"   Found {len(results.items)} species:")
    for i, sp in enumerate(results.items, 1):
        print(f"   {i}. {sp.atct_id}: {sp.name} ({sp.formula})")
    print()
    
    # 4. Get species by SMILES
    print("4. Get Species by SMILES:")
    smiles_results = get_species_by_smiles("CO")
    print(f"   Found {len(smiles_results.items)} species for SMILES 'CO':")
    for sp in smiles_results.items:
        print(f"   - {sp.atct_id}: {sp.name} ({sp.formula})")
    print()
    
    # 5. Get species by CAS Registry Number
    print("5. Get Species by CASRN:")
    casrn_results = get_species_by_casrn("67-56-1")
    print(f"   Found {len(casrn_results.items)} species for CASRN '67-56-1':")
    for sp in casrn_results.items:
        print(f"   - {sp.atct_id}: {sp.name} ({sp.formula})")
    print()
    
    # 5b. Get species by formula
    print("5b. Get Species by Formula:")
    formula_results = get_species_by_formula("H2O")
    print(f"   Found {len(formula_results.items)} species for formula 'H2O':")
    for sp in formula_results.items[:3]:  # Show first 3
        print(f"   - {sp.atct_id}: {sp.name} ({sp.formula})")
    print()
    
    # 5c. Get species by name
    print("5c. Get Species by Name:")
    name_results = get_species_by_name("water")
    print(f"   Found {len(name_results.items)} species for name 'water':")
    for sp in name_results.items[:3]:  # Show first 3
        print(f"   - {sp.atct_id}: {sp.name} ({sp.formula})")
    print()
    
    # 5d. Get species by InChI
    print("5d. Get Species by InChI:")
    inchi_results = get_species_by_inchi("InChI=1S/CH4O/c1-2/h2H,1H3")
    print(f"   Found {len(inchi_results.items)} species for InChI:")
    for sp in inchi_results.items:
        print(f"   - {sp.atct_id}: {sp.name} ({sp.formula})")
    print()
    
    # 5e. Get species by InChI Key
    print("5e. Get Species by InChI Key:")
    inchikey_results = get_species_by_inchikey("OKKJLVBELUTLKV-UHFFFAOYSA-N")
    print(f"   Found {len(inchikey_results.items)} species for InChI Key:")
    for sp in inchikey_results.items:
        print(f"   - {sp.atct_id}: {sp.name} ({sp.formula})")
    print()
    
    # 6. Get covariance between two species (legacy 2x2 method)
    print("6. Get Covariance Matrix (2x2):")
    try:
        cov = get_species_covariance_by_atctid("67-56-1*500", "67-56-1*0")
        print(f"   Covariance between methanol phases:")
        print(f"   Labels: {cov.labels}")
        print(f"   Units: {cov.units}")
        print(f"   Matrix: {cov.matrix}")
    except Exception as e:
        print(f"   Error getting covariance: {e}")
    print()
    
    # 6b. Get full covariance matrix for multiple species
    print("6b. Get Full Covariance Matrix:")
    try:
        atctids = ["67-56-1*0", "67-56-1*500", "7727-37-9*0"]  # Methanol gas, liquid, O2
        cov_matrix = get_species_covariance_matrix(atctids=atctids)
        print(f"   Full covariance matrix for {len(atctids)} species:")
        print(f"   Species ATcT IDs: {cov_matrix.species_atctids}")
        print(f"   Units: {cov_matrix.units}")
        print(f"   Matrix shape: {len(cov_matrix.matrix)}x{len(cov_matrix.matrix[0])}")
        print(f"   Matrix diagonal (variances): {[cov_matrix.matrix[i][i] for i in range(len(cov_matrix.matrix))]}")
    except Exception as e:
        print(f"   Error getting full covariance matrix: {e}")
    print()
    
    # 7. Convert to pandas DataFrame (if pandas is available)
    print("7. Convert to Pandas DataFrame:")
    try:
        df = as_dataframe(results.items)
        print(f"   DataFrame shape: {df.shape}")
        print(f"   Columns: {list(df.columns)}")
        print("   First few rows:")
        print(df[['ATcT_ID', 'Name', 'Formula', 'SMILES']].to_string(index=False))
    except ImportError:
        print("   Pandas not available. Install with: pip install atct[pandas]")
    except Exception as e:
        print(f"   Error creating DataFrame: {e}")

if __name__ == "__main__":
    main()
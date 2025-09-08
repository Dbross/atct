#!/usr/bin/env python3
"""Example usage of the pyATcT library."""

import os
from atct import (
    healthcheck,
    get_species,
    search_species,
    get_species_by_smiles,
    get_species_by_casrn,
    get_species_covariance_by_atctid,
    as_dataframe
)

def main():
    """Demonstrate pyATcT functionality."""
    print("=== pyATcT Example Usage ===\n")
    
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
    print(f"   ΔHf(298K): {species.delta_h_298k} {species.unit}")
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
    
    # 6. Get covariance between two species
    print("6. Get Covariance Matrix:")
    try:
        cov = get_species_covariance_by_atctid("67-56-1*500", "67-56-1*0")
        print(f"   Covariance between methanol phases:")
        print(f"   Labels: {cov.labels}")
        print(f"   Units: {cov.units}")
        print(f"   Matrix: {cov.matrix}")
    except Exception as e:
        print(f"   Error getting covariance: {e}")
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
        print("   Pandas not available. Install with: pip install pyATcT[pandas]")
    except Exception as e:
        print(f"   Error creating DataFrame: {e}")

if __name__ == "__main__":
    main()
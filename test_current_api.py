#!/usr/bin/env python3
"""Test the pyATcT package against the current API."""

import atct
from atct import get_species, search_species, get_covariance, ATCTError, NotFound

def test_current_api():
    """Test against the current API endpoints."""
    print("Testing atct package against current API")
    print("=" * 40)
    
    # Temporarily modify the base URL to use current API
    import os
    os.environ['ATCT_API_BASE_URL'] = 'https://atct.anl.gov/api'
    
    # Test 1: Get species by ATcT ID
    print("\n1. Get Species by ATcT ID:")
    try:
        species = get_species("67-56-1*0")  # Methanol
        print(f"   ✓ Name: {species.name}")
        print(f"   ✓ Formula: {species.formula}")
        print(f"   ✓ ATcT ID: {species.atct_id}")
        if species.delta_h_298k is not None:
            print(f"   ✓ ΔHf(298K): {species.delta_h_298k:.2f} kJ/mol")
            if species.delta_h_298k_unc:
                print(f"   ✓ Uncertainty: ±{species.delta_h_298k_unc.value:.2f} {species.delta_h_298k_unc.unit}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 2: Search species
    print("\n2. Search Species:")
    try:
        results = search_species("methanol", limit=3)
        print(f"   ✓ Found {results.total} results (showing {len(results.items)})")
        for i, species in enumerate(results.items, 1):
            print(f"   ✓ {i}. {species.name} ({species.atct_id})")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 3: Get covariance matrix
    print("\n3. Get Covariance Matrix:")
    try:
        cov_matrix = get_covariance("67-56-1*0", "7727-37-9*0")  # Methanol and N2
        print(f"   ✓ Species: {cov_matrix.atct_id1} and {cov_matrix.atct_id2}")
        print(f"   ✓ Matrix: {cov_matrix.matrix}")
        if cov_matrix.version:
            print(f"   ✓ ATcT Version: {cov_matrix.version}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    print("\n" + "=" * 40)
    print("Test completed!")

if __name__ == "__main__":
    test_current_api()

#!/usr/bin/env python3
"""Example usage of the pyATcT package."""

import atct
from atct import get_species, search_species, get_covariance, ATCTError, NotFound

def main():
    """Demonstrate basic usage of the atct package."""
    print("ATcT Python Client - Example Usage")
    print("=" * 40)
    
    # Test health check
    print("\n1. Health Check:")
    try:
        healthy = atct.healthcheck()
        print(f"   API is {'healthy' if healthy else 'unhealthy'}")
    except Exception as e:
        print(f"   Health check failed: {e}")
    
    # Example 1: Get species by ATcT ID
    print("\n2. Get Species by ATcT ID:")
    try:
        species = get_species("67-56-1*0")  # Methanol
        print(f"   Name: {species.name}")
        print(f"   Formula: {species.formula}")
        print(f"   ATcT ID: {species.atct_id}")
        if species.delta_h_298k is not None:
            print(f"   ΔHf(298K): {species.delta_h_298k:.2f} kJ/mol")
            if species.delta_h_298k_unc:
                print(f"   Uncertainty: ±{species.delta_h_298k_unc.value:.2f} {species.delta_h_298k_unc.unit}")
    except NotFound:
        print("   Species not found")
    except ATCTError as e:
        print(f"   Error: {e}")
    
    # Example 2: Search species
    print("\n3. Search Species:")
    try:
        results = search_species("methanol", limit=3)
        print(f"   Found {results.total} results (showing {len(results.items)})")
        for i, species in enumerate(results.items, 1):
            print(f"   {i}. {species.name} ({species.atct_id})")
    except ATCTError as e:
        print(f"   Error: {e}")
    
    # Example 3: Get covariance matrix
    print("\n4. Get Covariance Matrix:")
    try:
        cov_matrix = get_covariance("67-56-1*0", "7727-37-9*0")  # Methanol and N2
        print(f"   Species: {cov_matrix.atct_id1} and {cov_matrix.atct_id2}")
        print(f"   Matrix: {cov_matrix.matrix}")
        if cov_matrix.version:
            print(f"   ATcT Version: {cov_matrix.version}")
    except NotFound:
        print("   Covariance data not found")
    except ATCTError as e:
        print(f"   Error: {e}")
    
    print("\n" + "=" * 40)
    print("Example completed!")

if __name__ == "__main__":
    main()

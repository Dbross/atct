#!/usr/bin/env python3
"""Example usage of pyATcT reaction enthalpy calculations."""

import os
from atct.api import (
    get_species,
    create_reaction_calculator,
    calculate_reaction_enthalpy,
    get_species_covariance_by_atctid
)
from atct.api.models import ReactionSpecies, ReactionResult, ReactionCalculator

def main():
    """Demonstrate reaction enthalpy calculations."""
    print("=== pyATcT Reaction Enthalpy Calculations ===\n")
    
    # Set API base URL (optional - defaults to production)
    # os.environ["ATCT_API_BASE_URL"] = "https://atct.anl.gov/api/v1"
    
    # Example 1: Methanol Combustion
    print("1. Methanol Combustion Reaction:")
    print("   CH3OH + 1.5 O2 → CO2 + 2 H2O")
    print()
    
    species_data = {
        '67-56-1*0': -1.0,    # CH3OH (reactant)
        '7727-37-9*0': -1.5,  # O2 (reactant) 
        '124-38-9*0': 1.0,    # CO2 (product)
        '7732-18-5*500': 2.0    # H2O (product)
    }
    
    try:
        # Calculate using different methods
        print("   Calculating reaction enthalpy...")
        
        cov_result = calculate_reaction_enthalpy(species_data, 'covariance')
        ss_result = calculate_reaction_enthalpy(species_data, 'sum_squares')
        conv_result = calculate_reaction_enthalpy(species_data, 'conventional')
        
        print(f"   Covariance method:     {cov_result}")
        print(f"   Sum squares method:    {ss_result}")
        print(f"   Conventional method:   {conv_result}")
        print()
        
        # Compare methods
        calculator = create_reaction_calculator(species_data)
        comparison = calculator.compare_methods()
        
        print(f"   Difference in uncertainty: {comparison['difference']:.6f} kJ/mol")
        print(f"   Statistical significance: {comparison['significance']:.1%} of reference uncertainty")
        print()
        
    except Exception as e:
        print(f"   Error: {e}")
        print()
    
    # Example 2: Original example from extract_reaction_enthalpy_covariance.py
    print("2. Original Example (CH + N2 → HC(NN)):")
    print("   CH + N2 → HC(NN)")
    print()
    
    original_species_data = {
        '3315-37-5*1': -1.0,  # CH (reactant)
        '7727-37-9*0': -1.0,  # N2 (reactant)
        '69967-71-1*2': 1.0   # HC(NN) (product)
    }
    
    try:
        print("   Species data:")
        for atct_id, coeff in original_species_data.items():
            species = get_species(atct_id)
            print(f"   {species.name} ({atct_id}): {species.delta_h_298k} ± {species.delta_h_298k_uncertainty} kJ/mol")
        print()
        
        # Calculate reaction enthalpy
        result = calculate_reaction_enthalpy(original_species_data, 'covariance')
        print(f"   Reaction enthalpy: {result}")
        print()
        
        # Show method comparison
        calculator = create_reaction_calculator(original_species_data)
        comparison = calculator.compare_methods()
        
        print("   Method comparison:")
        print(f"   Covariance method:     {comparison['covariance']}")
        print(f"   Sum of squares method: {comparison['sum_squares']}")
        print(f"   Conventional method:   {comparison['conventional']}")
        print(f"   Difference in uncertainty: {comparison['difference']:.6f} kJ/mol")
        print(f"   Statistical significance: {comparison['significance']:.1%} of reference uncertainty")
        print()
        
    except Exception as e:
        print(f"   Error: {e}")
        print()
    
    # Example 3: Water Formation
    print("3. Water Formation Reaction:")
    print("   H2 + 0.5 O2 → H2O")
    print()
    
    water_formation_data = {
        '1333-74-0*0': -1.0,  # H2 (reactant)
        '7727-37-9*0': -0.5,  # O2 (reactant)
        '7732-18-5*0': 1.0    # H2O (product)
    }
    
    try:
        print("   Species data:")
        for atct_id, coeff in water_formation_data.items():
            species = get_species(atct_id)
            print(f"   {species.name} ({atct_id}): {species.delta_h_298k} ± {species.delta_h_298k_uncertainty} kJ/mol")
        print()
        
        result = calculate_reaction_enthalpy(water_formation_data, 'covariance')
        print(f"   Reaction enthalpy: {result}")
        print()
        
    except Exception as e:
        print(f"   Error: {e}")
        print()
    
    # Example 4: Demonstrate covariance matrix usage
    print("4. Covariance Matrix Example:")
    print("   Getting covariance between methanol phases...")
    
    try:
        cov = get_species_covariance_by_atctid('67-56-1*500', '67-56-1*0')
        print(f"   Covariance matrix: {cov.matrix}")
        print(f"   Labels: {cov.labels}")
        print(f"   Units: {cov.units}")
        print()
        
    except Exception as e:
        print(f"   Error: {e}")
        print()
    
    print("=== All examples completed! ===")

if __name__ == "__main__":
    main()

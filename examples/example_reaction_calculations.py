#!/usr/bin/env python3
"""Example usage of ATcT reaction enthalpy calculations."""

import asyncio
import os
from atct import (
    get_species,
    create_reaction_calculator,
    calculate_reaction_enthalpy,
    get_species_covariance_by_atctid,
    get_species_covariance_matrix
)
from atct import ReactionSpecies, ReactionResult, ReactionCalculator

async def main():
    """Demonstrate reaction enthalpy calculations."""
    print("=== ATcT Reaction Enthalpy Calculations ===\n")
    
    # Set API base URL (optional - defaults to production)
    # os.environ["ATCT_API_BASE_URL"] = "https://atct.anl.gov/api/v1"
    
    # Example 1: Methanol Combustion
    print("1. Methanol Combustion Reaction:")
    print("   CH3OH + 1.5 O2 → CO2 + 2 H2O")
    print()

    species_data = {
    '67-56-1*0':  -1.0,   # CH3OH (g)
    '7782-44-7*0': -1.5,  # O2 (g)  <-- fixed
    '124-38-9*0':  1.0,   # CO2 (g)
    '7732-18-5*500': 2.0  # H2O (cr,l)
    }
    
    
    try:
        # Calculate using different methods (298.15K values)
        print("   Calculating reaction enthalpy (298.15K values)...")
        
        cov_result_298k = await calculate_reaction_enthalpy(species_data, 'covariance', use_0k=False)
        conv_result_298k = await calculate_reaction_enthalpy(species_data, 'conventional', use_0k=False)
        
        print(f"   Covariance method (298.15K):     {cov_result_298k}")
        print(f"   Conventional method (298.15K):   {conv_result_298k}")
        print()
        
        # Try 0K values if available (conventional method only)
        print("   Calculating reaction enthalpy (0K values, conventional method only)...")
        try:
            conv_result_0k = await calculate_reaction_enthalpy(species_data, 'conventional', use_0k=True)
            print(f"   Conventional method (0K):       {conv_result_0k}")
            print()
            
        except ValueError as e:
            print(f"   Cannot use 0K values: {e}")
            print()
        
        # Compare methods (using 298.15K calculator)
        calculator = await create_reaction_calculator(species_data, use_0k=False)
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
        '3315-37-5*1': -1.0,  # CH (g, A2Δ)
        '7727-37-9*0': -1.0,  # N2 (g)
        '69967-71-1*2':  1.0  # HC(NN) (g, 2A')
    }
    
    try:
        print("   Species data:")
        for atct_id, coeff in original_species_data.items():
            species = await get_species(atct_id)
            print(f"   {species.name} ({atct_id}): {species.delta_h_298k} ± {species.delta_h_298k_uncertainty} kJ/mol")
        print()
        
        # Calculate reaction enthalpy (298.15K values)
        result_298k = await calculate_reaction_enthalpy(original_species_data, 'covariance', use_0k=False)
        print(f"   Reaction enthalpy (298.15K): {result_298k}")
        
        # Try 0K values if available (conventional method only)
        try:
            result_0k = await calculate_reaction_enthalpy(original_species_data, 'conventional', use_0k=True)
            print(f"   Reaction enthalpy (0K):      {result_0k}")
        except ValueError as e:
            print(f"   Cannot use 0K values: {e}")
        print()
        
        # Show method comparison
        calculator = await create_reaction_calculator(original_species_data, use_0k=False)
        comparison = calculator.compare_methods()
        
        print("   Method comparison:")
        if comparison['covariance']:
            print(f"   Covariance method:     {comparison['covariance']}")
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
        '7782-44-7*0': -0.5,  # O2 (reactant)
        '7732-18-5*0': 1.0    # H2O (product)
    }
    
    try:
        print("   Species data:")
        for atct_id, coeff in water_formation_data.items():
            species = await get_species(atct_id)
            print(f"   {species.name} ({atct_id}): {species.delta_h_298k} ± {species.delta_h_298k_uncertainty} kJ/mol")
        print()
        
        result_298k = await calculate_reaction_enthalpy(water_formation_data, 'covariance', use_0k=False)
        print(f"   Reaction enthalpy (298.15K): {result_298k}")
        
        # Try 0K values if available (conventional method only)
        try:
            result_0k = await calculate_reaction_enthalpy(water_formation_data, 'conventional', use_0k=True)
            print(f"   Reaction enthalpy (0K):      {result_0k}")
        except ValueError as e:
            print(f"   Cannot use 0K values: {e}")
        print()
        
    except Exception as e:
        print(f"   Error: {e}")
        print()
    
    # Example 4: Benzene Combustion
    print("4. Benzene Combustion Reaction:")
    print("   C6H6 + 7.5 O2 → 6 CO2 + 3 H2O")
    print()
    
    benzene_combustion_data = {
        '71-43-2*0': -1.0,     # C6H6 (reactant)
        '7782-44-7*0': -7.5,   # O2 (reactant)
        '124-38-9*0': 6.0,     # CO2 (product)
        '7732-18-5*500': 3.0     # H2O (product)
    }
    
    try:
        print("   Species data:")
        for atct_id, coeff in benzene_combustion_data.items():
            species = await get_species(atct_id)
            print(f"   {species.name} ({atct_id}): {species.delta_h_298k} ± {species.delta_h_298k_uncertainty} kJ/mol")
        print()
        
        # Calculate using different methods (298.15K values)
        cov_result_298k = await calculate_reaction_enthalpy(benzene_combustion_data, 'covariance', use_0k=False)
        conv_result_298k = await calculate_reaction_enthalpy(benzene_combustion_data, 'conventional', use_0k=False)
        
        print(f"   Covariance method (298.15K):     {cov_result_298k}")
        print(f"   Conventional method (298.15K):   {conv_result_298k}")
        
        # Try 0K values if available (conventional method only)
        try:
            conv_result_0k = await calculate_reaction_enthalpy(benzene_combustion_data, 'conventional', use_0k=True)
            print(f"   Conventional method (0K):       {conv_result_0k}")
        except ValueError as e:
            print(f"   Cannot use 0K values: {e}")
        print()
        
    except Exception as e:
        print(f"   Error: {e}")
        print()
    
    # Example 5: Benzene Ionization
    print("5. Benzene Ionization Reaction:")
    print("   C6H6 → C6H6+ + e-")
    print()
    
    benzene_ionization_data = {
        '71-43-2*0': -1.0,      # C6H6 (reactant)
        '34504-50-2*0': 1.0,    # C6H6+ (product)
    }
    
    try:
        print("   Species data:")
        for atct_id, coeff in benzene_ionization_data.items():
            species = await get_species(atct_id)
            print(f"   {species.name} ({atct_id}): {species.delta_h_298k} ± {species.delta_h_298k_uncertainty} kJ/mol")
        print()
        
        # Calculate using different methods (298.15K values)
        cov_result_298k = await calculate_reaction_enthalpy(benzene_ionization_data, 'covariance', use_0k=False)
        conv_result_298k = await calculate_reaction_enthalpy(benzene_ionization_data, 'conventional', use_0k=False)
        
        print(f"   Covariance method (298.15K):     {cov_result_298k}")
        print(f"   Conventional method (298.15K):   {conv_result_298k}")
        
        # Try 0K values if available (conventional method only)
        try:
            conv_result_0k = await calculate_reaction_enthalpy(benzene_ionization_data, 'conventional', use_0k=True)
            print(f"   Conventional method (0K):       {conv_result_0k}")
        except ValueError as e:
            print(f"   Cannot use 0K values: {e}")
        print()
        
    except Exception as e:
        print(f"   Error: {e}")
        print()
    
    # Example 6: Demonstrate covariance matrix usage
    print("6. Covariance Matrix Example:")
    print("   Getting covariance between methanol phases...")
    
    try:
        cov = await get_species_covariance_by_atctid('67-56-1*500', '67-56-1*0')
        print(f"   Covariance matrix: {cov.matrix}")
        print(f"   Labels: {cov.labels}")
        print(f"   Units: {cov.units}")
        print()
        
    except Exception as e:
        print(f"   Error: {e}")
        print()
    
    # Example 6b: Full covariance matrix for reaction species
    print("6b. Full Covariance Matrix for Reaction Species:")
    print("   Getting full covariance matrix for methanol combustion species...")
    
    try:
        # Get all species involved in methanol combustion
        reaction_atctids = ['67-56-1*0', '7727-37-9*0', '124-38-9*0', '7732-18-5*500']
        cov_matrix = await get_species_covariance_matrix(atctids=reaction_atctids)
        print(f"   Full covariance matrix shape: {len(cov_matrix.matrix)}x{len(cov_matrix.matrix[0])}")
        print(f"   Species ATcT IDs: {cov_matrix.species_atctids}")
        print(f"   Units: {cov_matrix.units}")
        print(f"   Diagonal (variances): {[cov_matrix.matrix[i][i] for i in range(len(cov_matrix.matrix))]}")
        print()
        
    except Exception as e:
        print(f"   Error: {e}")
        print()
    
    print("=== All examples completed! ===")

if __name__ == "__main__":
    asyncio.run(main())

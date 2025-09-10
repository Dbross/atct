#!/usr/bin/env python3
"""Example of fully async reaction enthalpy calculations using concurrent requests."""

import asyncio
import os
import time
from atct import (
    get_species,
    calculate_reaction_enthalpy,
    get_species_covariance_by_atctid,
    get_species_covariance_matrix
)
from atct import ReactionSpecies, ReactionResult

async def main():
    """Demonstrate fully async reaction enthalpy calculations with concurrent requests."""
    print("=== ATcT Async Reaction Enthalpy Calculations ===\n")
    
    # Set API base URL (optional - defaults to production)
    # os.environ["ATCT_API_BASE_URL"] = "https://atct.anl.gov/api/v1"
    
    # Define all reactions to calculate
    reactions = {
        "Methanol Combustion": {
            "equation": "CH3OH + 1.5 O2 → CO2 + 2 H2O",
            "species_data": {
                '67-56-1*0': -1.0,    # CH3OH (reactant)
                '7727-37-9*0': -1.5,  # O2 (reactant) 
                '124-38-9*0': 1.0,    # CO2 (product)
                '7732-18-5*500': 2.0    # H2O (product)
            }
        },
        "CH + N2 → HC(NN)": {
            "equation": "CH + N2 → HC(NN)",
            "species_data": {
                '3315-37-5*1': -1.0,  # CH (reactant)
                '7727-37-9*0': -1.0,  # N2 (reactant)
                '69967-71-1*2': 1.0   # HC(NN) (product)
            }
        },
        "Water Formation": {
            "equation": "H2 + 0.5 O2 → H2O",
            "species_data": {
                '1333-74-0*0': -1.0,  # H2 (reactant)
                '7727-37-9*0': -0.5,  # O2 (reactant)
                '7732-18-5*0': 1.0    # H2O (product)
            }
        },
        "Benzene Combustion": {
            "equation": "C6H6 + 7.5 O2 → 6 CO2 + 3 H2O",
            "species_data": {
                '71-43-2*0': -1.0,     # C6H6 (reactant)
                '7782-44-7*0': -7.5,   # O2 (reactant)
                '124-38-9*0': 6.0,     # CO2 (product)
                '7732-18-5*500': 3.0     # H2O (product)
            }
        },
        "Benzene Ionization": {
            "equation": "C6H6 → C6H6+ + e-",
            "species_data": {
                '71-43-2*0': -1.0,      # C6H6 (reactant)
                '34504-50-2*0': 1.0,    # C6H6+ (product)
                '183748-02-9*0': 1.0    # e- (product)
            }
        }
    }
    
    print("Firing off all requests concurrently...")
    print("(This may take a moment as we're making many API calls)\n")
    
    start_time = time.time()
    
    # Create all async tasks for species data retrieval
    species_tasks = {}
    for reaction_name, reaction_data in reactions.items():
        species_tasks[reaction_name] = {}
        for atct_id in reaction_data["species_data"].keys():
            species_tasks[reaction_name][atct_id] = get_species(atct_id)
    
    # Create all async tasks for reaction enthalpy calculations
    reaction_tasks = {}
    for reaction_name, reaction_data in reactions.items():
        # Create tasks for both covariance and conventional methods
        reaction_tasks[reaction_name] = {
            'covariance': calculate_reaction_enthalpy(reaction_data["species_data"], 'covariance'),
            'conventional': calculate_reaction_enthalpy(reaction_data["species_data"], 'conventional')
        }
    
    # Create tasks for covariance examples
    covariance_tasks = {
        'methanol_phases': get_species_covariance_by_atctid('67-56-1*500', '67-56-1*0'),
        'methanol_combustion_matrix': get_species_covariance_matrix(atctids=['67-56-1*0', '7727-37-9*0', '124-38-9*0', '7732-18-5*500'])
    }
    
    # Wait for ALL requests to complete
    print("Waiting for all API requests to complete...\n")
    
    # Gather all species data
    species_results = {}
    for reaction_name, tasks in species_tasks.items():
        species_results[reaction_name] = await asyncio.gather(*tasks.values(), return_exceptions=True)
        # Create a mapping from atct_id to species for easier access
        species_results[reaction_name] = dict(zip(tasks.keys(), species_results[reaction_name]))
    
    # Gather all reaction calculations
    reaction_results = {}
    for reaction_name, tasks in reaction_tasks.items():
        reaction_results[reaction_name] = await asyncio.gather(
            tasks['covariance'], 
            tasks['conventional'], 
            return_exceptions=True
        )
    
    # Gather covariance results
    covariance_results = await asyncio.gather(
        covariance_tasks['methanol_phases'],
        covariance_tasks['methanol_combustion_matrix'],
        return_exceptions=True
    )
    
    end_time = time.time()
    print(f"All API requests completed in {end_time - start_time:.2f} seconds!\n")
    
    print("All requests completed! Now processing results...\n")
    
    # Now process and display all results
    for i, (reaction_name, reaction_data) in enumerate(reactions.items(), 1):
        print(f"{i}. {reaction_name}:")
        print(f"   {reaction_data['equation']}")
        print()
        
        try:
            # Display species data
            print("   Species data:")
            species_data = species_results[reaction_name]
            for atct_id, coeff in reaction_data["species_data"].items():
                species = species_data[atct_id]
                if isinstance(species, Exception):
                    print(f"   Error getting {atct_id}: {species}")
                else:
                    print(f"   {species.name} ({atct_id}): {species.delta_h_298k} ± {species.delta_h_298k_uncertainty} kJ/mol")
            print()
            
            # Display reaction enthalpy results
            cov_result, conv_result = reaction_results[reaction_name]
            
            if isinstance(cov_result, Exception):
                print(f"   Covariance method error: {cov_result}")
            else:
                print(f"   Covariance method:     {cov_result}")
            
            if isinstance(conv_result, Exception):
                print(f"   Conventional method error: {conv_result}")
            else:
                print(f"   Conventional method:   {conv_result}")
            print()
            
            # Display method comparison (simplified)
            if not isinstance(cov_result, Exception) and not isinstance(conv_result, Exception):
                diff_unc = abs(cov_result.uncertainty - conv_result.uncertainty)
                ref_unc = max(cov_result.uncertainty, conv_result.uncertainty)
                significance_ratio = diff_unc / ref_unc if ref_unc > 0 else 0
                print(f"   Difference in uncertainty: {diff_unc:.6f} kJ/mol")
                print(f"   Statistical significance: {significance_ratio:.1%} of reference uncertainty")
            print()
            
        except Exception as e:
            print(f"   Error processing {reaction_name}: {e}")
            print()
    
    # Display covariance examples
    print("6. Covariance Matrix Examples:")
    print("   Getting covariance between methanol phases...")
    
    methanol_cov = covariance_results[0]
    if isinstance(methanol_cov, Exception):
        print(f"   Error: {methanol_cov}")
    else:
        print(f"   Covariance matrix: {methanol_cov.matrix}")
        print(f"   Labels: {methanol_cov.labels}")
        print(f"   Units: {methanol_cov.units}")
    print()
    
    print("6b. Full Covariance Matrix for Methanol Combustion Species:")
    print("   Getting full covariance matrix for methanol combustion species...")
    
    full_cov_matrix = covariance_results[1]
    if isinstance(full_cov_matrix, Exception):
        print(f"   Error: {full_cov_matrix}")
    else:
        print(f"   Full covariance matrix shape: {len(full_cov_matrix.matrix)}x{len(full_cov_matrix.matrix[0])}")
        print(f"   Species ATcT IDs: {full_cov_matrix.species_atctids}")
        print(f"   Units: {full_cov_matrix.units}")
        print(f"   Diagonal (variances): {[full_cov_matrix.matrix[i][i] for i in range(len(full_cov_matrix.matrix))]}")
    print()
    
    print("=== All async calculations completed! ===")
    print("Note: All API requests were fired concurrently and results were processed only after all completed.")
    print(f"Total execution time: {end_time - start_time:.2f} seconds")
    
    # Demonstrate the scaling benefit by showing how many requests we made
    total_requests = 0
    for reaction_data in reactions.values():
        total_requests += len(reaction_data["species_data"])  # Species requests
        total_requests += 2  # Two reaction calculation methods per reaction
    total_requests += 2  # Covariance requests
    
    print(f"Total API requests made: {total_requests}")
    print(f"Average time per request: {(end_time - start_time) / total_requests:.3f} seconds")
    print("\nTo see the real benefit of async, try running this with many more reactions!")
    print("The async version will scale much better as the number of requests increases.")

if __name__ == "__main__":
    asyncio.run(main())

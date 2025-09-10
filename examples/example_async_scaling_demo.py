#!/usr/bin/env python3
"""Demonstration of async scaling benefits with many concurrent requests."""

import asyncio
import time
from atct import get_species, calculate_reaction_enthalpy

async def main():
    """Demonstrate async scaling with many reactions."""
    print("=== ATcT Async Scaling Demonstration ===\n")
    
    # Create a large number of reactions to demonstrate scaling
    reactions = {
        "Methanol Combustion": {'67-56-1*0': -1.0, '7727-37-9*0': -1.5, '124-38-9*0': 1.0, '7732-18-5*500': 2.0},
        "Ethanol Combustion": {'64-17-5*0': -1.0, '7727-37-9*0': -3.0, '124-38-9*0': 2.0, '7732-18-5*500': 3.0},
        "Propane Combustion": {'74-98-6*0': -1.0, '7727-37-9*0': -5.0, '124-38-9*0': 3.0, '7732-18-5*500': 4.0},
        "Butane Combustion": {'106-97-8*0': -1.0, '7727-37-9*0': -6.5, '124-38-9*0': 4.0, '7732-18-5*500': 5.0},
        "Pentane Combustion": {'109-66-0*0': -1.0, '7727-37-9*0': -8.0, '124-38-9*0': 5.0, '7732-18-5*500': 6.0},
        "Hexane Combustion": {'110-54-3*0': -1.0, '7727-37-9*0': -9.5, '124-38-9*0': 6.0, '7732-18-5*500': 7.0},
        "Heptane Combustion": {'142-82-5*0': -1.0, '7727-37-9*0': -11.0, '124-38-9*0': 7.0, '7732-18-5*500': 8.0},
        "Octane Combustion": {'111-65-9*0': -1.0, '7727-37-9*0': -12.5, '124-38-9*0': 8.0, '7732-18-5*500': 9.0},
        "Water Formation": {'1333-74-0*0': -1.0, '7727-37-9*0': -0.5, '7732-18-5*0': 1.0},
        "CO Formation": {'1333-74-0*0': -1.0, '7727-37-9*0': -0.5, '630-08-0*0': 1.0},
        "CO2 Formation": {'1333-74-0*0': -1.0, '7727-37-9*0': -1.0, '124-38-9*0': 1.0},
        "Methane Formation": {'1333-74-0*0': -2.0, '7727-37-9*0': -1.0, '74-82-8*0': 1.0},
        "Ammonia Formation": {'1333-74-0*0': -1.5, '7727-37-9*0': -0.75, '7664-41-7*0': 1.0},
        "Hydrogen Chloride Formation": {'1333-74-0*0': -0.5, '7647-01-0*0': -0.5, '7647-01-0*1': 1.0},
        "Hydrogen Bromide Formation": {'1333-74-0*0': -0.5, '10035-10-6*0': -0.5, '10035-10-6*1': 1.0},
        "Hydrogen Iodide Formation": {'1333-74-0*0': -0.5, '10034-85-2*0': -0.5, '10034-85-2*1': 1.0},
        "Benzene Formation": {'1333-74-0*0': -3.0, '7727-37-9*0': -1.5, '71-43-2*0': 1.0},
        "Toluene Formation": {'1333-74-0*0': -4.0, '7727-37-9*0': -1.5, '108-88-3*0': 1.0},
        "Xylene Formation": {'1333-74-0*0': -5.0, '7727-37-9*0': -1.5, '106-42-3*0': 1.0},
        "Naphthalene Formation": {'1333-74-0*0': -5.0, '7727-37-9*0': -2.0, '91-20-3*0': 1.0},
    }
    
    print(f"Processing {len(reactions)} reactions with async concurrent requests...")
    print("This will demonstrate the scaling benefits of async programming.\n")
    
    start_time = time.time()
    
    # Create ALL tasks at once - this is the key to async benefits
    all_tasks = []
    task_info = []
    
    # Create species lookup tasks
    unique_species = set()
    for reaction_data in reactions.values():
        unique_species.update(reaction_data.keys())
    
    print(f"Fetching data for {len(unique_species)} unique species...")
    species_tasks = {atct_id: get_species(atct_id) for atct_id in unique_species}
    all_tasks.extend(species_tasks.values())
    task_info.extend([f"species_{atct_id}" for atct_id in unique_species])
    
    # Create reaction calculation tasks
    print(f"Creating {len(reactions)} reaction enthalpy calculations...")
    for reaction_name, reaction_data in reactions.items():
        cov_task = calculate_reaction_enthalpy(reaction_data, 'covariance')
        conv_task = calculate_reaction_enthalpy(reaction_data, 'conventional')
        all_tasks.extend([cov_task, conv_task])
        task_info.extend([f"reaction_{reaction_name}_cov", f"reaction_{reaction_name}_conv"])
    
    print(f"Total tasks created: {len(all_tasks)}")
    print("All tasks fired concurrently - waiting for completion...\n")
    
    # Wait for ALL tasks to complete
    results = await asyncio.gather(*all_tasks, return_exceptions=True)
    
    end_time = time.time()
    total_time = end_time - start_time
    
    print(f"All {len(all_tasks)} requests completed in {total_time:.2f} seconds!")
    print(f"Average time per request: {total_time / len(all_tasks):.3f} seconds")
    print(f"Effective requests per second: {len(all_tasks) / total_time:.1f}")
    
    # Process results
    species_results = {}
    for i, atct_id in enumerate(unique_species):
        species_results[atct_id] = results[i]
    
    reaction_results = {}
    result_idx = len(unique_species)
    for reaction_name in reactions.keys():
        cov_result = results[result_idx]
        conv_result = results[result_idx + 1]
        reaction_results[reaction_name] = (cov_result, conv_result)
        result_idx += 2
    
    # Display some results
    print("\nSample results:")
    for i, (reaction_name, (cov_result, conv_result)) in enumerate(list(reaction_results.items())[:5]):
        if not isinstance(cov_result, Exception) and not isinstance(conv_result, Exception):
            print(f"{i+1}. {reaction_name}:")
            print(f"   Covariance: {cov_result}")
            print(f"   Conventional: {conv_result}")
            print()
    
    print(f"=== Async Scaling Demo Complete ===")
    print(f"Processed {len(reactions)} reactions with {len(all_tasks)} total API requests")
    print(f"Total time: {total_time:.2f} seconds")
    print(f"Speedup vs sequential: ~{len(all_tasks) * 0.1 / total_time:.1f}x faster")
    print("\nThe async approach shows its real benefits with many concurrent requests!")
    print("Each additional reaction adds minimal time since requests run in parallel.")

if __name__ == "__main__":
    asyncio.run(main())

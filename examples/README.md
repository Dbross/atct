# atct Examples

This directory contains example scripts demonstrating how to use the atct library.

## Examples

### `example_usage.py`
Basic async usage examples of the atct library, including:
- Health check
- Species retrieval by ATcT ID
- Species search
- Species lookup by SMILES, CASRN, formula, name, InChI, and InChI Key
- Covariance matrix retrieval (both 2x2 and full N×N matrices)
- Pandas DataFrame conversion

### `example_blocking_usage.py`
Basic blocking usage examples of the atct library, demonstrating:
- How to use the `block=True` parameter for synchronous behavior
- Same functionality as `example_usage.py` but with blocking requests
- Useful for simple scripts that don't need async/await

### `example_reaction_calculations.py`
Advanced async reaction enthalpy calculation examples, including:
- Methanol combustion reaction
- CH + N2 → HC(NN) reaction example
- Water formation reaction
- **Benzene combustion reaction** (C6H6 + 7.5 O2 → 6 CO2 + 3 H2O)
- **Benzene ionization reaction** (C6H6 → C6H6+ + e-)
- Covariance matrix usage (both 2x2 and full N×N matrices)
- Comparison of covariance vs conventional uncertainty propagation methods

### `example_xyz_helpers.py`
XYZ coordinate helper examples for third-party library integration, including:
- Fetching species with XYZ data from the API
- Converting species XYZ data to ASE Atoms objects
- Converting to RDKit Mol objects
- Converting to pymatgen Molecule objects
- File I/O operations for XYZ data
- Error handling for missing dependencies
- Real-world usage patterns with actual API data

## Running the Examples

Make sure you have the atct package installed and the environment variable set:

```bash
# Install the package
pip install -e .

# Set the API base URL (optional - defaults to production)
export ATCT_API_BASE_URL="https://atct.anl.gov/api/v1"

# Run examples
python examples/example_usage.py              # Async examples
python examples/example_blocking_usage.py     # Blocking examples
python examples/example_reaction_calculations.py  # Async reaction calculations
python examples/example_xyz_helpers.py        # XYZ helper examples
```

## Async vs Blocking Usage

The atct library now supports both async and blocking usage patterns:

### Async Usage (Default)
```python
import asyncio
from atct import get_species, search_species

async def main():
    # These return asyncio.Task objects by default
    species_task = get_species("67-56-1*0")
    search_task = search_species("methanol")
    
    # Await the results
    species = await species_task
    results = await search_task
    
    print(f"Species: {species.name}")
    print(f"Found {len(results.items)} results")

asyncio.run(main())
```

### Blocking Usage
```python
from atct import get_species, search_species

def main():
    # Use block=True for synchronous behavior
    species = get_species("67-56-1*0", block=True)
    results = search_species("methanol", block=True)
    
    print(f"Species: {species.name}")
    print(f"Found {len(results.items)} results")

main()
```

## Available Functions

The atct library provides the following main functions. All functions support both async and blocking usage:

### Species Retrieval
- `get_species(atctid, block=False)` - Get species by ATcT ID
- `get_species_by_atctid(atctid, block=False)` - Alternative name for get_species
- `get_species_by_casrn(casrn, block=False)` - Get species by CAS Registry Number
- `get_species_by_smiles(smiles, block=False)` - Get species by SMILES notation
- `get_species_by_formula(formula, block=False)` - Get species by chemical formula
- `get_species_by_name(name, block=False)` - Get species by name
- `get_species_by_inchi(inchi, block=False)` - Get species by InChI
- `get_species_by_inchikey(inchikey, block=False)` - Get species by InChI Key
- `search_species(query, block=False)` - Search species by any identifier
- `get_all_species(block=False)` - Get all species from database

### Covariance Matrices
- `get_species_covariance_by_atctid(atctid1, atctid2, block=False)` - Get 2x2 covariance matrix
- `get_species_covariance_matrix(atctids=[], block=False)` - Get full N×N covariance matrix
- `get_species_covariance_by_ids(id1, id2, block=False)` - Get 2x2 covariance by internal IDs

### Reaction Calculations
- `create_reaction_calculator(species_data, block=False)` - Create reaction calculator
- `calculate_reaction_enthalpy(species_data, method, block=False)` - Calculate reaction enthalpy
- `healthcheck(block=False)` - Check API health
- `as_dataframe(species_list)` - Convert to pandas DataFrame

### XYZ Helper Functions
- `species.get_xyz_data()` - Get XYZ data as XYZData object
- `species.to_ase_atoms()` - Convert to ASE Atoms object
- `species.to_rdkit_mol()` - Convert to RDKit Mol object
- `species.to_pymatgen_molecule()` - Convert to pymatgen Molecule object
- `species.save_xyz_to_file(filename)` - Save XYZ data to file
- `species.get_xyz_string()` - Get XYZ data as formatted string

### Async-Only Functions
For advanced usage, you can also use the async-only versions:
- `*_async()` versions of all functions (e.g., `get_species_async()`, `search_species_async()`)

## Requirements

- Python 3.9+
- httpx (for async HTTP requests)
- numpy (for reaction calculations)
- pandas (optional, for DataFrame support)

### XYZ Helper Dependencies
The XYZ helper functions support integration with third-party libraries:
- **ASE** (optional): `pip install ase` - for atomic simulation
- **RDKit** (optional): `pip install rdkit` - for cheminformatics (may have limitations with some molecular structures, may show deprecation warnings)
- **pymatgen** (optional): `pip install pymatgen` - for materials science

Install with optional dependencies:
```bash
pip install atct[pandas]  # For DataFrame support
pip install ase           # For ASE integration
pip install rdkit         # For RDKit integration
pip install pymatgen      # For pymatgen integration
```

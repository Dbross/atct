# pyATcT Examples

This directory contains example scripts demonstrating how to use the pyATcT library.

## Examples

### `example_usage.py`
Basic usage examples of the pyATcT library, including:
- Health check
- Species retrieval by ATcT ID
- Species search
- Species lookup by SMILES, CASRN, formula, name, InChI, and InChI Key
- Covariance matrix retrieval (both 2x2 and full N×N matrices)
- Pandas DataFrame conversion

### `example_reaction_calculations.py`
Advanced reaction enthalpy calculation examples, including:
- Methanol combustion reaction
- CH + N2 → HC(NN) reaction example
- Water formation reaction
- **Benzene combustion reaction** (C6H6 + 7.5 O2 → 6 CO2 + 3 H2O)
- **Benzene ionization reaction** (C6H6 → C6H6+ + e-)
- Covariance matrix usage (both 2x2 and full N×N matrices)
- Comparison of covariance vs conventional uncertainty propagation methods

## Running the Examples

Make sure you have the pyATcT package installed and the environment variable set:

```bash
# Install the package
pip install -e .

# Set the API base URL (optional - defaults to production)
export ATCT_API_BASE_URL="https://atct.anl.gov/api/v1"

# Run examples
python examples/example_usage.py
python examples/example_reaction_calculations.py
```

## Available Functions

The pyATcT library provides the following main functions:

### Species Retrieval
- `get_species(atctid)` - Get species by ATcT ID
- `get_species_by_atctid(atctid)` - Alternative name for get_species
- `get_species_by_casrn(casrn)` - Get species by CAS Registry Number
- `get_species_by_smiles(smiles)` - Get species by SMILES notation
- `get_species_by_formula(formula)` - Get species by chemical formula
- `get_species_by_name(name)` - Get species by name
- `get_species_by_inchi(inchi)` - Get species by InChI
- `get_species_by_inchikey(inchikey)` - Get species by InChI Key
- `search_species(query)` - Search species by any identifier
- `get_all_species()` - Get all species from database

### Covariance Matrices
- `get_species_covariance_by_atctid(atctid1, atctid2)` - Get 2x2 covariance matrix
- `get_species_covariance_matrix(atctids=[])` - Get full N×N covariance matrix
- `get_species_covariance_by_ids(id1, id2)` - Get 2x2 covariance by internal IDs

### Reaction Calculations
- `create_reaction_calculator(species_data)` - Create reaction calculator
- `calculate_reaction_enthalpy(species_data, method)` - Calculate reaction enthalpy
- `healthcheck()` - Check API health
- `as_dataframe(species_list)` - Convert to pandas DataFrame

## Requirements

- Python 3.9+
- numpy (for reaction calculations)
- pandas (optional, for DataFrame support)

Install with optional dependencies:
```bash
pip install pyATcT[pandas]
```

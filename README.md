# pyATcT - Pythonic ATcT v1 API Client

A comprehensive Python client for the ATcT v1 API with advanced reaction enthalpy calculation capabilities.

## Features

- **Complete v1 API coverage** - all endpoints implemented with proper error handling
- **Reaction calculations** - advanced uncertainty propagation with multiple methods
- **Zero dependencies** - works with only Python standard library
- **Optional enhancements** - numpy for advanced calculations, pandas for DataFrames
- **Dual import support** - use `from atct.api import ...` or `from atct import ...`
- **Robust error handling** - comprehensive exception mapping and retry logic
- **Environment flexibility** - local/production URL switching

## Installation

```bash
# Basic installation (no external dependencies)
pip install atct

# With pandas support for DataFrames
pip install atct[pandas]

# With numpy support for advanced reaction calculations
pip install atct[numpy]

# With both pandas and numpy
pip install atct[all]
```

## Quickstart

```python
from atct.api import get_species, search_species, calculate_reaction_enthalpy

# Get species data
species = get_species("67-56-1*0")  # Methanol
print(f"{species.name}: {species.delta_h_298k} ± {species.delta_h_298k_uncertainty} kJ/mol")

# Search for species
results = search_species("methanol", limit=5)
for item in results.items:
    print(f"{item.name} ({item.atct_id})")

# Calculate reaction enthalpy
species_data = {
    '67-56-1*0': -1.0,    # CH3OH (reactant)
    '7727-37-9*0': -1.5,  # O2 (reactant) 
    '124-38-9*0': 1.0,    # CO2 (product)
    '7732-18-5*0': 2.0    # H2O (product)
}
result = calculate_reaction_enthalpy(species_data, 'covariance')
print(f"Reaction enthalpy: {result}")
```

## API Functions

### Species Data
- `get_species(atctid, expand_xyz=False)` - Get species by ATcT ID
- `search_species(q, limit=50, offset=0)` - Search species by name or formula
- `get_species_by_smiles(smiles)` - Get species by SMILES string
- `get_species_by_casrn(casrn)` - Get species by CAS Registry Number
- `get_species_by_inchi(inchi)` - Get species by InChI
- `get_species_by_formula(formula, phase=None, descriptor=None)` - Get species by formula
- `get_species_by_name(name)` - Get species by name

### Covariance Data
- `get_species_covariance_by_atctid(a_atctid, b_atctid)` - Get 2x2 covariance matrix
- `get_species_covariance_by_ids(a_id, b_id)` - Get covariance by numeric IDs

### Reaction Calculations
- `calculate_reaction_enthalpy(species_data, method)` - Calculate reaction enthalpy
- `create_reaction_calculator(species_data)` - Create reaction calculator

**Note:** The `covariance` method requires numpy (`pip install atct[numpy]`). 
The `conventional` method works without external dependencies.

### Utility
- `healthcheck()` - Check API health status
- `as_dataframe(items)` - Convert results to pandas DataFrame

## Data Models

### Species
```python
@dataclass
class Species:
    atct_tn_version: Optional[str]
    atct_id: str
    name: Optional[str]
    formula: Optional[str]
    delta_h_0k: Optional[str]
    delta_h_298k: Optional[str]
    delta_h_298k_uncertainty: Optional[str]
    unit: Optional[str]
    mass: Optional[str]
    mass_uncertainty: Optional[str]
    smiles: Optional[str]
    casrn: Optional[str]
    charge: Optional[int]
    xyz: Optional[Union[str, List[str]]]
```

### Covariance2x2
```python
@dataclass
class Covariance2x2:
    labels: List[str]
    units: str
    matrix: List[List[float]]  # 2x2 matrix
```

### Reaction Models
```python
@dataclass
class ReactionSpecies:
    species: Species
    stoichiometry: float

@dataclass
class ReactionResult:
    delta_h: float
    uncertainty: float
    method: str
    units: str = "kJ/mol"
```

## Error Handling

The package provides specific exceptions for different error conditions:

```python
from atct import ATCTError, NotFound, BadRequest, Unauthorized, ServerError, NetworkError

try:
    species = get_species("invalid-id")
except NotFound:
    print("Species not found")
except NetworkError:
    print("Network connection failed")
except ATCTError as e:
    print(f"API error: {e}")
```

## Configuration

Set environment variables to configure the client:

- `ATCT_API_BASE_URL` (default: `https://atct.anl.gov/api/v1`) - API endpoint
- `ATCT_API_KEY` (optional, for authenticated requests)
- `ATCT_TIMEOUT_S` (default: `10` seconds)
- `ATCT_RETRIES` (default: `3` retry attempts)
- `ATCT_USER_AGENT` (default: `atct/1.0.0`)

The client uses the v1 API exclusively.

## Development

```bash
# Install in editable mode
pip install -e .

# Run tests
pytest -q

# Run examples
python examples/example_usage.py
python examples/example_reaction_calculations.py
```

## Examples

See the `examples/` directory for comprehensive usage examples:
- `example_usage.py` - Basic API functionality
- `example_reaction_calculations.py` - Advanced reaction calculations

## v1 API Endpoints

This package targets the ATcT v1 API exclusively:

- `GET /api/v1/health` - Health check
- `GET /api/v1/species/get/by-atctid/` - Get species by ATcT ID
- `GET /api/v1/species/search/` - Search species
- `GET /api/v1/species/get/by-smiles/` - Get species by SMILES
- `GET /api/v1/species/get/by-casrn/` - Get species by CASRN
- `GET /api/v1/species/get/by-inchi/` - Get species by InChI
- `GET /api/v1/species/get/by-formula/` - Get species by formula
- `GET /api/v1/species/get/by-name/` - Get species by name
- `GET /api/v1/covariance/species/` - Get covariance matrices

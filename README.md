# atct - Active Thermochemical Tables Python Client

A minimal-dependency Python client for the ATcT v1 API, using only standard library components.

## Features

- **Zero external dependencies** - uses only Python standard library (urllib, json, dataclasses)
- **Pythonic API** - clean, simple functions like `atct.get_species("CH4")`
- **Robust error handling** - comprehensive exception mapping and retry logic
- **v1 API support** - targets the new stable v1 API at `https://atct.anl.gov/api/v1`
- **2x2 covariance matrices** - proper covariance matrix format for uncertainty propagation

## Installation

```bash
pip install atct
# or with pandas support for DataFrames
pip install atct[pandas]
```

## Quickstart

```python
import atct
from atct import get_species, search_species, get_covariance

# Get species data
species = get_species("67-56-1*0")  # Methanol
print(f"{species.name}: {species.delta_h_298k:.2f} ± {species.delta_h_298k_unc.value:.2f} kJ/mol")

# Search for species
results = search_species("methanol", limit=5)
for item in results.items:
    print(f"{item.name} ({item.atct_id})")

# Get covariance matrix between two species
cov_matrix = get_covariance("67-56-1*0", "7727-37-9*0")
print(f"2x2 covariance matrix: {cov_matrix.matrix}")
```

## API Functions

### Species Data
- `get_species(identifier)` - Get species by ATcT ID, formula, or name
- `search_species(query, limit=50, offset=0)` - Search species by name or formula
- `get_species_by_smiles(smiles)` - Get species by SMILES string
- `get_species_by_casrn(casrn)` - Get species by CAS Registry Number

### Covariance Data
- `get_covariance(atct_id1, atct_id2)` - Get 2x2 covariance matrix between two species

### Utility
- `healthcheck()` - Check API health status

## Data Models

### Species
```python
@dataclass
class Species:
    atct_id: str
    name: str
    formula: str
    phase: Optional[str] = None
    delta_h_298k: Optional[float] = None
    delta_h_298k_unc: Optional[Uncertainty] = None
    smiles: Optional[str] = None
    mass: Optional[float] = None
```

### CovarianceMatrix
```python
@dataclass
class CovarianceMatrix:
    atct_id1: str
    atct_id2: str
    matrix: List[List[float]]  # 2x2 matrix [[var1, cov], [cov, var2]]
    version: Optional[str] = None
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

- `ATCT_API_BASE_URL` (default: `https://atct.anl.gov/api/v1`)
- `ATCT_API_KEY` (optional, for authenticated requests)
- `ATCT_TIMEOUT_S` (default: `10` seconds)
- `ATCT_RETRIES` (default: `3` retry attempts)
- `ATCT_USER_AGENT` (default: `atct/1.0.0`)

## Development

```bash
# Install in editable mode
pip install -e .

# Run tests
pytest -q

# Run example
python example_usage.py
```

## v1 API Server

The package includes v1 API server endpoints (PHP) that can be deployed:

- `GET /api/v1/health` - Health check
- `GET /api/v1/species/{id}` - Get species by ATcT ID
- `GET /api/v1/species?q={query}` - Search species
- `GET /api/v1/covariance/{id1}/{id2}` - Get 2x2 covariance matrix

These endpoints provide a clean, RESTful interface with proper error handling and the new 2x2 covariance matrix format.

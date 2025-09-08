# pyATcT Examples

This directory contains example scripts demonstrating how to use the pyATcT library.

## Examples

### `example_usage.py`
Basic usage examples of the pyATcT library, including:
- Health check
- Species retrieval by ATcT ID
- Species search
- Species lookup by SMILES, CASRN, formula, name, and InChI
- Covariance matrix retrieval
- Pandas DataFrame conversion

### `example_reaction_calculations.py`
Advanced reaction enthalpy calculation examples, including:
- Methanol combustion reaction
- CH + N2 → HC(NN) reaction example
- Water formation reaction
- **Benzene combustion reaction** (C6H6 + 7.5 O2 → 6 CO2 + 3 H2O)
- **Benzene ionization reaction** (C6H6 → C6H6+ + e-)
- Covariance matrix usage
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

## Requirements

- Python 3.9+
- numpy (for reaction calculations)
- pandas (optional, for DataFrame support)

Install with optional dependencies:
```bash
pip install pyATcT[pandas]
```

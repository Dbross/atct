from __future__ import annotations
from typing import Optional, List
from ._config import settings
from ._http import request_json
from .models import Species, CovarianceMatrix, Page
from .exceptions import NotFound

def healthcheck() -> bool:
    """Check if the API is healthy."""
    url = f"{settings.base_url}/health"
    try:
        data = request_json("GET", url)
        return bool(data.get("ok", True))
    except Exception:
        return False

def get_species(identifier: str) -> Species:
    """Fetch a single species by ATcT ID, formula, or name.
    
    Args:
        identifier: ATcT ID (e.g., "67-56-1*0"), formula, or name
        
    Returns:
        Species object with thermodynamic data
        
    Raises:
        NotFound: If species not found
        ATCTError: For other API errors
    """
    # Use current API endpoint for now
    url = "https://atct.anl.gov/api/"
    params = {"atctid": identifier}
    data = request_json("GET", url, params=params)
    
    if not data or len(data) == 0:
        raise NotFound(f"No species found for identifier: {identifier}")
    
    return Species.from_dict(data[0])

def search_species(query: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Search for species by name or formula.
    
    Args:
        query: Search term (name or formula)
        limit: Maximum number of results (default 50)
        offset: Number of results to skip (default 0)
        
    Returns:
        Page object with search results
    """
    # Use current API endpoint for now
    url = "https://atct.anl.gov/api/"
    params = {"name": query}
    data = request_json("GET", url, params=params)
    
    if not data:
        data = []
    
    # Convert to Species objects
    items = [Species.from_dict(item) for item in data]
    
    # Apply pagination
    total = len(items)
    start = offset
    end = min(offset + limit, total)
    items = items[start:end]
    
    return Page(items=items, total=total, limit=limit, offset=offset)

def get_covariance(atct_id1: str, atct_id2: str) -> CovarianceMatrix:
    """Get covariance matrix between two species.
    
    Args:
        atct_id1: First species ATcT ID
        atct_id2: Second species ATcT ID
        
    Returns:
        CovarianceMatrix object with 2x2 matrix
        
    Raises:
        NotFound: If covariance data not found
        ATCTError: For other API errors
    """
    # Use current API endpoint for now
    url = "https://atct.anl.gov/api/covariance/"
    params = {"atctid1": atct_id1, "atctid2": atct_id2}
    data = request_json("GET", url, params=params)
    
    if not data or len(data) == 0:
        raise NotFound(f"No covariance data found for species {atct_id1} and {atct_id2}")
    
    return CovarianceMatrix.from_dict(data[0])

def get_species_by_smiles(smiles: str) -> Species:
    """Fetch species by SMILES string.
    
    Args:
        smiles: SMILES molecular representation
        
    Returns:
        Species object
        
    Raises:
        NotFound: If species not found
    """
    url = "https://atct.anl.gov/api/"
    params = {"smiles": smiles}
    data = request_json("GET", url, params=params)
    
    if not data or len(data) == 0:
        raise NotFound(f"No species found for SMILES: {smiles}")
    
    return Species.from_dict(data[0])

def get_species_by_casrn(casrn: str) -> Species:
    """Fetch species by CAS Registry Number.
    
    Args:
        casrn: CAS Registry Number
        
    Returns:
        Species object
        
    Raises:
        NotFound: If species not found
    """
    url = "https://atct.anl.gov/api/"
    params = {"casrn": casrn}
    data = request_json("GET", url, params=params)
    
    if not data or len(data) == 0:
        raise NotFound(f"No species found for CASRN: {casrn}")
    
    return Species.from_dict(data[0])

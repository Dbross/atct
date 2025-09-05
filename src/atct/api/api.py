from __future__ import annotations
from typing import Optional, List
from ._config import settings
from ._http import request_json
from .models import Species, CovarianceMatrix, Page
from .exceptions import NotFound, ATCTError

# API endpoint builder
def _api_url(endpoint: str) -> str:
    """Build API URL."""
    return f"{settings.base_url}/{endpoint.lstrip('/')}"

def healthcheck() -> bool:
    """Check if the API is healthy."""
    url = _api_url("health")
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
    url = _api_url(f"species/{identifier}")
    data = request_json("GET", url)
    return Species.from_dict(data)

def search_species(query: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Search for species by name or formula.
    
    Args:
        query: Search term (name or formula)
        limit: Maximum number of results (default 50)
        offset: Number of results to skip (default 0)
        
    Returns:
        Page object with search results
    """
    url = _api_url("species")
    params = {"q": query, "limit": limit, "offset": offset}
    data = request_json("GET", url, params=params)
    
    items = [Species.from_dict(item) for item in data.get("items", [])]
    return Page(
        items=items, 
        total=data.get("total", len(items)), 
        limit=limit, 
        offset=offset
    )

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
    url = _api_url(f"covariance/{atct_id1}/{atct_id2}")
    data = request_json("GET", url)
    return CovarianceMatrix.from_dict(data)

def get_species_by_smiles(smiles: str) -> Species:
    """Fetch species by SMILES string.
    
    Args:
        smiles: SMILES molecular representation
        
    Returns:
        Species object
        
    Raises:
        NotFound: If species not found
    """
    url = _api_url("species")
    params = {"smiles": smiles}
    data = request_json("GET", url, params=params)
    
    if data.get("items"):
        return Species.from_dict(data["items"][0])
    else:
        raise NotFound(f"No species found for SMILES: {smiles}")

def get_species_by_casrn(casrn: str) -> Species:
    """Fetch species by CAS Registry Number.
    
    Args:
        casrn: CAS Registry Number
        
    Returns:
        Species object
        
    Raises:
        NotFound: If species not found
    """
    url = _api_url("species")
    params = {"casrn": casrn}
    data = request_json("GET", url, params=params)
    
    if data.get("items"):
        return Species.from_dict(data["items"][0])
    else:
        raise NotFound(f"No species found for CASRN: {casrn}")

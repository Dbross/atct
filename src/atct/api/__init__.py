from __future__ import annotations
from typing import Optional, Dict, Any, List
from ._config import settings
from ._http import request_json
from .models import Species, Page, CovarianceMatrix, Covariance2x2, ReactionSpecies, ReactionResult, ReactionCalculator
from .pandas_io import as_dataframe

# -------- health --------
def healthcheck() -> bool:
    try:
        data = request_json("GET", f"{settings.base_url}/health/")
        return bool(data.get("ok", True))
    except Exception:
        return False

# -------- species: unique by ATcT_ID --------
def get_species_by_atctid(atctid: str, *, expand_xyz: bool = False) -> Species:
    """Get species by ATcT ID. Returns single species."""
    params: Dict[str, Any] = {"atctid": atctid}
    if expand_xyz:
        params["expand_xyz"] = "1"
    data = request_json("GET", f"{settings.base_url}/species/get/by-atctid/", params=params)
    # API returns list, take first item
    if isinstance(data, list) and len(data) > 0:
        return Species.from_dict(data[0])
    elif isinstance(data, dict):
        return Species.from_dict(data)
    else:
        raise ValueError("No species found with the given ATcT ID")

# Legacy alias
def get_species(atctid: str, *, expand_xyz: bool = False) -> Species:
    """Legacy alias for get_species_by_atctid."""
    return get_species_by_atctid(atctid, expand_xyz=expand_xyz)

# -------- species: multi "get by …" (paged) --------
def get_species_by_casrn(casrn: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Get species by CAS Registry Number. Returns paged results."""
    data = request_json("GET", f"{settings.base_url}/species/get/by-casrn/", params={"casrn": casrn, "limit": limit, "offset": offset})
    # Handle both list and paginated response formats
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_inchi(inchi: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Get species by InChI. Returns paged results."""
    data = request_json("GET", f"{settings.base_url}/species/get/by-inchi/", params={"inchi": inchi, "limit": limit, "offset": offset})
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_inchikey(inchikey: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Get species by InChI Key. Returns paged results."""
    data = request_json("GET", f"{settings.base_url}/species/get/by-inchikey/", params={"inchikey": inchikey, "limit": limit, "offset": offset})
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_smiles(smiles: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Get species by SMILES. Returns paged results."""
    data = request_json("GET", f"{settings.base_url}/species/get/by-smiles/", params={"smiles": smiles, "limit": limit, "offset": offset})
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_formula(formula: str, *, phase: Optional[str] = None, descriptor: Optional[str] = None, limit: int = 50, offset: int = 0) -> Page:
    """Get species by chemical formula. Returns paged results."""
    params: Dict[str, Any] = {"formula": formula, "limit": limit, "offset": offset}
    if phase is not None: params["phase"] = phase
    if descriptor is not None: params["descriptor"] = descriptor
    data = request_json("GET", f"{settings.base_url}/species/get/by-formula/", params=params)
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_name(name: str, *, default_permissive_search: bool = False, limit: int = 50, offset: int = 0) -> Page:
    """Get species by name. Returns paged results."""
    params: Dict[str, Any] = {"name": name, "limit": limit, "offset": offset}
    if default_permissive_search:
        params["default_permissive_search"] = "1"
    data = request_json("GET", f"{settings.base_url}/species/get/by-name/", params=params)
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

# -------- species: search (paged) --------
def search_species(q: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Search species by query (name, formula, SMILES, ATcT ID, or CASRN). Returns paged results."""
    data = request_json("GET", f"{settings.base_url}/species/search/", params={"q": q, "limit": limit, "offset": offset})
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

# -------- covariance (species only in v1) --------
def get_species_covariance_matrix(ids: Optional[List[str]] = None, atctids: Optional[List[str]] = None) -> CovarianceMatrix:
    """
    Get full N×N covariance matrix for multiple species.
    
    Args:
        ids: List of internal species IDs
        atctids: List of ATcT IDs
        
    Returns:
        CovarianceMatrix with full N×N matrix
        
    Note:
        Either ids or atctids must be provided (minimum 2 species required)
    """
    if not ids and not atctids:
        raise ValueError("Either ids or atctids must be provided")
    if ids and atctids:
        raise ValueError("Provide either ids or atctids, not both")
    
    params: Dict[str, Any] = {}
    if ids:
        params["ids"] = ",".join(ids)
    if atctids:
        params["atctids"] = ",".join(atctids)
    
    data = request_json("GET", f"{settings.base_url}/covariance/species/", params=params)
    return CovarianceMatrix.from_dict(data)

def get_species_covariance_by_ids(a_id: int, b_id: int) -> Covariance2x2:
    """Legacy function for 2x2 covariance matrix by internal IDs."""
    data = request_json("GET", f"{settings.base_url}/covariance/species/", params={"ids": f"{int(a_id)},{int(b_id)}"})
    return Covariance2x2.from_dict(data)

def get_species_covariance_by_atctid(a_atctid: str, b_atctid: str) -> Covariance2x2:
    """Legacy function for 2x2 covariance matrix by ATcT IDs."""
    data = request_json("GET", f"{settings.base_url}/covariance/species/", params={"atctids": f"{a_atctid},{b_atctid}"})
    return Covariance2x2.from_dict(data)

# -------- simple API endpoints --------
def get_species_simple(name: Optional[str] = None, smiles: Optional[str] = None, 
                      inchi: Optional[str] = None, inchikey: Optional[str] = None,
                      casrn: Optional[str] = None, atctid: Optional[str] = None) -> Species:
    """
    Simple API endpoint for quick species lookup.
    
    Args:
        name: Species name or chemical formula
        smiles: SMILES notation
        inchi: InChI identifier
        inchikey: InChI Key
        casrn: CAS Registry Number
        atctid: ATcT identifier
        
    Returns:
        Species object
        
    Note:
        Only one parameter should be provided
    """
    params: Dict[str, Any] = {}
    if name:
        params["name"] = name
    elif smiles:
        params["smiles"] = smiles
    elif inchi:
        params["inchi"] = inchi
    elif inchikey:
        params["inchikey"] = inchikey
    elif casrn:
        params["casrn"] = casrn
    elif atctid:
        params["atctid"] = atctid
    else:
        raise ValueError("At least one parameter must be provided")
    
    data = request_json("GET", f"{settings.base_url}/", params=params)
    # API returns list, take first item
    if isinstance(data, list) and len(data) > 0:
        return Species.from_dict(data[0])
    elif isinstance(data, dict):
        return Species.from_dict(data)
    else:
        raise ValueError("No species found with the given parameters")

def get_all_species() -> List[Species]:
    """Get all species from the database."""
    data = request_json("GET", f"{settings.base_url}/all/")
    if isinstance(data, list):
        return [Species.from_dict(x) for x in data]
    else:
        return [Species.from_dict(x) for x in data.get("items", [])]

# -------- reaction calculations --------
def create_reaction_calculator(species_data: Dict[str, float]) -> ReactionCalculator:
    """Create a ReactionCalculator from species ATcT IDs and stoichiometry.
    
    Args:
        species_data: Dictionary mapping ATcT IDs to stoichiometric coefficients
                     e.g., {"67-56-1*0": 1.0, "7727-37-9*0": -1.0}
    
    Returns:
        ReactionCalculator instance
    """
    reaction_species = []
    
    for atct_id, stoichiometry in species_data.items():
        species = get_species(atct_id)
        reaction_species.append(ReactionSpecies(species, stoichiometry))
    
    return ReactionCalculator(reaction_species)

def calculate_reaction_enthalpy(species_data: Dict[str, float], method: str = "conventional") -> ReactionResult:
    """Calculate reaction enthalpy for a given reaction.
    
    Args:
        species_data: Dictionary mapping ATcT IDs to stoichiometric coefficients
        method: Calculation method ("covariance" or "conventional")
                Defaults to "conventional" for compatibility without numpy
    
    Returns:
        ReactionResult with enthalpy and uncertainty
    """
    calculator = create_reaction_calculator(species_data)
    
    if method == "covariance":
        return calculator.calculate_covariance_method()
    elif method == "conventional":
        return calculator.calculate_conventional_method()
    else:
        raise ValueError(f"Unknown method: {method}. Use 'covariance' or 'conventional'")

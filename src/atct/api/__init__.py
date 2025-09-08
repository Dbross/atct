from __future__ import annotations
from typing import Optional, Dict, Any, List
from ._config import settings
from ._http import request_json
from .models import Species, Page, Covariance2x2, ReactionSpecies, ReactionResult, ReactionCalculator
from .pandas_io import as_dataframe

# -------- health --------
def healthcheck() -> bool:
    try:
        data = request_json("GET", f"{settings.base_url}/health")
        return bool(data.get("ok", True))
    except Exception:
        return False

# -------- species: unique by ATcT_ID --------
def get_species(atctid: str, *, expand_xyz: bool = False) -> Species:
    params: Dict[str, Any] = {"atctid": atctid}
    if expand_xyz:
        params["expand_xyz"] = "1"
    data = request_json("GET", f"{settings.base_url}/species/get/by-atctid/", params=params)
    return Species.from_dict(data)

# -------- species: multi "get by …" (paged) --------
def get_species_by_casrn(casrn: str, *, limit: int = 50, offset: int = 0) -> Page:
    data = request_json("GET", f"{settings.base_url}/species/get/by-casrn/", params={"casrn": casrn, "limit": limit, "offset": offset})
    return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_inchi(inchi: str, *, limit: int = 50, offset: int = 0) -> Page:
    data = request_json("GET", f"{settings.base_url}/species/get/by-inchi/", params={"inchi": inchi, "limit": limit, "offset": offset})
    return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_smiles(smiles: str, *, limit: int = 50, offset: int = 0) -> Page:
    data = request_json("GET", f"{settings.base_url}/species/get/by-smiles/", params={"smiles": smiles, "limit": limit, "offset": offset})
    return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_formula(formula: str, *, phase: Optional[str] = None, descriptor: Optional[str] = None, limit: int = 50, offset: int = 0) -> Page:
    params: Dict[str, Any] = {"formula": formula, "limit": limit, "offset": offset}
    if phase is not None: params["phase"] = phase
    if descriptor is not None: params["descriptor"] = descriptor
    data = request_json("GET", f"{settings.base_url}/species/get/by-formula/", params=params)
    return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_name(name: str, *, limit: int = 50, offset: int = 0) -> Page:
    data = request_json("GET", f"{settings.base_url}/species/get/by-name/", params={"name": name, "limit": limit, "offset": offset})
    return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

# -------- species: search (paged) --------
def search_species(q: str, *, limit: int = 50, offset: int = 0) -> Page:
    data = request_json("GET", f"{settings.base_url}/species/search/", params={"q": q, "limit": limit, "offset": offset})
    return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

# -------- covariance (species only in v1) --------
def get_species_covariance_by_ids(a_id: int, b_id: int) -> Covariance2x2:
    data = request_json("GET", f"{settings.base_url}/covariance/species/", params={"a_id": int(a_id), "b_id": int(b_id)})
    return Covariance2x2.from_dict(data)

def get_species_covariance_by_atctid(a_atctid: str, b_atctid: str) -> Covariance2x2:
    data = request_json("GET", f"{settings.base_url}/covariance/species/", params={"a_atctid": a_atctid, "b_atctid": b_atctid})
    return Covariance2x2.from_dict(data)

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

def calculate_reaction_enthalpy(species_data: Dict[str, float], method: str = "sum_squares") -> ReactionResult:
    """Calculate reaction enthalpy for a given reaction.
    
    Args:
        species_data: Dictionary mapping ATcT IDs to stoichiometric coefficients
        method: Calculation method ("covariance", "sum_squares", or "conventional")
                Defaults to "sum_squares" for compatibility without numpy
    
    Returns:
        ReactionResult with enthalpy and uncertainty
    """
    calculator = create_reaction_calculator(species_data)
    
    if method == "covariance":
        return calculator.calculate_covariance_method()
    elif method == "sum_squares":
        return calculator.calculate_sum_squares_method()
    elif method == "conventional":
        return calculator.calculate_conventional_method()
    else:
        raise ValueError(f"Unknown method: {method}. Use 'covariance', 'sum_squares', or 'conventional'")

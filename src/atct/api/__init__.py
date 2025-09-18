from __future__ import annotations
import asyncio
from typing import Optional, Dict, Any, List, Union
from ._config import settings
from ._http import request_json, request_json_async
from .models import Species, Page, CovarianceMatrix, Covariance2x2, ReactionSpecies, ReactionResult, ReactionCalculator
from .pandas_io import as_dataframe

# -------- health --------
async def healthcheck_async() -> bool:
    """Async health check."""
    try:
        data = await request_json_async("GET", settings.url("health/"))
        return bool(data.get("ok", True))
    except Exception:
        return False

def healthcheck(block: bool = False) -> Union[bool, asyncio.Task]:
    """Check API health.
    
    Args:
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Boolean health status
        If block=False: asyncio.Task that resolves to boolean health status
    """
    task = healthcheck_async()
    if block:
        return asyncio.run(task)
    return task

# -------- species: unique by ATcT_ID --------
async def get_species_by_atctid_async(atctid: str, *, expand_xyz: bool = False) -> Species:
    """Async get species by ATcT ID. Returns single species."""
    params: Dict[str, Any] = {"atctid": atctid}
    if expand_xyz:
        params["expand_xyz"] = "1"
    data = await request_json_async("GET", settings.url("species/get/by-atctid/"), params=params)
    # API returns list, take first item
    if isinstance(data, list) and len(data) > 0:
        return Species.from_dict(data[0])
    elif isinstance(data, dict):
        return Species.from_dict(data)
    else:
        raise ValueError("No species found with the given ATcT ID")

def get_species_by_atctid(atctid: str, *, expand_xyz: bool = False, block: bool = False) -> Union[Species, asyncio.Task]:
    """Get species by ATcT ID. Returns single species.
    
    Args:
        atctid: ATcT ID of the species
        expand_xyz: Whether to expand XYZ coordinates
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Species object
        If block=False: asyncio.Task that resolves to Species object
    """
    task = get_species_by_atctid_async(atctid, expand_xyz=expand_xyz)
    if block:
        return asyncio.run(task)
    return task

# Legacy alias
async def get_species_async(atctid: str, *, expand_xyz: bool = False) -> Species:
    """Async legacy alias for get_species_by_atctid."""
    return await get_species_by_atctid_async(atctid, expand_xyz=expand_xyz)

def get_species(atctid: str, *, expand_xyz: bool = False, block: bool = False) -> Union[Species, asyncio.Task]:
    """Legacy alias for get_species_by_atctid.
    
    Args:
        atctid: ATcT ID of the species
        expand_xyz: Whether to expand XYZ coordinates
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Species object
        If block=False: asyncio.Task that resolves to Species object
    """
    task = get_species_async(atctid, expand_xyz=expand_xyz)
    if block:
        return asyncio.run(task)
    return task

# -------- species: multi "get by …" (paged) --------
async def get_species_by_casrn_async(casrn: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Async get species by CAS Registry Number. Returns paged results."""
    data = await request_json_async("GET", f"{settings.base_url}/species/get/by-casrn/", params={"casrn": casrn, "limit": limit, "offset": offset})
    # Handle both list and paginated response formats
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_casrn(casrn: str, *, limit: int = 50, offset: int = 0, block: bool = False) -> Union[Page, asyncio.Task]:
    """Get species by CAS Registry Number. Returns paged results.
    
    Args:
        casrn: CAS Registry Number
        limit: Maximum number of results
        offset: Number of results to skip
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Page object with species
        If block=False: asyncio.Task that resolves to Page object
    """
    task = get_species_by_casrn_async(casrn, limit=limit, offset=offset)
    if block:
        return asyncio.run(task)
    return task

async def get_species_by_inchi_async(inchi: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Async get species by InChI. Returns paged results."""
    data = await request_json_async("GET", f"{settings.base_url}/species/get/by-inchi/", params={"inchi": inchi, "limit": limit, "offset": offset})
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_inchi(inchi: str, *, limit: int = 50, offset: int = 0, block: bool = False) -> Union[Page, asyncio.Task]:
    """Get species by InChI. Returns paged results.
    
    Args:
        inchi: InChI identifier
        limit: Maximum number of results
        offset: Number of results to skip
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Page object with species
        If block=False: asyncio.Task that resolves to Page object
    """
    task = get_species_by_inchi_async(inchi, limit=limit, offset=offset)
    if block:
        return asyncio.run(task)
    return task

async def get_species_by_inchikey_async(inchikey: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Async get species by InChI Key. Returns paged results."""
    data = await request_json_async("GET", f"{settings.base_url}/species/get/by-inchikey/", params={"inchikey": inchikey, "limit": limit, "offset": offset})
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_inchikey(inchikey: str, *, limit: int = 50, offset: int = 0, block: bool = False) -> Union[Page, asyncio.Task]:
    """Get species by InChI Key. Returns paged results.
    
    Args:
        inchikey: InChI Key identifier
        limit: Maximum number of results
        offset: Number of results to skip
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Page object with species
        If block=False: asyncio.Task that resolves to Page object
    """
    task = get_species_by_inchikey_async(inchikey, limit=limit, offset=offset)
    if block:
        return asyncio.run(task)
    return task

async def get_species_by_smiles_async(smiles: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Async get species by SMILES. Returns paged results."""
    data = await request_json_async("GET", f"{settings.base_url}/species/get/by-smiles/", params={"smiles": smiles, "limit": limit, "offset": offset})
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_smiles(smiles: str, *, limit: int = 50, offset: int = 0, block: bool = False) -> Union[Page, asyncio.Task]:
    """Get species by SMILES. Returns paged results.
    
    Args:
        smiles: SMILES notation
        limit: Maximum number of results
        offset: Number of results to skip
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Page object with species
        If block=False: asyncio.Task that resolves to Page object
    """
    task = get_species_by_smiles_async(smiles, limit=limit, offset=offset)
    if block:
        return asyncio.run(task)
    return task

async def get_species_by_formula_async(formula: str, *, phase: Optional[str] = None, descriptor: Optional[str] = None, limit: int = 50, offset: int = 0) -> Page:
    """Async get species by chemical formula. Returns paged results."""
    params: Dict[str, Any] = {"formula": formula, "limit": limit, "offset": offset}
    if phase is not None: params["phase"] = phase
    if descriptor is not None: params["descriptor"] = descriptor
    data = await request_json_async("GET", f"{settings.base_url}/species/get/by-formula/", params=params)
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_formula(formula: str, *, phase: Optional[str] = None, descriptor: Optional[str] = None, limit: int = 50, offset: int = 0, block: bool = False) -> Union[Page, asyncio.Task]:
    """Get species by chemical formula. Returns paged results.
    
    Args:
        formula: Chemical formula
        phase: Phase filter (optional)
        descriptor: Descriptor filter (optional)
        limit: Maximum number of results
        offset: Number of results to skip
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Page object with species
        If block=False: asyncio.Task that resolves to Page object
    """
    task = get_species_by_formula_async(formula, phase=phase, descriptor=descriptor, limit=limit, offset=offset)
    if block:
        return asyncio.run(task)
    return task

async def get_species_by_name_async(name: str, *, default_permissive_search: bool = False, limit: int = 50, offset: int = 0) -> Page:
    """Async get species by name. Returns paged results."""
    params: Dict[str, Any] = {"name": name, "limit": limit, "offset": offset}
    if default_permissive_search:
        params["default_permissive_search"] = "1"
    data = await request_json_async("GET", f"{settings.base_url}/species/get/by-name/", params=params)
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def get_species_by_name(name: str, *, default_permissive_search: bool = False, limit: int = 50, offset: int = 0, block: bool = False) -> Union[Page, asyncio.Task]:
    """Get species by name. Returns paged results.
    
    Args:
        name: Species name
        default_permissive_search: Whether to use permissive search
        limit: Maximum number of results
        offset: Number of results to skip
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Page object with species
        If block=False: asyncio.Task that resolves to Page object
    """
    task = get_species_by_name_async(name, default_permissive_search=default_permissive_search, limit=limit, offset=offset)
    if block:
        return asyncio.run(task)
    return task

# -------- species: search (paged) --------
async def search_species_async(q: str, *, limit: int = 50, offset: int = 0) -> Page:
    """Async search species by query (name, formula, SMILES, ATcT ID, or CASRN). Returns paged results."""
    data = await request_json_async("GET", settings.url("species/search/"), params={"q": q, "limit": limit, "offset": offset})
    if isinstance(data, list):
        return Page(items=[Species.from_dict(x) for x in data], total=len(data), limit=limit, offset=offset)
    else:
        return Page(items=[Species.from_dict(x) for x in data.get("items", [])], total=int(data.get("total", 0)), limit=limit, offset=offset)

def search_species(q: str, *, limit: int = 50, offset: int = 0, block: bool = False) -> Union[Page, asyncio.Task]:
    """Search species by query (name, formula, SMILES, ATcT ID, or CASRN). Returns paged results.
    
    Args:
        q: Search query
        limit: Maximum number of results
        offset: Number of results to skip
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Page object with species
        If block=False: asyncio.Task that resolves to Page object
    """
    task = search_species_async(q, limit=limit, offset=offset)
    if block:
        return asyncio.run(task)
    return task

# -------- covariance (species only in v1) --------
async def get_species_covariance_matrix_async(ids: Optional[List[str]] = None, atctids: Optional[List[str]] = None) -> CovarianceMatrix:
    """
    Async get full N×N covariance matrix for multiple species.
    
    Args:
        ids: List of internal species IDs
        atctids: List of ATcT IDs
        
    Returns:
        CovarianceMatrix with full N×N matrix
        
    Note:
        Either ids or atctids must be provided (minimum 2 species required)
        For large requests (>25 species), automatically uses batch processing
    """
    if not ids and not atctids:
        raise ValueError("Either ids or atctids must be provided")
    if ids and atctids:
        raise ValueError("Provide either ids or atctids, not both")
    
    # Determine the list of identifiers and count
    identifier_list = ids if ids else atctids
    count = len(identifier_list)
    
    # For large requests, use manual batching with regular endpoint
    if count > 25:
        if atctids:
            # Process in batches of 25 using the regular endpoint
            batch_size = 25
            all_covariance_data = []
            
            for i in range(0, count, batch_size):
                batch = atctids[i:i + batch_size]
                batch_num = i // batch_size + 1
                
                try:
                    batch_data = await request_json_async("GET", settings.url("covariance/species/"), 
                                                        params={"atctids": ",".join(batch)})
                    if batch_data:
                        all_covariance_data.append(batch_data)
                        print(f"✓ Batch {batch_num}: {len(batch)} species")
                except Exception as e:
                    # If batch fails, try smaller sub-batches
                    print(f"✗ Batch {batch_num} failed: {e}")
                    print(f"  Trying smaller sub-batches for batch {batch_num}...")
                    
                    # Try sub-batches of 10 species each
                    sub_batch_size = 10
                    for j in range(0, len(batch), sub_batch_size):
                        sub_batch = batch[j:j + sub_batch_size]
                        sub_batch_num = j // sub_batch_size + 1
                        
                        try:
                            sub_batch_data = await request_json_async("GET", settings.url("covariance/species/"), 
                                                                    params={"atctids": ",".join(sub_batch)})
                            if sub_batch_data:
                                all_covariance_data.append(sub_batch_data)
                                print(f"  ✓ Sub-batch {batch_num}.{sub_batch_num}: {len(sub_batch)} species")
                        except Exception as sub_e:
                            print(f"  ✗ Sub-batch {batch_num}.{sub_batch_num} failed: {sub_e}")
                            # Try individual species as last resort
                            for k, single_atctid in enumerate(sub_batch):
                                try:
                                    single_data = await request_json_async("GET", settings.url("covariance/species/"), 
                                                                          params={"atctids": single_atctid})
                                    if single_data:
                                        all_covariance_data.append(single_data)
                                        print(f"    ✓ Individual {batch_num}.{sub_batch_num}.{k+1}: {single_atctid}")
                                except Exception as single_e:
                                    print(f"    ✗ Individual {batch_num}.{sub_batch_num}.{k+1} failed: {single_atctid}")
                                    continue
            
            if not all_covariance_data:
                raise ValueError("No covariance data retrieved from any batch")
            
            # Combine all successful batches into a single covariance matrix
            print(f"Combining covariance data from {len(all_covariance_data)} successful batches...")
            
            # Collect all species and their data from all batches
            all_species = []
            all_species_atctids = []
            species_to_batch = {}  # Maps ATcT ID to batch index
            
            for batch_idx, batch_data in enumerate(all_covariance_data):
                if 'species' in batch_data and 'species_atctids' in batch_data:
                    for i, atctid in enumerate(batch_data['species_atctids']):
                        if atctid not in species_to_batch:
                            all_species.append(batch_data['species'][i])
                            all_species_atctids.append(atctid)
                            species_to_batch[atctid] = (batch_idx, i)
            
            # Create combined matrix
            n_total = len(all_species_atctids)
            combined_matrix = [[0.0 for _ in range(n_total)] for _ in range(n_total)]
            
            # Fill in the diagonal blocks from each batch
            for batch_idx, batch_data in enumerate(all_covariance_data):
                if 'matrix' in batch_data and 'species_atctids' in batch_data:
                    batch_matrix = batch_data['matrix']
                    batch_atctids = batch_data['species_atctids']
                    
                    for i, atctid_i in enumerate(batch_atctids):
                        if atctid_i in species_to_batch:
                            global_i = all_species_atctids.index(atctid_i)
                            for j, atctid_j in enumerate(batch_atctids):
                                if atctid_j in species_to_batch:
                                    global_j = all_species_atctids.index(atctid_j)
                                    combined_matrix[global_i][global_j] = batch_matrix[i][j]
            
            # Create combined CovarianceMatrix object
            combined_data = {
                'species': all_species,
                'species_atctids': all_species_atctids,
                'units': all_covariance_data[0].get('units', 'kJ/mol'),
                'matrix': combined_matrix
            }
            
            print(f"Combined matrix size: {n_total}x{n_total}")
            return CovarianceMatrix.from_dict(combined_data)
        else:
            # For large ID requests, we need to convert IDs to ATcT IDs first
            # This is a limitation - batch endpoint only supports ATcT IDs
            raise ValueError("Large requests with internal IDs not supported. Use ATcT IDs instead.")
    
    # For smaller requests, use the regular endpoint
    params: Dict[str, Any] = {}
    if ids:
        params["ids"] = ",".join(ids)
    if atctids:
        params["atctids"] = ",".join(atctids)
    
    data = await request_json_async("GET", settings.url("covariance/species/"), params=params)
    return CovarianceMatrix.from_dict(data)

def get_species_covariance_matrix(ids: Optional[List[str]] = None, atctids: Optional[List[str]] = None, block: bool = False) -> Union[CovarianceMatrix, asyncio.Task]:
    """
    Get full N×N covariance matrix for multiple species.
    
    Args:
        ids: List of internal species IDs
        atctids: List of ATcT IDs
        block: If True, blocks and returns result. If False, returns async task.
        
    Returns:
        If block=True: CovarianceMatrix with full N×N matrix
        If block=False: asyncio.Task that resolves to CovarianceMatrix
        
    Note:
        Either ids or atctids must be provided (minimum 2 species required)
    """
    task = get_species_covariance_matrix_async(ids=ids, atctids=atctids)
    if block:
        return asyncio.run(task)
    return task

async def get_species_covariance_by_ids_async(a_id: int, b_id: int) -> Covariance2x2:
    """Async legacy function for 2x2 covariance matrix by internal IDs."""
    data = await request_json_async("GET", f"{settings.base_url}/covariance/species/", params={"ids": f"{int(a_id)},{int(b_id)}"})
    return Covariance2x2.from_dict(data)

def get_species_covariance_by_ids(a_id: int, b_id: int, block: bool = False) -> Union[Covariance2x2, asyncio.Task]:
    """Legacy function for 2x2 covariance matrix by internal IDs.
    
    Args:
        a_id: First species internal ID
        b_id: Second species internal ID
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Covariance2x2 object
        If block=False: asyncio.Task that resolves to Covariance2x2 object
    """
    task = get_species_covariance_by_ids_async(a_id, b_id)
    if block:
        return asyncio.run(task)
    return task

async def get_species_covariance_by_atctid_async(a_atctid: str, b_atctid: str) -> Covariance2x2:
    """Async legacy function for 2x2 covariance matrix by ATcT IDs."""
    data = await request_json_async("GET", settings.url("covariance/species/"), params={"atctids": f"{a_atctid},{b_atctid}"})
    return Covariance2x2.from_dict(data)

def get_species_covariance_by_atctid(a_atctid: str, b_atctid: str, block: bool = False) -> Union[Covariance2x2, asyncio.Task]:
    """Legacy function for 2x2 covariance matrix by ATcT IDs.
    
    Args:
        a_atctid: First species ATcT ID
        b_atctid: Second species ATcT ID
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: Covariance2x2 object
        If block=False: asyncio.Task that resolves to Covariance2x2 object
    """
    task = get_species_covariance_by_atctid_async(a_atctid, b_atctid)
    if block:
        return asyncio.run(task)
    return task

# -------- simple API endpoints --------
async def get_species_simple_async(name: Optional[str] = None, smiles: Optional[str] = None, 
                      inchi: Optional[str] = None, inchikey: Optional[str] = None,
                      casrn: Optional[str] = None, atctid: Optional[str] = None) -> Species:
    """
    Async simple API endpoint for quick species lookup.
    
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
    
    data = await request_json_async("GET", f"{settings.base_url}/", params=params)
    # API returns list, take first item
    if isinstance(data, list) and len(data) > 0:
        return Species.from_dict(data[0])
    elif isinstance(data, dict):
        return Species.from_dict(data)
    else:
        raise ValueError("No species found with the given parameters")

def get_species_simple(name: Optional[str] = None, smiles: Optional[str] = None, 
                      inchi: Optional[str] = None, inchikey: Optional[str] = None,
                      casrn: Optional[str] = None, atctid: Optional[str] = None, block: bool = False) -> Union[Species, asyncio.Task]:
    """
    Simple API endpoint for quick species lookup.
    
    Args:
        name: Species name or chemical formula
        smiles: SMILES notation
        inchi: InChI identifier
        inchikey: InChI Key
        casrn: CAS Registry Number
        atctid: ATcT identifier
        block: If True, blocks and returns result. If False, returns async task.
        
    Returns:
        If block=True: Species object
        If block=False: asyncio.Task that resolves to Species object
        
    Note:
        Only one parameter should be provided
    """
    task = get_species_simple_async(name=name, smiles=smiles, inchi=inchi, inchikey=inchikey, casrn=casrn, atctid=atctid)
    if block:
        return asyncio.run(task)
    return task

async def get_all_species_async() -> List[Species]:
    """Async get all species from the database."""
    data = await request_json_async("GET", f"{settings.base_url}/all/")
    if isinstance(data, list):
        return [Species.from_dict(x) for x in data]
    else:
        return [Species.from_dict(x) for x in data.get("items", [])]

def get_all_species(block: bool = False) -> Union[List[Species], asyncio.Task]:
    """Get all species from the database.
    
    Args:
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: List of Species objects
        If block=False: asyncio.Task that resolves to List of Species objects
    """
    task = get_all_species_async()
    if block:
        return asyncio.run(task)
    return task

# -------- reaction calculations --------
async def create_reaction_calculator_async(species_data: Dict[str, float]) -> ReactionCalculator:
    """Async create a ReactionCalculator from species ATcT IDs and stoichiometry.
    
    Args:
        species_data: Dictionary mapping ATcT IDs to stoichiometric coefficients
                     e.g., {"67-56-1*0": 1.0, "7727-37-9*0": -1.0}
    
    Returns:
        ReactionCalculator instance
    """
    reaction_species = []
    
    # Collect all species requests
    species_tasks = []
    for atct_id in species_data.keys():
        species_tasks.append(get_species_async(atct_id))
    
    # Wait for all species to be fetched
    species_list = await asyncio.gather(*species_tasks)
    
    # Create reaction species
    for species, (atct_id, stoichiometry) in zip(species_list, species_data.items()):
        reaction_species.append(ReactionSpecies(species, stoichiometry))
    
    # Create calculator without building covariance matrix initially
    calculator = ReactionCalculator(reaction_species, cov_matrix=None, build_cov_matrix=False)
    try:
        import numpy as np
        # If numpy is available, build async covariance matrix
        calculator.cov_matrix = await calculator._build_covariance_matrix_async()
    except ImportError:
        pass  # numpy not available, leave cov_matrix as None
    
    return calculator

def create_reaction_calculator(species_data: Dict[str, float], block: bool = False) -> Union[ReactionCalculator, asyncio.Task]:
    """Create a ReactionCalculator from species ATcT IDs and stoichiometry.
    
    Args:
        species_data: Dictionary mapping ATcT IDs to stoichiometric coefficients
                     e.g., {"67-56-1*0": 1.0, "7727-37-9*0": -1.0}
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: ReactionCalculator instance
        If block=False: asyncio.Task that resolves to ReactionCalculator instance
    """
    task = create_reaction_calculator_async(species_data)
    if block:
        return asyncio.run(task)
    return task

async def calculate_reaction_enthalpy_async(species_data: Dict[str, float], method: str = "conventional") -> ReactionResult:
    """Async calculate reaction enthalpy for a given reaction.
    
    Args:
        species_data: Dictionary mapping ATcT IDs to stoichiometric coefficients
        method: Calculation method ("covariance" or "conventional")
                Defaults to "conventional" for compatibility without numpy
    
    Returns:
        ReactionResult with enthalpy and uncertainty
    """
    calculator = await create_reaction_calculator_async(species_data)
    
    if method == "covariance":
        return calculator.calculate_covariance_method()
    elif method == "conventional":
        return calculator.calculate_conventional_method()
    else:
        raise ValueError(f"Unknown method: {method}. Use 'covariance' or 'conventional'")

def calculate_reaction_enthalpy(species_data: Dict[str, float], method: str = "conventional", block: bool = False) -> Union[ReactionResult, asyncio.Task]:
    """Calculate reaction enthalpy for a given reaction.
    
    Args:
        species_data: Dictionary mapping ATcT IDs to stoichiometric coefficients
        method: Calculation method ("covariance" or "conventional")
                Defaults to "conventional" for compatibility without numpy
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: ReactionResult with enthalpy and uncertainty
        If block=False: asyncio.Task that resolves to ReactionResult
    """
    task = calculate_reaction_enthalpy_async(species_data, method=method)
    if block:
        return asyncio.run(task)
    return task

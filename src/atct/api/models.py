from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional, Union
from math import sqrt

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False
    np = None

# ---------- Shared pagination ----------
@dataclass
class Page:
    items: List[Any]
    total: int
    limit: int
    offset: int
    def to_dict(self) -> Dict[str, Any]:
        return {
            "items": [getattr(x, "to_dict", lambda: x)() for x in self.items],
            "total": self.total,
            "limit": self.limit,
            "offset": self.offset,
        }

# ---------- Species (original API fields only) ----------
# NOTE: XYZ may be a filename (str), expanded list[str], or None
XYZType = Optional[Union[str, List[str]]]

@dataclass
class Species:
    atct_tn_version: Optional[str]
    atct_id: str
    name: Optional[str]
    formula: Optional[str]
    delta_h_0k: Optional[str]
    delta_h_298k: Optional[str]
    delta_h_298k_uncertainty: Optional[str]
    smiles: Optional[str]
    casrn: Optional[str]
    inchi: Optional[str]
    inchi_key: Optional[str]
    charge: Optional[int]
    xyz: XYZType

    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "Species":
        return Species(
            atct_tn_version=d.get("ATcT_TN_Version"),
            atct_id=d.get("ATcT_ID") or "",
            name=d.get("Name"),
            formula=d.get("Formula"),
            delta_h_0k=d.get("∆fH_0K"),
            delta_h_298k=d.get("∆fH_298K"),
            delta_h_298k_uncertainty=d.get("∆fH_298K_uncertainty"),
            smiles=d.get("SMILES"),
            casrn=d.get("CASRN"),
            inchi=d.get("InChI"),
            inchi_key=d.get("InChI_Key"),
            charge=(int(d["charge"]) if d.get("charge") is not None else None),
            xyz=d.get("XYZ"),
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ATcT_TN_Version": self.atct_tn_version,
            "ATcT_ID": self.atct_id,
            "Name": self.name,
            "Formula": self.formula,
            "∆fH_0K": self.delta_h_0k,
            "∆fH_298K": self.delta_h_298k,
            "∆fH_298K_uncertainty": self.delta_h_298k_uncertainty,
            "SMILES": self.smiles,
            "CASRN": self.casrn,
            "InChI": self.inchi,
            "InChI_Key": self.inchi_key,
            "charge": self.charge,
            "XYZ": self.xyz,
        }


# ---------- Covariance ----------
@dataclass
class CovarianceMatrix:
    species_ids: List[str]
    species_atctids: List[str]
    matrix: List[List[float]]
    units: str = "kJ/mol"

    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "CovarianceMatrix":
        # Extract species info from the species array if available
        species_ids = []
        species_atctids = []
        
        if "species" in d and isinstance(d["species"], list):
            # API returns species in an array with id and atctid fields
            species_ids = [str(species.get("id", "")) for species in d["species"]]
            species_atctids = [species.get("atctid", "") for species in d["species"]]
        else:
            # Fallback to old format
            species_ids = list(d.get("species_ids", []))
            species_atctids = list(d.get("species_atctids", []))
        
        return CovarianceMatrix(
            species_ids=species_ids,
            species_atctids=species_atctids,
            matrix=[[float(x) for x in row] for row in d.get("matrix", [])],
            units=str(d.get("units", "kJ/mol")),
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "species_ids": self.species_ids,
            "species_atctids": self.species_atctids,
            "matrix": self.matrix,
            "units": self.units
        }

# Legacy support for 2x2 covariance
@dataclass
class Covariance2x2:
    labels: List[str]
    units: str
    matrix: List[List[float]]

    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "Covariance2x2":
        return Covariance2x2(
            labels=list(d.get("labels") or ["ΔH(A)", "ΔH(B)"]),
            units=str(d.get("units") or ""),
            matrix=[[float(x) for x in row] for row in d.get("matrix", [])],
        )

    def to_dict(self) -> Dict[str, Any]:
        return {"labels": self.labels, "units": self.units, "matrix": self.matrix}


# ---------- Reaction Calculation ----------
@dataclass
class ReactionSpecies:
    """A species in a reaction with stoichiometric coefficient."""
    species: Species
    stoichiometry: float
    
    def __post_init__(self):
        if not isinstance(self.species, Species):
            raise ValueError("species must be a Species object")
        if not isinstance(self.stoichiometry, (int, float)):
            raise ValueError("stoichiometry must be a number")


@dataclass
class ReactionResult:
    """Result of a reaction enthalpy calculation."""
    delta_h: float
    uncertainty: float
    method: str
    units: str = "kJ/mol"
    
    def __str__(self) -> str:
        return f"{self.delta_h:.2f} ± {self.uncertainty:.6f} {self.units} ({self.method})"


class ReactionCalculator:
    """Calculate reaction enthalpies with proper uncertainty propagation."""
    
    def __init__(self, reaction_species: List[ReactionSpecies], cov_matrix: Optional[np.ndarray] = None, build_cov_matrix: bool = True):
        """Initialize with list of ReactionSpecies objects."""
        self.reaction_species = reaction_species
        
        if NUMPY_AVAILABLE:
            self.stoichiometry = np.array([rs.stoichiometry for rs in reaction_species])
            self.enthalpies = np.array([rs.species.delta_h_298k for rs in reaction_species])
            self.uncertainties = np.array([rs.species.delta_h_298k_uncertainty for rs in reaction_species])
            
            # Convert string values to floats, handling 'exact' uncertainties
            self.enthalpies = np.array([float(h) if h is not None else 0.0 for h in self.enthalpies])
            self.uncertainties = np.array([
                0.0 if unc == 'exact' or unc is None else float(unc) 
                for unc in self.uncertainties
            ])
            
            # Use provided covariance matrix or build it (if requested)
            if cov_matrix is not None:
                self.cov_matrix = cov_matrix
            elif build_cov_matrix:
                self.cov_matrix = self._build_covariance_matrix()
            else:
                self.cov_matrix = None
        else:
            # Fallback to lists when numpy is not available
            self.stoichiometry = [rs.stoichiometry for rs in reaction_species]
            self.enthalpies = [float(rs.species.delta_h_298k) if rs.species.delta_h_298k is not None else 0.0 for rs in reaction_species]
            self.uncertainties = [
                0.0 if rs.species.delta_h_298k_uncertainty == 'exact' or rs.species.delta_h_298k_uncertainty is None 
                else float(rs.species.delta_h_298k_uncertainty) 
                for rs in reaction_species
            ]
            self.cov_matrix = None
    
    def _build_covariance_matrix(self):
        """Build covariance matrix from species data."""
        if not NUMPY_AVAILABLE:
            return None
            
        n = len(self.reaction_species)
        
        # Get all species ATcT IDs
        atct_ids = [rs.species.atct_id for rs in self.reaction_species]
        
        # Import here to avoid circular imports
        from . import get_species_covariance_matrix
        import asyncio
        
        try:
            # Check if we're in an async context
            try:
                loop = asyncio.get_running_loop()
                # We're in an async context, can't use asyncio.run()
                # Return diagonal matrix to avoid the warning
                matrix = np.zeros((n, n))
                for i in range(n):
                    matrix[i, i] = self.uncertainties[i] ** 2
                return matrix
            except RuntimeError:
                # No running loop, safe to use asyncio.run()
                cov_data = get_species_covariance_matrix(atctids=atct_ids, block=True)
                
                # The API returns the matrix in the order of cov_data.species_atctids
                # We need to reorder it to match our atct_ids order
                api_atct_ids = cov_data.species_atctids
                if not api_atct_ids:
                    raise Exception("API returned empty species_atctids")
                
                # Create mapping from ATcT ID to index in our order
                id_to_index = {atct_id: i for i, atct_id in enumerate(atct_ids)}
                
                # Reorder the matrix to match our species order
                reordered_matrix = np.zeros((n, n))
                for i, api_atct_id in enumerate(api_atct_ids):
                    if api_atct_id in id_to_index:
                        our_index = id_to_index[api_atct_id]
                        for j, api_atct_id2 in enumerate(api_atct_ids):
                            if api_atct_id2 in id_to_index:
                                our_index2 = id_to_index[api_atct_id2]
                                reordered_matrix[our_index, our_index2] = cov_data.matrix[i][j]
                
                return reordered_matrix
                
        except Exception:
            # If we can't get the full covariance matrix, fall back to diagonal only
            matrix = np.zeros((n, n))
            for i in range(n):
                matrix[i, i] = self.uncertainties[i] ** 2
            return matrix
    
    async def _build_covariance_matrix_async(self):
        """Build covariance matrix from species data asynchronously."""
        if not NUMPY_AVAILABLE:
            return None
            
        n = len(self.reaction_species)
        
        # Get all species ATcT IDs
        atct_ids = [rs.species.atct_id for rs in self.reaction_species]
        
        # Import here to avoid circular imports
        from . import get_species_covariance_matrix_async
        
        try:
            # Get the full covariance matrix for all species at once
            cov_data = await get_species_covariance_matrix_async(atctids=atct_ids)
            
            # The API returns the matrix in the order of cov_data.species_atctids
            # We need to reorder it to match our atct_ids order
            api_atct_ids = cov_data.species_atctids
            if not api_atct_ids:
                raise Exception("API returned empty species_atctids")
            
            # Create mapping from ATcT ID to index in our order
            id_to_index = {atct_id: i for i, atct_id in enumerate(atct_ids)}
            
            # Reorder the matrix to match our species order
            reordered_matrix = np.zeros((n, n))
            for i, api_atct_id in enumerate(api_atct_ids):
                if api_atct_id in id_to_index:
                    our_index = id_to_index[api_atct_id]
                    for j, api_atct_id2 in enumerate(api_atct_ids):
                        if api_atct_id2 in id_to_index:
                            our_index2 = id_to_index[api_atct_id2]
                            reordered_matrix[our_index, our_index2] = cov_data.matrix[i][j]
            
            return reordered_matrix
            
        except Exception:
            # If we can't get the full covariance matrix, fall back to diagonal only
            matrix = np.zeros((n, n))
            for i in range(n):
                matrix[i, i] = self.uncertainties[i] ** 2
            return matrix
    
    def calculate_covariance_method(self) -> ReactionResult:
        """Calculate reaction enthalpy using full covariance matrix."""
        if not NUMPY_AVAILABLE:
            raise ImportError("Install with `atct[numpy]` for covariance method support.")
        
        variance = np.dot(self.stoichiometry, np.dot(self.cov_matrix, self.stoichiometry))
        uncertainty = sqrt(variance)
        delta_h = np.dot(self.stoichiometry, self.enthalpies)
        return ReactionResult(delta_h, uncertainty, "covariance")
    
    
    def calculate_conventional_method(self) -> ReactionResult:
        """Calculate reaction enthalpy using conventional uncertainty propagation."""
        if NUMPY_AVAILABLE:
            uncertainty = sqrt(sum((self.uncertainties * self.stoichiometry) ** 2))
            delta_h = np.dot(self.stoichiometry, self.enthalpies)
        else:
            # Fallback calculation without numpy
            uncertainty = sqrt(sum((u * s) ** 2 for u, s in zip(self.uncertainties, self.stoichiometry)))
            delta_h = sum(h * s for h, s in zip(self.enthalpies, self.stoichiometry))
        return ReactionResult(delta_h, uncertainty, "conventional")
    
    def compare_methods(self) -> Dict[str, Any]:
        """Compare calculation methods."""
        results = {}
        
        # Always available method
        conv_result = self.calculate_conventional_method()
        results['conventional'] = conv_result
        
        # Covariance method only if numpy is available
        if NUMPY_AVAILABLE:
            cov_result = self.calculate_covariance_method()
            results['covariance'] = cov_result
            
            diff_unc = abs(cov_result.uncertainty - conv_result.uncertainty)
            ref_unc = max(cov_result.uncertainty, conv_result.uncertainty)
            significance_ratio = diff_unc / ref_unc if ref_unc > 0 else 0
            results['difference'] = diff_unc
            results['significance'] = significance_ratio
        else:
            results['covariance'] = None
            results['difference'] = 0.0
            results['significance'] = 0.0
        
        return results
    
    def get_uncertainties_from_diagonal(self):
        """Get uncertainties from covariance matrix diagonal."""
        if not NUMPY_AVAILABLE:
            raise ImportError("Install with `atct[numpy]` for covariance matrix support.")
        return np.sqrt(np.diag(self.cov_matrix))

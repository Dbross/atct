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
        return CovarianceMatrix(
            species_ids=list(d.get("species_ids", [])),
            species_atctids=list(d.get("species_atctids", [])),
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
    
    def __init__(self, reaction_species: List[ReactionSpecies]):
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
            
            # Build covariance matrix
            self.cov_matrix = self._build_covariance_matrix()
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
        matrix = np.zeros((n, n))
        
        # Fill diagonal elements with variances
        for i in range(n):
            matrix[i, i] = self.uncertainties[i] ** 2
        
        # Fill off-diagonal elements with actual covariances
        # Import here to avoid circular imports
        from . import get_species_covariance_by_atctid
        
        for i in range(n):
            for j in range(i + 1, n):
                try:
                    # Get covariance between species i and j
                    cov_data = get_species_covariance_by_atctid(
                        self.reaction_species[i].species.atct_id,
                        self.reaction_species[j].species.atct_id,
                        block=True
                    )
                    # The covariance is the off-diagonal element of the 2x2 matrix
                    covariance = cov_data.matrix[0][1]  # or [1][0], they should be the same
                    matrix[i, j] = covariance
                    matrix[j, i] = covariance  # Symmetric matrix
                    matrix[i, i] = cov_data.matrix[0][0]
                    matrix[j, j] = cov_data.matrix[1][1] 
                except Exception:
                    # If covariance data is not available, assume zero correlation
                    matrix[i, j] = 0.0
                    matrix[j, i] = 0.0
        
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

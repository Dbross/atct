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
    unit: Optional[str]
    mass: Optional[str]
    mass_uncertainty: Optional[str]
    smiles: Optional[str]
    casrn: Optional[str]
    charge: Optional[int]
    xyz: XYZType

    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "Species":
        return Species(
            atct_tn_version=d.get("ATcT_TN_Version"),
            atct_id=d.get("ATcT_ID") or "",
            name=d.get("Name"),
            formula=d.get("Formula"),
            delta_h_0k=d.get("Delta_Hf_0K"),
            delta_h_298k=d.get("Delta_Hf_298K"),
            delta_h_298k_uncertainty=d.get("Delta_Hf298K_uncertainty"),
            unit=d.get("unit"),
            mass=d.get("mass"),
            mass_uncertainty=d.get("mass_uncertainty"),
            smiles=d.get("SMILES"),
            casrn=d.get("CASRN"),
            charge=(int(d["charge"]) if d.get("charge") is not None else None),
            xyz=d.get("XYZ"),
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ATcT_TN_Version": self.atct_tn_version,
            "ATcT_ID": self.atct_id,
            "Name": self.name,
            "Formula": self.formula,
            "Delta_Hf_0K": self.delta_h_0k,
            "Delta_Hf_298K": self.delta_h_298k,
            "Delta_Hf298K_uncertainty": self.delta_h_298k_uncertainty,
            "unit": self.unit,
            "mass": self.mass,
            "mass_uncertainty": self.mass_uncertainty,
            "SMILES": self.smiles,
            "CASRN": self.casrn,
            "charge": self.charge,
            "XYZ": self.xyz,
        }


# ---------- Covariance ----------
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
        
        # For now, we'll use a simplified approach
        # In a full implementation, you'd fetch actual covariance data
        for i in range(n):
            for j in range(n):
                if i == j:
                    # Diagonal elements are variances (uncertainty squared)
                    matrix[i, j] = self.uncertainties[i] ** 2
                else:
                    # Off-diagonal elements would be actual covariances
                    # For now, assume zero correlation
                    matrix[i, j] = 0.0
        
        return matrix
    
    def calculate_covariance_method(self) -> ReactionResult:
        """Calculate reaction enthalpy using full covariance matrix."""
        if not NUMPY_AVAILABLE:
            raise ImportError("Install with `atct[numpy]` for covariance method support.")
        
        variance = np.dot(self.stoichiometry, np.dot(self.cov_matrix, self.stoichiometry))
        uncertainty = sqrt(variance)
        delta_h = np.dot(self.stoichiometry, self.enthalpies)
        return ReactionResult(delta_h, uncertainty, "covariance")
    
    def calculate_sum_squares_method(self) -> ReactionResult:
        """Calculate reaction enthalpy using sum of squares method."""
        if NUMPY_AVAILABLE:
            uncertainty = sqrt(sum((self.uncertainties * self.stoichiometry) ** 2))
            delta_h = np.dot(self.stoichiometry, self.enthalpies)
        else:
            # Fallback calculation without numpy
            uncertainty = sqrt(sum((u * s) ** 2 for u, s in zip(self.uncertainties, self.stoichiometry)))
            delta_h = sum(h * s for h, s in zip(self.enthalpies, self.stoichiometry))
        return ReactionResult(delta_h, uncertainty, "sum_squares")
    
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
        """Compare all three calculation methods."""
        results = {}
        
        # Always available methods
        ss_result = self.calculate_sum_squares_method()
        conv_result = self.calculate_conventional_method()
        results['sum_squares'] = ss_result
        results['conventional'] = conv_result
        
        # Covariance method only if numpy is available
        if NUMPY_AVAILABLE:
            cov_result = self.calculate_covariance_method()
            results['covariance'] = cov_result
            
            diff_unc = abs(cov_result.uncertainty - ss_result.uncertainty)
            ref_unc = max(cov_result.uncertainty, ss_result.uncertainty)
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

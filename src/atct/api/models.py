from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, List, Union, Dict, Any

@dataclass
class Uncertainty:
    value: float
    unit: str
    method: Optional[str] = None

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
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> Species:
        """Create Species from API response dict."""
        # Handle both old and new API field names
        atct_id = data.get("ATcT_ID") or data.get("atct_id", "")
        name = data.get("Name") or data.get("name", "")
        formula = data.get("Formula") or data.get("formula", "")
        phase = data.get("phase")
        
        # Handle enthalpy data
        delta_h_298k = data.get("Delta_Hf_298K") or data.get("delta_h_298k")
        if delta_h_298k is not None:
            delta_h_298k = float(delta_h_298k)
        
        delta_h_298k_unc = None
        if delta_h_298k is not None:
            unc_value = data.get("Delta_Hf298K_uncertainty") or data.get("delta_h_298k_unc")
            if unc_value is not None and unc_value != "exact":
                delta_h_298k_unc = Uncertainty(
                    value=float(unc_value),
                    unit="kJ/mol",
                    method="ATcT"
                )
        
        smiles = data.get("SMILES") or data.get("smiles")
        mass = data.get("mass")
        if mass is not None:
            mass = float(mass)
        
        return cls(
            atct_id=atct_id,
            name=name,
            formula=formula,
            phase=phase,
            delta_h_298k=delta_h_298k,
            delta_h_298k_unc=delta_h_298k_unc,
            smiles=smiles,
            mass=mass
        )

@dataclass
class CovarianceMatrix:
    """2x2 covariance matrix for two species."""
    atct_id1: str
    atct_id2: str
    matrix: List[List[float]]  # 2x2 matrix
    version: Optional[str] = None
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> CovarianceMatrix:
        """Create CovarianceMatrix from API response dict."""
        atct_id1 = data.get("atctid1", "")
        atct_id2 = data.get("atctid2", "")
        version = data.get("version") or data.get("ATcT_TN_Version")
        
        # Handle v1 API format with 2x2 matrix
        if "matrix" in data:
            matrix = data["matrix"]
        else:
            # Fallback to old format
            covariance = float(data.get("covariance", 0.0))
            matrix = [[0.0, covariance], [covariance, 0.0]]
        
        return cls(
            atct_id1=atct_id1,
            atct_id2=atct_id2,
            matrix=matrix,
            version=version
        )

@dataclass
class Page:
    items: List[Species]
    total: int
    limit: int
    offset: int

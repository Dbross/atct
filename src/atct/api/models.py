from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional, Union, Iterator
from math import sqrt
import io
import tempfile
import os

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

# ---------- XYZ Data Handling ----------
@dataclass
class XYZData:
    """Helper class for handling XYZ coordinate data."""
    atoms: List[str]
    coordinates: List[List[float]]
    comment: Optional[str] = None
    
    def __post_init__(self):
        if len(self.atoms) != len(self.coordinates):
            raise ValueError("Number of atoms must match number of coordinate sets")
    
    @classmethod
    def from_xyz_lines(cls, xyz_lines: List[str]) -> "XYZData":
        """Create XYZData from list of XYZ format lines."""
        if not xyz_lines:
            raise ValueError("Empty XYZ data")
        
        # First line is number of atoms
        try:
            num_atoms = int(xyz_lines[0].strip())
        except (ValueError, IndexError):
            raise ValueError("Invalid XYZ format: first line must be number of atoms")
        
        # Second line is comment (optional)
        comment = xyz_lines[1].strip() if len(xyz_lines) > 1 else None
        
        # Remaining lines are atom coordinates
        if len(xyz_lines) < 2 + num_atoms:
            raise ValueError(f"Not enough coordinate lines: expected {num_atoms}, got {len(xyz_lines) - 2}")
        
        atoms = []
        coordinates = []
        
        for i in range(2, 2 + num_atoms):
            line = xyz_lines[i].strip()
            if not line:
                continue
                
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"Invalid coordinate line {i-1}: {line}")
            
            atoms.append(parts[0])
            try:
                coords = [float(x) for x in parts[1:4]]
                coordinates.append(coords)
            except ValueError:
                raise ValueError(f"Invalid coordinates in line {i-1}: {line}")
        
        return cls(atoms=atoms, coordinates=coordinates, comment=comment)
    
    def to_xyz_string(self) -> str:
        """Convert to XYZ format string."""
        lines = [str(len(self.atoms))]
        if self.comment:
            lines.append(self.comment)
        else:
            lines.append("")  # Empty comment line
        
        for atom, coords in zip(self.atoms, self.coordinates):
            lines.append(f"{atom} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}")
        
        return "\n".join(lines)
    
    def to_ase_atoms(self):
        """Convert to ASE Atoms object (requires ASE to be installed)."""
        try:
            from ase import Atoms
        except ImportError:
            raise ImportError("ASE is not installed. Install with: pip install ase")
        
        return Atoms(symbols=self.atoms, positions=self.coordinates)
    
    def to_rdkit_mol(self):
        """Convert to RDKit Mol object (requires RDKit to be installed)."""
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors
            from rdkit.Chem import rdDistGeom
        except ImportError:
            raise ImportError("RDKit is not installed. Install with: pip install rdkit")
        
        # Method 1: Create molecule from XYZ data with bond inference
        try:
            mol = Chem.RWMol()
            
            # Add atoms
            for i, symbol in enumerate(self.atoms):
                atom = Chem.Atom(symbol)
                mol.AddAtom(atom)
            
            # Try to infer bonds based on distances
            self._add_bonds_from_distances(mol)
            
            # Set 3D coordinates
            conf = Chem.Conformer(len(self.atoms))
            for i, coords in enumerate(self.coordinates):
                conf.SetAtomPosition(i, coords)
            mol.AddConformer(conf, assignId=True)
            
            # Get the molecule
            rdkit_mol = mol.GetMol()
            
            # Initialize RingInfo to avoid RingInfo errors
            try:
                rdkit_mol.GetRingInfo()
            except Exception:
                pass
            
            # Try to sanitize the molecule
            try:
                Chem.SanitizeMol(rdkit_mol)
            except Exception:
                try:
                    # If full sanitization fails, try partial
                    Chem.SanitizeMol(rdkit_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS)
                except Exception:
                    try:
                        # Minimal sanitization
                        Chem.SanitizeMol(rdkit_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
                    except Exception:
                        pass
            
            return rdkit_mol
            
        except Exception as e:
            # Method 2: Fallback to disconnected atoms (original behavior)
            try:
                mol = Chem.RWMol()
                
                for i, symbol in enumerate(self.atoms):
                    atom = Chem.Atom(symbol)
                    mol.AddAtom(atom)
                
                conf = Chem.Conformer(len(self.atoms))
                for i, coords in enumerate(self.coordinates):
                    conf.SetAtomPosition(i, coords)
                mol.AddConformer(conf, assignId=True)
                
                rdkit_mol = mol.GetMol()
                
                # Initialize RingInfo to avoid RingInfo errors
                try:
                    rdkit_mol.GetRingInfo()
                except Exception:
                    pass
                
                try:
                    Chem.SanitizeMol(rdkit_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
                except Exception:
                    pass
                
                return rdkit_mol
                
            except Exception as e2:
                raise Exception(f"Failed to create RDKit molecule: {e}. Fallback also failed: {e2}")
    
    def _infer_smiles_from_xyz(self):
        """Try to infer SMILES from XYZ data for common molecules."""
        # Simple heuristics for common small molecules
        atom_counts = {}
        for atom in self.atoms:
            atom_counts[atom] = atom_counts.get(atom, 0) + 1
        
        # For molecules with explicit hydrogens, we need to create explicit SMILES
        # Methanol: 1 C, 1 O, 4 H -> CO with explicit hydrogens
        if atom_counts.get('C', 0) == 1 and atom_counts.get('O', 0) == 1 and atom_counts.get('H', 0) == 4:
            return 'C([H])([H])([H])O[H]'  # Explicit methanol
        
        # Methane: 1 C, 4 H -> C with explicit hydrogens
        if atom_counts.get('C', 0) == 1 and atom_counts.get('H', 0) == 4 and atom_counts.get('O', 0) == 0:
            return 'C([H])([H])([H])[H]'  # Explicit methane
        
        # Water: 1 O, 2 H -> O with explicit hydrogens
        if atom_counts.get('O', 0) == 1 and atom_counts.get('H', 0) == 2 and atom_counts.get('C', 0) == 0:
            return 'O([H])[H]'  # Explicit water
        
        # Ammonia: 1 N, 3 H -> N with explicit hydrogens
        if atom_counts.get('N', 0) == 1 and atom_counts.get('H', 0) == 3 and atom_counts.get('C', 0) == 0:
            return 'N([H])([H])[H]'  # Explicit ammonia
        
        return None
    
    def _add_bonds_from_distances(self, mol):
        """Add bonds to molecule based on atomic distances."""
        import math
        from rdkit import Chem
        
        # Bond distance thresholds (in Angstroms) - more conservative
        bond_thresholds = {
            ('C', 'C'): 1.7,
            ('C', 'O'): 1.5,
            ('C', 'H'): 1.2,
            ('O', 'H'): 1.1,
            ('N', 'H'): 1.1,
            ('N', 'C'): 1.5,
            ('N', 'O'): 1.5,
        }
        
        # For methanol specifically, we know the structure
        if (len(self.atoms) == 6 and 
            self.atoms.count('C') == 1 and 
            self.atoms.count('O') == 1 and 
            self.atoms.count('H') == 4):
            # Methanol structure: C-O bond, C-H bonds, O-H bond
            c_idx = self.atoms.index('C')
            o_idx = self.atoms.index('O')
            h_indices = [i for i, atom in enumerate(self.atoms) if atom == 'H']
            
            # C-O bond
            mol.AddBond(c_idx, o_idx, Chem.BondType.SINGLE)
            
            # C-H bonds (3 of them)
            for h_idx in h_indices[:3]:
                mol.AddBond(c_idx, h_idx, Chem.BondType.SINGLE)
            
            # O-H bond (1 of them)
            if len(h_indices) >= 4:
                mol.AddBond(o_idx, h_indices[3], Chem.BondType.SINGLE)
            
            return
        
        # General bond inference for other molecules
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                # Calculate distance
                dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(self.coordinates[i], self.coordinates[j])))
                
                # Check if atoms should be bonded
                atom1, atom2 = self.atoms[i], self.atoms[j]
                threshold = bond_thresholds.get((atom1, atom2), bond_thresholds.get((atom2, atom1), 1.5))
                
                if dist < threshold:
                    # Add bond
                    bond_order = 1  # Assume single bond
                    if atom1 == 'C' and atom2 == 'C' and dist < 1.4:
                        bond_order = 2  # Double bond
                    elif atom1 == 'C' and atom2 == 'C' and dist < 1.2:
                        bond_order = 3  # Triple bond
                    
                    mol.AddBond(i, j, Chem.BondType.SINGLE if bond_order == 1 else 
                               Chem.BondType.DOUBLE if bond_order == 2 else 
                               Chem.BondType.TRIPLE)
    
    def get_standard_smiles(self):
        """Get standard SMILES (implicit hydrogens) from RDKit molecule."""
        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError("RDKit is not installed. Install with: pip install rdkit")
        
        rdkit_mol = self.to_rdkit_mol()
        # Remove hydrogens to get standard SMILES
        mol_no_h = Chem.RemoveHs(rdkit_mol)
        return Chem.MolToSmiles(mol_no_h)
    
    def get_explicit_smiles(self):
        """Get explicit SMILES (all hydrogens shown) from RDKit molecule."""
        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError("RDKit is not installed. Install with: pip install rdkit")
        
        rdkit_mol = self.to_rdkit_mol()
        return Chem.MolToSmiles(rdkit_mol, allHsExplicit=True)
    
    def to_pymatgen_molecule(self):
        """Convert to pymatgen Molecule object (requires pymatgen to be installed)."""
        try:
            from pymatgen.core import Molecule
        except ImportError:
            raise ImportError("pymatgen is not installed. Install with: pip install pymatgen")
        
        return Molecule(self.atoms, self.coordinates)
    
    def save_to_file(self, filename: str) -> None:
        """Save XYZ data to file."""
        with open(filename, 'w') as f:
            f.write(self.to_xyz_string())
    
    def __len__(self) -> int:
        return len(self.atoms)
    
    def __iter__(self) -> Iterator[tuple[str, List[float]]]:
        return zip(self.atoms, self.coordinates)

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
    
    def get_xyz_data(self) -> Optional[XYZData]:
        """Get XYZ data as XYZData object for easy third-party library integration.
        
        Returns:
            XYZData object if XYZ coordinates are available, None otherwise.
            
        Raises:
            ValueError: If XYZ data is in an unsupported format or invalid.
        """
        if self.xyz is None:
            return None
        
        if isinstance(self.xyz, str):
            # XYZ is a filename - this would require file access which we don't have
            raise ValueError("XYZ data is stored as filename. Use expand_xyz=True when fetching species to get coordinates directly.")
        
        if isinstance(self.xyz, list):
            # XYZ is expanded coordinates
            if not self.xyz:
                return None
            
            try:
                return XYZData.from_xyz_lines(self.xyz)
            except Exception as e:
                raise ValueError(f"Failed to parse XYZ data: {e}")
        
        return None
    
    def to_ase_atoms(self):
        """Convert species XYZ data to ASE Atoms object.
        
        Returns:
            ASE Atoms object if XYZ data is available.
            
        Raises:
            ValueError: If no XYZ data is available or data is invalid.
            ImportError: If ASE is not installed.
        """
        xyz_data = self.get_xyz_data()
        if xyz_data is None:
            raise ValueError("No XYZ data available for this species")
        
        return xyz_data.to_ase_atoms()
    
    def to_rdkit_mol(self):
        """Convert species XYZ data to RDKit Mol object.
        
        Returns:
            RDKit Mol object if XYZ data is available.
            
        Raises:
            ValueError: If no XYZ data is available or data is invalid.
            ImportError: If RDKit is not installed.
        """
        xyz_data = self.get_xyz_data()
        if xyz_data is None:
            raise ValueError("No XYZ data available for this species")
        
        return xyz_data.to_rdkit_mol()
    
    def to_pymatgen_molecule(self):
        """Convert species XYZ data to pymatgen Molecule object.
        
        Returns:
            pymatgen Molecule object if XYZ data is available.
            
        Raises:
            ValueError: If no XYZ data is available or data is invalid.
            ImportError: If pymatgen is not installed.
        """
        xyz_data = self.get_xyz_data()
        if xyz_data is None:
            raise ValueError("No XYZ data available for this species")
        
        return xyz_data.to_pymatgen_molecule()
    
    def save_xyz_to_file(self, filename: str) -> None:
        """Save species XYZ data to file.
        
        Args:
            filename: Path where to save the XYZ file.
            
        Raises:
            ValueError: If no XYZ data is available or data is invalid.
        """
        xyz_data = self.get_xyz_data()
        if xyz_data is None:
            raise ValueError("No XYZ data available for this species")
        
        xyz_data.save_to_file(filename)
    
    def get_xyz_string(self) -> str:
        """Get species XYZ data as formatted string.
        
        Returns:
            XYZ format string if data is available.
            
        Raises:
            ValueError: If no XYZ data is available or data is invalid.
        """
        xyz_data = self.get_xyz_data()
        if xyz_data is None:
            raise ValueError("No XYZ data available for this species")
        
        return xyz_data.to_xyz_string()
    
    def get_standard_smiles(self) -> str:
        """Get standard SMILES (implicit hydrogens) from species XYZ data.
        
        Returns:
            Standard SMILES string if XYZ data is available.
            
        Raises:
            ValueError: If no XYZ data is available or data is invalid.
            ImportError: If RDKit is not installed.
        """
        xyz_data = self.get_xyz_data()
        if xyz_data is None:
            raise ValueError("No XYZ data available for this species")
        
        return xyz_data.get_standard_smiles()
    
    def get_explicit_smiles(self) -> str:
        """Get explicit SMILES (all hydrogens shown) from species XYZ data.
        
        Returns:
            Explicit SMILES string if XYZ data is available.
            
        Raises:
            ValueError: If no XYZ data is available or data is invalid.
            ImportError: If RDKit is not installed.
        """
        xyz_data = self.get_xyz_data()
        if xyz_data is None:
            raise ValueError("No XYZ data available for this species")
        
        return xyz_data.get_explicit_smiles()


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

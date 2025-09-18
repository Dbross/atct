#!/usr/bin/env python3
"""Tests for XYZ helper functionality without external dependencies."""

import pytest
import tempfile
import os
from atct.api.models import Species, XYZData


class TestXYZData:
    """Test XYZData class functionality."""
    
    def test_xyz_data_creation(self):
        """Test basic XYZData creation."""
        atoms = ["C", "O", "H", "H", "H", "H"]
        coordinates = [
            [0.0, 0.0, 0.0],
            [1.4, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0]
        ]
        
        xyz_data = XYZData(atoms=atoms, coordinates=coordinates, comment="Methanol")
        
        assert len(xyz_data) == 6
        assert xyz_data.atoms == atoms
        assert xyz_data.coordinates == coordinates
        assert xyz_data.comment == "Methanol"
    
    def test_xyz_data_creation_mismatch_error(self):
        """Test XYZData creation with mismatched atoms and coordinates."""
        atoms = ["C", "O", "H"]
        coordinates = [
            [0.0, 0.0, 0.0],
            [1.4, 0.0, 0.0]
        ]
        
        with pytest.raises(ValueError, match="Number of atoms must match number of coordinate sets"):
            XYZData(atoms=atoms, coordinates=coordinates)
    
    def test_from_xyz_lines(self):
        """Test creating XYZData from XYZ format lines."""
        xyz_lines = [
            "6",
            "Methanol",
            "C 0.0 0.0 0.0",
            "O 1.4 0.0 0.0",
            "H 0.0 1.0 0.0",
            "H 0.0 -1.0 0.0",
            "H 0.0 0.0 1.0",
            "H 0.0 0.0 -1.0"
        ]
        
        xyz_data = XYZData.from_xyz_lines(xyz_lines)
        
        assert len(xyz_data) == 6
        assert xyz_data.atoms == ["C", "O", "H", "H", "H", "H"]
        assert xyz_data.comment == "Methanol"
        assert len(xyz_data.coordinates) == 6
        assert xyz_data.coordinates[0] == [0.0, 0.0, 0.0]
        assert xyz_data.coordinates[1] == [1.4, 0.0, 0.0]
    
    def test_from_xyz_lines_empty(self):
        """Test creating XYZData from empty lines."""
        with pytest.raises(ValueError, match="Empty XYZ data"):
            XYZData.from_xyz_lines([])
    
    def test_from_xyz_lines_invalid_format(self):
        """Test creating XYZData from invalid format."""
        # Missing number of atoms
        with pytest.raises(ValueError, match="Invalid XYZ format"):
            XYZData.from_xyz_lines(["invalid"])
        
        # Not enough coordinate lines
        with pytest.raises(ValueError, match="Not enough coordinate lines"):
            XYZData.from_xyz_lines(["2", "Comment", "C 0 0 0"])
        
        # Invalid coordinate line
        with pytest.raises(ValueError, match="Invalid coordinate line"):
            XYZData.from_xyz_lines(["1", "Comment", "C 0 0"])  # Missing coordinate
    
    def test_to_xyz_string(self):
        """Test converting XYZData to XYZ format string."""
        atoms = ["C", "O", "H", "H", "H", "H"]
        coordinates = [
            [0.0, 0.0, 0.0],
            [1.4, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0]
        ]
        
        xyz_data = XYZData(atoms=atoms, coordinates=coordinates, comment="Methanol")
        xyz_string = xyz_data.to_xyz_string()
        
        lines = xyz_string.split('\n')
        assert lines[0] == "6"
        assert lines[1] == "Methanol"
        assert "C 0.000000 0.000000 0.000000" in lines
        assert "O 1.400000 0.000000 0.000000" in lines
    
    def test_to_xyz_string_no_comment(self):
        """Test converting XYZData to string without comment."""
        atoms = ["C", "O"]
        coordinates = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]
        
        xyz_data = XYZData(atoms=atoms, coordinates=coordinates)
        xyz_string = xyz_data.to_xyz_string()
        
        lines = xyz_string.split('\n')
        assert lines[0] == "2"
        assert lines[1] == ""  # Empty comment line
    
    def test_save_to_file(self):
        """Test saving XYZData to file."""
        atoms = ["C", "O", "H", "H", "H", "H"]
        coordinates = [
            [0.0, 0.0, 0.0],
            [1.4, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0]
        ]
        
        xyz_data = XYZData(atoms=atoms, coordinates=coordinates, comment="Methanol")
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            temp_filename = f.name
        
        try:
            xyz_data.save_to_file(temp_filename)
            
            with open(temp_filename, 'r') as f:
                content = f.read()
            
            lines = content.split('\n')
            assert lines[0] == "6"
            assert lines[1] == "Methanol"
            assert "C 0.000000 0.000000 0.000000" in content
        finally:
            if os.path.exists(temp_filename):
                os.unlink(temp_filename)
    
    def test_iteration(self):
        """Test XYZData iteration."""
        atoms = ["C", "O", "H"]
        coordinates = [
            [0.0, 0.0, 0.0],
            [1.4, 0.0, 0.0],
            [0.0, 1.0, 0.0]
        ]
        
        xyz_data = XYZData(atoms=atoms, coordinates=coordinates)
        
        atom_coord_pairs = list(xyz_data)
        assert len(atom_coord_pairs) == 3
        assert atom_coord_pairs[0] == ("C", [0.0, 0.0, 0.0])
        assert atom_coord_pairs[1] == ("O", [1.4, 0.0, 0.0])
        assert atom_coord_pairs[2] == ("H", [0.0, 1.0, 0.0])
    
    def test_external_dependencies_not_imported(self):
        """Test that external dependencies are not imported at module level."""
        # This test ensures that ASE, RDKit, and pymatgen are not imported
        # when the module is loaded, only when the methods are called
        
        atoms = ["C", "O"]
        coordinates = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]
        xyz_data = XYZData(atoms=atoms, coordinates=coordinates)
        
        # Check which dependencies are available
        ase_available = False
        rdkit_available = False
        pymatgen_available = False
        
        try:
            import ase
            ase_available = True
        except ImportError:
            pass
        
        try:
            import rdkit
            rdkit_available = True
        except ImportError:
            pass
        
        try:
            import pymatgen
            pymatgen_available = True
        except ImportError:
            pass
        
        # Test ASE
        if ase_available:
            # If ASE is available, the method should work
            ase_atoms = xyz_data.to_ase_atoms()
            assert ase_atoms is not None
            assert len(ase_atoms) == 2
        else:
            with pytest.raises(ImportError, match="ASE is not installed"):
                xyz_data.to_ase_atoms()
        
        # Test RDKit
        if rdkit_available:
            # If RDKit is available, the method should work
            rdkit_mol = xyz_data.to_rdkit_mol()
            assert rdkit_mol is not None
        else:
            with pytest.raises(ImportError, match="RDKit is not installed"):
                xyz_data.to_rdkit_mol()
        
        # Test pymatgen
        if pymatgen_available:
            # If pymatgen is available, the method should work
            pymatgen_mol = xyz_data.to_pymatgen_molecule()
            assert pymatgen_mol is not None
        else:
            with pytest.raises(ImportError, match="pymatgen is not installed"):
                xyz_data.to_pymatgen_molecule()


class TestSpeciesXYZHelpers:
    """Test Species XYZ helper methods."""
    
    def test_species_with_expanded_xyz(self):
        """Test Species with expanded XYZ data."""
        species = Species(
            atct_tn_version="1.0",
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            delta_h_0k="-201.0",
            delta_h_298k="-200.7",
            delta_h_298k_uncertainty="0.1",
            smiles="CO",
            casrn="67-56-1",
            inchi="InChI=1S/CH4O/c1-2/h2H,1H3",
            inchi_key="OKKJLVBELUTLKV-UHFFFAOYSA-N",
            charge=0,
            xyz=[
                "6",
                "Methanol",
                "C 0.0 0.0 0.0",
                "O 1.4 0.0 0.0",
                "H 0.0 1.0 0.0",
                "H 0.0 -1.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 0.0 -1.0"
            ]
        )
        
        # Test get_xyz_data
        xyz_data = species.get_xyz_data()
        assert xyz_data is not None
        assert len(xyz_data) == 6
        assert xyz_data.atoms == ["C", "O", "H", "H", "H", "H"]
        assert xyz_data.comment == "Methanol"
        
        # Test get_xyz_string
        xyz_string = species.get_xyz_string()
        assert "6" in xyz_string
        assert "Methanol" in xyz_string
        assert "C 0.000000 0.000000 0.000000" in xyz_string
    
    def test_species_with_no_xyz(self):
        """Test Species with no XYZ data."""
        species = Species(
            atct_tn_version="1.0",
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            delta_h_0k="-201.0",
            delta_h_298k="-200.7",
            delta_h_298k_uncertainty="0.1",
            smiles="CO",
            casrn="67-56-1",
            inchi="InChI=1S/CH4O/c1-2/h2H,1H3",
            inchi_key="OKKJLVBELUTLKV-UHFFFAOYSA-N",
            charge=0,
            xyz=None
        )
        
        # Test get_xyz_data
        xyz_data = species.get_xyz_data()
        assert xyz_data is None
        
        # Test methods that require XYZ data
        with pytest.raises(ValueError, match="No XYZ data available"):
            species.get_xyz_string()
        
        with pytest.raises(ValueError, match="No XYZ data available"):
            species.to_ase_atoms()
        
        with pytest.raises(ValueError, match="No XYZ data available"):
            species.to_rdkit_mol()
        
        with pytest.raises(ValueError, match="No XYZ data available"):
            species.to_pymatgen_molecule()
    
    def test_species_with_filename_xyz(self):
        """Test Species with XYZ filename (not expanded)."""
        species = Species(
            atct_tn_version="1.0",
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            delta_h_0k="-201.0",
            delta_h_298k="-200.7",
            delta_h_298k_uncertainty="0.1",
            smiles="CO",
            casrn="67-56-1",
            inchi="InChI=1S/CH4O/c1-2/h2H,1H3",
            inchi_key="OKKJLVBELUTLKV-UHFFFAOYSA-N",
            charge=0,
            xyz="methanol.xyz"
        )
        
        # Test get_xyz_data with filename
        with pytest.raises(ValueError, match="XYZ data is stored as filename"):
            species.get_xyz_data()
    
    def test_species_with_empty_xyz_list(self):
        """Test Species with empty XYZ list."""
        species = Species(
            atct_tn_version="1.0",
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            delta_h_0k="-201.0",
            delta_h_298k="-200.7",
            delta_h_298k_uncertainty="0.1",
            smiles="CO",
            casrn="67-56-1",
            inchi="InChI=1S/CH4O/c1-2/h2H,1H3",
            inchi_key="OKKJLVBELUTLKV-UHFFFAOYSA-N",
            charge=0,
            xyz=[]
        )
        
        # Test get_xyz_data with empty list
        xyz_data = species.get_xyz_data()
        assert xyz_data is None
    
    def test_species_save_xyz_to_file(self):
        """Test saving Species XYZ data to file."""
        species = Species(
            atct_tn_version="1.0",
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            delta_h_0k="-201.0",
            delta_h_298k="-200.7",
            delta_h_298k_uncertainty="0.1",
            smiles="CO",
            casrn="67-56-1",
            inchi="InChI=1S/CH4O/c1-2/h2H,1H3",
            inchi_key="OKKJLVBELUTLKV-UHFFFAOYSA-N",
            charge=0,
            xyz=[
                "6",
                "Methanol",
                "C 0.0 0.0 0.0",
                "O 1.4 0.0 0.0",
                "H 0.0 1.0 0.0",
                "H 0.0 -1.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 0.0 -1.0"
            ]
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            temp_filename = f.name
        
        try:
            species.save_xyz_to_file(temp_filename)
            
            with open(temp_filename, 'r') as f:
                content = f.read()
            
            assert "6" in content
            assert "Methanol" in content
            assert "C 0.000000 0.000000 0.000000" in content
        finally:
            if os.path.exists(temp_filename):
                os.unlink(temp_filename)
    
    def test_species_external_dependencies_not_imported(self):
        """Test that external dependencies are not imported when calling conversion methods."""
        species = Species(
            atct_tn_version="1.0",
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            delta_h_0k="-201.0",
            delta_h_298k="-200.7",
            delta_h_298k_uncertainty="0.1",
            smiles="CO",
            casrn="67-56-1",
            inchi="InChI=1S/CH4O/c1-2/h2H,1H3",
            inchi_key="OKKJLVBELUTLKV-UHFFFAOYSA-N",
            charge=0,
            xyz=[
                "6",
                "Methanol",
                "C 0.0 0.0 0.0",
                "O 1.4 0.0 0.0",
                "H 0.0 1.0 0.0",
                "H 0.0 -1.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 0.0 -1.0"
            ]
        )
        
        # Check which dependencies are available
        ase_available = False
        rdkit_available = False
        pymatgen_available = False
        
        try:
            import ase
            ase_available = True
        except ImportError:
            pass
        
        try:
            import rdkit
            rdkit_available = True
        except ImportError:
            pass
        
        try:
            import pymatgen
            pymatgen_available = True
        except ImportError:
            pass
        
        # Test ASE
        if ase_available:
            # If ASE is available, the method should work
            ase_atoms = species.to_ase_atoms()
            assert ase_atoms is not None
            assert len(ase_atoms) == 6
        else:
            with pytest.raises(ImportError, match="ASE is not installed"):
                species.to_ase_atoms()
        
        # Test RDKit
        if rdkit_available:
            # If RDKit is available, the method should work
            rdkit_mol = species.to_rdkit_mol()
            assert rdkit_mol is not None
        else:
            with pytest.raises(ImportError, match="RDKit is not installed"):
                species.to_rdkit_mol()
        
        # Test pymatgen
        if pymatgen_available:
            # If pymatgen is available, the method should work
            pymatgen_mol = species.to_pymatgen_molecule()
            assert pymatgen_mol is not None
        else:
            with pytest.raises(ImportError, match="pymatgen is not installed"):
                species.to_pymatgen_molecule()
    
    def test_species_invalid_xyz_data(self):
        """Test Species with invalid XYZ data."""
        species = Species(
            atct_tn_version="1.0",
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            delta_h_0k="-201.0",
            delta_h_298k="-200.7",
            delta_h_298k_uncertainty="0.1",
            smiles="CO",
            casrn="67-56-1",
            inchi="InChI=1S/CH4O/c1-2/h2H,1H3",
            inchi_key="OKKJLVBELUTLKV-UHFFFAOYSA-N",
            charge=0,
            xyz=[
                "invalid",
                "Methanol",
                "C 0.0 0.0 0.0"
            ]
        )
        
        # Test get_xyz_data with invalid data
        with pytest.raises(ValueError, match="Failed to parse XYZ data"):
            species.get_xyz_data()


if __name__ == "__main__":
    pytest.main([__file__])

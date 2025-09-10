"""Comprehensive tests for ATcT v1 API functionality."""

import pytest
import os
from unittest.mock import patch, MagicMock
from atct.api import (
    healthcheck,
    get_species,
    get_species_by_atctid,
    get_species_by_casrn,
    get_species_by_inchi,
    get_species_by_inchikey,
    get_species_by_smiles,
    get_species_by_formula,
    get_species_by_name,
    search_species,
    get_species_covariance_matrix,
    get_species_covariance_by_ids,
    get_species_covariance_by_atctid,
    get_species_simple,
    get_all_species,
    create_reaction_calculator,
    calculate_reaction_enthalpy,
)
from atct.api.models import Species, CovarianceMatrix, Covariance2x2, Page, ReactionSpecies, ReactionResult, ReactionCalculator
from atct.api.pandas_io import as_dataframe
from atct.api.exceptions import NotFound, ATCTError


class TestHealthCheck:
    """Test health check functionality."""
    
    def test_healthcheck_success(self):
        """Test successful health check."""
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = {"ok": True, "service": "ATcT API"}
            result = healthcheck()
            assert result is True
            mock_request.assert_called_once()
    
    def test_healthcheck_failure(self):
        """Test health check failure."""
        with patch('atct.api.request_json') as mock_request:
            mock_request.side_effect = Exception("Network error")
            result = healthcheck()
            assert result is False


class TestGetSpecies:
    """Test get_species functionality."""
    
    def test_get_species_basic(self):
        """Test basic species retrieval by ATcT ID."""
        mock_data = {
            "ATcT_TN_Version": "1.123",
            "ATcT_ID": "67-56-1*0",
            "Name": "Methanol",
            "Formula": "CH4O",
            "∆fH_0K": "-201.0",
            "∆fH_298K": "-201.0",
            "∆fH_298K_uncertainty": "0.1",
            "SMILES": "CO",
            "CASRN": "67-56-1",
            "InChI": "InChI=1S/CH4O/c1-2/h2H,1H3",
            "InChI_Key": "OKKJLVBELUTLKV-UHFFFAOYSA-N",
            "charge": 0,
            "XYZ": None
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            species = get_species("67-56-1*0")
            
            assert isinstance(species, Species)
            assert species.atct_id == "67-56-1*0"
            assert species.name == "Methanol"
            assert species.formula == "CH4O"
            assert species.smiles == "CO"
            assert species.casrn == "67-56-1"
            assert species.inchi == "InChI=1S/CH4O/c1-2/h2H,1H3"
            assert species.inchi_key == "OKKJLVBELUTLKV-UHFFFAOYSA-N"
            assert species.charge == 0
            mock_request.assert_called_once()
    
    def test_get_species_with_expand_xyz(self):
        """Test species retrieval with XYZ expansion."""
        mock_data = {
            "ATcT_ID": "67-56-1*0",
            "Name": "Methanol",
            "Formula": "CH4O",
            "XYZ": ["C 0.0 0.0 0.0", "O 1.4 0.0 0.0", "H 0.0 1.0 0.0", "H 0.0 -1.0 0.0", "H 0.0 0.0 1.0", "H 0.0 0.0 -1.0"]
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            species = get_species("67-56-1*0", expand_xyz=True)
            
            assert species.atct_id == "67-56-1*0"
            assert isinstance(species.xyz, list)
            assert len(species.xyz) == 6
            mock_request.assert_called_once()
    
    def test_get_species_by_atctid_list_response(self):
        """Test get_species_by_atctid with list response."""
        mock_data = [{
            "ATcT_ID": "67-56-1*0",
            "Name": "Methanol",
            "Formula": "CH4O",
            "∆fH_298K": "-201.0",
            "SMILES": "CO",
            "CASRN": "67-56-1",
            "InChI": "InChI=1S/CH4O/c1-2/h2H,1H3",
            "InChI_Key": "OKKJLVBELUTLKV-UHFFFAOYSA-N",
            "charge": 0,
            "XYZ": None
        }]
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            species = get_species_by_atctid("67-56-1*0")
            
            assert isinstance(species, Species)
            assert species.atct_id == "67-56-1*0"
            assert species.name == "Methanol"
            mock_request.assert_called_once()


class TestGetSpeciesBy:
    """Test species retrieval by various identifiers."""
    
    def test_get_species_by_casrn(self):
        """Test species retrieval by CAS Registry Number."""
        mock_data = {
            "items": [{
                "ATcT_ID": "67-56-1*0",
                "Name": "Methanol",
                "Formula": "CH4O",
                "CASRN": "67-56-1"
            }],
            "total": 1,
            "limit": 50,
            "offset": 0
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            page = get_species_by_casrn("67-56-1")
            
            assert isinstance(page, Page)
            assert len(page.items) == 1
            assert page.total == 1
            assert isinstance(page.items[0], Species)
            assert page.items[0].casrn == "67-56-1"
    
    def test_get_species_by_inchi(self):
        """Test species retrieval by InChI."""
        mock_data = {
            "items": [{
                "ATcT_ID": "67-56-1*0",
                "Name": "Methanol",
                "Formula": "CH4O",
                "InChI": "1S/CH4O/c1-2/h2H,1H3"
            }],
            "total": 1,
            "limit": 50,
            "offset": 0
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            page = get_species_by_inchi("1S/CH4O/c1-2/h2H,1H3")
            
            assert isinstance(page, Page)
            assert len(page.items) == 1
            assert page.total == 1
    
    def test_get_species_by_smiles(self):
        """Test species retrieval by SMILES."""
        mock_data = {
            "items": [{
                "ATcT_ID": "67-56-1*0",
                "Name": "Methanol",
                "Formula": "CH4O",
                "SMILES": "CO"
            }],
            "total": 1,
            "limit": 50,
            "offset": 0
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            page = get_species_by_smiles("CO")
            
            assert isinstance(page, Page)
            assert len(page.items) == 1
            assert page.items[0].smiles == "CO"
    
    def test_get_species_by_inchikey(self):
        """Test species retrieval by InChI Key."""
        mock_data = {
            "items": [{
                "ATcT_ID": "67-56-1*0",
                "Name": "Methanol",
                "Formula": "CH4O",
                "InChI_Key": "OKKJLVBELUTLKV-UHFFFAOYSA-N"
            }],
            "total": 1,
            "limit": 50,
            "offset": 0
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            page = get_species_by_inchikey("OKKJLVBELUTLKV-UHFFFAOYSA-N")
            
            assert isinstance(page, Page)
            assert len(page.items) == 1
            assert page.items[0].inchi_key == "OKKJLVBELUTLKV-UHFFFAOYSA-N"
    
    def test_get_species_by_formula(self):
        """Test species retrieval by formula."""
        mock_data = {
            "items": [{
                "ATcT_ID": "67-56-1*0",
                "Name": "Methanol",
                "Formula": "CH4O"
            }],
            "total": 1,
            "limit": 50,
            "offset": 0
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            page = get_species_by_formula("CH4O")
            
            assert isinstance(page, Page)
            assert len(page.items) == 1
            assert page.items[0].formula == "CH4O"
    
    def test_get_species_by_formula_with_phase(self):
        """Test species retrieval by formula with phase filter."""
        mock_data = {
            "items": [{
                "ATcT_ID": "67-56-1*0",
                "Name": "Methanol",
                "Formula": "CH4O"
            }],
            "total": 1,
            "limit": 50,
            "offset": 0
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            page = get_species_by_formula("CH4O", phase="g")
            
            assert isinstance(page, Page)
            # Verify phase parameter was passed
            call_args = mock_request.call_args
            assert "phase" in call_args[1]["params"]
            assert call_args[1]["params"]["phase"] == "g"
    
    def test_get_species_by_name(self):
        """Test species retrieval by name."""
        mock_data = {
            "items": [{
                "ATcT_ID": "67-56-1*0",
                "Name": "Methanol",
                "Formula": "CH4O"
            }],
            "total": 1,
            "limit": 50,
            "offset": 0
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            page = get_species_by_name("Methanol")
            
            assert isinstance(page, Page)
            assert len(page.items) == 1
            assert page.items[0].name == "Methanol"


class TestSimpleAPI:
    """Test simple API endpoints."""
    
    def test_get_species_simple_by_smiles(self):
        """Test simple API with SMILES."""
        mock_data = [{
            "ATcT_ID": "67-56-1*0",
            "Name": "Methanol",
            "Formula": "CH4O",
            "SMILES": "CO"
        }]
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            species = get_species_simple(smiles="CO")
            
            assert isinstance(species, Species)
            assert species.atct_id == "67-56-1*0"
            assert species.name == "Methanol"
            mock_request.assert_called_once()
    
    def test_get_species_simple_by_name(self):
        """Test simple API with name."""
        mock_data = [{
            "ATcT_ID": "67-56-1*0",
            "Name": "Methanol",
            "Formula": "CH4O"
        }]
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            species = get_species_simple(name="methanol")
            
            assert isinstance(species, Species)
            assert species.name == "Methanol"
            mock_request.assert_called_once()
    
    def test_get_species_simple_no_params(self):
        """Test simple API with no parameters."""
        with pytest.raises(ValueError, match="At least one parameter must be provided"):
            get_species_simple()
    
    def test_get_all_species(self):
        """Test get all species."""
        mock_data = [{
            "ATcT_ID": "67-56-1*0",
            "Name": "Methanol",
            "Formula": "CH4O"
        }]
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            species_list = get_all_species()
            
            assert isinstance(species_list, list)
            assert len(species_list) == 1
            assert isinstance(species_list[0], Species)
            assert species_list[0].name == "Methanol"
            mock_request.assert_called_once()


class TestSearchSpecies:
    """Test species search functionality."""
    
    def test_search_species(self):
        """Test general species search."""
        mock_data = {
            "items": [
                {
                    "ATcT_ID": "67-56-1*0",
                    "Name": "Methanol",
                    "Formula": "CH4O"
                },
                {
                    "ATcT_ID": "67-56-1*500",
                    "Name": "Methanol (different phase)",
                    "Formula": "CH4O"
                }
            ],
            "total": 2,
            "limit": 50,
            "offset": 0
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            page = search_species("methanol")
            
            assert isinstance(page, Page)
            assert len(page.items) == 2
            assert page.total == 2
            assert all(isinstance(item, Species) for item in page.items)
    
    def test_search_species_with_pagination(self):
        """Test species search with custom pagination."""
        mock_data = {
            "items": [{
                "ATcT_ID": "67-56-1*0",
                "Name": "Methanol",
                "Formula": "CH4O"
            }],
            "total": 10,
            "limit": 5,
            "offset": 5
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            page = search_species("methanol", limit=5, offset=5)
            
            assert page.limit == 5
            assert page.offset == 5
            assert page.total == 10


class TestCovariance:
    """Test covariance functionality."""
    
    def test_get_species_covariance_by_ids(self):
        """Test covariance retrieval by numeric IDs."""
        mock_data = {
            "labels": ["ΔH(A)", "ΔH(B)"],
            "units": "kJ/mol",
            "matrix": [[0.01, 0.005], [0.005, 0.02]]
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            cov = get_species_covariance_by_ids(500, 0)
            
            assert isinstance(cov, Covariance2x2)
            assert cov.labels == ["ΔH(A)", "ΔH(B)"]
            assert cov.units == "kJ/mol"
            assert cov.matrix == [[0.01, 0.005], [0.005, 0.02]]
    
    def test_get_species_covariance_by_atctid(self):
        """Test covariance retrieval by ATcT IDs."""
        mock_data = {
            "labels": ["ΔH(A)", "ΔH(B)"],
            "units": "kJ/mol",
            "matrix": [[0.01, 0.005], [0.005, 0.02]]
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            cov = get_species_covariance_by_atctid("67-56-1*500", "67-56-1*0")
            
            assert isinstance(cov, Covariance2x2)
            assert cov.labels == ["ΔH(A)", "ΔH(B)"]
            assert cov.units == "kJ/mol"
            assert cov.matrix == [[0.01, 0.005], [0.005, 0.02]]
    
    def test_get_species_covariance_matrix(self):
        """Test covariance matrix retrieval for multiple species."""
        mock_data = {
            "species_ids": ["1", "2"],
            "species_atctids": ["67-56-1*0", "67-56-1*500"],
            "matrix": [[0.01, 0.005], [0.005, 0.02]],
            "units": "kJ/mol"
        }
        
        with patch('atct.api.request_json') as mock_request:
            mock_request.return_value = mock_data
            cov_matrix = get_species_covariance_matrix(atctids=["67-56-1*0", "67-56-1*500"])
            
            assert isinstance(cov_matrix, CovarianceMatrix)
            assert cov_matrix.species_atctids == ["67-56-1*0", "67-56-1*500"]
            assert cov_matrix.units == "kJ/mol"
            assert cov_matrix.matrix == [[0.01, 0.005], [0.005, 0.02]]
    
    def test_get_species_covariance_matrix_no_params(self):
        """Test covariance matrix with no parameters."""
        with pytest.raises(ValueError, match="Either ids or atctids must be provided"):
            get_species_covariance_matrix()
    
    def test_get_species_covariance_matrix_both_params(self):
        """Test covariance matrix with both parameters."""
        with pytest.raises(ValueError, match="Provide either ids or atctids, not both"):
            get_species_covariance_matrix(ids=["1", "2"], atctids=["67-56-1*0", "67-56-1*500"])


class TestModels:
    """Test data models."""
    
    def test_species_from_dict(self):
        """Test Species creation from API response."""
        data = {
            "ATcT_TN_Version": "1.123",
            "ATcT_ID": "67-56-1*0",
            "Name": "Methanol",
            "Formula": "CH4O",
            "∆fH_0K": "-201.0",
            "∆fH_298K": "-201.0",
            "∆fH_298K_uncertainty": "0.1",
            "SMILES": "CO",
            "CASRN": "67-56-1",
            "InChI": "InChI=1S/CH4O/c1-2/h2H,1H3",
            "InChI_Key": "OKKJLVBELUTLKV-UHFFFAOYSA-N",
            "charge": 0,
            "XYZ": None
        }
        
        species = Species.from_dict(data)
        
        assert species.atct_tn_version == "1.123"
        assert species.atct_id == "67-56-1*0"
        assert species.name == "Methanol"
        assert species.formula == "CH4O"
        assert species.delta_h_0k == "-201.0"
        assert species.delta_h_298k == "-201.0"
        assert species.delta_h_298k_uncertainty == "0.1"
        assert species.smiles == "CO"
        assert species.casrn == "67-56-1"
        assert species.inchi == "InChI=1S/CH4O/c1-2/h2H,1H3"
        assert species.inchi_key == "OKKJLVBELUTLKV-UHFFFAOYSA-N"
        assert species.charge == 0
        assert species.xyz is None
    
    def test_species_to_dict(self):
        """Test Species serialization to dictionary."""
        species = Species(
            atct_tn_version="1.123",
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            delta_h_0k="-201.0",
            delta_h_298k="-201.0",
            delta_h_298k_uncertainty="0.1",
            smiles="CO",
            casrn="67-56-1",
            inchi="InChI=1S/CH4O/c1-2/h2H,1H3",
            inchi_key="OKKJLVBELUTLKV-UHFFFAOYSA-N",
            charge=0,
            xyz=None
        )
        
        data = species.to_dict()
        
        assert data["ATcT_TN_Version"] == "1.123"
        assert data["ATcT_ID"] == "67-56-1*0"
        assert data["Name"] == "Methanol"
        assert data["Formula"] == "CH4O"
        assert data["SMILES"] == "CO"
        assert data["CASRN"] == "67-56-1"
        assert data["charge"] == 0
    
    def test_covariance2x2_from_dict(self):
        """Test Covariance2x2 creation from API response."""
        data = {
            "labels": ["ΔH(A)", "ΔH(B)"],
            "units": "kJ/mol",
            "matrix": [[0.01, 0.005], [0.005, 0.02]]
        }
        
        cov = Covariance2x2.from_dict(data)
        
        assert cov.labels == ["ΔH(A)", "ΔH(B)"]
        assert cov.units == "kJ/mol"
        assert cov.matrix == [[0.01, 0.005], [0.005, 0.02]]
    
    def test_covariance2x2_defaults(self):
        """Test Covariance2x2 with default values."""
        data = {
            "matrix": [[0.01, 0.005], [0.005, 0.02]]
        }
        
        cov = Covariance2x2.from_dict(data)
        
        assert cov.labels == ["ΔH(A)", "ΔH(B)"]
        assert cov.units == ""
        assert cov.matrix == [[0.01, 0.005], [0.005, 0.02]]
    
    def test_page_creation(self):
        """Test Page creation and serialization."""
        species1 = Species(
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            atct_tn_version=None,
            delta_h_0k=None,
            delta_h_298k=None,
            delta_h_298k_uncertainty=None,
            smiles=None,
            casrn=None,
            inchi=None,
            inchi_key=None,
            charge=None,
            xyz=None
        )
        
        page = Page(items=[species1], total=1, limit=50, offset=0)
        
        assert len(page.items) == 1
        assert page.total == 1
        assert page.limit == 50
        assert page.offset == 0
        
        # Test to_dict method
        page_dict = page.to_dict()
        assert page_dict["total"] == 1
        assert page_dict["limit"] == 50
        assert page_dict["offset"] == 0
        assert len(page_dict["items"]) == 1


class TestReactionCalculations:
    """Test reaction calculation functionality."""
    
    def test_reaction_species_creation(self):
        """Test ReactionSpecies creation."""
        species = Species(
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            atct_tn_version=None,
            delta_h_0k=None,
            delta_h_298k="-200.85",
            delta_h_298k_uncertainty="0.1",
            smiles=None,
            casrn=None,
            inchi=None,
            inchi_key=None,
            charge=None,
            xyz=None
        )
        
        reaction_species = ReactionSpecies(species, -1.0)
        assert reaction_species.species == species
        assert reaction_species.stoichiometry == -1.0
    
    def test_reaction_result_creation(self):
        """Test ReactionResult creation."""
        result = ReactionResult(-200.85, 0.1, "covariance")
        assert result.delta_h == -200.85
        assert result.uncertainty == 0.1
        assert result.method == "covariance"
        assert result.units == "kJ/mol"
        assert "covariance" in str(result)
    
    def test_reaction_calculator_creation(self):
        """Test ReactionCalculator creation."""
        species1 = Species(
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            atct_tn_version=None,
            delta_h_0k=None,
            delta_h_298k="-200.85",
            delta_h_298k_uncertainty="0.1",
            inchi=None,
            inchi_key=None,
            smiles=None,
            casrn=None,
            charge=None,
            xyz=None
        )
        
        species2 = Species(
            atct_id="7727-37-9*0",
            name="Dinitrogen",
            formula="N2",
            atct_tn_version=None,
            delta_h_0k=None,
            delta_h_298k="0",
            delta_h_298k_uncertainty="exact",
            inchi=None,
            inchi_key=None,
            smiles=None,
            casrn=None,
            charge=None,
            xyz=None
        )
        
        reaction_species = [
            ReactionSpecies(species1, -1.0),
            ReactionSpecies(species2, -1.0)
        ]
        
        calculator = ReactionCalculator(reaction_species)
        assert len(calculator.reaction_species) == 2
        assert calculator.stoichiometry[0] == -1.0
        assert calculator.stoichiometry[1] == -1.0
    
    def test_reaction_calculation_methods(self):
        """Test different reaction calculation methods."""
        species1 = Species(
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            atct_tn_version=None,
            delta_h_0k=None,
            delta_h_298k="-200.85",
            delta_h_298k_uncertainty="0.1",
            inchi=None,
            inchi_key=None,
            smiles=None,
            casrn=None,
            charge=None,
            xyz=None
        )
        
        species2 = Species(
            atct_id="7727-37-9*0",
            name="Dinitrogen",
            formula="N2",
            atct_tn_version=None,
            delta_h_0k=None,
            delta_h_298k="0",
            delta_h_298k_uncertainty="exact",
            inchi=None,
            inchi_key=None,
            smiles=None,
            casrn=None,
            charge=None,
            xyz=None
        )
        
        reaction_species = [
            ReactionSpecies(species1, -1.0),
            ReactionSpecies(species2, -1.0)
        ]
        
        calculator = ReactionCalculator(reaction_species)
        
        cov_result = calculator.calculate_covariance_method()
        conv_result = calculator.calculate_conventional_method()
        
        assert isinstance(cov_result, ReactionResult)
        assert isinstance(conv_result, ReactionResult)
        
        assert cov_result.method == "covariance"
        assert conv_result.method == "conventional"
    
    def test_reaction_calculator_comparison(self):
        """Test reaction calculator method comparison."""
        species1 = Species(
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            atct_tn_version=None,
            delta_h_0k=None,
            delta_h_298k="-200.85",
            delta_h_298k_uncertainty="0.1",
            inchi=None,
            inchi_key=None,
            smiles=None,
            casrn=None,
            charge=None,
            xyz=None
        )
        
        reaction_species = [ReactionSpecies(species1, -1.0)]
        calculator = ReactionCalculator(reaction_species)
        
        comparison = calculator.compare_methods()
        
        assert 'covariance' in comparison
        assert 'conventional' in comparison
        assert 'difference' in comparison
        assert 'significance' in comparison
        
        assert isinstance(comparison['covariance'], ReactionResult)
        assert isinstance(comparison['conventional'], ReactionResult)


class TestPandasIO:
    """Test pandas integration."""
    
    def test_as_dataframe_with_species(self):
        """Test DataFrame creation from Species objects."""
        species = Species(
            atct_id="67-56-1*0",
            name="Methanol",
            formula="CH4O",
            atct_tn_version=None,
            delta_h_0k=None,
            delta_h_298k=None,
            delta_h_298k_uncertainty=None,
            smiles=None,
            casrn=None,
            inchi=None,
            inchi_key=None,
            charge=None,
            xyz=None
        )
        
        with patch('atct.api.pandas_io.pd') as mock_pd:
            mock_df = MagicMock()
            mock_pd.DataFrame.return_value = mock_df
            
            df = as_dataframe([species])
            
            mock_pd.DataFrame.assert_called_once()
            call_args = mock_pd.DataFrame.call_args[0][0]
            assert len(call_args) == 1
            assert call_args[0]["ATcT_ID"] == "67-56-1*0"
            assert call_args[0]["Name"] == "Methanol"
    
    def test_as_dataframe_with_dicts(self):
        """Test DataFrame creation from dictionaries."""
        data = [{"name": "test", "value": 123}]
        
        with patch('atct.api.pandas_io.pd') as mock_pd:
            mock_df = MagicMock()
            mock_pd.DataFrame.return_value = mock_df
            
            df = as_dataframe(data)
            
            mock_pd.DataFrame.assert_called_once_with(data)
    
    def test_as_dataframe_no_pandas(self):
        """Test DataFrame creation when pandas is not available."""
        with patch('atct.api.pandas_io.PANDAS_AVAILABLE', False):
            with pytest.raises(ImportError, match="pandas is required"):
                as_dataframe([{"test": "data"}])


class TestIntegration:
    """Integration tests with real API (if environment variable is set)."""
    
    @pytest.mark.skipif(
        not os.environ.get("ATCT_API_BASE_URL"),
        reason="ATCT_API_BASE_URL not set - skipping integration tests"
    )
    def test_real_healthcheck(self):
        """Test real health check against API."""
        result = healthcheck()
        assert isinstance(result, bool)
    
    @pytest.mark.skipif(
        not os.environ.get("ATCT_API_BASE_URL"),
        reason="ATCT_API_BASE_URL not set - skipping integration tests"
    )
    def test_real_get_species(self):
        """Test real species retrieval."""
        species = get_species("67-56-1*0")
        assert isinstance(species, Species)
        assert species.atct_id == "67-56-1*0"
        assert species.name == "Methanol"
        assert "CH3OH" in species.formula
        assert species.smiles == "CO"
        assert species.casrn == "67-56-1"
    
    @pytest.mark.skipif(
        not os.environ.get("ATCT_API_BASE_URL"),
        reason="ATCT_API_BASE_URL not set - skipping integration tests"
    )
    def test_real_search_species(self):
        """Test real species search."""
        page = search_species("methanol")
        assert isinstance(page, Page)
        assert len(page.items) > 0
        assert all(isinstance(item, Species) for item in page.items)
    
    @pytest.mark.skipif(
        not os.environ.get("ATCT_API_BASE_URL"),
        reason="ATCT_API_BASE_URL not set - skipping integration tests"
    )
    def test_real_covariance(self):
        """Test real covariance retrieval."""
        cov = get_species_covariance_by_atctid("67-56-1*500", "67-56-1*0")
        assert isinstance(cov, Covariance2x2)
        assert len(cov.labels) == 2
        assert len(cov.matrix) == 2
        assert len(cov.matrix[0]) == 2

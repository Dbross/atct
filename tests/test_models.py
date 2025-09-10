"""Test the data models."""

import pytest
from atct.api.models import Species, CovarianceMatrix, Covariance2x2, Page


def test_species_from_dict():
    """Test Species creation from API response."""
    data = {
        "ATcT_ID": "67-56-1*0",
        "Name": "Methanol",
        "Formula": "CH4O (g)",
        "∆fH_298K": "-201.0",
        "∆fH_298K_uncertainty": "0.1",
        "SMILES": "CO",
        "CASRN": "67-56-1",
        "InChI": "InChI=1S/CH4O/c1-2/h2H,1H3",
        "InChI_Key": "OKKJLVBELUTLKV-UHFFFAOYSA-N"
    }
    
    species = Species.from_dict(data)
    
    assert species.atct_id == "67-56-1*0"
    assert species.name == "Methanol"
    assert species.formula == "CH4O (g)"
    assert species.delta_h_298k == "-201.0"
    assert species.delta_h_298k_uncertainty == "0.1"
    assert species.smiles == "CO"
    assert species.casrn == "67-56-1"
    assert species.inchi == "InChI=1S/CH4O/c1-2/h2H,1H3"
    assert species.inchi_key == "OKKJLVBELUTLKV-UHFFFAOYSA-N"


def test_species_from_dict_exact_uncertainty():
    """Test Species with exact uncertainty."""
    data = {
        "ATcT_ID": "67-56-1*0",
        "Name": "Methanol",
        "Formula": "CH4O (g)",
        "∆fH_298K": "-201.0",
        "∆fH_298K_uncertainty": "exact"
    }
    
    species = Species.from_dict(data)
    
    assert species.delta_h_298k == "-201.0"
    assert species.delta_h_298k_uncertainty == "exact"


def test_covariance2x2_from_dict():
    """Test Covariance2x2 creation from API response."""
    data = {
        "labels": ["ΔH(A)", "ΔH(B)"],
        "units": "kJ/mol",
        "matrix": [[0.01, 0.005], [0.005, 0.02]]
    }
    
    cov_matrix = Covariance2x2.from_dict(data)
    
    assert cov_matrix.labels == ["ΔH(A)", "ΔH(B)"]
    assert cov_matrix.units == "kJ/mol"
    assert cov_matrix.matrix == [[0.01, 0.005], [0.005, 0.02]]


def test_covariance2x2_defaults():
    """Test Covariance2x2 with default values."""
    data = {
        "matrix": [[0.01, 0.005], [0.005, 0.02]]
    }
    
    cov_matrix = Covariance2x2.from_dict(data)
    
    assert cov_matrix.labels == ["ΔH(A)", "ΔH(B)"]
    assert cov_matrix.units == ""
    assert cov_matrix.matrix == [[0.01, 0.005], [0.005, 0.02]]


def test_covariance_matrix_from_dict():
    """Test CovarianceMatrix creation from API response."""
    data = {
        "species_ids": ["1", "2"],
        "species_atctids": ["67-56-1*0", "67-56-1*500"],
        "matrix": [[0.01, 0.005], [0.005, 0.02]],
        "units": "kJ/mol"
    }
    
    cov_matrix = CovarianceMatrix.from_dict(data)
    
    assert cov_matrix.species_ids == ["1", "2"]
    assert cov_matrix.species_atctids == ["67-56-1*0", "67-56-1*500"]
    assert cov_matrix.units == "kJ/mol"
    assert cov_matrix.matrix == [[0.01, 0.005], [0.005, 0.02]]


def test_covariance_matrix_to_dict():
    """Test CovarianceMatrix serialization to dictionary."""
    cov_matrix = CovarianceMatrix(
        species_ids=["1", "2"],
        species_atctids=["67-56-1*0", "67-56-1*500"],
        matrix=[[0.01, 0.005], [0.005, 0.02]],
        units="kJ/mol"
    )
    
    data = cov_matrix.to_dict()
    
    assert data["species_ids"] == ["1", "2"]
    assert data["species_atctids"] == ["67-56-1*0", "67-56-1*500"]
    assert data["units"] == "kJ/mol"
    assert data["matrix"] == [[0.01, 0.005], [0.005, 0.02]]


def test_page_creation():
    """Test Page dataclass."""
    species1 = Species(
        atct_id="id1",
        name="name1",
        formula="formula1",
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
    species2 = Species(
        atct_id="id2",
        name="name2",
        formula="formula2",
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
    
    page = Page(items=[species1, species2], total=2, limit=50, offset=0)
    
    assert len(page.items) == 2
    assert page.total == 2
    assert page.limit == 50
    assert page.offset == 0
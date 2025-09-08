"""Test the data models."""

import pytest
from atct.api.models import Species, Covariance2x2, Page


def test_species_from_dict():
    """Test Species creation from API response."""
    data = {
        "ATcT_ID": "67-56-1*0",
        "Name": "Methanol",
        "Formula": "CH4O (g)",
        "Delta_Hf_298K": "-201.0",
        "Delta_Hf298K_uncertainty": "0.1",
        "SMILES": "CO",
        "mass": "32.04"
    }
    
    species = Species.from_dict(data)
    
    assert species.atct_id == "67-56-1*0"
    assert species.name == "Methanol"
    assert species.formula == "CH4O (g)"
    assert species.delta_h_298k == "-201.0"
    assert species.delta_h_298k_uncertainty == "0.1"
    assert species.smiles == "CO"
    assert species.mass == "32.04"


def test_species_from_dict_exact_uncertainty():
    """Test Species with exact uncertainty."""
    data = {
        "ATcT_ID": "67-56-1*0",
        "Name": "Methanol",
        "Formula": "CH4O (g)",
        "Delta_Hf_298K": "-201.0",
        "Delta_Hf298K_uncertainty": "exact"
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
        unit=None,
        mass=None,
        mass_uncertainty=None,
        smiles=None,
        casrn=None,
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
        unit=None,
        mass=None,
        mass_uncertainty=None,
        smiles=None,
        casrn=None,
        charge=None,
        xyz=None
    )
    
    page = Page(items=[species1, species2], total=2, limit=50, offset=0)
    
    assert len(page.items) == 2
    assert page.total == 2
    assert page.limit == 50
    assert page.offset == 0
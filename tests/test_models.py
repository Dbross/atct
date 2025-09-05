"""Test the data models."""

import pytest
from atct.api.models import Species, CovarianceMatrix, Uncertainty, Page


def test_species_from_dict():
    """Test Species creation from API response."""
    data = {
        "ATcT_ID": "67-56-1*0",
        "Name": "Methanol",
        "Formula": "CH4O (g)",
        "phase": "g",
        "Delta_Hf_298K": -201.0,
        "Delta_Hf298K_uncertainty": 0.1,
        "SMILES": "CO",
        "mass": 32.04
    }
    
    species = Species.from_dict(data)
    
    assert species.atct_id == "67-56-1*0"
    assert species.name == "Methanol"
    assert species.formula == "CH4O (g)"
    assert species.phase == "g"
    assert species.delta_h_298k == -201.0
    assert species.delta_h_298k_unc is not None
    assert species.delta_h_298k_unc.value == 0.1
    assert species.delta_h_298k_unc.unit == "kJ/mol"
    assert species.smiles == "CO"
    assert species.mass == 32.04


def test_species_from_dict_exact_uncertainty():
    """Test Species with exact uncertainty."""
    data = {
        "ATcT_ID": "67-56-1*0",
        "Name": "Methanol",
        "Formula": "CH4O (g)",
        "Delta_Hf_298K": -201.0,
        "Delta_Hf298K_uncertainty": "exact"
    }
    
    species = Species.from_dict(data)
    
    assert species.delta_h_298k == -201.0
    assert species.delta_h_298k_unc is None


def test_covariance_matrix_from_dict_v1():
    """Test CovarianceMatrix creation from v1 API response."""
    data = {
        "atctid1": "67-56-1*0",
        "atctid2": "7727-37-9*0",
        "matrix": [[0.01, 0.005], [0.005, 0.02]],
        "version": "1.123"
    }
    
    cov_matrix = CovarianceMatrix.from_dict(data)
    
    assert cov_matrix.atct_id1 == "67-56-1*0"
    assert cov_matrix.atct_id2 == "7727-37-9*0"
    assert cov_matrix.matrix == [[0.01, 0.005], [0.005, 0.02]]
    assert cov_matrix.version == "1.123"


def test_covariance_matrix_from_dict_legacy():
    """Test CovarianceMatrix creation from legacy API response."""
    data = {
        "atctid1": "67-56-1*0",
        "atctid2": "7727-37-9*0",
        "covariance": 0.005,
        "ATcT_TN_Version": "1.123"
    }
    
    cov_matrix = CovarianceMatrix.from_dict(data)
    
    assert cov_matrix.atct_id1 == "67-56-1*0"
    assert cov_matrix.atct_id2 == "7727-37-9*0"
    assert cov_matrix.matrix == [[0.0, 0.005], [0.005, 0.0]]
    assert cov_matrix.version == "1.123"


def test_uncertainty_creation():
    """Test Uncertainty dataclass."""
    unc = Uncertainty(value=0.1, unit="kJ/mol", method="ATcT")
    
    assert unc.value == 0.1
    assert unc.unit == "kJ/mol"
    assert unc.method == "ATcT"


def test_page_creation():
    """Test Page dataclass."""
    species1 = Species("id1", "name1", "formula1")
    species2 = Species("id2", "name2", "formula2")
    
    page = Page(items=[species1, species2], total=2, limit=50, offset=0)
    
    assert len(page.items) == 2
    assert page.total == 2
    assert page.limit == 50
    assert page.offset == 0

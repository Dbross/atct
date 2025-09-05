"""Public API for the ATcT package.

Usage:
    from atct import get_species, search_species, Species
"""
from .api import (
    get_species,
    search_species,
    get_covariance,
    get_species_by_smiles,
    get_species_by_casrn,
    healthcheck,
)
from .models import Species, CovarianceMatrix, Page, Uncertainty
from .exceptions import (
    ATCTError,
    NotFound,
    BadRequest,
    Unauthorized,
    ServerError,
    NetworkError,
)

__all__ = [
    "get_species",
    "search_species",
    "get_covariance",
    "get_species_by_smiles",
    "get_species_by_casrn",
    "healthcheck",
    "Species",
    "CovarianceMatrix",
    "Page",
    "Uncertainty",
    "ATCTError",
    "NotFound",
    "BadRequest",
    "Unauthorized",
    "ServerError",
    "NetworkError",
]

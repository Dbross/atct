"""ATcT: Active Thermochemical Tables - Python client for v1 API."""

from .api import (
    get_species,
    search_species,
    get_covariance,
    get_species_by_smiles,
    get_species_by_casrn,
    healthcheck,
)
from .api.models import Species, CovarianceMatrix, Page, Uncertainty
from .api.exceptions import ATCTError, NotFound, BadRequest, Unauthorized, ServerError, NetworkError

__version__ = "1.0.0"
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

def test_import_style():
    # Ensure the requested import style works
    import atct
    from atct import get_species, search_species, Species, get_species_covariance_by_atctid
    from atct import get_species_by_atctid, get_species_by_inchikey, get_species_covariance_matrix
    from atct import get_species_simple, get_all_species
    from atct import ATCTError, NotFound, BadRequest, Unauthorized, ServerError, NetworkError, TooManyRequests
    from atct.api import create_reaction_calculator, calculate_reaction_enthalpy
    from atct.api.models import ReactionSpecies, ReactionResult, ReactionCalculator, CovarianceMatrix
    
    # Test that we can create objects
    assert hasattr(atct, 'get_species')
    assert hasattr(atct, 'get_species_by_atctid')
    assert hasattr(atct, 'get_species_by_inchikey')
    assert hasattr(atct, 'get_species_covariance_matrix')
    assert hasattr(atct, 'get_species_simple')
    assert hasattr(atct, 'get_all_species')
    assert hasattr(atct, 'search_species')
    assert hasattr(atct, 'get_species_covariance_by_atctid')
    assert hasattr(atct, 'Species')
    assert hasattr(atct, 'CovarianceMatrix')
    assert hasattr(atct, 'Covariance2x2')
    assert hasattr(atct, 'ReactionSpecies')
    assert hasattr(atct, 'ReactionResult')
    assert hasattr(atct, 'ReactionCalculator')
    
    # Test that exceptions are available
    assert hasattr(atct, 'ATCTError')
    assert hasattr(atct, 'NotFound')
    assert hasattr(atct, 'TooManyRequests')

def test_import_style():
    # Ensure the requested import style works
    import atct
    from atct import get_species, search_species, Species, get_covariance
    from atct import ATCTError, NotFound, BadRequest, Unauthorized, ServerError, NetworkError
    
    # Test that we can create objects
    assert hasattr(atct, 'get_species')
    assert hasattr(atct, 'search_species')
    assert hasattr(atct, 'get_covariance')
    assert hasattr(atct, 'Species')
    assert hasattr(atct, 'CovarianceMatrix')
    
    # Test that exceptions are available
    assert hasattr(atct, 'ATCTError')
    assert hasattr(atct, 'NotFound')

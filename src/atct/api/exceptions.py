class ATCTError(Exception):
    """Base error for the ATcT library."""

class NotFound(ATCTError):
    pass

class BadRequest(ATCTError):
    pass

class Unauthorized(ATCTError):
    pass

class ServerError(ATCTError):
    pass

class NetworkError(ATCTError):
    pass

class TooManyRequests(ATCTError):
    pass

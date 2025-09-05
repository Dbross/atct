from __future__ import annotations
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Dict, Any, Optional
from ._config import settings
from .exceptions import ATCTError, BadRequest, NotFound, Unauthorized, ServerError, NetworkError

def _raise_for_status(response: urllib.request.addinfourl, content: bytes) -> None:
    """Check HTTP status and raise appropriate exception."""
    status_code = response.getcode()
    if 200 <= status_code < 300:
        return
    if status_code == 400:
        raise BadRequest(content.decode('utf-8', errors='replace'))
    if status_code in (401, 403):
        raise Unauthorized(content.decode('utf-8', errors='replace'))
    if status_code == 404:
        raise NotFound(content.decode('utf-8', errors='replace'))
    if 500 <= status_code < 600:
        raise ServerError(content.decode('utf-8', errors='replace'))
    raise ATCTError(f"Unexpected status {status_code}: {content.decode('utf-8', errors='replace')}")

def _headers() -> Dict[str, str]:
    """Get default headers for requests."""
    h = {"User-Agent": settings.user_agent}
    if settings.api_key:
        h["Authorization"] = f"Bearer {settings.api_key}"
    return h

def _build_url(url: str, params: Optional[Dict[str, Any]] = None) -> str:
    """Build URL with query parameters."""
    if not params:
        return url
    query_string = urllib.parse.urlencode(params, doseq=True)
    separator = "&" if "?" in url else "?"
    return f"{url}{separator}{query_string}"

def request_json(method: str, url: str, *, params: Optional[Dict[str, Any]] = None, 
                json_data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Make HTTP request and return JSON response with retries."""
    delay = 0.5
    headers = _headers()
    
    # Add JSON content type if sending JSON data
    if json_data is not None:
        headers["Content-Type"] = "application/json"
    
    for attempt in range(settings.retries + 1):
        try:
            # Build URL with parameters
            full_url = _build_url(url, params)
            
            # Prepare request data
            data = None
            if json_data is not None:
                data = json.dumps(json_data).encode('utf-8')
            
            # Create request
            req = urllib.request.Request(full_url, data=data, headers=headers, method=method)
            
            # Make request with timeout
            with urllib.request.urlopen(req, timeout=settings.timeout_s) as response:
                content = response.read()
                _raise_for_status(response, content)
                return json.loads(content.decode('utf-8'))
                
        except urllib.error.HTTPError as e:
            # Handle HTTP errors
            if attempt >= settings.retries:
                content = e.read() if hasattr(e, 'read') else b''
                _raise_for_status(e, content)
            time.sleep(delay)
            delay *= 2
        except (urllib.error.URLError, OSError) as e:
            # Handle network errors
            if attempt >= settings.retries:
                raise NetworkError(str(e)) from e
            time.sleep(delay)
            delay *= 2
        except json.JSONDecodeError as e:
            # Handle JSON decode errors
            if attempt >= settings.retries:
                raise ATCTError(f"Invalid JSON response: {e}") from e
            time.sleep(delay)
            delay *= 2

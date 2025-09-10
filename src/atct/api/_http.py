from __future__ import annotations
import asyncio
import json
import time
from typing import Dict, Any, Optional, Union
import httpx
from ._config import settings
from .exceptions import ATCTError, BadRequest, NotFound, Unauthorized, ServerError, NetworkError

def _raise_for_status(response: httpx.Response) -> None:
    """Check HTTP status and raise appropriate exception."""
    status_code = response.status_code
    if 200 <= status_code < 300:
        return
    if status_code == 400:
        raise BadRequest(response.text)
    if status_code in (401, 403):
        raise Unauthorized(response.text)
    if status_code == 404:
        raise NotFound(response.text)
    if 500 <= status_code < 600:
        raise ServerError(response.text)
    raise ATCTError(f"Unexpected status {status_code}: {response.text}")

def _headers() -> Dict[str, str]:
    """Get default headers for requests."""
    h = {"User-Agent": settings.user_agent}
    if settings.api_key:
        h["Authorization"] = f"Bearer {settings.api_key}"
    return h

def _get_client() -> httpx.AsyncClient:
    """Get configured httpx client."""
    return httpx.AsyncClient(
        headers=_headers(),
        timeout=settings.timeout_s,
        follow_redirects=True
    )

async def request_json_async(method: str, url: str, *, params: Optional[Dict[str, Any]] = None, 
                           json_data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Make async HTTP request and return JSON response with retries."""
    delay = 0.5
    
    for attempt in range(settings.retries + 1):
        try:
            async with _get_client() as client:
                response = await client.request(
                    method=method,
                    url=url,
                    params=params,
                    json=json_data
                )
                _raise_for_status(response)
                return response.json()
                
        except httpx.HTTPStatusError as e:
            # Handle HTTP errors
            if attempt >= settings.retries:
                _raise_for_status(e.response)
            await asyncio.sleep(delay)
            delay *= 2
        except (httpx.ConnectError, httpx.TimeoutException, httpx.RequestError) as e:
            # Handle network errors
            if attempt >= settings.retries:
                raise NetworkError(str(e)) from e
            await asyncio.sleep(delay)
            delay *= 2
        except json.JSONDecodeError as e:
            # Handle JSON decode errors
            if attempt >= settings.retries:
                raise ATCTError(f"Invalid JSON response: {e}") from e
            await asyncio.sleep(delay)
            delay *= 2

def request_json(method: str, url: str, *, params: Optional[Dict[str, Any]] = None, 
                json_data: Optional[Dict[str, Any]] = None, block: bool = False) -> Union[Dict[str, Any], asyncio.Task]:
    """Make HTTP request and return JSON response.
    
    Args:
        method: HTTP method
        url: Request URL
        params: Query parameters
        json_data: JSON data to send
        block: If True, blocks and returns result. If False, returns async task.
    
    Returns:
        If block=True: JSON response data
        If block=False: asyncio.Task that resolves to JSON response data
    """
    task = request_json_async(method, url, params=params, json_data=json_data)
    if block:
        return asyncio.run(task)
    return task

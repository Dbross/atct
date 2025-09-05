from __future__ import annotations
import os
from dataclasses import dataclass

@dataclass(frozen=True)
class Settings:
    base_url: str = os.environ.get("ATCT_API_BASE_URL", "https://atct.anl.gov/api/v1")
    api_key: str | None = os.environ.get("ATCT_API_KEY")
    timeout_s: float = float(os.environ.get("ATCT_TIMEOUT_S", "10"))
    retries: int = int(os.environ.get("ATCT_RETRIES", "3"))
    user_agent: str = os.environ.get("ATCT_USER_AGENT", "atct/1.0.0")

settings = Settings()

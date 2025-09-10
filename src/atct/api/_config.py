from __future__ import annotations
import os
from dataclasses import dataclass
from urllib.parse import urljoin

@dataclass(frozen=True)
class Settings:
    # API endpoint
    base_url: str = os.environ.get("ATCT_API_BASE_URL", "https://atct.anl.gov/api/v1/")
    
    # Client configuration
    api_key: str | None = os.environ.get("ATCT_API_KEY")
    timeout_s: float = float(os.environ.get("ATCT_TIMEOUT_S", "10"))
    retries: int = int(os.environ.get("ATCT_RETRIES", "3"))
    user_agent: str = os.environ.get("ATCT_USER_AGENT", "Mozilla/5.0 (compatible; ATcT-Python/1.0.0)")
    
    def url(self, path: str) -> str:
        """Join base_url with a path, handling trailing slashes properly."""
        return urljoin(self.base_url.rstrip('/') + '/', path.lstrip('/'))

settings = Settings()

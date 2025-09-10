from __future__ import annotations
from typing import Iterable, Any
try:
    import pandas as pd  # type: ignore
except Exception:
    pd = None

def as_dataframe(items: Iterable[Any]):
    if pd is None:
        raise ImportError("Install with `atct[pandas]` for DataFrame support.")
    rows = []
    for x in items:
        if hasattr(x, "to_dict"):
            rows.append(x.to_dict())
        elif isinstance(x, dict):
            rows.append(x)
        else:
            rows.append(vars(x))
    return pd.DataFrame(rows)

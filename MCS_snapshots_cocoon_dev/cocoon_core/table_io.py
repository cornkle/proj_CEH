"""Safe development-output functions."""
from pathlib import Path

def save_parquet(df, path, overwrite=False):
    p = Path(path)
    if p.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite {p}")
    if 'MCS_5000km2_tables_v3' in str(p) and '_DEV' not in str(p):
        raise ValueError("Refusing write to observational production directory")
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(p, index=False)
    return p

"""Compare parquet frames with strict key matching and diagnostics."""
import numpy as np
import pandas as pd

def compare_frames(new, old, keys, atol=1e-8, rtol=1e-8):
    if len(new)==0: raise AssertionError('No reproduced storm rows')
    for key in keys:
        if key not in new or key not in old:
            raise AssertionError(f'Missing key {key}')
    if new.duplicated(keys).any() or old.duplicated(keys).any():
        raise AssertionError('Duplicate comparison keys')
    shared_cols=[c for c in new.columns if c in old.columns and c not in keys]
    x=new.set_index(keys).sort_index()
    old_indexed=old.set_index(keys)
    if not x.index.isin(old_indexed.index).all():
        raise AssertionError('Some reproduced keys absent from reference')
    y=old_indexed.reindex(x.index)
    failures=[]
    for c in shared_cols:
        a=x[c]; b=y[c]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            valid=np.isclose(a.to_numpy(dtype=float),b.to_numpy(dtype=float),atol=atol,rtol=rtol,equal_nan=True)
        else:
            valid=(a.fillna('<NA>').astype(str).to_numpy()==b.fillna('<NA>').astype(str).to_numpy())
        if not valid.all(): failures.append((c,int((~valid).sum())))
    if failures: raise AssertionError(f'Mismatched columns (counts): {failures[:25]}')
    return len(new),len(shared_cols)

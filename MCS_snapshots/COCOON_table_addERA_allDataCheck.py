import pandas as pd
import glob
import os
from shared.paths import load_paths

paths = load_paths()
ROOT = paths.datasets["obs_mcs_5000km2_tables_ERA"]

records = []

for f in glob.glob(
    f"{ROOT}/*/region=*/year=*/part.parquet"
):
    print('Doing', f)
    try:
        df = pd.read_parquet(f)

        records.append({
            "file": f,
            "nrows": len(df),
            "nan_fraction": df["value"].isna().mean(),
        })

    except Exception:

        records.append({
            "file": f,
            "nrows": -1,
            "nan_fraction": np.nan,
        })

audit = pd.DataFrame(records)

print(audit[audit["nan_fraction"] == 1])
print(audit.sort_values("nan_fraction", ascending=False).head(50))
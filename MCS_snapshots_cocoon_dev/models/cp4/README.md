# COCOON CP4-Africa base-table adapter

Self-contained standalone script for converting HadREM-CP4A 4.5-km `rlut`
(OLR, W m-2) and `pr` (kg m-2 s-1) into *observational-style* hourly
MCS base-table Parquets. No changes or writes to the observational COCOON
production scripts, datasets, or ERA5 Parquets.

## Default input and output locations

**Input**: `/gws/pw/j25/cp4a/<GCM>-<HIST|FUT>/{rlut,pr}/*.nc`

**Output**: `/gws/ssde/j25b/cocoon/MCS_5000km2_tables_v3/CP4_base_tables_tir-prcp/<RUN>/year_<YYYY>/month_<MM>.parquet`

For example:
`.../CP4_base_tables_tir-prcp/HadGEM3-GC31-MM-FUT/year_2080/month_08.parquet`

The output covers the **whole CP4 domain** (not clipped to WAf, EAf, etc.),
so region and local-time filtering remain analysis decisions. Monthly Parquets
avoid a year-long wait for an initial comparison. The file schema preserves
all OBS tabular storm and rainfall columns and adds model/run metadata. It
intentionally **does not add compactness**.

## Install/use within your existing Python workspace

You can copy this whole directory into `MCS_snapshots_cocoon_dev/models/cp4/`
(or run it directly from here) without changing the original `MCS_snapshots`
package. Dependencies: `numpy pandas xarray scipy scikit-image shapely pyproj
cftime netCDF4 pyarrow` (or replace `netCDF4` with `h5netcdf`).

Run from a terminal in the directory containing `build_cp4_base_tables.py`:

```bash
# Check which files will be used for the earliest available August year
# in each of the five HIST and five FUT runs; no output written.
python build_cp4_base_tables.py --month 8 --dry-run

# Test first 24 hours of one month in just one run.
# These files live under CP4_base_tables_tir-prcp/_smoke_tests/,
# so they cannot be confused with complete monthly files.
python build_cp4_base_tables.py --run HadGEM3-GC31-MM-FUT --month 8 --max-hours 24

# Produce a complete August for one run.
python build_cp4_base_tables.py --run HadGEM3-GC31-MM-FUT --month 8

# Produce one whole August month for all 10 available runs.
python build_cp4_base_tables.py --month 8

# Choose a specific year, e.g. HadGEM3 future August 2089.
python build_cp4_base_tables.py --run HadGEM3-GC31-MM-FUT --year 2089 --month 8

# Overwrite a completed output after changing code/thresholds.
python build_cp4_base_tables.py --run HadGEM3-GC31-MM-FUT --year 2089 --month 8 --overwrite
```

When `--year` is not given, the **first year with a uniquely matching pair**
of `rlut`/`pr` files for the requested month is used independently in each run.
The source filename need not contain `-FUT` or `-HIST`; the parent simulation
directory identifies the run. Missing runs are reported and skipped. A failed
run does not stop the next, but the script exits nonzero after any failures.

## Exact storm definitions from observational COCOON

- OLR -> brightness temperature through the empirical Yang/Slingo relation
  used in the existing legacy CP4 scripts; cloud threshold in °C.
- 4-neighbour labelled cold cloud at Tb < -40°C.
- Physical area 5,000–300,000 km²; Tmin <= -50°C.
- Maximum **raw** precipitation under the cold cloud >= 1 mm/h.
- `pr × 3600` -> rainfall rate in mm/h.
- `area`, `area50`, `area60`, `area70`, cloud temperature statistics,
  minimum-temperature location, bounding box, `shape_complexity`,
  `major_axis_km`, `minor_axis_km`, `aspect_ratio`.
- All observational rainfall columns: masked mean, 80/90/95/99th percentiles,
  area >/= 1/3/8 mm/h, and observationally defined `precipitation_max`:
  mean within a +/-0.075° lat/lon window around the maximum-rain grid cell.
- `cloudMask` and full `tir` arrays are not saved, matching the final OBS
  Parquet schema.

Geometry contours are computed on padded object subarrays for efficiency;
objects touching the outer model-domain boundary can have slightly different
contour handling than the full-domain OBS function. A one-storm synthetic
parity check of the regular case is provided under `tests/`.

## Timestamp convention

Model OLR and precipitation are hourly means stamped at `hh:30`. For example
17:30 represents the hour ending at 18:00 *if* the `time_bnds` for the source
confirm the expected centered one-hour interval. The Parquet therefore
retains OBS-compatible `date`, `hour` and `minute` fields for the source center
time and adds `storm_hour_end_utc` for identifying the intended 18 UTC storms.
The time-bounds convention should be checked once on the actual inputs.

Each row also has `storm_id` (globally unique by model run, time and label),
`storm_time`, `calendar`, `local_solar_hour`, `pixel_count`, and source file
names. The original OBS-like `stormID` remains present.

No regional filter, solar-time filter, shape-complexity filter, resampling or
coarsening is applied during detection. All-hour MCS catalogues can later feed
the native 12/18/00 environmental sampling without another TIR/pr pass.

## Read outputs

```python
from pathlib import Path
import pandas as pd

root = Path('/gws/ssde/j25b/cocoon/MCS_5000km2_tables_v3/CP4_base_tables_tir-prcp')
run_dir = root / 'HadGEM3-GC31-MM-FUT'
# One month
mcs = pd.read_parquet(run_dir / 'year_2080' / 'month_08.parquet')
# 18 UTC hour-ending storms (timestamps centered at 17:30)
mcs18 = mcs.loc[mcs.storm_hour_end_utc == 18]
print(len(mcs18), mcs18[['area', 'precipitation_max', 'shape_complexity']].describe())
```

## Validation performed

The core `extract_storms()` routine has been compared in this environment
against the supplied `COCOON_table_create_01deg.py` on a synthetic storm, and
all OBS metric columns matched within floating-point tolerance. **It has not
been run against the JASMIN `/gws/pw/j25/cp4a` data in this environment.**
Do the 24-hour smoke test and inspect a few storm distributions/images before
launching the full ensemble.

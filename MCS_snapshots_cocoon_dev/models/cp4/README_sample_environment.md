# CP4 environment sampling — mandatory 0.7° spatial footprint (v3)

Copy `sample_cp4_environments.py`, `submit_cp4_august_env.sbatch`, and `tests/test_env_sampling.py` to `MCS_snapshots_cocoon_dev/models/cp4/` (with the test under `tests/`). This does **not** change the observational scripts, Parquets, or CP4 base tables.

## Critical scientific choices

- OBS production samples the **nearest ERA5 0.7° grid point** at `(tminlat, tminlon)`. This CP4 sampler therefore averages its ~4.5 km native values in a **mandatory 0.7° latitude × 0.7° longitude window** (±0.35° from the storm's Tmin coordinate) to produce a comparable environmental footprint. These are matched *spatial scales*, not numerically identical methods of coarse gridding.
- No native 4.5 km nearest-pixel sampling option exists in this script. Results record `spatial_sampling=box_mean_0p7deg`, `sampling_box_radius_deg=0.35`, and `sampling_box_width_deg=0.7`.
- Only native 6-hourly synoptic atmosphere timestamps are sampled: for an 18 UTC storm, `rel_hour=-6,0,+6` corresponds to **12:00, 18:00 and next-day 00:00 UTC**.
- On hourly *time-mean* surface fields, when no exact synoptic timestamp exists but a mean stamped at the preceding :30 does, that mean is used **only** if explicitly marked `time: mean`. If `time_bnds` are present, the end bound must equal the requested synoptic time. This is not interpolation; all other time mismatches are reported and skipped.
- Wind components are averaged over the box first and then used to form shear magnitude. Different native variable grids (including the extra 2101st latitude for pressure-level winds) are sampled independently.
- A 0.7° box which extends outside the model grid is `NaN` rather than silently shrunk to native-pixel sampling.
- The default output is deliberately **separate** from any older native-pixel results: `/gws/ssde/j25b/cocoon/MCS_5000km2_tables_v3/CP4_storm_timeseries_box07/`.

## Start with one run, August

From your `models/cp4/` directory:

```bash
python -m unittest discover -s tests -p 'test_env_sampling.py' -v

python sample_cp4_environments.py \
  --run CanESM5-HIST --month 8 --all-years --dry-run

# Creates ONLY smoke-test parts, not production parts:
python sample_cp4_environments.py \
  --run CanESM5-HIST --month 8 --year 1998 --max-storms 30

# Complete finished August base tables for one run:
python sample_cp4_environments.py \
  --run CanESM5-HIST --month 8 --all-years
```

The current stage uses `--storm-hours 18 --rel-hours -6 0 6 --profile core` by default, and restricts storm centres to the union of the existing African OBS boxes: WAf, CAf, EAf and SAf.

## SLURM parallel run

Copy the SLURM `--account=...` line from your working `submit_cp4_full.sbatch`. The account associated with the storage path is **not necessarily** the job account. Also check partition/QoS. From the directory holding the script, with `logs/` created:

```bash
mkdir -p logs
sbatch submit_cp4_august_env.sbatch
```

The job array has ten independent run tasks, at most two concurrent, uses `SLURM_SUBMIT_DIR`, and operates only on completed monthly base Parquets. Reading shared NetCDF data is safe; parts for different runs never overlap.

## Rerun, fill months, or add locations later

Existing monthly, run, variable, storm-hour and scope-specific parts are skipped. Failed/incomplete variable parts are not written and can be retried.

```bash
# When more monthly CP4 base tables appear:
python sample_cp4_environments.py --month 8 --all-years

# After all base tables finish:
python sample_cp4_environments.py --all

# Sample ONLY storms outside African OBS boxes without repeating inside storms:
python sample_cp4_environments.py --all --outside-boxes

# Or process both scopes, skipping completed in-box parts:
python sample_cp4_environments.py --all --full-domain
```

In-box part: `<output>/<RUN>/<VARIABLE>/year_YYYY/month_MM/part-18utc.parquet`.
Outside part: `.../part-18utc-outside.parquet`. Both represent the **same logical per-variable dataset** with different storm IDs; the region list and storm ID are recorded in every row.

If you later request further storm hours or relative-hour sets, output filenames distinguish storm hour and offsets; no need to change or overwrite the first production.

## About `tests/test_env_sampling.py`

This **does not produce any storm or environment data**. It is a fast collection of synthetic regression tests verifying mandatory 0.7° averaging (not native-pixel extraction), protection against incomplete spatial footprints, 0–360 longitude conventions, disjoint box/outside parts, vector shear, 360-day month rollover and timed sampling from staggered model grids. It creates short-lived synthetic NetCDF files in temporary folders. It never reads or writes your JASMIN production input/output trees. Some tests require `cftime`/`pyarrow` and run on JASMIN's Python environment.

## Validation status

The local non-`cftime` tests run successfully. A real CP4 NetCDF + Parquet smoke test on JASMIN is still required. Inspect its `SAMPLED` logs and final `spatial_sampling` metadata before large production.

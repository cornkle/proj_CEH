# COCOON development workspace (isolated)

The uploaded observation production files are stored **unchanged** in `reference/obs_production/`. No new scripts modify `MCS_snapshots/` or write into the production OBS or ERA5 output roots. New model readers will go under `models/adapters/` and `models/configs/`.

## Install on JASMIN / in the existing Python environment

1. Extract `MCS_snapshots_cocoon_dev.zip` next to `MCS_snapshots/` (or another location).
2. Activate the environment you use to run COCOON, with the same `MCS_snapshots`, `GLOBAL`, `utils`, `pyarrow`, etc. available.
3. `cd MCS_snapshots_cocoon_dev && bash install.sh`

`pip install --no-deps -e .` makes the new modules importable from that environment. The editable install does not install or upgrade scientific dependencies. The original production package is never patched.

## Tests

Start in the new root:

```bash
python -m tests.test_obs_base --region GPlains --year 2008 --max-files 25
python -m tests.test_era5 --region GPlains --year 2008 --storms 2 --variables tcwv q650 ushear100_650
```

Or `bash tests/run_checks.sh GPlains 2008`.

`test_obs_base` recreates the **same observation base-table calculations** from a few input scenes and compares matching `(stormID, date)` rows to the existing production base parquet. The sample is not a whole year; if no storms occur early in the year, raise `--max-files`.

`test_era5` uses the **frozen ERA5 sampling functions**, selects a small number of storms after reproducing static and land/sea filtering, writes test parquets, and compares rows keyed by `(storm_id, rel_hour)` against the existing `ERA5_storm_timeseries` files. It retains the original ±hour window and current variable definitions. Requires JASMIN/ERA5 data and original project imports. Exact mismatches may indicate that your live production version differs from the uploaded snapshot or that reference outputs were generated with a different configuration.

All outputs default to `outputs/tests/`. Tests abort on errors; a PASS is printed only after comparisons succeed. The ERA sampler may report missing input data and use NaNs, which will count as mismatches if reference had data.

The scripts are not run against the JASMIN `/prj` filesystem in this environment; the tests must be executed where those files exist.

## Full-year obs reproduction (optional, DEV output only)

```bash
python -m scripts.build_obs_base_tables_dev --region GPlains --year 2008
```

The archived driver, patched in memory to use the archived detector, runs the frozen algorithm and writes to `outputs/obs_base`, never to production. It can be expensive, unlike the targeted tests.

## Model additions

After discussing the model data format, put its adapter under `models/adapters/`, run-specific configuration under `models/configs/`, and reuse the common detection/rainfall processing functions. The legacy model scripts are archived as references only; do not execute them for production.

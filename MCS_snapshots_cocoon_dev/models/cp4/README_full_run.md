# CP4 COCOON: complete Africa-domain storm catalogue

This is v2 of the base-table script. It **does not change any COCOON storm
threshold, rain metric, shape metric, or Parquet column** relative to v1.
It adds a full-production loop and retains original single-month behaviour.

Files:

- `build_cp4_base_tables.py` — replace the v1 script under
  `MCS_snapshots_cocoon_dev/models/cp4/` with this script.
- `submit_cp4_full.sbatch` — save alongside it.

Output stays under:

`/gws/ssde/j25b/cocoon/MCS_5000km2_tables_v3/CP4_base_tables_tir-prcp/`

Input stays under `/gws/pw/j25/cp4a/`.

## CLI examples

```bash
# Previous behaviour remains unchanged: first available August year per run
python build_cp4_base_tables.py --month 8

# All available August years for all ten runs
python build_cp4_base_tables.py --month 8 --all-years

# Full production (ALL available months and years across all ten runs)
python build_cp4_base_tables.py --all

# One simulation only, all available months and years
python build_cp4_base_tables.py --run CanESM5-HIST --all

# Dry run: inspect all inputs and planned outputs without writing
python build_cp4_base_tables.py --run CanESM5-HIST --all --dry-run
```

`--all` processes **every paired** `rlut` + `pr` year/month in each run's
input directories. August is deliberately processed first, across all years.
An existing completed monthly output is skipped unless `--overwrite` is used.
If a monthly processing error occurs, the script logs it, processes subsequent
months, and exits nonzero at the end.

## Recommended LOTUS job array

From the `models/cp4/` folder, create the log folder before submitting:

```bash
mkdir -p logs
sbatch submit_cp4_full.sbatch
squeue --me
```

The array includes 10 runs, capped at 3 simultaneous jobs (`--array=0-9%3`).
Each job works on a different run and different output subdirectory, so there
are no overlapping writes. All jobs read source model NetCDFs read-only.

**Before submission:** confirm that `j25b` is a valid Slurm account for you;
otherwise change `#SBATCH --account=...`. The wrapper assumes the notebook's
Miniforge environment at `~/miniforge3/envs/work/bin/python`; edit `PYTHON`
if your environment differs. LOTUS's `standard` partition/QoS supports up to
24 hours per job, as documented on JASMIN's help site.

Logs are `logs/cp4_<jobid>_<arrayindex>.out` and `.err`. Use:

```bash
squeue --me
# To monitor completed/failed array tasks, substitute the returned job ID:
sacct -j JOBID --format=JobID,JobName,State,Elapsed,MaxRSS
```

If a task times out or fails, simply resubmit the same array. Completed
months remain intact and are skipped on the next attempt. Avoid launching two
instances of the **same run** simultaneously, because both would write the
same monthly output files.

## Compatibility notes

- Your observational Parquet outputs and OBS scripts are not touched.
- The MCS detection/metrics computation in v2 is byte-for-byte unchanged from v1.
- The new `--all` mode is strictly orchestration over existing `process_month`.
- All raw hourly `rlut` and `pr` objects are sampled, for the full model domain.
- No environmental sampling is performed by this script.
- If model input is incomplete (missing `rlut` or `pr` for a month), that month
  is not processed. Inspect the log and input inventory for completeness.

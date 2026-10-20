#!/usr/bin/env python3
"""Incremental CP4 environmental sampling for COCOON storm snapshots.

Reads ONLY completed CP4 monthly base Parquets. Independently samples each
model fields over a mandatory 0.7 x 0.7 degree box centred on the
coldest-cloud pixel; creates OBS-like per-variable Parquets (storm_id, rel_hour, value).

Initial pass: --month 8 --all-years --storm-hours 18 --rel-hours -6 0 6.
Rerunning skips complete per-variable/month/storm-hour parts.

No existing OBS scripts/outputs or CP4 base tables are modified.
"""
from __future__ import annotations

import argparse
from datetime import timedelta
from functools import lru_cache
import logging
import os
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
import xarray as xr

LOG = logging.getLogger('cp4_env')
BASE_ROOT = Path('/gws/ssde/j25b/cocoon/MCS_5000km2_tables_v3')
INPUT_ROOT = Path('/gws/pw/j25/cp4a')
BASE_TABLES = BASE_ROOT / 'CP4_base_tables_tir-prcp'
OUTPUT_ROOT = BASE_ROOT / 'CP4_storm_timeseries_box07'
# This is deliberately separate from old native-point CP4_storm_timeseries.
BOX_RADIUS_DEG = 0.35   # mandatory half-width (0.7 x 0.7 degrees)
MODELS = ('CanESM5', 'HadGEM3-GC31-MM', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'TaiESM1')

# The exact African OBS comparison boxes used in the plotting notebook.
# Boxes (west, east, south, north); no duplicate sampling of overlaps.
AFRICA_BOXES = {
    'WAf': (-18., 25., 4., 25.),
    'SAf': (17., 33., -32., -14.),
    'CAf': (13.5, 29., -14., 4.),
    'EAf': (29., 44., -2.5, 15.5),
}

# Input native variable -> output field name (a single output per input).
# These are raw values in the NetCDF units, not interpolated vertically.
DIRECT_INPUTS = {
    'huss': 'q_sfc', 'hus850': 'q850', 'hus700': 'q700',
    'hus600': 'q600', 'hus500': 'q500', 'hus200': 'q200',
    'tas': 't2m', 'ta850': 't850', 'ta700': 't700',
    'ta600': 't600', 'ta500': 't500', 'ta200': 't200',
    'uas': 'u_sfc', 'vas': 'v_sfc',
    'ua850': 'u850', 'va850': 'v850',
    'ua700': 'u700', 'va700': 'v700',
    'ua600': 'u600', 'va600': 'v600',
    'ua500': 'u500', 'va500': 'v500',
    'ua200': 'u200', 'va200': 'v200',
    'mlcape': 'cape', 'mlcin': 'cin',
    'ps': 'p_srfc', 'hurs': 'rh_sfc',
    'wa850': 'w850', 'wa700': 'w700', 'wa500': 'w500',
    'zg850': 'zg850', 'zg500': 'zg500', 'zg200': 'zg200',
}

# Same sample coordinates on independent native grids; shear derived after
# extracting both components (no grid broadcasting/interpolation).
SHEAR_LAYERS = {
    'sfc_600': ('uas', 'vas', 'ua600', 'va600'),
    '850_600': ('ua850', 'va850', 'ua600', 'va600'),
}
SHEAR_OUTPUTS = {
    'sfc_600': ('ushear_sfc_600', 'vshear_sfc_600', 'shear_sfc_600'),
    '850_600': ('ushear850_600', 'vshear850_600', 'shear850_600'),
}

CORE_OUTPUTS = {
    'q_sfc', 'q850', 'q700', 'q600', 't2m', 't850', 't600',
    'u_sfc', 'v_sfc', 'u850', 'v850', 'u600', 'v600',
    'cape', 'cin',
    *SHEAR_OUTPUTS['sfc_600'],
    *SHEAR_OUTPUTS['850_600'],
}
ALL_OUTPUTS = set(DIRECT_INPUTS.values()) | {
    v for names in SHEAR_OUTPUTS.values() for v in names
}

# Matches both names with and without '-FUT' etc.
NC_RANGE = re.compile(r'_(\d{12})-(\d{12})\.nc$')
BASE_FILE = re.compile(r'year_(\d{4})/month_(\d{2})\.parquet$')


def native_key(t):
    """Calendar-independent, sortable YYYYMMDDHHMM key."""
    if isinstance(t, np.datetime64):
        t = pd.Timestamp(t)
    return (f'{int(t.year):04d}{int(t.month):02d}{int(t.day):02d}'
            f'{int(t.hour):02d}{int(t.minute):02d}')


def pretty_key(k):
    return f'{k[:4]}-{k[4:6]}-{k[6:8]} {k[8:10]}:{k[10:12]}:00'


def half_hour_before_key(key, calendar):
    """Previous :30 mean timestamp, including month/year in CF calendars."""
    import cftime
    t = cftime.datetime(int(key[:4]), int(key[4:6]), int(key[6:8]),
                        int(key[8:10]), int(key[10:12]),
                        calendar=str(calendar))
    return native_key(t - timedelta(minutes=30))


def add_hours(y, m, d, h, minute, calendar, hours):
    """Use model-specific CF calendar, never convert to Gregorian pandas."""
    import cftime
    t = cftime.datetime(int(y), int(m), int(d), int(h), int(minute),
                        calendar=str(calendar))
    return t + timedelta(hours=int(hours))


def locate_regions(lon, lat):
    """Return a list of region labels for the storm centre, no duplicates."""
    lon = float((float(lon) + 180) % 360 - 180)
    lat = float(lat)
    return [name for name, (w, e, s, n) in AFRICA_BOXES.items()
            if w <= lon <= e and s <= lat <= n]


def _round_nominal_storm_time(row):
    """One-hour average at :30 is represented by the hour ending at hh+1."""
    # CP4 base tables include true calendar, hour and minute.
    return add_hours(row.year, row.month, row.day, row.hour, row.minute,
                     row.calendar, 0) + timedelta(minutes=30)


def make_requests(storms, hours, rel_hours, restrict_boxes=True):
    needed = {'run', 'storm_id', 'year', 'month', 'day', 'hour', 'minute',
              'calendar', 'clon', 'clat', 'tminlat', 'tminlon',
              'storm_hour_end_utc'}
    missing = sorted(needed - set(storms.columns))
    if missing:
        raise KeyError('Base parquet missing: ' + ', '.join(missing))

    # The parquet may store year/month integer, but the storm hour is :30.
    hour_values = np.mod(storms['storm_hour_end_utc'].to_numpy(dtype=float), 24)
    use = storms.loc[np.isclose(hour_values[:,None], np.asarray(hours)[None,:]).any(axis=1)].copy()
    # Note: pass specific hours, e.g. [18]; passing [0] catches 23:30 storms.
    rows = []
    for s in use.itertuples(index=False):
        regions = locate_regions(s.clon, s.clat)
        if restrict_boxes and not regions:
            continue
        anchor = _round_nominal_storm_time(s)
        for rel in rel_hours:
            etime = anchor + timedelta(hours=int(rel))
            rows.append({
                'run': s.run, 'storm_id': s.storm_id,
                'region': regions[0] if regions else 'Outside',
                'regions': ','.join(regions) if regions else '',
                'storm_utc_time': pretty_key(native_key(anchor)),
                'storm_snapshot_time': pretty_key(native_key(add_hours(
                    s.year, s.month, s.day, s.hour, s.minute, s.calendar, 0))),
                'era_utc_time': pretty_key(native_key(etime)),
                'env_time_key': native_key(etime),
                'storm_utc_hour': int(anchor.hour),
                'storm_snapshot_hour_utc': float(s.hour) + float(s.minute)/60,
                'era_utc_hour': int(etime.hour),
                'rel_hour': int(rel),
                'calendar': s.calendar,
                'year': int(s.year), 'month': int(s.month), 'day': int(s.day),
                'hour': int(s.hour), 'minute': int(s.minute),
                'tminlat': float(s.tminlat), 'tminlon': float(s.tminlon),
                'clat': float(s.clat), 'clon': float(s.clon),
                'storm_local_solar_hour': (int(anchor.hour)+float(s.tminlon)/15.) % 24,
            })
    df = pd.DataFrame(rows)
    if not df.empty:
        if df.duplicated(['storm_id', 'rel_hour']).any():
            raise ValueError('Duplicate storm_id/rel_hour in requests')
    return df


@lru_cache(maxsize=512)
def source_files(input_root, run, var):
    p = Path(input_root) / run / var
    out = []
    for f in sorted(p.glob(f'{var}_*.nc')):
        mat = NC_RANGE.search(f.name)
        if mat:
            out.append((mat.group(1), mat.group(2), f))
    return tuple(out)


def file_for_time(files, key):
    choices = [f for start, end, f in files if start <= key <= end]
    # Normally only one; if a boundary occurs in two files, prefer last start.
    return choices[-1] if choices else None


def _box_index_ranges(axis, targets, halfwidth):
    """Inclusive coordinate ranges as [start, stop) for independent 1-D grid.

    Handles ascending or descending lat/lon axes. Never falls back to a
    single native pixel. The user-facing width is fixed at 0.7 degrees.
    """
    axis = np.asarray(axis, dtype=float)
    targets = np.asarray(targets, dtype=float)
    if axis.ndim != 1 or axis.size < 2:
        raise ValueError('Expected 1D grid coordinate')
    ascending = axis[0] < axis[-1]
    aa = axis if ascending else axis[::-1]
    if not np.all(np.diff(aa) > 0):
        raise ValueError('Non-monotonic lat/lon grid')
    left = np.searchsorted(aa, targets - halfwidth, side='left')
    right = np.searchsorted(aa, targets + halfwidth, side='right')
    if ascending:
        return left, right
    return len(aa) - right, len(aa) - left


def _normalize_lon_targets(lon_grid, storm_lon):
    """Map longitude onto the continuous coordinate convention of this grid.

    In particular, a native 340..410E African grid is supported without
    wrongly converting 5E to 5 instead of 365E.
    """
    lon_grid = np.asarray(lon_grid, dtype=float)
    out = np.asarray(storm_lon, dtype=float)
    centre = (np.min(lon_grid) + np.max(lon_grid)) / 2
    return centre + (out - centre + 180) % 360 - 180


def sample_box_means(da_2d, lat_grid, lon_grid, storm_lat, storm_lon):
    """Unweighted average of native cells in each mandatory 0.7-degree box.

    Nearest grid cell is NOT supported. A partly out-of-domain footprint is
    flagged NaN instead of silently using a smaller spatial averaging scale.
    Missing native pixels are excluded from the average.
    """
    storm_lat = np.asarray(storm_lat, dtype=float)
    storm_lon = _normalize_lon_targets(lon_grid, storm_lon)
    lat_grid = np.asarray(lat_grid, dtype=float)
    lon_grid = np.asarray(lon_grid, dtype=float)
    y0, y1 = _box_index_ranges(lat_grid, storm_lat, BOX_RADIUS_DEG)
    x0, x1 = _box_index_ranges(lon_grid, storm_lon, BOX_RADIUS_DEG)
    result = np.full(len(storm_lat), np.nan, dtype='float64')
    eps = 1e-7
    for j in range(len(result)):
        if not (np.isfinite(storm_lat[j]) and np.isfinite(storm_lon[j])):
            continue
        # Never replace 0.7 degree sampling with a smaller clipped box.
        if (storm_lat[j] - BOX_RADIUS_DEG < lat_grid.min()-eps or
            storm_lat[j] + BOX_RADIUS_DEG > lat_grid.max()+eps or
            storm_lon[j] - BOX_RADIUS_DEG < lon_grid.min()-eps or
            storm_lon[j] + BOX_RADIUS_DEG > lon_grid.max()+eps):
            continue
        patch = np.asarray(da_2d.isel(
            lat=slice(int(y0[j]), int(y1[j])),
            lon=slice(int(x0[j]), int(x1[j]))
        ).values, dtype='float64')
        patch = patch[np.isfinite(patch) & (np.abs(patch) < 1e19)]
        if patch.size:
            result[j] = float(np.mean(patch))
    return result


def time_dict(times):
    """Source times, not filename or possibly stale global metadata, are truth."""
    return {native_key(t): i for i, t in enumerate(times)}


def sample_one_input(run, var, requests, input_root):
    """Native-time, 0.7-degree box means; return N values or None if incomplete.

    Filename ranges are used only to locate files. We verify requested
    times against decoded CF coordinates; no temporal interpolation is used.
    """
    files = source_files(str(input_root), run, var)
    if not files:
        LOG.warning('%s: no input NetCDF files for %s', run, var)
        return None

    times = requests['env_time_key'].unique()
    calendar = str(requests['calendar'].iloc[0])
    # Prefer exact synoptic times. Some surface time-mean fields are stamped
    # at :30; their mean over 11:00-12:00 represents the 12:00 environment.
    # Only use the preceding midpoint if its field is explicitly time: mean.
    preceding = {str(k): half_hour_before_key(str(k), calendar) for k in times}
    assigned = {str(k): (file_for_time(files, str(k)) or
                         file_for_time(files, preceding[str(k)])) for k in times}
    missing_files = [k for k, f in assigned.items() if f is None]
    if missing_files:
        LOG.warning('%s/%s: no files for %d/%d native times, first=%s',
                    run, var, len(missing_files), len(times), missing_files[:3])
        return None  # do not write incomplete outputs which would be skipped later

    out = np.full(len(requests), np.nan, dtype='float64')
    # Group by path: open each monthly source at most once.
    by_file = {}
    for k, path in assigned.items():
        by_file.setdefault(path, []).append(k)

    # Requests have RangeIndex by construction; use positional indices.
    for path, keys in by_file.items():
        try:
            with xr.open_dataset(path, decode_times=True, use_cftime=True) as ds:
                if var not in ds:
                    raise KeyError(f'{var} not in {path.name}')
                da = ds[var]
                if not all(d in da.dims for d in ('time', 'lat', 'lon')):
                    raise ValueError(f'{var} has unexpected dimensions {da.dims}')
                # Use lat/lon of this exact variable, including staggered UV grid.
                lat_grid = ds['lat'].values
                lon_grid = ds['lon'].values
                decoded = time_dict(ds['time'].values)
                actual_times = {}
                for k in keys:
                    if k in decoded:
                        actual_times[k] = k
                        continue
                    earlier = preceding[k]
                    if earlier not in decoded:
                        LOG.warning('%s/%s: no exact or preceding :30 time for %s in %s',
                                    run, var, k, path.name)
                        return None
                    cell_methods = str(da.attrs.get('cell_methods', '')).lower()
                    if 'time: mean' not in cell_methods:
                        LOG.warning('%s/%s: refusing :30 substitution for non-mean field (%s)',
                                    run, var, path.name)
                        return None
                    if 'time_bnds' in ds:
                        bound = ds['time_bnds'].isel(time=decoded[earlier], bnds=1).item()
                        if native_key(bound) != k:
                            LOG.warning('%s/%s: bound end %s does not match target %s',
                                        run, var, native_key(bound), k)
                            return None
                    actual_times[k] = earlier
                    LOG.info('%s/%s: hourly mean %s (ending %s)',
                             run, var, pretty_key(earlier), pretty_key(k))
                for k in keys:
                    pos = np.flatnonzero(requests['env_time_key'].to_numpy() == k)
                    lats = requests['tminlat'].to_numpy(dtype=float)[pos]
                    lons = requests['tminlon'].to_numpy(dtype=float)[pos]
                    # Never sample a native CP4 pixel. All fields are averaged
                    # across the same 0.7 x 0.7 degree footprint.
                    selected = sample_box_means(
                        da.isel(time=decoded[actual_times[k]]), lat_grid, lon_grid,
                        lats, lons
                    )
                    if len(selected) != len(pos):
                        raise ValueError(f'Point indexing produced {selected.shape}, expected {len(pos)}')
                    out[pos] = np.where(np.abs(selected) < 1e19, selected, np.nan)
        except Exception as exc:
            LOG.warning('Cannot sample %s/%s from %s: %s', run, var, path, exc)
            return None
    return out


def compute_output(var, sampled):
    for input_var, out in DIRECT_INPUTS.items():
        if var == out:
            return sampled.get(input_var)
    for layer, (u0, v0, u1, v1) in SHEAR_LAYERS.items():
        outputs = SHEAR_OUTPUTS[layer]
        if var not in outputs:
            continue
        if not all(k in sampled and sampled[k] is not None for k in (u0,v0,u1,v1)):
            return None
        du, dv = sampled[u1]-sampled[u0], sampled[v1]-sampled[v0]
        return (du, dv, np.hypot(du,dv))[outputs.index(var)]
    raise KeyError(var)


def needed_inputs(outputs):
    needed = set()
    for out in outputs:
        match = [raw for raw, name in DIRECT_INPUTS.items() if name == out]
        if match:
            needed.add(match[0]); continue
        for layer, names in SHEAR_OUTPUTS.items():
            if out in names:
                needed.update(SHEAR_LAYERS[layer])
                break
    return sorted(needed)


def output_path(root, run, var, year, month, storm_hour, rel_hours=(-6,0,6), scope="boxes"):
    # No hive partitions: source rows already have year/month, and this keeps
    # reading individual files from pandas free of schema name collisions.
    suffix = '' if list(rel_hours) == [-6,0,6] else '_rel_' + '_'.join(
        ('m'+str(abs(r)) if r<0 else 'p'+str(r)) for r in rel_hours)
    scope_suffix = '' if scope == 'boxes' else '-outside'
    return (root / run / var / f'year_{year:04d}' / f'month_{month:02d}' /
            f'part-{storm_hour:02d}utc{suffix}{scope_suffix}.parquet')


def save_one(var, values, requests, dest):
    out = requests.drop(columns=['env_time_key']).copy()
    out['variable'] = var
    out['value'] = values
    out['spatial_sampling'] = 'box_mean_0p7deg'
    out['sampling_box_radius_deg'] = BOX_RADIUS_DEG
    out['sampling_box_width_deg'] = 2 * BOX_RADIUS_DEG
    dest.parent.mkdir(parents=True, exist_ok=True)
    temp = dest.with_name(f'.{dest.name}.{os.getpid()}.tmp')
    try:
        out.to_parquet(temp, index=False)
        os.replace(temp, dest)
    finally:
        if temp.exists():
            temp.unlink()
    LOG.info('SAVED %-22s rows=%d %s', var, len(out), dest)


def process_base(base, args):
    match = BASE_FILE.search(base.as_posix())
    if not match:
        raise ValueError(f'Unexpected CP4 base parquet path: {base}')
    year, month = map(int, match.groups())
    run = base.parent.parent.name
    LOG.info('BASE %s %04d-%02d',run,year,month)
    # Read only metadata required for environmental requests.
    cols = ['run','storm_id','year','month','day','hour','minute','calendar',
            'clon','clat','tminlat','tminlon','storm_hour_end_utc']
    storms = pd.read_parquet(base, columns=cols)
    # Build all requests once; choose disjoint box/outside subsets below.
    # This means --full-domain can later fill missing outside rows without
    # recomputing or duplicating any already-completed boxed data.
    requests = make_requests(storms, args.storm_hours, args.rel_hours,
                             restrict_boxes=False)
    if requests.empty:
        LOG.info('NO MATCHING STORMS %s %04d-%02d',run,year,month)
        return

    root = (args.output_root / '_smoke_tests' / f'first_{args.max_storms:04d}_storms'
            if args.max_storms else args.output_root)
    requested_outputs = CORE_OUTPUTS if args.profile == 'core' else ALL_OUTPUTS
    scopes = (('boxes', 'outside') if args.full_domain else
              ('outside',) if args.outside_boxes else ('boxes',))

    # Each (hour, scope, variable) owns exactly one Parquet part.
    for shour in args.storm_hours:
        this_hour = requests.loc[requests['storm_utc_hour'] == shour]
        for scope in scopes:
            in_boxes = this_hour['regions'].ne('')
            req = this_hour.loc[in_boxes if scope == 'boxes' else ~in_boxes].reset_index(drop=True)
            if req.empty:
                continue
            if args.max_storms:
                keep = set(req['storm_id'].drop_duplicates().head(args.max_storms))
                req = req.loc[req['storm_id'].isin(keep)].reset_index(drop=True)
            dirs = {var: output_path(root,run,var,year,month,shour,args.rel_hours,scope=scope)
                    for var in requested_outputs}
            todo = sorted(v for v, p in dirs.items() if args.overwrite or not p.exists())
            if not todo:
                LOG.info('ALL EXISTS %s %04d-%02d %02d UTC %s',run,year,month,shour,scope)
                continue
            LOG.info('REQUEST %s %04d-%02d %02dUTC scope=%s storms=%d outputs=%d',
                     run,year,month,shour,scope,req['storm_id'].nunique(),len(todo))
            if args.dry_run:
                continue
            sampled = {}
            for var in needed_inputs(todo):
                sampled[var] = sample_one_input(run,var,req,args.input_root)
                if sampled[var] is not None:
                    LOG.info('SAMPLED %-10s %5d values (%5d non-NaN)',
                             var,len(sampled[var]),np.isfinite(sampled[var]).sum())
            for out in todo:
                values = compute_output(out,sampled)
                if values is None:
                    LOG.warning('SKIP output %s: required native model data missing',out)
                    continue
                if not np.isfinite(values).any():
                    LOG.warning('SKIP output %s: all sampled values are non-finite; will retry on rerun',out)
                    continue
                save_one(out,values,req,dirs[out])


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', action='append', help='Run directory, repeatable; defaults to all')
    p.add_argument('--month', type=int, default=8,
                   help='Month to process (default August); use --all for all completed months')
    p.add_argument('--year', type=int, help='Single year, otherwise first available unless --all-years')
    p.add_argument('--all-years', action='store_true',help='All completed years for selected month')
    p.add_argument('--all', action='store_true',help='All completed base-table months/years')
    p.add_argument('--storm-hours', type=int, nargs='+', default=[18],
                   help='Nominal UTC hours ending hourly storm snapshots (default 18)')
    p.add_argument('--rel-hours', type=int, nargs='+', default=[-6,0,6],
                   help='Native UTC environmental offsets, 6-hourly (default -6 0 6)')
    p.add_argument('--profile', choices=['core','all'], default='core')
    p.add_argument('--full-domain', action='store_true',
                   help='Produce both in-box and outside-box parts, skipping already completed parts')
    p.add_argument('--outside-boxes', action='store_true',
                   help='Sample ONLY storms outside OBS African boxes; no recomputation inside boxes')
    p.add_argument('--input-root',type=Path,default=INPUT_ROOT)
    p.add_argument('--base-root',type=Path,default=BASE_TABLES)
    p.add_argument('--output-root',type=Path,default=None,
                   help='Separate output root; default always uses dedicated 0.7-degree directory')
    p.add_argument('--max-storms',type=int,default=None,
                   help='Smoke test first N matching storms per base file, separate output')
    p.add_argument('--overwrite',action='store_true')
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--log-level',default='INFO')
    args=p.parse_args(argv)
    logging.basicConfig(level=args.log_level,
                        format='%(asctime)s %(levelname)s %(message)s')
    if args.output_root is None:
        args.output_root = OUTPUT_ROOT
    if args.full_domain and args.outside_boxes:
        p.error('--outside-boxes and --full-domain are mutually exclusive')
    if args.all and (args.year is not None or args.all_years):
        p.error('--all cannot be combined with --year or --all-years')
    if args.all_years and args.year is not None:
        p.error('--all-years cannot be combined with --year')
    if not args.storm_hours or any(h not in (0,6,12,18) for h in args.storm_hours):
        p.error('For this fixed-time version, --storm-hours must be among 0, 6, 12, 18')
    if not args.rel_hours or any(h % 6 != 0 for h in args.rel_hours):
        p.error('--rel-hours must be native six-hourly offsets')
    if args.max_storms is not None and args.max_storms < 1:
        p.error('--max-storms must be positive')
    if args.output_root.resolve()==args.base_root.resolve() or 'base_tables_tir-prcp' in str(args.output_root):
        p.error('Output must be separate from observational and CP4 base tables')
    # Prevent writing to OBS production trees even if an incorrect override is supplied.
    if 'ERA5_storm_timeseries' in str(args.output_root):
        p.error('Refusing to write into OBS ERA5 outputs')
    if args.output_root.resolve() == (BASE_ROOT/'CP4_storm_timeseries').resolve():
        p.error('Refusing to mix 0.7-degree boxes with older native-pixel CP4 outputs')

    runs=args.run or [f'{m}-{e}' for m in MODELS for e in ('HIST','FUT')]
    failures=[]
    processed=0
    for run in runs:
        run_dir=args.base_root/run
        if not run_dir.exists():
            LOG.warning('NO CP4 BASE TABLE DIRECTORY: %s',run_dir)
            continue
        paths=sorted(run_dir.glob('year_*/month_??.parquet'))
        if not args.all:
            paths=[f for f in paths if int(BASE_FILE.search(f.as_posix()).group(2))==args.month]
        if args.year is not None:
            paths=[f for f in paths if int(BASE_FILE.search(f.as_posix()).group(1))==args.year]
        elif not args.all and not args.all_years:
            paths=paths[:1]
        if args.all:
            # August first (then other months), as in full CP4 detection production.
            paths.sort(key=lambda f:(int(BASE_FILE.search(f.as_posix()).group(2))!=8,
                                     int(BASE_FILE.search(f.as_posix()).group(1)),
                                     int(BASE_FILE.search(f.as_posix()).group(2))))
        if not paths:
            LOG.warning('No completed %s base Parquets for requested months/years',run)
        for path in paths:
            try:
                process_base(path,args)
                processed+=1
            except Exception:
                LOG.exception('FAILED %s',path)
                failures.append(str(path))
    LOG.info('FINISHED %d monthly base files, %d failures',processed,len(failures))
    if failures:
        return 1
    return 0

if __name__=='__main__':
    sys.exit(main())

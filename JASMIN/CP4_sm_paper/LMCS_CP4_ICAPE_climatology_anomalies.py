#!/usr/bin/env python3
"""True CP4 ICAPE climatologies and anomalies for the existing MCS tables.

Defaults: +/-7 CP4 360-day-calendar days, including day zero, for each year
1998--2006. Compute ICAPE for EACH day before taking within-year means, then
weight the yearly means equally. Hist/future raw data are kept separate.

Each case keeps its original location, four boxes, and storm-day lsRain_noon
validity mask throughout the reference period. Surface parcels use that
BACKGROUND DAY's t2, q2 and p_srfc. Pressure profiles are at 12 UTC; surface
fields default to instantaneous 12 UTC, matching the original ICAPE script.
Use --surface-time-mode 3hmean ONLY if the absolute ICAPE tables were produced
with that option too. No new daily rain masks or storm-area (Smean) are used.

Inputs (defaults in mean3h_v2/tables):
  {tag}_mean_table_JASMIN_3hmeansVersion_rainMask_fullBox_ICAPE.csv
  {tag}_anom_table_JASMIN_3hmeansVersion_rainMask_fullBox.csv
Output: the latter filename with _ICAPE appended, containing the eight ICAPE
ANOMALY columns and all existing anomaly-table fields. A separate case-level
_ICAPE_climatology.csv records absolute ICAPE, background means, anomalies and
coverage; _ICAPE_coverage_by_year.csv records all nine annual means/counts.
All ICAPE quantities are J m-2. Input files are never overwritten.
The two tables need not contain identical case lists. Unmatched anomaly rows
are retained with NaN ICAPE; no absolute value is borrowed from another storm.
Matching reports list matched, anomaly-only and mean-only cases. A tiny
1e-5-degree coordinate tolerance accommodates numeric precision differences.
Use --audit-only to inspect matching without reading raw fields or storm masks.

Date checkpoints are grouped by RAW BACKGROUND DATE: each of the five raw
fields is opened once per worker task, shared across all relevant cases.
Cached storm masks are memory mapped. Restart the same command to resume.
Do not run two copies against the same work directory simultaneously.

Examples (run from proj_CEH, in the environment used for absolute ICAPE):
  python JASMIN/CP4_sm_paper/LMCS_CP4_ICAPE_climatology_anomalies.py \
      --tags hist --processes 1 --case-limit 3 --benchmark-dates 3
  python JASMIN/CP4_sm_paper/LMCS_CP4_ICAPE_climatology_anomalies.py \
      --tags hist --processes 1 --case-limit 3
  python JASMIN/CP4_sm_paper/LMCS_CP4_ICAPE_climatology_anomalies.py --processes 4

The first command tests only three background dates; it deliberately defers
merging. The second finishes all 135 reference dates for those three cases
(more dates if their calendar days differ) and produces SAMPLE outputs. Omit
both limits for production. Sample and production checkpoints are separate.
Before background processing, the first three usable storm cases are also
recalculated and compared against their absolute table values. A mismatch
stops the run (e.g. a surface-time setting inconsistent with the input table).

Like the original climatology, available valid days are averaged within each
year (--min-valid-days defaults to 1). By default ALL requested years must
contribute; use --min-valid-days 12 for stricter 12-of-15-day coverage. Missing
data are never replaced by zero. All-zero physical CAPE remains valid zero.

The calendar intentionally fixes the old Gregorian/day-31-clipping behavior:
there are 15 distinct model-calendar dates per year, not duplicated day 30.
The event itself remains in its year's reference window, as in the original.

Dependencies: numpy, pandas, xarray, netCDF4 or another xarray NetCDF backend,
cftime, scipy, MetPy, and the existing JASMIN.MetUM_variables and
shared.utils.u_interpolate modules. Numerical CAPE settings match
LMCS_CP4_add_ICAPE_to_tables.py, including source top 800 hPa and CAPE top
500 hPa/all available levels. See:
https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.cape_cin.html
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd


RAW_ROOT = Path('/home/users/cornkle/linked_CP4')
BOX_ROOT = Path('/gws/ssde/j25b/lmcs/cklein/CP_models/MCS_files/WAf/CP4_box_JASMIN/mean3h_v2')
TABLE_STEM = 'table_JASMIN_3hmeansVersion_rainMask_fullBox'
DOMAIN = (-18., 25., 5., 25.)
BOXES = {'0.25deg': 3, '1deg': 11, '2deg': 22, 'fullBox': None}
COLUMNS = [f'{kind}_{box}' for kind in ('ICAPE500', 'ICAPEfull') for box in BOXES]
KEY_FIELDS = ('year', 'month', 'day', 'hour', 'lon', 'lat')
HALF_BOX, ENV_HOUR, SOURCE_TOP, PARTIAL_TOP, G = 57, 12, 800., 500., 9.80665
MIN_AREA_FRACTION = 0.5
ALGORITHM = 'icape-daily-climatology-v1'
# np.round(one_element_array, 1) produces e.g. [11.], [-0.] or [ 7.2].
# Allow trailing decimal points, whitespace and scientific notation; require
# the complete latitude token rather than accepting a numeric prefix.
FILENAME_NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
LONLAT_RE = re.compile(rf'_lonXlat_\[?\s*({FILENAME_NUMBER})\s*\]?_\[?\s*({FILENAME_NUMBER})\s*\]?(?=_|$)')
# Status is per case and spatial scale (each scale produces two diagnostics).
INACTIVE, OK, PROFILE_MISSING, ERROR = 0, 1, 2, 3
WORKER = {}


def scalar(value, field):
    if pd.isna(value):
        raise ValueError(f'Missing {field}')
    if isinstance(value, (int, float, np.number)):
        result = float(value)
    else:
        text = str(value).strip().strip('[]').strip()
        result = float(text)
    if not np.isfinite(result):
        raise ValueError(f'Nonfinite {field}: {value!r}')
    return result


def read_table(path):
    # low_memory=False avoids the pandas mixed-chunk/index parser error.
    raw = pd.read_csv(path, index_col=0, low_memory=False)
    needed = set(KEY_FIELDS)
    if needed.issubset(raw.index.astype(str)):
        frame, transposed = raw.T.copy(), True
    elif needed.issubset(raw.columns.astype(str)):
        frame, transposed = raw.copy(), False
    else:
        raise ValueError(f'{path}: cannot identify case metadata/orientation')
    if not frame.columns.is_unique or not frame.index.is_unique:
        raise ValueError(f'{path}: duplicate field names or row labels')
    return frame, transposed


def records_from_table(frame):
    records = []
    for i, (_, row) in enumerate(frame.iterrows()):
        rec = {key: scalar(row[key], key) for key in KEY_FIELDS}
        for key in KEY_FIELDS[:4]:
            if rec[key] != int(rec[key]):
                raise ValueError(f'Noninteger {key} for case {i}')
            rec[key] = int(rec[key])
        if not (1 <= rec['month'] <= 12 and 1 <= rec['day'] <= 30 and 0 <= rec['hour'] <= 23):
            raise ValueError(f'Invalid 360-day date/hour for case {i}: {rec}')
        for key in ('index', 'indx', 'indy'):
            if key in row and not pd.isna(row[key]):
                rec[key] = scalar(row[key], key)
        rec['case_id'] = i
        records.append(rec)
    return records


def align_cases(mean_records, anomaly_records, tolerance=1e-5, strict=False):
    """Left join onto anomaly rows; -1 means genuinely absent from mean table.

    Date/hour must match exactly. The small absolute coordinate tolerance is
    only for stored-number precision, NOT nearest-storm matching. Ambiguous
    or many-to-one matches still raise rather than silently choosing a storm.
    """
    def stamp(rec):
        return tuple(rec[k] for k in KEY_FIELDS[:4])
    lookup = defaultdict(list)
    for rec in mean_records:
        lookup[stamp(rec)].append(rec)
    positions = []
    for rec in anomaly_records:
        candidates = [r for r in lookup.get(stamp(rec), [])
                      if abs(r['lon'] - rec['lon']) <= tolerance and abs(r['lat'] - rec['lat']) <= tolerance]
        if len(candidates) > 1 and 'index' in rec:
            preferred = [r for r in candidates if r.get('index') == rec['index']]
            if len(preferred) == 1:
                candidates = preferred
        if len(candidates) > 1:
            raise ValueError(f'Anomaly case {rec["case_id"]}: ambiguous mean-table matches '
                             f'{[r["case_id"] for r in candidates]} at {stamp(rec)}, {rec["lon"]}, {rec["lat"]}')
        positions.append(candidates[0]['case_id'] if candidates else -1)
    matched = [i for i in positions if i >= 0]
    if len(set(matched)) != len(matched):
        raise ValueError('Two anomaly-table rows match the same mean-table case')
    if strict and (len(matched) != len(positions) or set(matched) != set(range(len(mean_records)))):
        raise ValueError('Mean and anomaly tables do not contain exactly the same cases (--strict-case-match)')
    return np.asarray(positions, dtype=int)


def matching_report(mean_records, anomaly_records, positions):
    means = {r['case_id']: r for r in mean_records}
    rows = []
    for rec, mid in zip(anomaly_records, positions):
        mean = means.get(int(mid))
        row = {key: rec[key] for key in KEY_FIELDS[:4]}
        row.update(status='matched' if mean is not None else 'anomaly_only',
                   anomaly_case_id=rec['case_id'], mean_case_id=int(mid) if mean is not None else np.nan,
                   anomaly_lon=rec['lon'], anomaly_lat=rec['lat'], anomaly_index=rec.get('index', np.nan),
                   mean_lon=mean['lon'] if mean is not None else np.nan,
                   mean_lat=mean['lat'] if mean is not None else np.nan,
                   mean_index=mean.get('index', np.nan) if mean is not None else np.nan,
                   delta_lon_deg=mean['lon'] - rec['lon'] if mean is not None else np.nan,
                   delta_lat_deg=mean['lat'] - rec['lat'] if mean is not None else np.nan)
        rows.append(row)
    used = set(positions[positions >= 0])
    for mid, rec in means.items():
        if mid not in used:
            row = {key: rec[key] for key in KEY_FIELDS[:4]}
            row.update(status='mean_only', mean_case_id=mid, mean_lon=rec['lon'], mean_lat=rec['lat'],
                       mean_index=rec.get('index', np.nan), anomaly_case_id=np.nan)
            rows.append(row)
    return pd.DataFrame(rows)


def shift_360(year, month, day, offset):
    """Calendar arithmetic without Gregorian February or day-31 clipping."""
    ordinal = year * 360 + (month - 1) * 30 + day - 1 + offset
    yy, rem = divmod(ordinal, 360)
    mm, dd = divmod(rem, 30)
    return f'{yy:04d}{mm + 1:02d}{dd + 1:02d}'


def make_schedule(records, active, years, half_window):
    # Each row in a job has a case ID and the CENTRAL reference-year index.
    # The raw date's year can differ at a December/January window boundary.
    grouped = defaultdict(lambda: [[], []])
    for rec in records:
        cid = rec['case_id']
        if not active[cid].any():
            continue
        for yi, year in enumerate(years):
            for offset in range(-half_window, half_window + 1):
                date = shift_360(year, rec['month'], rec['day'], offset)
                grouped[date][0].append(cid)
                grouped[date][1].append(yi)
    return [(date, np.asarray(ids, dtype=np.int32), np.asarray(yis, dtype=np.int16))
            for date, (ids, yis) in sorted(grouped.items())]


def file_hash(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + f'.tmp.{os.getpid()}')
    with open(temporary, 'w') as handle:
        json.dump(data, handle, indent=2, sort_keys=True, allow_nan=False)
    os.replace(temporary, path)


def atomic_csv(path, frame, index=True):
    path = Path(path)
    temporary = path.with_name(path.name + f'.tmp.{os.getpid()}')
    frame.to_csv(temporary, index=index)
    os.replace(temporary, path)


def atomic_npz(path, **arrays):
    path = Path(path)
    temporary = path.with_name(path.name + f'.tmp.{os.getpid()}')
    with open(temporary, 'wb') as handle:
        np.savez_compressed(handle, **arrays)
    os.replace(temporary, path)


def find_2d_file(row, box_root, tag):
    stamp = f'{row["year"]:04d}-{row["month"]:02d}-{row["day"]:02d}_{row["hour"]:02d}:00:00'
    folder = Path(box_root) / f'mean_{tag}'
    candidates = sorted(folder.glob(f'{stamp}_*_lonXlat_*.nc'))
    target_lon = row.get('indx', round(row['lon'], 1))
    target_lat = row.get('indy', round(row['lat'], 1))
    matching, parsed, unparsed = [], [], []
    for path in candidates:
        match = LONLAT_RE.search(path.name)
        if match:
            lon, lat = float(match[1]), float(match[2])
            distance = (lon - target_lon) ** 2 + (lat - target_lat) ** 2
            parsed.append((distance, path.name, lon, lat))
            if distance <= 2 * 0.05 ** 2:
                matching.append((distance, path))
        else:
            unparsed.append(path.name)
    if not matching:
        nearest = [(name, lon, lat) for _, name, lon, lat in sorted(parsed)[:3]]
        raise FileNotFoundError(
            f'No location-matched mean 2-D file for case {row["case_id"]} ({stamp}). '
            f'Folder={folder}, exists={folder.is_dir()}; index={row.get("index")}; '
            f'target filename coordinates=({target_lon}, {target_lat}); timestamp candidates={len(candidates)}; '
            f'nearest parsed (filename, lon, lat)={nearest}; unparsed filenames={unparsed[:5]}')
    # Prefer the original storm index only among geographically valid candidates.
    if 'index' in row:
        preferred = [item for item in matching if item[1].name.startswith(f'{stamp}_{int(row["index"])}_')]
        matching = preferred or matching
    matching.sort()
    if len(matching) > 1 and np.isclose(matching[0][0], matching[1][0]):
        raise ValueError(f'Ambiguous 2-D location match for case {row["case_id"]}')
    return matching[0][1]


def scale_slice(distance):
    return slice(None) if distance is None else slice(HALF_BOX - distance, HALF_BOX + distance + 1)


def prepare_masks(run_dir, records, original, box_root, tag):
    """Read each storm file once; processes share the disk-backed bool array."""
    mask_path, active_path = run_dir / 'masks.npy', run_dir / 'active.npy'
    count = len(records)
    if mask_path.exists() and active_path.exists():
        masks = np.load(mask_path, mmap_mode='r', allow_pickle=False)
        active = np.load(active_path, allow_pickle=False)
        if masks.shape != (count, 115, 115) or active.shape != (count, 4):
            raise ValueError('Invalid mask cache shape; use a new --work-dir')
        return active
    import xarray as xr
    temporary = mask_path.with_name(f'masks.tmp.{os.getpid()}.npy')
    masks = np.lib.format.open_memmap(temporary, mode='w+', dtype=bool, shape=(count, 115, 115))
    masks[:] = False
    active = np.isfinite(original[:, :4]) | np.isfinite(original[:, 4:])
    sources = []
    for rec in records:
        cid = rec['case_id']
        if active[cid].any():
            path = find_2d_file(rec, box_root, tag)
            with xr.open_dataset(path) as ds:
                rain = np.asarray(ds['lsRain_noon'].squeeze().values, dtype=float)
            if rain.shape != (115, 115):
                raise ValueError(f'{path}: expected a 115x115 noon-rain mask, found {rain.shape}')
            masks[cid] = np.isfinite(rain) & (rain <= 0)
            for bi, distance in enumerate(BOXES.values()):
                sl = scale_slice(distance)
                if active[cid, bi] and masks[cid, sl, sl].mean() < MIN_AREA_FRACTION:
                    raise ValueError(f'Case {cid}, {list(BOXES)[bi]}: finite table ICAPE but mask fails 50%; check inputs')
            sources.append({'case_id': cid, 'file': str(path), 'mtime_ns': path.stat().st_mtime_ns})
        if (cid + 1) % 250 == 0 or cid + 1 == count:
            print(f'{tag}: prepared masks {cid + 1}/{count}', flush=True)
    masks.flush()
    del masks
    os.replace(temporary, mask_path)
    temp_active = active_path.with_name(f'active.tmp.{os.getpid()}.npy')
    with open(temp_active, 'wb') as handle:
        np.save(handle, active, allow_pickle=False)
    os.replace(temp_active, active_path)
    atomic_csv(run_dir / 'mask_sources.csv', pd.DataFrame(sources), index=False)
    return active


def runtime_dependencies():
    import xarray as xr
    from metpy import calc
    from metpy.units import units
    from JASMIN import MetUM_variables as mu
    from shared.utils import u_interpolate as u_int
    return xr, calc, units, mu, u_int


def load_raw(raw_dir, variable, date, mode):
    xr, _, _, mu, _ = WORKER['deps']
    paths = sorted((raw_dir / variable).glob(f'*_{date}*.nc'))
    if not paths:
        raise FileNotFoundError(f'No {variable} raw file on {date} in {raw_dir}')
    if len(paths) != 1:
        raise ValueError(f'Multiple {variable} files on {date}: {paths}')
    with xr.open_dataset(paths[0]) as ds:
        da = ds[mu.create_CP4_filename(variable)]
        if float(da.longitude.max()) > 180:
            da = da.assign_coords(longitude=da.longitude - 360)
        # Do not silently reverse the storm-mask spatial ordering.
        for coordinate in ('latitude', 'longitude'):
            if not np.all(np.diff(da[coordinate].values) > 0):
                raise ValueError(f'{variable}: nonascending {coordinate}; verify alignment with original ICAPE loader')
        da = da.sel(longitude=slice(DOMAIN[0], DOMAIN[1]), latitude=slice(DOMAIN[2], DOMAIN[3]))
        if 'time' in da.dims:
            hours = da.time.dt.hour
            if mode == '3hmean':
                da = da.where((hours >= ENV_HOUR - 2) & (hours <= ENV_HOUR), drop=True)
                if da.sizes['time'] == 0:
                    raise ValueError(f'{variable}: no 10--12 UTC data on {date}')
                da = da.mean('time')
            else:
                da = da.where(hours == ENV_HOUR, drop=True)
                if da.sizes['time'] != 1:
                    raise ValueError(f'{variable}: expected one 12 UTC sample on {date}')
                da = da.squeeze('time', drop=True)
        da = da.squeeze().load()
    return da


def normalise(values, kind):
    out = np.asarray(values, dtype=float)
    finite = out[np.isfinite(out)]
    if kind == 'q':
        finite = finite[finite > 0]
    if finite.size == 0:
        raise ValueError(f'No finite {kind} values')
    median = float(np.median(finite))
    if kind == 'T':
        if median > 1000:
            out = out / 100.
        valid = 150 < float(np.median(out[np.isfinite(out)])) < 350
    elif kind == 'q':
        if median > .2:
            out = out / 1000.
            median /= 1000.
        valid = 0 < median < .1
    else:
        if median > 2000:
            out = out / 100.
            median /= 100.
        valid = 700 < median < 1100
    if not valid:
        raise ValueError(f'Implausible {kind} median/units: {median}')
    return out


def regrid_surface(da, target):
    da = da.transpose('latitude', 'longitude')
    if (np.array_equal(da.longitude.values, target.longitude.values)
            and np.array_equal(da.latitude.values, target.latitude.values)):
        return np.asarray(da.values, dtype=float)
    # Only one coordinate/weight cache per worker to bound memory use.
    signature = tuple(np.asarray(a).tobytes() for a in
                      (da.longitude.values, da.latitude.values, target.longitude.values, target.latitude.values))
    if WORKER.get('grid_signature') != signature:
        WORKER['grid_weights'] = WORKER['deps'][4].interpolation_weights(
            da.longitude, da.latitude, target.longitude, target.latitude)
        WORKER['grid_signature'] = signature
    inds, weights, shape = WORKER['grid_weights']
    return np.asarray(WORKER['deps'][4].interpolate_data(da.values, inds, weights, shape), dtype=float)


def load_day_fields(date):
    raw_dir = Path(WORKER['config']['raw_root']) / ('hist' if WORKER['config']['tag'] == 'hist' else 'future')
    t_da = load_raw(raw_dir, 't_pl', date, 'instant').sortby('pressure', ascending=False)
    q_da = load_raw(raw_dir, 'q_pl', date, 'instant').sortby('pressure', ascending=False)
    t_da, q_da = (da.transpose('pressure', 'latitude', 'longitude') for da in (t_da, q_da))
    for coord in ('pressure', 'latitude', 'longitude'):
        if not np.array_equal(t_da[coord].values, q_da[coord].values):
            raise ValueError(f'T/q {coord} grids differ on {date}')
    p = np.asarray(t_da.pressure.values, dtype=float)
    if np.nanmedian(p) > 2000:
        p = p / 100.
    t, q = normalise(t_da.values, 'T'), normalise(q_da.values, 'q')
    t[t == 0], q[q == 0] = np.nan, np.nan
    result = {'pressure': p, 'temperature': t, 'humidity': q,
              'longitude': np.asarray(t_da.longitude.values), 'latitude': np.asarray(t_da.latitude.values)}
    for variable, name, kind in (('t2', 't2', 'T'), ('q2', 'q2', 'q'), ('p_srfc', 'ps', 'p')):
        da = load_raw(raw_dir, variable, date, WORKER['config']['surface_time_mode'])
        result[name] = normalise(regrid_surface(da, t_da), kind)
    return result


def centred_cut(fields, lon, lat):
    xi = int(np.nanargmin(np.abs(fields['longitude'] - lon)))
    yi = int(np.nanargmin(np.abs(fields['latitude'] - lat)))
    if min(xi, yi) < HALF_BOX:
        raise ValueError(f'Incomplete fullBox at {lon}, {lat}')
    xs, ys = slice(xi - HALF_BOX, xi + HALF_BOX + 1), slice(yi - HALF_BOX, yi + HALF_BOX + 1)
    cut = {'pressure': fields['pressure']}
    for name in ('temperature', 'humidity'):
        cut[name] = fields[name][:, ys, xs]
    for name in ('t2', 'q2', 'ps'):
        cut[name] = fields[name][ys, xs]
    if cut['temperature'].shape[1:] != (115, 115):
        raise ValueError(f'Incomplete fullBox at {lon}, {lat}')
    return cut


def finite_mean(values, axis=None):
    """nanmean-equivalent for finite inputs, without empty-level warnings."""
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    count = valid.sum(axis=axis)
    total = np.where(valid, values, 0.).sum(axis=axis)
    result = np.full(np.shape(total), np.nan, dtype=float)
    np.divide(total, count, out=result, where=count > 0)
    return result


def area_mean_profile(cut, mask, distance):
    sl = scale_slice(distance)
    selected = mask[sl, sl]
    if selected.mean() < MIN_AREA_FRACTION:
        raise ValueError('Cached active box unexpectedly fails the fixed rain mask')
    t = finite_mean(np.where(selected[None], cut['temperature'][:, sl, sl], np.nan), axis=(1, 2))
    q = finite_mean(np.where(selected[None], cut['humidity'][:, sl, sl], np.nan), axis=(1, 2))
    t2, q2, ps = [float(finite_mean(np.where(selected, cut[name][sl, sl], np.nan))) for name in ('t2', 'q2', 'ps')]
    return cut['pressure'], t, q, ps, (t2, q2)


def cape_pair_for_sources(pressure, temperature, humidity, ps, surface_state):
    _, calc, units, _, _ = WORKER['deps']
    if not np.all(np.isfinite([ps, *surface_state])):
        return np.nan, np.nan
    p, t, q = (np.asarray(a, dtype=float) for a in (pressure, temperature, humidity))
    good = np.isfinite(p) & np.isfinite(t) & np.isfinite(q) & (p < ps - .05)
    p, t, q = p[good], t[good], q[good]
    if p.size < 5:
        return np.nan, np.nan
    order = np.argsort(p)[::-1]
    p, t, q = p[order], t[order], q[order]
    p = np.concatenate(([ps], p))
    t = np.concatenate(([surface_state[0]], t))
    q = np.concatenate(([surface_state[1]], q))
    if np.any(np.diff(p) >= 0):
        raise ValueError('Sounding pressure is not strictly decreasing')
    if p.min() > PARTIAL_TOP or not np.any(np.isclose(p, SOURCE_TOP, atol=.6)):
        return np.nan, np.nan
    # Raw CP4 must include 500 hPa; do not silently terminate at 550 hPa.
    if not np.any(np.isclose(p, PARTIAL_TOP, atol=.01)):
        return np.nan, np.nan
    pu, tu = p * units.hPa, t * units.K
    td = calc.dewpoint_from_specific_humidity(pu, tu, q * units('kg/kg')).to('K')
    td = np.minimum(td.magnitude, t) * units.K
    sources = np.where(p >= SOURCE_TOP - .01)[0]
    if sources.size < 2:
        return np.nan, np.nan
    values = []
    for idx in sources:
        pp, tt, dd = pu[idx:], tu[idx:], td[idx:]
        partial = pp >= PARTIAL_TOP * units.hPa
        if pp.size < 3 or partial.sum() < 3:
            return np.nan, np.nan
        parcel = calc.parcel_profile(pp, tu[idx], td[idx]).to('K')
        full, _ = calc.cape_cin(pp, tt, dd, parcel, which_lfc='bottom', which_el='top')
        low, _ = calc.cape_cin(pp[partial], tt[partial], dd[partial], parcel[partial],
                               which_lfc='bottom', which_el='top')
        pair = np.array([float(low.to('J/kg').magnitude), float(full.to('J/kg').magnitude)])
        if not np.isfinite(pair).all():
            return np.nan, np.nan  # NEVER let max(0, nan) turn a failed result into zero.
        values.append(np.maximum(0., pair))
    integrate = getattr(np, 'trapezoid', None) or np.trapz
    result = integrate(np.asarray(values)[::-1], p[sources][::-1] * 100., axis=0) / G
    return float(result[0]), float(result[1])


def init_worker(config, run_dir, records):
    WORKER.clear()
    WORKER.update(config=config, run_dir=Path(run_dir), records=records,
                  deps=runtime_dependencies(), masks=np.load(Path(run_dir) / 'masks.npy', mmap_mode='r'),
                  active=np.load(Path(run_dir) / 'active.npy'))


def validate_absolute(records, original, limit, run_dir):
    """Check raw loader/masks/numerics against existing absolute table values."""
    ids = np.flatnonzero(np.isfinite(original).any(axis=1))[:limit]
    grouped = defaultdict(list)
    for cid in ids:
        rec = records[cid]
        grouped[f'{rec["year"]:04d}{rec["month"]:02d}{rec["day"]:02d}'].append(cid)
    checked, max_difference = 0, 0.
    for date, case_ids in grouped.items():
        fields = load_day_fields(date)
        for cid in case_ids:
            rec = records[cid]
            cut = centred_cut(fields, rec['lon'], rec['lat'])
            for bi, distance in enumerate(BOXES.values()):
                expected = original[cid, [bi, bi + 4]]
                finite = np.isfinite(expected)
                if not finite.any():
                    continue
                actual = np.asarray(cape_pair_for_sources(*area_mean_profile(cut, WORKER['masks'][cid], distance)))
                if not np.allclose(actual[finite], expected[finite], rtol=1e-5, atol=.1):
                    raise ValueError(f'Absolute validation failed: case {cid}, {list(BOXES)[bi]}, '
                                     f'table={expected}, recalculated={actual}. Check surface-time-mode, '
                                     'raw inputs, mask matching and the absolute ICAPE calculation before continuing.')
                max_difference = max(max_difference, float(np.abs(actual[finite] - expected[finite]).max()))
                checked += int(finite.sum())
        del fields
    result = {'cases': len(ids), 'values_checked': checked, 'max_absolute_difference_J_m2': max_difference}
    atomic_json(run_dir / 'absolute_validation.json', result)
    print(f'Absolute-table validation passed: {result}', flush=True)


def process_day(job):
    date, ids, year_indices = job
    started = time.monotonic()
    values = np.full((len(ids), 8), np.nan)
    active = WORKER['active'][ids]
    status = np.zeros((len(ids), 4), dtype=np.uint8)
    errors, day_error = [], ''
    try:
        fields = load_day_fields(date)
    except Exception as exc:
        status[active] = ERROR
        day_error = f'{type(exc).__name__}: {exc}'
        fields = None
    last_log = started
    if fields is not None:
        for ri, cid in enumerate(ids):
            rec = WORKER['records'][cid]
            try:
                cut = centred_cut(fields, rec['lon'], rec['lat'])
            except Exception as exc:
                status[ri, active[ri]] = ERROR
                errors.append(f'case {cid}: {type(exc).__name__}: {exc}')
                continue
            for bi, distance in enumerate(BOXES.values()):
                if not active[ri, bi]:
                    continue
                try:
                    pair = cape_pair_for_sources(*area_mean_profile(cut, WORKER['masks'][cid], distance))
                    values[ri, [bi, bi + 4]] = pair
                    status[ri, bi] = OK if np.isfinite(pair).all() else PROFILE_MISSING
                except Exception as exc:
                    status[ri, bi] = ERROR
                    errors.append(f'case {cid}, {list(BOXES)[bi]}: {type(exc).__name__}: {exc}')
            now = time.monotonic()
            if now - last_log >= 60:
                print(f'{WORKER["config"]["tag"]} background {date}: {ri + 1}/{len(ids)} cases', flush=True)
                last_log = now
    path = WORKER['run_dir'] / 'daily' / f'{date}.npz'
    atomic_npz(path, case_id=ids, year_index=year_indices, values=values, status=status,
               date=np.array(date), errors=np.asarray(errors, dtype=str), day_error=np.array(day_error),
               run_signature=np.array(WORKER['config']['run_signature']))
    summary = {'tag': WORKER['config']['tag'], 'background_date': date, 'cases': len(ids),
               'valid_boxes': int((status == OK).sum()), 'missing_profiles': int((status == PROFILE_MISSING).sum()),
               'error_boxes': int((status == ERROR).sum()), 'seconds': round(time.monotonic() - started, 1)}
    if day_error:
        summary['error'] = day_error
    elif errors:
        summary['first_error'] = errors[0]
    return summary


def read_checkpoint(path, job, signature):
    with np.load(path, allow_pickle=False) as checkpoint:
        data = {key: checkpoint[key] for key in checkpoint.files}
    date, ids, years = job
    if (str(data['run_signature']) != signature or str(data['date']) != date
            or not np.array_equal(data['case_id'], ids) or not np.array_equal(data['year_index'], years)
            or data['values'].shape != (len(ids), 8) or data['status'].shape != (len(ids), 4)):
        raise ValueError(f'Checkpoint does not match this run: {path}')
    return data


def finish_statistics(totals, counts, min_days, min_years):
    annual = np.full(totals.shape, np.nan)
    np.divide(totals, counts, out=annual, where=counts >= min_days)
    n_years = np.isfinite(annual).sum(axis=0)
    background = finite_mean(annual, axis=0)
    background[n_years < min_years] = np.nan
    return annual, background, n_years


def merge_results(args, config, run_dir, jobs, records, original, anomaly, anomaly_order, transposed, output, matches):
    n_years, n_cases = len(args.years), len(records)
    totals = np.zeros((n_years, n_cases, 8), dtype=float)
    counts = np.zeros((n_years, n_cases, 8), dtype=np.int16)
    issues = []
    for job in jobs:
        date, ids, yis = job
        data = read_checkpoint(run_dir / 'daily' / f'{date}.npz', job, config['run_signature'])
        vals = data['values']
        good = np.isfinite(vals)
        totals[yis, ids] += np.where(good, vals, 0.)
        counts[yis, ids] += good.astype(np.int16)
        if str(data['day_error']):
            issues.append({'background_date': date, 'error': str(data['day_error'])})
        issues.extend({'background_date': date, 'error': str(error)} for error in data['errors'])
    annual, background, n_year_valid = finish_statistics(totals, counts, args.min_valid_days, args.min_valid_years)
    deviations = original - background
    # Never use -1 to index deviations: it denotes an unmatched row, not the last mean case.
    matched = anomaly_order >= 0
    output_values = np.full((len(anomaly_order), 8), np.nan)
    output_values[matched] = deviations[anomaly_order[matched]]
    enriched = anomaly.copy()
    for ci, column in enumerate(COLUMNS):
        enriched[column] = output_values[:, ci]
    metadata = pd.DataFrame(records).set_index('case_id')
    audit = metadata.copy()
    for ci, column in enumerate(COLUMNS):
        audit[f'{column}_absolute'] = original[:, ci]
        audit[f'{column}_clim'] = background[:, ci]
        audit[f'{column}_anom'] = deviations[:, ci]
        audit[f'{column}_n_days'] = counts[:, :, ci].sum(axis=0)
        audit[f'{column}_n_years'] = n_year_valid[:, ci]
    by_year = []
    for yi, year in enumerate(args.years):
        block = metadata.copy()
        block['reference_year'] = year
        for ci, column in enumerate(COLUMNS):
            block[f'{column}_mean'] = annual[yi, :, ci]
            block[f'{column}_n_days'] = counts[yi, :, ci]
        by_year.append(block)
    expected = 2 * args.half_window + 1
    report = {
        'settings': config, 'min_valid_days_per_year': args.min_valid_days,
        'min_valid_years': args.min_valid_years, 'expected_days_per_year': expected,
        'original_missing_values': int((~np.isfinite(original)).sum()),
        'additional_missing_anomalies': int((np.isfinite(original) & ~np.isfinite(deviations)).sum()),
        'valid_anomaly_values': int(np.isfinite(output_values).sum()),
        'unmatched_anomaly_rows_in_output': int((~matched).sum()),
        'full_input_case_matching': {str(k): int(v) for k, v in matches['status'].value_counts().items()},
        'exception_records': len(issues), 'units': 'J m-2',
        'definition': 'absolute table ICAPE minus equally weighted yearly means of daily ICAPE',
        'output': str(output),
    }
    # Preserve all original files, even if a custom output was accidentally set to an input.
    stem = output.with_suffix('')
    destinations = [output, Path(str(stem) + '_climatology.csv'), Path(str(stem) + '_coverage_by_year.csv'),
                    Path(str(stem) + '_background_errors.csv'), Path(str(stem) + '_metadata.json'),
                    Path(str(stem) + '_case_matching.csv')]
    sources = {Path(config[key]).resolve() for key in ('mean_table', 'anomaly_table')}
    for path in destinations:
        if path.resolve() in sources:
            raise ValueError(f'Refusing to overwrite input table: {path}')
        if path.exists() and not args.overwrite_output:
            raise FileExistsError(f'{path} exists; use --overwrite-output to replace generated outputs')
    output.parent.mkdir(parents=True, exist_ok=True)
    atomic_csv(destinations[1], audit)
    atomic_csv(destinations[2], pd.concat(by_year))
    atomic_csv(destinations[3], pd.DataFrame(issues, columns=['background_date', 'error']), index=False)
    atomic_json(destinations[4], report)
    atomic_csv(destinations[5], matches, index=False)
    atomic_csv(output, enriched.T if transposed else enriched)
    print(f'Saved {output}\nValid anomaly values: {report["valid_anomaly_values"]}; '
          f'additional missing: {report["additional_missing_anomalies"]}; '
          f'unmatched anomaly rows: {int((~matched).sum())}; exception records: {len(issues)}', flush=True)


def run_tag(args, tag):
    table_dir = args.table_dir or args.box_root / 'tables'
    mean_path = getattr(args, f'{tag}_mean_table') or table_dir / f'{tag}_mean_{TABLE_STEM}_ICAPE.csv'
    anomaly_path = getattr(args, f'{tag}_anom_table') or table_dir / f'{tag}_anom_{TABLE_STEM}.csv'
    frame, _ = read_table(mean_path)
    anomaly, transposed = read_table(anomaly_path)
    mean_records = records_from_table(frame)
    anomaly_records = records_from_table(anomaly)
    anomaly_order = align_cases(mean_records, anomaly_records, args.match_tolerance, args.strict_case_match)
    matches = matching_report(mean_records, anomaly_records, anomaly_order)
    for column in COLUMNS:
        if column not in frame:
            raise ValueError(f'{mean_path}: missing {column}')
    original = frame[COLUMNS].apply(pd.to_numeric, errors='raise').to_numpy(dtype=float)
    sample_suffix = ''
    if args.case_limit:
        limit = min(args.case_limit, len(mean_records))
        selected = (anomaly_order >= 0) & (anomaly_order < limit)
        anomaly, anomaly_order = anomaly.iloc[np.flatnonzero(selected)].copy(), anomaly_order[selected]
        mean_records, original = mean_records[:limit], original[:limit]
        sample_suffix = f'_sample{limit}'
    if not mean_records:
        raise ValueError('No cases in input table')
    output = getattr(args, f'{tag}_output') or anomaly_path.with_name(anomaly_path.stem + '_ICAPE' + sample_suffix + '.csv')
    if output.resolve() in (mean_path.resolve(), anomaly_path.resolve()):
        raise ValueError('Output must not overwrite an input table')
    config = {'algorithm': ALGORITHM, 'tag': tag, 'years': args.years, 'half_window': args.half_window,
              'calendar': '360_day', 'include_central_day': True, 'hour': ENV_HOUR,
              'surface_time_mode': args.surface_time_mode, 'raw_root': str(args.raw_root.resolve()),
              'box_root': str(args.box_root.resolve()), 'mean_table': str(mean_path.resolve()),
              'anomaly_table': str(anomaly_path.resolve()), 'mean_sha256': file_hash(mean_path),
              'anomaly_sha256': file_hash(anomaly_path), 'n_cases': len(mean_records), 'boxes': BOXES,
              'source_top_hPa': SOURCE_TOP, 'partial_top_hPa': PARTIAL_TOP,
              'minimum_valid_area': MIN_AREA_FRACTION, 'matching_policy': 'anomaly-left-join-v2',
              'match_tolerance_deg': args.match_tolerance}
    signature = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    config['run_signature'] = signature
    root = args.work_dir or table_dir / 'ICAPE_climatology_work'
    run_dir = root / f'{tag}_{signature[:12]}'
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = run_dir / 'manifest.json'
    if manifest_path.exists():
        with open(manifest_path) as handle:
            if json.load(handle) != config:
                raise ValueError('Manifest mismatch; use a new --work-dir')
    else:
        atomic_json(manifest_path, config)
    # A nonblocking advisory lock prevents duplicate runs corrupting checkpoints.
    import fcntl
    with open(run_dir / 'run.lock', 'a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f'Another process is using {run_dir}') from exc
        print(f'{tag}: {len(mean_records)} cases; work={run_dir}', flush=True)
        summary = {str(k): int(v) for k, v in matches['status'].value_counts().items()}
        atomic_csv(run_dir / 'case_matching.csv', matches, index=False)
        print(f'{tag}: full-input case matching: {summary}; report={run_dir / "case_matching.csv"}', flush=True)
        if args.audit_only:
            return
        if not np.any(anomaly_order >= 0):
            raise ValueError('No matching anomaly rows in the selected sample; inspect case_matching.csv or increase --case-limit')
        print(f'Climatology: +/-{args.half_window} days, {args.years}, equal year weights; '
              f'require >= {args.min_valid_days} valid days/year and >= {args.min_valid_years} years', flush=True)
        if args.merge_only:
            active = np.load(run_dir / 'active.npy', allow_pickle=False)
        else:
            active = prepare_masks(run_dir, mean_records, original, args.box_root, tag)
        jobs = make_schedule(mean_records, active, args.years, args.half_window)
        (run_dir / 'daily').mkdir(exist_ok=True)
        pending = []
        for job in jobs:
            path = run_dir / 'daily' / f'{job[0]}.npz'
            if not path.exists() or args.overwrite_dates:
                pending.append(job)
            else:
                data = read_checkpoint(path, job, signature)
                if args.retry_failed_dates and np.any(data['status'] == ERROR):
                    pending.append(job)
        evaluations = sum(len(job[1]) for job in jobs)
        print(f'{tag}: {len(jobs)} unique background dates; {len(pending)} to process; '
              f'{evaluations:,} case-background-day evaluations in the full schedule', flush=True)
        if args.prepare_only:
            return
        selected_jobs = pending[:args.benchmark_dates] if args.benchmark_dates else pending
        if selected_jobs and not args.merge_only:
            # Import failures are fatal before launching workers, not thousands of failed checkpoints.
            init_worker(config, run_dir, mean_records)
            if args.validate_cases:
                validate_absolute(mean_records, original, args.validate_cases, run_dir)
            started = time.monotonic()
            if args.processes == 1:
                for job in selected_jobs:
                    print(process_day(job), flush=True)
            else:
                WORKER.clear()  # Release the parent's interpolation cache before spawning workers.
                with mp.get_context('spawn').Pool(args.processes, initializer=init_worker,
                                                 initargs=(config, run_dir, mean_records)) as pool:
                    for summary in pool.imap_unordered(process_day, selected_jobs, chunksize=1):
                        print(summary, flush=True)
            elapsed = time.monotonic() - started
            evaluations = sum(len(job[1]) for job in selected_jobs)
            print(f'{tag}: {evaluations} case-background-day evaluations in {elapsed / 60:.1f} min', flush=True)
        missing = [job[0] for job in jobs if not (run_dir / 'daily' / f'{job[0]}.npz').exists()]
        if missing:
            print(f'{tag}: {len(missing)} background dates still missing; anomaly-table merge deferred.', flush=True)
            if args.merge_only:
                raise RuntimeError('Cannot merge incomplete climatology checkpoints')
            return
        merge_results(args, config, run_dir, jobs, mean_records, original, anomaly, anomaly_order, transposed, output, matches)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--tags', nargs='+', choices=('hist', 'fut'), default=['hist', 'fut'])
    p.add_argument('--processes', type=int, default=4)
    p.add_argument('--raw-root', type=Path, default=RAW_ROOT)
    p.add_argument('--box-root', type=Path, default=BOX_ROOT)
    p.add_argument('--table-dir', type=Path)
    p.add_argument('--work-dir', type=Path)
    for tag in ('hist', 'fut'):
        p.add_argument(f'--{tag}-mean-table', type=Path, help='Absolute ICAPE-enriched mean table')
        p.add_argument(f'--{tag}-anom-table', type=Path, help='Existing environmental anomaly table')
        p.add_argument(f'--{tag}-output', type=Path)
    p.add_argument('--years', nargs='+', type=int, default=list(range(1998, 2007)))
    p.add_argument('--half-window', type=int, default=7)
    p.add_argument('--surface-time-mode', choices=('instant', '3hmean'), default='instant')
    p.add_argument('--min-valid-days', type=int, default=1, help='Required valid days per yearly window; legacy available-day mean=1')
    p.add_argument('--min-valid-years', type=int, help='Default: all requested years must contribute')
    p.add_argument('--case-limit', type=int, default=0, help='Test first N mean-table cases; uses separate sample outputs/checkpoints')
    p.add_argument('--validate-cases', type=int, default=3, help='Recompute first N usable storm cases to validate against absolute table (default 3)')
    p.add_argument('--match-tolerance', type=float, default=1e-5, help='Absolute tolerance in degrees per coordinate for stored numeric precision')
    p.add_argument('--strict-case-match', action='store_true', help='Require identical case lists; default retains unmatched anomaly rows as NaN')
    p.add_argument('--audit-only', action='store_true', help='Write case-matching report and stop; no raw or 2-D storm files opened')
    p.add_argument('--benchmark-dates', type=int, default=0, help='Process only N pending RAW background dates, not storm dates')
    p.add_argument('--prepare-only', action='store_true', help='Cache storm masks and print work estimate, without raw calculations')
    p.add_argument('--merge-only', action='store_true', help='Recreate outputs from a complete set of checkpoints')
    p.add_argument('--retry-failed-dates', action='store_true', help='Retry checkpoints with caught exceptions, not ordinary missing profiles')
    p.add_argument('--overwrite-dates', action='store_true', help='Recompute all selected background dates')
    p.add_argument('--overwrite-output', action='store_true', help='Replace generated outputs, never inputs')
    args = p.parse_args(argv)
    args.tags = list(dict.fromkeys(args.tags))
    args.years = sorted(set(args.years))
    if args.min_valid_years is None:
        args.min_valid_years = len(args.years)
    if args.processes < 1 or not 0 <= args.half_window <= 30 or min(args.case_limit, args.benchmark_dates, args.validate_cases) < 0:
        p.error('Need processes >= 1, half-window in 0..30 and nonnegative benchmark limits')
    if not 1 <= args.min_valid_days <= 2 * args.half_window + 1 or not 1 <= args.min_valid_years <= len(args.years):
        p.error('Invalid minimum background-day/year coverage')
    if not all(1 <= year <= 9999 for year in args.years):
        p.error('Reference years must be in 1..9999')
    if not np.isfinite(args.match_tolerance) or not 0 <= args.match_tolerance <= 1e-3:
        p.error('--match-tolerance must be in 0..0.001 degrees; do not use it to match different storms')
    if args.merge_only and (args.overwrite_dates or args.retry_failed_dates or args.prepare_only or args.benchmark_dates or args.audit_only):
        p.error('--merge-only cannot be combined with calculation/retry/prepare flags')
    return args


def main():
    args = parse_args()
    # Catch a missing future/historical input before spending days on the other tag.
    table_dir = args.table_dir or args.box_root / 'tables'
    for tag in args.tags:
        for kind, suffix in (('mean', '_ICAPE'), ('anom', '')):
            path = getattr(args, f'{tag}_{kind}_table') or table_dir / f'{tag}_{kind}_{TABLE_STEM}{suffix}.csv'
            if not path.is_file():
                raise FileNotFoundError(f'Required input not found: {path}; set --{tag}-{kind}-table if needed')
    for tag in args.tags:
        run_tag(args, tag)


if __name__ == '__main__':
    mp.freeze_support()
    main()

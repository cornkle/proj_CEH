"""Sample storm-centred CP4 time series for SM, sensible heat flux and rainfall.

Examples
--------
python LMCS_CP4_surfaceScales_timeseries_v2.py 17 hist
python LMCS_CP4_surfaceScales_timeseries_v2.py 17 future
python LMCS_CP4_surfaceScales_timeseries_v2.py 17 both
python LMCS_CP4_surfaceScales_timeseries_v2.py 17 both --max-storms 5
"""

from __future__ import annotations

import argparse
import glob
import os
import re
from collections import defaultdict
from datetime import timedelta
from pathlib import Path

import cftime
import numpy as np
import xarray as xr

from JASMIN import MetUM_variables as mu
time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)


CP4_ROOT = Path('/home/users/cornkle/linked_CP4')
STORM_ROOT = Path('/gws/ssde/j25b/lmcs/cklein/CP_models/MCS_files/WAf/CP4_box_JASMIN/mean3h_v2')
OUTPUT_ROOT = Path('/gws/ssde/j25b/lmcs/cklein/CP_models/MCS_files/WAf/CP4_box_JASMIN/mean3h_v2/saves/')

VARIABLES = {'SM': 'SM', 'sh': 'sh', 'lsRain': 'lsRain'}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('hour', type=int, help='Storm occurrence hour, e.g. 17')
    parser.add_argument('experiment', choices=['hist', 'future', 'both'], help='Historical, future, or both experiments')
    parser.add_argument('--days', type=int, default=10, help='Days before and after storm occurrence')
    parser.add_argument('--box-pixels', type=int, default=23, help='Odd number of native pixels per side; 23 is about 101 km')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing output files')
    parser.add_argument('--max-storms', type=int, default=None, help='Maximum number of storms per experiment for testing')
    return parser.parse_args()


def date_from_native_filename(path: str) -> cftime.Datetime360Day | None:
    """Read the first YYYYMMDDHHMM timestamp from a native CP4 filename."""
    match = re.search(r'_(\d{4})(\d{2})(\d{2})\d{4}-', os.path.basename(path))
    if match is None:
        return None

    year, month, day = map(int, match.groups())
    return cftime.Datetime360Day(year, month, day)


def event_from_storm_filename(path: str) -> cftime.Datetime360Day:
    """Read model-calendar storm time from YYYY-MM-DD_HH:MM:SS at filename start."""
    match = re.match(r'(\d{4})-(\d{2})-(\d{2})_(\d{2}):(\d{2}):(\d{2})', os.path.basename(path))
    if match is None:
        raise ValueError(f'Cannot parse storm time from {path}')

    year, month, day, hour, minute, second = map(int, match.groups())
    return cftime.Datetime360Day(year, month, day, hour, minute, second)


def date_key(date) -> tuple[int, int, int]:
    return date.year, date.month, date.day


def scalar_attr(ds: xr.Dataset, name: str) -> float:
    if name not in ds.attrs:
        raise KeyError(f'{name!r} missing from {ds.encoding.get("source", "storm file")}')
    return float(np.asarray(ds.attrs[name]).squeeze())


def index_native_files(cp4_dir: Path, var: str) -> dict[tuple[int, int, int], str]:
    """Build date-to-file lookup from filenames without opening the files."""
    lookup = {}

    for path in sorted(glob.glob(str(cp4_dir / var / '*.nc'))):
        date = date_from_native_filename(path)
        if date is not None:
            lookup[date_key(date)] = path

    if not lookup:
        raise FileNotFoundError(f'No dated files found in {cp4_dir / var}')

    return lookup


def shift_longitude(obj: xr.Dataset | xr.DataArray) -> xr.Dataset | xr.DataArray:
    if 'longitude' in obj.coords and float(obj.longitude.max()) > 180:
        obj = obj.assign_coords(longitude=obj.longitude - 360).sortby('longitude')
    return obj


def prepare_variable(da: xr.DataArray, var: str) -> xr.DataArray:
    da = da.rename(var)

    if var == 'SM':
        if 'depth' in da.coords:
            da = da.sel(depth=0.05, method='nearest').squeeze(drop=True)
        da = da.where(da < 500)

    elif var == 'lsRain':
        attrs = dict(da.attrs)
        da = da * 3600.0
        da.attrs = attrs
        da.attrs['units'] = 'mm h-1'
        da.attrs['conversion'] = 'Native precipitation rate multiplied by 3600'

    return da


def open_native_window(paths: list[str], var: str) -> xr.DataArray:
    native_name = mu.create_CP4_filename(var)

    def preprocess(ds: xr.Dataset) -> xr.Dataset:
        ds = shift_longitude(ds)
        if native_name not in ds:
            raise KeyError(f'{native_name!r} not found in {ds.encoding.get("source", "native file")}')
        return ds[[native_name]]

    ds = xr.open_mfdataset(paths, combine='by_coords', preprocess=preprocess, data_vars='minimal', coords='minimal',
                           compat='override', parallel=False, chunks={'time': 24}, decode_times=time_coder)

    return prepare_variable(ds[native_name], var)


def spatial_box(da: xr.DataArray, lon: float, lat: float, box_pixels: int) -> xr.DataArray:
    if box_pixels < 1 or box_pixels % 2 == 0:
        raise ValueError('box_pixels must be a positive odd integer')

    point = da.sel(longitude=lon, latitude=lat, method='nearest')
    xpos = int(np.abs(da.longitude.values - point.longitude.item()).argmin())
    ypos = int(np.abs(da.latitude.values - point.latitude.item()).argmin())
    half = box_pixels // 2

    box = da.isel(latitude=slice(ypos - half, ypos + half + 1),
                  longitude=slice(xpos - half, xpos + half + 1))

    if box.sizes.get('latitude') != box_pixels or box.sizes.get('longitude') != box_pixels:
        raise ValueError(f'{box_pixels}x{box_pixels} box crosses domain edge at lon={lon:.2f}, lat={lat:.2f}')

    return box


def sample_one_variable(da: xr.DataArray, event: cftime.Datetime360Day, lon: float, lat: float,
                        days: int, box_pixels: int, var: str) -> xr.Dataset:
    start, end = event - timedelta(days=days), event + timedelta(days=days)
    window = da.sel(time=slice(start, end))
    box = spatial_box(window, lon, lat, box_pixels)

    if box.sizes.get('time', 0) == 0:
        raise ValueError(f'No {var} times within ±{days} days of {event}')

    # One value per native model time step, averaged over the local ~100 x 100 km box.
    series = box.mean(('latitude', 'longitude'), skipna=True).load()
    
    hourly_mean = series.groupby('time.hour').mean('time', skipna=True)
    anom = series.groupby('time.hour') - hourly_mean
    
    lag = np.array([(time - event).total_seconds() / 3600 for time in series.time.values])
    
    time_dim, lag_name = f'time_{var}', f'lag_hour_{var}'
    series = series.rename({'time': time_dim}).assign_coords({lag_name: (time_dim, lag)})
    anom = anom.rename({'time': time_dim}).assign_coords({lag_name: (time_dim, lag)})

    series.name = var
    anom.name = f'{var}_anom'

    return xr.merge([series.to_dataset(), anom.to_dataset()], compat='override')


def run_experiment(experiment: str, hour: int, days: int, box_pixels: int,
                   overwrite: bool, max_storms: int | None = None) -> None:
    ftag = 'hist' if experiment == 'hist' else 'fut'
    cp4_dir = CP4_ROOT / experiment
    storm_dir = STORM_ROOT / f'mean_{ftag}'
    out_dir = OUTPUT_ROOT / f'SM_sh_lsRain_{ftag}_{hour:02d}UTC'
    out_dir.mkdir(parents=True, exist_ok=True)

    storm_files = sorted(glob.glob(str(storm_dir / f'*_{hour:02d}:*.nc')))

    if max_storms is not None:
        storm_files = storm_files[:max_storms]

    if not storm_files:
        raise FileNotFoundError(f'No {hour:02d} UTC storm files in {storm_dir}')

    print(f'{experiment}: processing {len(storm_files)} storms')

    storms_by_date = defaultdict(list)

    for path in storm_files:
        event = event_from_storm_filename(path)
        storms_by_date[date_key(event)].append((path, event))

    file_indexes = {outvar: index_native_files(cp4_dir, native_var)
                    for outvar, native_var in VARIABLES.items()}

    for storm_key, storms in sorted(storms_by_date.items()):
        year, month, day = storm_key
        centre_date = cftime.Datetime360Day(year, month, day)
        dates = [centre_date + timedelta(days=offset) for offset in range(-days, days + 1)]
        date_keys = [date_key(date) for date in dates]
        date_label = f'{year:04d}-{month:02d}-{day:02d}'

        native_data = {}

        for outvar, native_var in VARIABLES.items():
            paths = [file_indexes[outvar][key] for key in date_keys if key in file_indexes[outvar]]

            if not paths:
                print(f'No native {native_var} files for {date_label}')
                native_data[outvar] = None
                continue

            print(f'Opening {len(paths)} {native_var} files for {date_label}')
            native_data[outvar] = open_native_window(paths, native_var)

        for storm_file, event in storms:
            out_file = out_dir / f'{Path(storm_file).stem}_timeseries.nc'

            if out_file.exists() and not overwrite:
                print(f'Exists, skipping: {out_file.name}')
                continue

            try:
                with xr.open_dataset(storm_file, decode_times=time_coder) as storm_ds:
                    lon = scalar_attr(storm_ds, 'minlon')
                    lat = scalar_attr(storm_ds, 'minlat')
            except Exception as exc:
                print(f'Could not read storm location from {storm_file}: {exc}')
                continue

            output_parts = []

            for var, da in native_data.items():
                if da is None:
                    continue

                try:
                    output_parts.append(sample_one_variable(da, event, lon, lat, days, box_pixels, var))
                except Exception as exc:
                    print(f'Could not sample {var} for {Path(storm_file).name}: {exc}')

            if not output_parts:
                print(f'No variables sampled for {Path(storm_file).name}')
                continue

            output = xr.merge(output_parts, compat='override')
            output.attrs = {'storm_time': str(event), 'storm_lon': lon, 'storm_lat': lat,
                            'calendar': '360_day', 'box_pixels': box_pixels,
                            'approx_box_km': box_pixels * 4.4, 'window_days': days,
                            'experiment': experiment}

            encoding = {name: {'zlib': True, 'complevel': 5} for name in output.data_vars}
            output.to_netcdf(out_file, encoding=encoding)
            print(f'Saved {out_file}')

        for da in native_data.values():
            if da is not None:
                da.close()


def main() -> None:
    args = parse_args()
    experiments = ['hist', 'future'] if args.experiment == 'both' else [args.experiment]

    for experiment in experiments:
        run_experiment(experiment, args.hour, args.days, args.box_pixels, args.overwrite, args.max_storms)


if __name__ == '__main__':
    main()
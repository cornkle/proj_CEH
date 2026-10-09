#!/usr/bin/env python3
"""
Build a common 25-km daily JAS soil-moisture cache for CP4 and P25, then
calculate within-year temporal SM standard deviation.

Design choices
--------------
* Daily fixed-hour snapshots (default 12 UTC), matching the old analysis and
  avoiding a diurnal-sampling difference between models.
* July-September processed separately for each model year. Temporal std is
  calculated within each year, then years are averaged. No differences are
  ever taken across the discontinuous Sep -> next-Jul boundary.
* CP4 is area-averaged onto the native P25 grid using a sparse lookup matrix.
  The lookup is built once, saved, and reused.
* Daily cached 25-km fields are saved and reused by the wavelet/persistence
  production script, so CP4 never needs to be regridded again.

Outputs are written below:
    paths.cp4_2d / 'CP4_box_JASMIN/SM_model_comparison25'
"""

from __future__ import annotations

import argparse
import glob
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy import sparse

REPO_ROOT = Path('/home/users/cornkle/pythonWorkspace/proj_CEH')
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from shared.paths import load_paths

paths = load_paths()
time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

CP4_DIR = {
    'hist': Path('/home/users/cornkle/linked_CP4/hist/SM'),
    'future': Path('/home/users/cornkle/linked_CP4/future/SM'),
}
P25_DIR = {
    'hist': Path('/home/users/cornkle/linked_CP4/hist25/SM'),
    'future': Path('/home/users/cornkle/linked_CP4/future25/SM'),
}

OUT_ROOT = Path(paths.cp4_2d) / 'CP4_box_JASMIN/SM_model_comparison25'
CACHE_DIR = OUT_ROOT / 'daily_cache'
LOOKUP_DIR = OUT_ROOT / 'lookup'

MONTHS = (7, 8, 9)
HOUR = 12

# Cache a buffer larger than the final analysis belt so that the spatial
# wavelet calculation can later be done on the buffer and cropped inward.
BUFFER = {'lat0': 5.0, 'lat1': 19.0, 'lon0': -15.0, 'lon1': 25.0}
ANALYSIS = {'lat0': 8.0, 'lat1': 16.0, 'lon0': -12.0, 'lon1': 20.0}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hour', type=int, default=HOUR)
    p.add_argument('--overwrite-cache', action='store_true')
    p.add_argument('--rebuild-lookup', action='store_true')
    p.add_argument('--max-years', type=int, default=None, help='Testing only')
    return p.parse_args()


def file_date(path: str | Path):
    m = re.search(r'_(\d{4})(\d{2})(\d{2})\d{4}-', Path(path).name)
    if m is None:
        return None
    return tuple(map(int, m.groups()))


def files_by_year(root: Path):
    files = sorted(glob.glob(str(root / '**' / 'e08223*.nc'), recursive=True))
    out = {}
    for f in files:
        d = file_date(f)
        if d is None:
            continue
        year, month, _ = d
        if month in MONTHS:
            out.setdefault(year, []).append(f)
    return out


def find_sm_var(ds: xr.Dataset) -> str:
    for name in ('SM', 'sm', 'e08223'):
        if name in ds.data_vars:
            return name
    for name, da in ds.data_vars.items():
        text = ' '.join(str(da.attrs.get(k, '')) for k in ('standard_name', 'long_name', 'title', 'name')).lower()
        if 'soil' in text and 'moist' in text:
            return name
    raise KeyError(f'Could not identify soil-moisture variable. Available: {list(ds.data_vars)}')


def standardise_da(ds: xr.Dataset) -> xr.DataArray:
    da = ds[find_sm_var(ds)]

    # Top soil layer in the two file conventions.
    for z in ('depth', 'level6', 'level'):
        if z in da.dims:
            da = da.isel({z: 0}, drop=True)
            break

    if 't' in da.dims and 'time' not in da.dims:
        da = da.rename({'t': 'time'})
    if 't' in da.coords and 'time' not in da.coords:
        da = da.rename({'t': 'time'})

    if 'latitude' not in da.coords or 'longitude' not in da.coords:
        raise KeyError('Expected latitude/longitude coordinates')

    # Put longitudes on the familiar -180..180 convention and sort both axes.
    lon = ((da.longitude.astype(float) + 180.0) % 360.0) - 180.0
    da = da.assign_coords(longitude=lon).sortby('longitude').sortby('latitude')
    da = da.where(da < 500)  # same sanity filter used in the old CP4 analysis
    return da


def sel_box(da: xr.DataArray, box: dict) -> xr.DataArray:
    return da.sel(
        latitude=slice(box['lat0'], box['lat1']),
        longitude=slice(box['lon0'], box['lon1']),
    )


def read_grid(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    with xr.open_dataset(path, decode_times=time_coder) as ds:
        da = standardise_da(ds)
        return da.latitude.values.astype(float), da.longitude.values.astype(float)


def cell_edges(c: np.ndarray) -> np.ndarray:
    c = np.asarray(c, float)
    mid = 0.5 * (c[:-1] + c[1:])
    return np.r_[c[0] - (mid[0] - c[0]), mid, c[-1] + (c[-1] - mid[-1])]


def prepare_target_grid(sample_p25: str | Path):
    lat, lon = read_grid(sample_p25)
    ilat = (lat >= BUFFER['lat0']) & (lat <= BUFFER['lat1'])
    ilon = (lon >= BUFFER['lon0']) & (lon <= BUFFER['lon1'])
    return lat[ilat], lon[ilon]


def build_lookup(sample_cp4: str | Path, tgt_lat: np.ndarray, tgt_lon: np.ndarray):
    """Build area-weighted CP4 -> P25 cell-average sparse matrix."""
    src_lat_full, src_lon_full = read_grid(sample_cp4)
    lat_edges, lon_edges = cell_edges(tgt_lat), cell_edges(tgt_lon)

    iy = np.where((src_lat_full >= lat_edges[0]) & (src_lat_full < lat_edges[-1]))[0]
    ix = np.where((src_lon_full >= lon_edges[0]) & (src_lon_full < lon_edges[-1]))[0]
    src_lat, src_lon = src_lat_full[iy], src_lon_full[ix]

    ybin = np.searchsorted(lat_edges, src_lat, side='right') - 1
    xbin = np.searchsorted(lon_edges, src_lon, side='right') - 1

    YY = np.broadcast_to(ybin[:, None], (len(src_lat), len(src_lon)))
    XX = np.broadcast_to(xbin[None, :], (len(src_lat), len(src_lon)))
    valid = (YY >= 0) & (YY < len(tgt_lat)) & (XX >= 0) & (XX < len(tgt_lon))

    rows = (YY * len(tgt_lon) + XX)[valid].ravel()
    cols = np.arange(len(src_lat) * len(src_lon)).reshape(len(src_lat), len(src_lon))[valid].ravel()
    # Approximate grid-cell area weighting on a lat/lon grid.
    w2d = np.broadcast_to(np.cos(np.deg2rad(src_lat))[:, None], (len(src_lat), len(src_lon)))
    weights = w2d[valid].ravel()

    W = sparse.coo_matrix(
        (weights, (rows, cols)),
        shape=(len(tgt_lat) * len(tgt_lon), len(src_lat) * len(src_lon)),
    ).tocsr()
    rowsum = np.asarray(W.sum(axis=1)).ravel()
    if np.any(rowsum == 0):
        bad = np.where(rowsum == 0)[0]
        raise RuntimeError(f'{len(bad)} target cells have no CP4 source pixels')
    W = sparse.diags(1.0 / rowsum) @ W
    return W.tocsr(), src_lat, src_lon


def lookup_paths():
    return LOOKUP_DIR / 'CP4_to_P25_JAS_buffer_weights.npz', LOOKUP_DIR / 'CP4_to_P25_JAS_buffer_coords.npz'


def get_lookup(sample_cp4, tgt_lat, tgt_lon, rebuild=False):
    LOOKUP_DIR.mkdir(parents=True, exist_ok=True)
    wfile, cfile = lookup_paths()

    if wfile.exists() and cfile.exists() and not rebuild:
        W = sparse.load_npz(wfile)
        meta = np.load(cfile)
        if np.allclose(meta['tgt_lat'], tgt_lat) and np.allclose(meta['tgt_lon'], tgt_lon):
            print('Using saved CP4 -> P25 lookup:', wfile)
            return W, meta['src_lat'], meta['src_lon']
        print('Saved lookup target grid differs; rebuilding.')

    print('Building CP4 -> P25 lookup once...')
    W, src_lat, src_lon = build_lookup(sample_cp4, tgt_lat, tgt_lon)
    sparse.save_npz(wfile, W)
    np.savez_compressed(cfile, src_lat=src_lat, src_lon=src_lon, tgt_lat=tgt_lat, tgt_lon=tgt_lon)
    print('Saved lookup:', wfile)
    return W, src_lat, src_lon


def select_hour(da: xr.DataArray, hour: int) -> xr.DataArray:
    if 'time' not in da.coords:
        raise KeyError('No decoded time coordinate')
    return da.where(da.time.dt.hour == hour, drop=True)


def read_year_native(files: list[str], hour: int, lat=None, lon=None, box=None) -> xr.DataArray:
    """Read only the chosen hour from each file; concatenate within one year."""
    parts = []
    for f in files:
        with xr.open_dataset(f, decode_times=time_coder) as ds:
            da = standardise_da(ds)
            if lat is not None and lon is not None:
                da = da.sel(latitude=lat, longitude=lon)
            elif box is not None:
                da = sel_box(da, box)
            da = select_hour(da, hour).load().astype('float32')
            if da.sizes.get('time', 0):
                parts.append(da)

    if not parts:
        raise RuntimeError('No selected-hour data found')

    out = xr.concat(parts, dim='time').sortby('time')
    # Remove any duplicated boundary timestamp safely.
    labels = np.array([str(t) for t in out.time.values])
    _, keep = np.unique(labels, return_index=True)
    return out.isel(time=np.sort(keep))


def regrid_cp4(da: xr.DataArray, W, tgt_lat, tgt_lon) -> xr.DataArray:
    arr = np.asarray(da.values, dtype='float32')
    nt = arr.shape[0]
    flat = arr.reshape(nt, -1)
    good = np.isfinite(flat)
    num = (W @ np.nan_to_num(flat, nan=0.0).T).T
    den = (W @ good.astype('float32').T).T
    with np.errstate(invalid='ignore', divide='ignore'):
        out = num / den
    out[den == 0] = np.nan
    out = out.reshape(nt, len(tgt_lat), len(tgt_lon)).astype('float32')
    return xr.DataArray(
        out,
        coords={'time': da.time.values, 'latitude': tgt_lat, 'longitude': tgt_lon},
        dims=('time', 'latitude', 'longitude'),
        name='SM',
        attrs={'units': 'kg m-2', 'grid': 'P25 native grid; CP4 area-averaged where applicable'},
    )


def direct_p25_to_target(da, tgt_lat, tgt_lon):
    """
    Put native P25 data onto the common P25 target grid.

    HIST/FUT P25 files are nominally on the same ~25-km grid, but their
    coordinate values need not be bit-for-bit identical. If the coordinates
    coincide to numerical precision, use direct nearest indexing. If the
    grids are genuinely offset, interpolate onto the common target grid.
    """

    src_lat = np.asarray(da.latitude.values, dtype=float)
    src_lon = np.asarray(da.longitude.values, dtype=float)

    # Nearest native coordinate corresponding to every target coordinate.
    iy = np.abs(src_lat[:, None] - tgt_lat[None, :]).argmin(axis=0)
    ix = np.abs(src_lon[:, None] - tgt_lon[None, :]).argmin(axis=0)

    lat_err = np.abs(src_lat[iy] - tgt_lat)
    lon_err = np.abs(src_lon[ix] - tgt_lon)

    dlat = np.nanmedian(np.abs(np.diff(tgt_lat)))
    dlon = np.nanmedian(np.abs(np.diff(tgt_lon)))

    print(
        f"P25 grid check: max latitude mismatch={lat_err.max():.6f}°, "
        f"max longitude mismatch={lon_err.max():.6f}°"
    )

    # Tiny differences are just coordinate precision: no interpolation needed.
    if lat_err.max() < 0.01 * dlat and lon_err.max() < 0.01 * dlon:
        out = da.isel(
            latitude=xr.DataArray(iy, dims="latitude"),
            longitude=xr.DataArray(ix, dims="longitude"),
        )

        # Make coordinates exactly identical across all four experiments.
        out = out.assign_coords(
            latitude=tgt_lat,
            longitude=tgt_lon,
        )

    else:
        # The two P25 integrations genuinely use slightly different grids.
        # Interpolation is cheap here because these are already 25-km fields.
        print("P25 grid is offset: interpolating onto common target grid.")

        out = da.interp(
            latitude=tgt_lat,
            longitude=tgt_lon,
            method="linear",
        )

    return out.rename("SM").astype("float32")


def cache_path(tag: str, year: int, hour: int):
    return CACHE_DIR / f'SM_daily{hour:02d}_JAS_{tag}_{year}_onP25.nc'


def build_cache(tag, model, climate, year, files, hour, W, src_lat, src_lon, tgt_lat, tgt_lon, overwrite=False):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    outfile = cache_path(tag, year, hour)
    if outfile.exists() and not overwrite:
        print('Cache exists:', outfile.name)
        return outfile

    print(f'Building {tag} {year}: {len(files)} input files')
    if model == 'CP4':
        native = read_year_native(files, hour, lat=src_lat, lon=src_lon)
        da = regrid_cp4(native, W, tgt_lat, tgt_lon)
    else:
        native = read_year_native(files, hour, box=BUFFER)
        da = direct_p25_to_target(native, tgt_lat, tgt_lon)

    da.attrs.update({'model': model, 'climate': climate, 'year': year, 'months': '7,8,9', 'hour_utc': hour})
    da.to_dataset().to_netcdf(outfile, encoding={'SM': {'zlib': True, 'complevel': 4, 'dtype': 'float32'}})
    print('Saved:', outfile)
    return outfile


def main():
    args = parse_args()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    cp4_files = {c: files_by_year(CP4_DIR[c]) for c in ('hist', 'future')}
    p25_files = {c: files_by_year(P25_DIR[c]) for c in ('hist', 'future')}

    sample_p25 = p25_files['hist'][sorted(p25_files['hist'])[0]][0]
    sample_cp4 = cp4_files['hist'][sorted(cp4_files['hist'])[0]][0]
    tgt_lat, tgt_lon = prepare_target_grid(sample_p25)
    W, src_lat, src_lon = get_lookup(sample_cp4, tgt_lat, tgt_lon, rebuild=args.rebuild_lookup)

    configs = [
        ('CP4hist', 'CP4', 'hist', cp4_files['hist']),
        ('CP4fut', 'CP4', 'future', cp4_files['future']),
        ('P25hist', 'P25', 'hist', p25_files['hist']),
        ('P25fut', 'P25', 'future', p25_files['future']),
    ]

    std_by_tag = []
    summaries = []

    for tag, model, climate, grouped in configs:
        years = sorted(grouped)
        if args.max_years is not None:
            years = years[:args.max_years]
        yearly = []

        for year in years:
            cfile = build_cache(
                tag, model, climate, year, grouped[year], args.hour,
                W, src_lat, src_lon, tgt_lat, tgt_lon,
                overwrite=args.overwrite_cache,
            )
            with xr.open_dataset(cfile) as ds:
                da = sel_box(ds['SM'], ANALYSIS)
                # Crucial: std only within the single JAS season/year.
                st = da.std('time', skipna=True).load().astype('float32')
            st = st.expand_dims(year=[year])
            yearly.append(st)

        if not yearly:
            continue

        ys = xr.concat(yearly, dim='year').rename('SM_std')
        ys = ys.expand_dims(experiment=[tag])
        std_by_tag.append(ys)

        vals = ys.mean('year', skipna=True).values.ravel()
        vals = vals[np.isfinite(vals)]
        summaries.append({
            'experiment': tag,
            'n_years': len(years),
            'mean_of_yearmean_std_map': float(np.mean(vals)),
            'median_of_yearmean_std_map': float(np.median(vals)),
            'p25': float(np.percentile(vals, 25)),
            'p75': float(np.percentile(vals, 75)),
        })

    std_all = xr.concat(std_by_tag, dim='experiment', join='outer')
    dsout = xr.Dataset({
        'SM_std_yearly': std_all,
        'SM_std_mean_years': std_all.mean('year', skipna=True),
    })
    dsout.attrs.update({
        'description': 'JAS daily fixed-hour temporal SM std; calculated within each year before multi-year averaging',
        'hour_utc': args.hour,
        'analysis_box': str(ANALYSIS),
        'grid': 'common native P25 grid; CP4 area-averaged with saved sparse lookup',
    })

    outfile = OUT_ROOT / f'SM_temporal_std_JAS{args.hour:02d}_onP25.nc'
    dsout.to_netcdf(outfile, encoding={k: {'zlib': True, 'complevel': 4} for k in dsout.data_vars})
    pd.DataFrame(summaries).to_csv(OUT_ROOT / f'SM_temporal_std_JAS{args.hour:02d}_summary.csv', index=False)

    print('\nWrote:')
    print(outfile)
    print(OUT_ROOT / f'SM_temporal_std_JAS{args.hour:02d}_summary.csv')
    print('\nThe daily common-grid caches are in:')
    print(CACHE_DIR)


if __name__ == '__main__':
    main()

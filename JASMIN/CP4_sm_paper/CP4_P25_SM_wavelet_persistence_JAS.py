#!/usr/bin/env python3
"""
Spatial wavelet power + scale composition + temporal persistence of SM patterns
for the common 25-km JAS daily cache created by CP4_P25_SM_stddev_JAS.py.

The analysis deliberately stays simple:
  1) same 2-D Mexican-hat wavelet package already used in the paper;
  2) FUT-HIST change in mean power at each spatial scale = pattern-amplitude change;
  3) FUT-HIST change in power fraction = redistribution among spatial scales;
  4) autocorrelation of spatial wavelet coefficient maps = persistence of the
     actual spatial pattern at each scale (not merely persistence of total power).

Persistence is calculated within each JAS/year only, so no correlation crosses
Sep -> next-Jul discontinuities.
"""

from __future__ import annotations

import argparse
import glob
import re
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.ndimage import distance_transform_edt

REPO_ROOT = Path('/home/users/cornkle/pythonWorkspace/proj_CEH')
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from shared.paths import load_paths
from land_wavelet import wav

paths = load_paths()

OUT_ROOT = Path(paths.cp4_2d) / 'CP4_box_JASMIN/SM_model_comparison25'
CACHE_DIR = OUT_ROOT / 'daily_cache'

ANALYSIS = {'lat0': 8.0, 'lat1': 16.0, 'lon0': -12.0, 'lon1': 20.0}
DX_KM = 25.0

# At 25-km resolution, do not interpret sub-50-km scales.  These settings
# produce approximately 50, 71, 100, 141, 200, 283, 400 km.
WAVELET_START_KM = 50.0
WAVELET_DJ = 0.5
WAVELET_NB = 6
MAX_LAG_DAYS = 10


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hour', type=int, default=12)
    p.add_argument('--max-lag-days', type=int, default=MAX_LAG_DAYS)
    p.add_argument('--max-years', type=int, default=None, help='Testing only')
    return p.parse_args()


def cache_files(tag, hour):
    return sorted(glob.glob(str(CACHE_DIR / f'SM_daily{hour:02d}_JAS_{tag}_*_onP25.nc')))


def year_from_cache(path):
    m = re.search(r'_([0-9]{4})_onP25\.nc$', Path(path).name)
    if m is None:
        raise ValueError(f'Cannot parse year from {path}')
    return int(m.group(1))


def analysis_mask(lat, lon):
    yy = (lat >= ANALYSIS['lat0']) & (lat <= ANALYSIS['lat1'])
    xx = (lon >= ANALYSIS['lon0']) & (lon <= ANALYSIS['lon1'])
    return yy[:, None] & xx[None, :]


def common_land_mask(files_by_tag):
    """Intersection of cells that contain SM in every model/climate sample."""
    mask = None
    for tag, files in files_by_tag.items():
        if not files:
            raise FileNotFoundError(f'No cache files for {tag}')
        with xr.open_dataset(files[0]) as ds:
            m = np.isfinite(ds['SM'].isel(time=0).values)
        mask = m if mask is None else (mask & m)
    return mask


def scale_masks(land, base_analysis, scales):
    """
    Avoid averaging wavelet power right against missing/ocean boundaries.
    The transform itself uses zero anomaly over missing cells; each scale is
    averaged only where the centre is at least ~half that real scale inland.
    """
    dist_pix = distance_transform_edt(land)
    out = []
    for s in scales:
        radius_pix = max(1.0, (s / 2.0) / DX_KM)
        out.append(base_analysis & land & (dist_pix >= radius_pix))
    return out


def spatial_corr(a, b, mask):
    x = a[mask].astype(float)
    y = b[mask].astype(float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 20:
        return np.nan
    x, y = x[ok], y[ok]
    x -= x.mean(); y -= y.mean()
    den = np.sqrt(np.sum(x*x) * np.sum(y*y))
    return float(np.sum(x*y) / den) if den > 0 else np.nan


def acf_coeff_maps(coeff, masks, max_lag):
    """Mean spatial-pattern correlation for coefficient maps separated by lag days."""
    nt, ns = coeff.shape[:2]
    lags = np.arange(max_lag + 1)
    out = np.full((ns, len(lags)), np.nan, dtype=float)
    out[:, 0] = 1.0

    for s in range(ns):
        mask = masks[s]
        for lag in range(1, max_lag + 1):
            if lag >= nt:
                continue
            vals = [spatial_corr(coeff[t, s], coeff[t + lag, s], mask) for t in range(nt - lag)]
            out[s, lag] = np.nanmean(vals) if np.any(np.isfinite(vals)) else np.nan
    return out


def decorrelation_time(acf, lags):
    target = 1.0 / np.e
    out = np.full(acf.shape[0], np.nan)
    for s in range(acf.shape[0]):
        y = acf[s]
        good = np.isfinite(y)
        x, yy = lags[good], y[good]
        if len(x) < 2:
            continue
        hit = np.where(yy <= target)[0]
        if len(hit) == 0:
            continue
        j = hit[0]
        if j == 0:
            out[s] = x[0]
        else:
            x0, x1 = x[j-1], x[j]
            y0, y1 = yy[j-1], yy[j]
            out[s] = x0 + (target - y0) * (x1 - x0) / (y1 - y0) if y1 != y0 else x1
    return out


def process_year(path, wobj, masks, max_lag):
    with xr.open_dataset(path) as ds:
        da = ds['SM'].load()
        lat = da.latitude.values.astype(float)
        lon = da.longitude.values.astype(float)
        arr = da.values.astype('float32')

    nt = arr.shape[0]
    ns = len(wobj.scales)
    coeff_store = np.empty((nt, ns, arr.shape[1], arr.shape[2]), dtype='float32')
    power_ts = np.full((nt, ns), np.nan, dtype='float64')

    for t in range(nt):
        field = arr[t].astype(float)
        good = np.isfinite(field)
        if good.sum() == 0:
            continue

        # Important correction for SM: subtract the valid-domain mean first,
        # then set missing/ocean anomalies to zero.  This avoids turning ocean
        # points into -mean anomalies.
        anom = field - np.nanmean(field[good])
        anom[~good] = 0.0

        coeff, power = wobj.calc_coeffs(anom, power_normed='scale')
        coeff_store[t] = coeff.astype('float32')

        for s, mask in enumerate(masks):
            power_ts[t, s] = np.nanmean(power[s][mask]) if np.any(mask) else np.nan

    annual_power = np.nanmean(power_ts, axis=0)
    annual_fraction = annual_power / np.nansum(annual_power)
    acf = acf_coeff_maps(coeff_store, masks, max_lag)
    lags = np.arange(max_lag + 1, dtype=float)
    tau = decorrelation_time(acf, lags)

    return annual_power, annual_fraction, acf, tau


def main():
    args = parse_args()

    tags = ['CP4hist', 'CP4fut', 'P25hist', 'P25fut']
    files_by_tag = {tag: cache_files(tag, args.hour) for tag in tags}

    for tag, fs in files_by_tag.items():
        print(tag, len(fs), 'cached years')
        if not fs:
            raise FileNotFoundError(
                f'No {tag} daily cache. Run CP4_P25_SM_stddev_JAS.py first.'
            )

    with xr.open_dataset(files_by_tag['P25hist'][0]) as ds:
        lat = ds.latitude.values.astype(float)
        lon = ds.longitude.values.astype(float)

    land = common_land_mask(files_by_tag)
    amask = analysis_mask(lat, lon)

    # Use the existing 2-D Mexican-hat implementation, but start at 2*dx.
    wobj = wav.wavelet(
        res=DX_KM,
        dist=WAVELET_DJ,
        nb=WAVELET_NB,
        start=WAVELET_START_KM,
    )
    scales = np.asarray(wobj.scales, float)
    masks = scale_masks(land, amask, scales)
    lags = np.arange(args.max_lag_days + 1, dtype=float)

    print('Wavelet scales (km):', np.round(scales, 2))
    print('Valid analysis pixels by scale:', [int(m.sum()) for m in masks])

    outputs = []

    for tag in tags:
        fs = files_by_tag[tag]
        if args.max_years is not None:
            fs = fs[:args.max_years]

        years, powers, fractions, acfs, taus = [], [], [], [], []
        for f in fs:
            year = year_from_cache(f)
            print('Processing', tag, year)
            power, frac, acf, tau = process_year(f, wobj, masks, args.max_lag_days)
            years.append(year)
            powers.append(power)
            fractions.append(frac)
            acfs.append(acf)
            taus.append(tau)

        ds_tag = xr.Dataset(
            {
                'power_yearly': (('year', 'scale'), np.asarray(powers)),
                'power_fraction_yearly': (('year', 'scale'), np.asarray(fractions)),
                'pattern_acf_yearly': (('year', 'scale', 'lag_day'), np.asarray(acfs)),
                'decorrelation_days_yearly': (('year', 'scale'), np.asarray(taus)),
            },
            coords={
                'year': years,
                'scale': scales,
                'lag_day': lags,
            },
        ).expand_dims(experiment=[tag])
        outputs.append(ds_tag)

    out = xr.concat(outputs, dim='experiment', join='outer')
    out['power_mean_years'] = out['power_yearly'].mean('year', skipna=True)
    out['power_fraction_mean_years'] = out['power_fraction_yearly'].mean('year', skipna=True)
    out['pattern_acf_mean_years'] = out['pattern_acf_yearly'].mean('year', skipna=True)
    out['decorrelation_days_mean_years'] = out['decorrelation_days_yearly'].mean('year', skipna=True)

    out.attrs.update({
        'description': 'JAS daily fixed-hour spatial Mexican-hat SM wavelet power and coefficient-map persistence',
        'hour_utc': args.hour,
        'grid_km': DX_KM,
        'wavelet_start_km': WAVELET_START_KM,
        'wavelet_scale_spacing_log2': WAVELET_DJ,
        'analysis_box': str(ANALYSIS),
        'power_normalisation': 'wavelet power divided by wavelet scale (power_normed=scale)',
        'missing_data_handling': 'subtract valid spatial mean, then set missing anomaly to zero; scale-dependent inland mask used for spatial averages/correlations',
    })

    outfile = OUT_ROOT / f'SM_wavelet_persistence_JAS{args.hour:02d}_onP25.nc'
    out.to_netcdf(outfile, encoding={k: {'zlib': True, 'complevel': 4} for k in out.data_vars})
    print('\nWrote:', outfile)


if __name__ == '__main__':
    main()

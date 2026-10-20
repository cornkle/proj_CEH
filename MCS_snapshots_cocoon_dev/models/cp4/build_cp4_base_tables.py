#!/usr/bin/env python3
"""CP4 hourly COCOON storm base tables; no writes to OBS paths.

Read paired monthly HadREM-CP4A rlut/pr files and generate incremental,
all-Africa storm catalogues with OBS-compatible metric names.

Uses the physical thresholds and geometric/precipitation statistic definitions
of COCOON_table_create_01deg.py but processes *cropped storm objects* to avoid
rescanning a 2100x2000 array for every object. Extra compactness is intentionally
NOT calculated because OBS does not calculate it.

Dependencies: numpy, pandas, xarray, netCDF4 or h5netcdf, cftime, scipy,
scikit-image, shapely, pyproj, pyarrow.
"""
from __future__ import annotations
import argparse
import logging
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.ndimage import label, find_objects
from shapely.geometry import Polygon
from skimage import measure
import pyproj

LOG = logging.getLogger('cp4_cocoon')
INPUT_ROOT = Path('/gws/pw/j25/cp4a')
OUTPUT_ROOT = Path('/gws/ssde/j25b/cocoon/MCS_5000km2_tables_v3/CP4_base_tables_tir-prcp')
MODELS = ('CanESM5', 'HadGEM3-GC31-MM', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'TaiESM1')
OBS_COLUMNS = [
    'date','month','hour','minute','year','day','area','area70','area60','area50',
    'minlon','minlat','maxlon','maxlat','clon','clat','tminlon','tminlat',
    'tmin','tmean','tp1','tp99','stormID','shape_complexity',
    'major_axis_km','minor_axis_km','aspect_ratio',
    'precipitation_mean','precipitation_max','precipitation_p80',
    'precipitation_p90','precipitation_p95','precipitation_p99',
    'precipitation_area1mm','precipitation_area3mm','precipitation_area8mm',
]
EXTRA_COLUMNS = [
    'storm_id','run','model','experiment','calendar',
    'storm_time','storm_hour_utc','storm_hour_end_utc','local_solar_hour',
    'pixel_count','pr_source','rlut_source',
]
SCHEMA = OBS_COLUMNS + EXTRA_COLUMNS
SIGMA=5.670373e-8

def olr_to_bt_kelvin(olr):
    """Use Yang & Slingo OLR -> brightness temperature, as in legacy CP4.

    Input W/m2; output kelvin, then convert to Celsius for cloud objects.
    """
    tf = (np.maximum(olr, 0) / SIGMA) ** 0.25
    a,b = 1.228,-1.106e-3
    disc = a*a+4*b*tf
    out=np.full(tf.shape, np.nan, dtype='float32')
    valid=np.isfinite(olr) & (olr>0) & (disc>=0)
    out[valid] = ((-a+np.sqrt(disc[valid]))/(2*b)).astype('float32')
    return out

def pixel_area_latlon_km2(lat, lon):
    """Match observational COCOON_table_main_01deg.py cell-area formula."""
    R=6371.0
    dlat=np.abs(np.gradient(lat))
    dlon=np.abs(np.gradient(lon))
    return (R**2 * np.deg2rad(dlat)[:,None]
            * np.deg2rad(dlon)[None,:]
            * np.cos(np.deg2rad(lat))[:,None]).astype('float32')

def polygon_from_cropped_mask(mask, lat, lon, transformer):
    """OBS definition: largest valid find_contours polygon, equal area CRS.

    Crop is padded by a false-pixel border when possible. For objects on the
    physical model domain edge, this closes the contour, unlike the OBS
    full-array contour behaviour; mark such edge cases separately if relevant.
    """
    contours=measure.find_contours(mask.astype(float), 0.5)
    polygons=[]
    for contour in contours:
        rows,cols=contour[:,0],contour[:,1]
        clat=np.interp(rows,np.arange(len(lat)),lat)
        clon=np.interp(cols,np.arange(len(lon)),lon)
        x,y=transformer.transform(clon,clat)
        poly=Polygon(zip(x,y))
        if poly.is_valid and poly.area>0:
            polygons.append(poly)
    return max(polygons,key=lambda p:p.area) if polygons else None

def geometry_metrics(poly):
    if poly is None or poly.is_empty:
        return (np.nan,)*4
    hull=poly.convex_hull.area
    complexity=(1-poly.area/hull) if hull>0 else np.nan
    try:
        rect=poly.minimum_rotated_rectangle
        xy=np.asarray(rect.exterior.coords)
        d=np.sqrt(np.sum(np.diff(xy,axis=0)**2,axis=1))
        major,minor=float(np.max(d))/1000,float(np.min(d))/1000
        aspect=major/minor if minor>0 else np.nan
    except (AttributeError, ValueError):
        major=minor=aspect=np.nan
    return complexity,major,minor,aspect

def _timestamp_info(t):
    """Handle cftime 360_day and 365_day without pandas datetime coercion."""
    y,m,d,h,mi = (int(t.year),int(t.month),int(t.day),int(t.hour),int(t.minute))
    return y,m,d,h,mi,f'{y:04d}-{m:02d}-{d:02d}_{h:02d}:{mi:02d}'

def _slice_for_interval(axis, centre, radius):
    # Existing OBS uses coordinate slices rather than a nearest-pixel average.
    return (axis >= centre-radius) & (axis <= centre+radius)

def extract_storms(tb_c, rain_mmh, lat, lon, area, time, run,
                   transformer, pr_file, rlut_file,
                   min_area=5000,max_area=300000,min_t=-50,min_rain=1,
                   max_t=-40):
    """One time slice; all obs names and metrics; return list of dicts."""
    cold = np.isfinite(tb_c) & (tb_c < max_t)
    labels,nlabel = label(cold) # scipy.ndimage default 4-neighbour connectivity
    if nlabel==0:
        return []
    # Remove candidate objects using true physical area, before expensive metrics.
    areas=np.bincount(labels.ravel(),weights=area.ravel(),minlength=nlabel+1)
    keep=np.nonzero((areas>=min_area)&(areas<=max_area))[0]
    keep=keep[keep>0]
    if keep.size==0:
        return []
    bounds=find_objects(labels)
    yy,mm,dd,hh,mn,datestr=_timestamp_info(time)
    year_model=run.rsplit('-',1)[0]
    exp=run.rsplit('-',1)[-1]
    cal=getattr(time,'calendar','standard')
    rows=[]
    # Data may have ascending or descending axis, but no axis flipping required
    # for the counting/statistical metrics.
    for idx in keep:
        region=bounds[int(idx)-1]
        if region is None:
            continue
        ys,xs=region
        mask=(labels[region]==idx)
        ctt=tb_c[region]
        rain=rain_mmh[region]
        pixel_area=area[region]
        rvals=rain[mask]
        rvals=rvals[np.isfinite(rvals)]
        if not rvals.size or float(np.max(rvals))<min_rain:
            continue
        tvals=ctt[mask]
        tmin=float(np.min(tvals))
        if tmin>min_t:
            continue
        cell_lat=lat[ys]
        cell_lon=lon[xs]
        rr,cc=np.where(mask)
        ir,ic=np.unravel_index(int(np.argmin(np.where(mask,ctt,np.inf))),mask.shape)
        tmin_lat=float(cell_lat[ir]);tmin_lon=float(cell_lon[ic])
        rain_2d=np.where(mask,rain,np.nan)
        # OBS: find maximum precipitation pixel inside cold cloud, then average
        # original masked rain field in a +/- 0.075 degree box around it.
        max_rain_index=int(np.nanargmax(rain_2d))
        rmax_y,rmax_x=np.unravel_index(max_rain_index,mask.shape)
        lon_r=float(cell_lon[rmax_x]);lat_r=float(cell_lat[rmax_y])
        p_lat=_slice_for_interval(cell_lat,lat_r,0.075)
        p_lon=_slice_for_interval(cell_lon,lon_r,0.075)
        selected=rain_2d[np.ix_(p_lat,p_lon)]
        pmax=float(np.nanmean(selected)) if np.isfinite(selected).any() else np.nan
        # The reference OBS contour is generated using the full native mask.
        # To reduce full-domain scans, make a padded local bounding mask.
        y0=max(0,ys.start-1);y1=min(len(lat),ys.stop+1)
        x0=max(0,xs.start-1);x1=min(len(lon),xs.stop+1)
        localmask=(labels[y0:y1,x0:x1]==idx)
        # Match OBS: failures in contour geometry should not discard a storm.
        try:
            poly=polygon_from_cropped_mask(localmask,lat[y0:y1],lon[x0:x1],transformer)
            sci,major,minor,aspect=geometry_metrics(poly)
        except Exception as exc:
            LOG.warning('Shape metric failed for %s %s label %d: %s',run,datestr,int(idx),exc)
            sci=major=minor=aspect=np.nan
        row={
          'date':datestr,'year':yy,'month':mm,'day':dd,'hour':hh,'minute':mn,
          'area':float(np.sum(pixel_area[mask])),
          'area50':float(np.sum(pixel_area[mask & (ctt<=-50)])),
          'area60':float(np.sum(pixel_area[mask & (ctt<=-60)])),
          'area70':float(np.sum(pixel_area[mask & (ctt<=-70)])),
          'minlon':float(cell_lon[cc].min()),'maxlon':float(cell_lon[cc].max()),
          'minlat':float(cell_lat[rr].min()),'maxlat':float(cell_lat[rr].max()),
          'tminlat':tmin_lat,'tminlon':tmin_lon,'tmin':tmin,
          'tmean':float(np.mean(tvals)), 'tp1':float(np.percentile(tvals,1)),
          'tp99':float(np.percentile(tvals,99)),
          'stormID':f'{datestr}_{idx}',
          'shape_complexity':sci,'major_axis_km':major,
          'minor_axis_km':minor,'aspect_ratio':aspect,
          'precipitation_mean':float(np.mean(rvals)),
          'precipitation_max':pmax,
          'precipitation_p80':float(np.percentile(rvals,80)),
          'precipitation_p90':float(np.percentile(rvals,90)),
          'precipitation_p95':float(np.percentile(rvals,95)),
          'precipitation_p99':float(np.percentile(rvals,99)),
          'precipitation_area1mm':float(np.sum(pixel_area[mask & (rain>=1)])),
          'precipitation_area3mm':float(np.sum(pixel_area[mask & (rain>=3)])),
          'precipitation_area8mm':float(np.sum(pixel_area[mask & (rain>=8)])),
          'run':run,'model':year_model,'experiment':exp,'calendar':cal,
          'storm_time':str(time),'storm_hour_utc':hh+mn/60,
          'storm_hour_end_utc':(hh+mn/60+0.5)%24,
          'local_solar_hour':(hh+mn/60+tmin_lon/15)%24,
          'pixel_count':int(mask.sum()),
          'pr_source':Path(pr_file).name,'rlut_source':Path(rlut_file).name,
        }
        row['clon']=(row['minlon']+row['maxlon'])/2
        row['clat']=(row['minlat']+row['maxlat'])/2
        row['storm_id']=f'{run}_{datestr}_{idx}'
        rows.append(row)
    return rows

def find_month_file(run_dir, var, year, month):
    path=run_dir/var
    # Some files omit -FUT/-HIST from basename. Use parent run as identity.
    candidates=sorted(path.glob(f'{var}_*_{year:04d}{month:02d}*.nc'))
    if len(candidates)!=1:
        raise FileNotFoundError(
            f'Expected one {var} NetCDF in {path} for {year:04d}-{month:02d}; '
            f'found {len(candidates)}: {[p.name for p in candidates[:4]]}')
    return candidates[0]

def available_pairs(run_dir, month):
    def dates(v):
        found={}
        for p in (run_dir/v).glob(f'{v}_*.nc'):
            m=re.search(r'_(\d{4})(\d{2})\d{2}\d{4}-\d{12}\.nc$',p.name)
            if m and int(m.group(2))==month:
                found.setdefault(int(m.group(1)),[]).append(p)
        return found
    a,b=dates('rlut'),dates('pr')
    return sorted(y for y in a.keys()&b.keys() if len(a[y])==1 and len(b[y])==1)

def _require_safe_output(path):
    resolved=path.resolve()
    # Explicit never-write barriers for observational production.
    if str(resolved).startswith('/prj/global_water/MCS_5000km2_tables_v3'):
        raise ValueError('Refusing to write into OBS /prj production tree')
    if resolved.name!='CP4_base_tables_tir-prcp':
        raise ValueError('Output root directory must be named CP4_base_tables_tir-prcp')

def process_month(run,year,month,args):
    run_dir=args.input_root/run
    rlut_file=find_month_file(run_dir,'rlut',year,month)
    pr_file=find_month_file(run_dir,'pr',year,month)
    # Smoke tests always write to a dedicated directory so a partial month
    # can never be mistaken for a completed monthly production part.
    root = args.output_root / ('_smoke_tests' if args.max_hours is not None else '')
    # Avoid hive partition names (year=/month=) that duplicate the year/month
    # columns in OBS-compatible storm rows when reading a whole run.
    dest = (root / run / f'year_{year:04d}' /
            (f'month_{month:02d}_first_{args.max_hours:03d}h.parquet'
             if args.max_hours is not None else f'month_{month:02d}.parquet'))
    if dest.exists() and not args.overwrite:
        LOG.info('EXISTS - skipping %s',dest)
        return dest
    LOG.info('PROCESSING %s %04d-%02d',run,year,month)
    if args.dry_run:
        LOG.info('RLUT=%s PR=%s -> %s',rlut_file,pr_file,dest)
        return dest
    transformer=pyproj.Transformer.from_crs('EPSG:4326','EPSG:6933',always_xy=True)
    with xr.open_dataset(rlut_file,decode_times=True,use_cftime=True) as olrds, \
         xr.open_dataset(pr_file,decode_times=True,use_cftime=True) as prds:
        olr=olrds['rlut'];pr=prds['pr']
        if olr.shape!=pr.shape or olr.ndim!=3:
            raise ValueError('OLR and PR must have matching (time,lat,lon) dimensions')
        if not olr.time.equals(pr.time):
            raise ValueError('RLUT and PR timestamps differ; not using nearest time')
        if not np.allclose(olr.lat.values,pr.lat.values) or not np.allclose(olr.lon.values,pr.lon.values):
            raise ValueError('RLUT and PR grids do not match; interpolation not automatic')
        lat=olr.lat.values;lon=olr.lon.values
        if not (np.all(np.diff(lat)>0) or np.all(np.diff(lat)<0)):
            raise ValueError('Non-monotonic latitude; cannot calculate pixel areas')
        if not (np.all(np.diff(lon)>0) or np.all(np.diff(lon)<0)):
            raise ValueError('Non-monotonic longitude; check seam/grid')
        area=pixel_area_latlon_km2(lat,lon)
        all_rows=[]
        count=int(olr.sizes['time'])
        if args.max_hours is not None:
            count=min(count,args.max_hours)
        for it in range(count):
            t=olr.time.values[it]
            olr2=np.asarray(olr.isel(time=it).values,dtype=np.float32)
            pr2=np.asarray(pr.isel(time=it).values,dtype=np.float32)
            # CF _FillValue decoded by xarray, additional guard for bad raw data.
            olr2=np.where(np.abs(olr2)<1e19,olr2,np.nan)
            pr2=np.where(np.abs(pr2)<1e19,pr2,np.nan)
            tb_c=olr_to_bt_kelvin(olr2)-273.15
            rain=pr2*3600.0
            batch=extract_storms(tb_c,rain,lat,lon,area,t,run,transformer,
                       pr_file,rlut_file)
            all_rows.extend(batch)
            if (it+1)%24==0 or it+1==count:
                LOG.info('%s %04d-%02d: %d/%d hours, %d storms',run,year,month,it+1,count,len(all_rows))
    df=pd.DataFrame(all_rows,columns=SCHEMA)
    dest.parent.mkdir(parents=True,exist_ok=True)
    tmp=dest.with_name('part.tmp.parquet')
    df.to_parquet(tmp,index=False)
    os.replace(tmp,dest)
    LOG.info('SAVED %s (%d storm-hour rows)',dest,len(df))
    return dest

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root', type=Path, default=INPUT_ROOT)
    p.add_argument('--output-root', type=Path, default=OUTPUT_ROOT)
    p.add_argument('--run', action='append',
                   help='Exact simulation directory; repeatable. Defaults to all 10 runs')
    p.add_argument('--month', type=int, default=None,
                   help='Single month (default: August, unless --all)')
    p.add_argument('--year', type=int,
                   help='Explicit year; default is first common year unless --all/--all-years')
    p.add_argument('--all', action='store_true',
                   help='Process ALL available paired months and years for each run')
    p.add_argument('--all-years', action='store_true',
                   help='For the selected month, process every available paired year')
    p.add_argument('--max-hours', type=int, default=None,
                   help='SMOKE TEST: first N hours, saved separately under _smoke_tests')
    p.add_argument('--overwrite', action='store_true')
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--log-level', default='INFO')
    args = p.parse_args()

    logging.basicConfig(level=args.log_level,
                        format='%(asctime)s %(levelname)s %(message)s')
    if args.all and (args.year is not None or args.month is not None or args.all_years):
        p.error('--all cannot be combined with --year, --month or --all-years')
    if args.all_years and args.year is not None:
        p.error('--all-years cannot be combined with --year')
    if args.month is not None and not (1 <= args.month <= 12):
        p.error('--month must be 1..12')
    if args.year is not None and args.year < 1850:
        p.error('--year should be a model output year')
    if args.max_hours is not None and args.max_hours < 1:
        p.error('--max-hours must be positive')
    if args.all and args.max_hours is not None:
        p.error('--all cannot be combined with --max-hours (smoke tests)')
    _require_safe_output(args.output_root)

    runs = args.run or [f'{m}-{e}' for m in MODELS for e in ('HIST', 'FUT')]
    # In full production, do all available Augusts first so seasonal analysis
    # can begin before the remaining months have completed.
    months = ([8] + [m for m in range(1, 13) if m != 8]) if args.all \
             else [args.month if args.month is not None else 8]

    completed = []
    failures = []
    missing = []
    valid_runs = {f'{m}-{e}' for m in MODELS for e in ('HIST', 'FUT')}

    for run in runs:
        if run not in valid_runs:
            LOG.warning('Non-standard run name %s', run)
        run_dir = args.input_root / run
        if not run_dir.is_dir():
            LOG.warning('MISSING RUN DIRECTORY %s', run_dir)
            missing.append(run)
            continue
        for month in months:
            years = ([args.year] if args.year is not None
                     else available_pairs(run_dir, month))
            if not args.all and not args.all_years:
                years = years[:1]
            if not years:
                LOG.warning('No matching rlut/pr files for %s month %02d',
                            run, month)
                continue
            for year in years:
                try:
                    completed.append(process_month(run, year, month, args))
                except Exception as exc:
                    LOG.exception('FAILED %s %04d-%02d', run, year, month)
                    failures.append((run, year, month, str(exc)))
                    # Other months/runs remain independent. Completed month
                    # files are safe and automatically skipped on resubmission.
                    continue

    LOG.info('FINISHED: %d completed/skipped monthly files; '
             '%d failures; %d missing run directories',
             len(completed), len(failures), len(missing))
    if failures:
        for run, year, month, msg in failures:
            LOG.error('FAILED %s %04d-%02d: %s', run, year, month, msg)
        sys.exit(1)
    if missing:
        sys.exit(1)

if __name__ == '__main__':
    main()

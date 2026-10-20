#!/usr/bin/env python3
"""Extract small ERA5 subset using frozen source logic and compare to production ERA5 parquet.
WARNING: requires original project's GLOBAL/utils modules and ERA5 archive mounted.
"""
import argparse,sys
from pathlib import Path
import numpy as np
import pandas as pd
from tests.compare import compare_frames

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--region',default='GPlains')
    ap.add_argument('--year',type=int,default=2008)
    ap.add_argument('--storms',type=int,default=2)
    ap.add_argument('--variables',nargs='+',default=['tcwv','q650','ushear100_650'])
    ap.add_argument('--reference-root',default='/prj/global_water/MCS_5000km2_tables_v3/ERA5_storm_timeseries')
    ap.add_argument('--out-dir',default=str(Path(__file__).resolve().parents[1]/'outputs'/'tests'/'era5'))
    a=ap.parse_args()
    from reference.obs_production import COCOON_table_addERA_toMain as era
    # No writes to production: explicitly override output and never call production driver.
    out=Path(a.out_dir).resolve()
    refroot=Path(a.reference_root).resolve()
    if out==refroot or refroot in out.parents: raise ValueError('Unsafe output root')
    for var in a.variables:
        if var not in era.OUTPUT_VARIABLES: raise ValueError(f'Unknown variable: {var}')
    base=era.read_base_table(a.region,a.year)
    if base is None: raise FileNotFoundError('Missing observation base table')
    # Recreate ordering of the original ERA5 pipeline: static fields, land/sea filtering, then subset.
    base=era.add_static_fields_to_storms(base)
    base=era.filter_land_ocean_storms(base)
    base=base.head(a.storms)
    if base.empty: raise RuntimeError('No eligible storms')
    req=era.make_request_table(base,a.region)
    produced={v:[] for v in a.variables}
    for _,group in req.groupby('era_time'):
        result=era.extract_for_one_era_time(group.copy())
        if result is None: continue
        for v in a.variables: produced[v].append(result[v])
    out.mkdir(parents=True,exist_ok=True)
    for v in a.variables:
        if not produced[v]: raise RuntimeError(f'No ERA results for {v}')
        new=pd.concat(produced[v],ignore_index=True)
        if new.empty: raise RuntimeError(f'Empty ERA sample for {v}')
        newpath=out/f'{v}_{a.region}_{a.year}_subset.parquet'
        new.to_parquet(newpath,index=False)
        reference=refroot/v/f'region={a.region}'/f'year={a.year}'/'part.parquet'
        if not reference.is_file(): raise FileNotFoundError(f'Production ERA parquet not found: {reference}')
        old=pd.read_parquet(reference)
        n,nc=compare_frames(new,old,['storm_id','rel_hour'])
        print(f'PASS ERA {v}: {n} rows, {nc} columns. Test parquet: {newpath}')
if __name__=='__main__': main()

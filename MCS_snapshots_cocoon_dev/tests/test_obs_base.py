#!/usr/bin/env python3
"""Recreate a minimal observation scene and verify its rows against production base parquet.
The source scene is selected by checking the first input files for the year; no production write.
"""
import argparse,glob,sys
from pathlib import Path
import pandas as pd
import xarray as xr
from datasets.obs_feng_01deg import process_scene
from reference.obs_production import COCOON_table_main_01deg as obs
from tests.compare import compare_frames

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--region',default='GPlains')
    ap.add_argument('--year',type=int,default=2008)
    ap.add_argument('--max-files',type=int,default=8)
    ap.add_argument('--input-root',default=obs.LMCS)
    ap.add_argument('--reference-root',default=obs.OUTROOT)
    ap.add_argument('--out-dir',default=str(Path(__file__).resolve().parents[1]/'outputs'/'tests'))
    a=ap.parse_args()
    ref=Path(a.reference_root)/a.region/f'{a.year}_{obs.OUT_SUFFIX}'
    if not ref.is_file(): raise FileNotFoundError(f'Production base parquet not found: {ref}')
    candidates=sorted(glob.glob(str(Path(a.input_root)/f'{a.year}*'/'*.nc')))
    if not candidates: raise FileNotFoundError('No observation scene files found')
    produced=[]
    for f in candidates[:a.max_files]:
        print('Try',f,flush=True)
        with xr.open_dataset(f) as ds:
            df=process_scene(ds,a.region)
        if not df.empty: produced.append(df)
    if not produced: raise RuntimeError('No MCS found in selected scenes. Increase --max-files')
    new=pd.concat(produced,ignore_index=True)
    out=Path(a.out_dir).resolve()/(f'base_{a.region}_{a.year}_subset.parquet')
    out.parent.mkdir(parents=True,exist_ok=True)
    new.to_parquet(out,index=False)
    old=pd.read_parquet(ref)
    # stormID is per-image label; date ensures same comparison key.
    n,nc=compare_frames(new,old,['stormID','date'])
    print(f'PASS base: {n} storms, {nc} columns compared. Test parquet: {out}')
if __name__=='__main__': main()

#!/usr/bin/env python3
"""Run frozen observation production logic with forced non-production OUTROOT."""
import argparse
from pathlib import Path
from reference.obs_production import COCOON_table_main_01deg as obs
from reference.obs_production import COCOON_table_create_01deg as frozen_detector

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--region',required=True)
    ap.add_argument('--year',required=True,type=int)
    ap.add_argument('--out-root',default=str(Path(__file__).resolve().parents[1]/'outputs'/'obs_base'))
    a=ap.parse_args()
    out=Path(a.out_root).resolve()
    original=Path(obs.OUTROOT).resolve()
    if out == original or original in out.parents:
        raise ValueError('DEV script refuses to write under observation production root')
    obs.COCOON_table_create_01deg=frozen_detector  # force archived logic, not live production module
    obs.OUTROOT=str(out)
    obs.YEARS=[a.year]
    obs.make_table(a.region)
if __name__=='__main__': main()

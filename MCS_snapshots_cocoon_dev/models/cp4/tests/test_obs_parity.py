#!/usr/bin/env python3
"""Synthetic metric parity check against the actual OBS COCOON module.

Run from within user's existing workspace or set PYTHONPATH so Python imports
COCOON_table_create_01deg (or MCS_snapshots.COCOON_table_create_01deg).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import pyproj

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_cp4_base_tables import extract_storms, pixel_area_latlon_km2, OBS_COLUMNS
try:
    from MCS_snapshots import COCOON_table_create_01deg as obs
except ImportError:
    import COCOON_table_create_01deg as obs

lat=np.linspace(-2,3,130)
lon=np.linspace(-1,4,130)
y,x=np.mgrid[:130,:130]
cloud=(((y-70)/20)**2+((x-60)/29)**2)<1
cloud[65:77,85:97]=False
TB=np.full((130,130),285.0)
TB[cloud]=215.0
TB[70,60]=190.0
rain=np.zeros((130,130))
rain[cloud]=3
rain[68:75,60:65]=15
area=pixel_area_latlon_km2(lat,lon)
trans=pyproj.Transformer.from_crs('EPSG:4326','EPSG:6933',always_xy=True)
t=pd.Timestamp('2089-08-10 17:30')
ds=xr.Dataset({'tb':(('lat','lon'),TB),
               'precipitation':(('lat','lon'),rain)},
              coords={'lat':lat,'lon':lon,'time':t})
obs_dic=obs.process_tir_image(ds,4.5,area_grid_km2=area,
                              transformer=trans,rainvar_name='precipitation',
                              min_rain_pmax=1)
obs_res=obs.add_environment_toTable(obs_dic,ds.expand_dims('time'),4.5,
                                   rainvar_name='precipitation',area_grid_km2=area)
new=extract_storms(TB-273.15,rain,lat,lon,area,t,
                   'HadGEM3-GC31-MM-FUT',trans,'pr.nc','rlut.nc')
assert len(new)==len(obs_res['date'])==1
for key in OBS_COLUMNS:
    observed=obs_res[key][0]
    actual=new[0][key]
    if isinstance(observed,(float,int,np.number)):
        assert np.isclose(observed,actual,rtol=1e-5,atol=1e-5,equal_nan=True),(key,observed,actual)
    else:
        assert observed==actual,(key,observed,actual)
assert 'compactness' not in new[0]
print('PASS: all observational storm and precipitation metrics match on test storm.')

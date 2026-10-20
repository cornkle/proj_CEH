"""Observation adapter matching frozen COCOON processing behaviour."""
from reference.obs_production import COCOON_table_create_01deg as detector
from reference.obs_production import COCOON_table_main_01deg as prod

def process_scene(ds, region):
    import numpy as np
    import pandas as pd
    import pyproj
    box = prod.REGIONS[region][0]
    ds = ds.sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3])).load()
    if not ds.lat.size or not ds.lon.size:
        return pd.DataFrame()
    area = prod.pixel_area_latlon_km2(ds.lat.values, ds.lon.values)
    transformer = pyproj.Transformer.from_crs('EPSG:4326','EPSG:6933',always_xy=True)
    basic = detector.process_tir_image(ds, prod.DATA_RES_KM, t_thresh=prod.T_THRESH,
        min_mcs_size=prod.MIN_MCS_SIZE, area_grid_km2=area, transformer=transformer,
        rainvar_name='precipitation', min_rain_pmax=1)
    tab = detector.add_environment_toTable(basic,ds,prod.DATA_RES_KM,
        envvar_take=[],rainvar_name='precipitation',area_grid_km2=area)
    tab.pop('cloudMask',None)
    tab.pop('tir',None)
    if not tab['date']:
        return pd.DataFrame()
    df=pd.DataFrame.from_dict(tab)
    df['region']=region
    return df

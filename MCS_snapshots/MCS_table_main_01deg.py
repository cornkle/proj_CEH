from MCS_snapshots import MCS_table_create_01deg
import xarray as xr
import pandas as pd
import glob
import os
from shared.utils import constants as cnst
from shared.LMCS import glob_util
import numpy as np
import multiprocessing
import pyproj


REGIONS = glob_util.REGIONS

def pixel_area_latlon_km2(lat, lon):
    """
    Approximate pixel area for regular lat/lon grid.
    lat, lon are 1D arrays in degrees.
    Returns 2D area array in km2.
    """
    R = 6371.0  # km

    dlat = np.abs(np.gradient(lat))
    dlon = np.abs(np.gradient(lon))

    lat_rad = np.deg2rad(lat)
    dlat_rad = np.deg2rad(dlat)
    dlon_rad = np.deg2rad(dlon)

    area = (
        R**2
        * dlat_rad[:, None]
        * dlon_rad[None, :]
        * np.cos(lat_rad[:, None])
    )

    return area


def make_table(reg):
    """
    Start with scanning image for MCSs as defined in MCS_table_create.process_tir_image.
    :return:
    """
    lmcs = cnst.lmcs_drive + '/MCS_Feng/global_v2/2d_fields/'
    out = cnst.lmcs_drive + '/MCS_5000km2_tables_v2/'+reg+'/'
    box = REGIONS[reg][0] 

    for yy in range(2000,2021): # 2021
        infiles = sorted(glob.glob(lmcs + str(yy) + '*/*.nc'))
        full_year = []
        out_dic = {}
        print('Doing', yy)

        os.makedirs(out, exist_ok=True)
        outfile = out + str(yy)+'_MCS_5000km2_-40C_0.1degTIR-IMERG_hourly.csv'
        if os.path.isfile(outfile):
           print(outfile, ' exists, continue') 
           continue

        for infile in infiles:
            print('Doing', infile)
            try:
                da = (xr.open_dataset(infile))
            except:
                print('2d file open error, continue')
                continue
            da = da.sel(lon=slice(box[0],box[1]), lat=slice(box[2], box[3])) 

            area_grid_km2 = pixel_area_latlon_km2(da.lat.values, da.lon.values)

            transformer = pyproj.Transformer.from_crs(
                "EPSG:4326",
                "EPSG:6933",   # equal-area, metres
                always_xy=True)

            basic_tab = MCS_table_create_01deg.process_tir_image(da, 10, t_thresh=-40, min_mcs_size=5000, area_grid_km2=area_grid_km2, transformer=transformer)
            merge_tab = MCS_table_create_01deg.add_environment_toTable(basic_tab, da,  10, envvar_take=[], rainvar_name='precipitation', area_grid_km2=area_grid_km2)

            merge_tab.pop('cloudMask')
            merge_tab.pop('tir')
            if len(merge_tab['date']) ==0:
               continue

            full_year.append(merge_tab)
            print('Did' , infile)
        
        for key in full_year[0].keys():
            out_dic[key] = []
        for single_tab in full_year:
            for key in single_tab.keys():
                out_dic[key].extend(single_tab[key])
        
        pd_out = pd.DataFrame.from_dict(out_dic)
        pd_out.to_csv(outfile)
        del out_dic
        del pd_out
        del merge_tab
        del basic_tab
        del da


# for reg in REGIONS.keys():
#     print('Doing region', reg)
#     make_table(reg)

if __name__ == "__main__":
    with multiprocessing.Pool(processes=4) as pool:
        res = pool.map(make_table, list(REGIONS.keys()))
